"""Stage 1 loaders: native deposited formats -> one raw-count AnnData per specimen and cohort.

No cells are filtered here except that cohort A's raw droplet matrices are restricted to the
author-annotated nuclei, as the assignment instructs. Cohorts B and C keep every deposited cell;
B cells carry an `in_author_annotation` flag. Counts are stored as float32 after checking they are
non-negative integers small enough to be exact in float32.

Cell key: ``<sample_id>_<source barcode>``. Donor key: ``<cohort>:<source donor id>``.
"""

from __future__ import annotations

import gzip
import shutil
import tarfile
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from opengo import config

FLOAT32_EXACT = 2**24  # integers up to this are exactly representable in float32


@dataclass
class SpecimenRecord:
    """One row of the Stage 1 reconciliation table."""

    cohort: str
    sample_id: str
    donor_key: str
    condition: str
    assay: str
    manifest_cells: int  # `paper_annotated_cells` in sample_manifest.tsv
    raw_barcodes: int  # barcodes in the deposited matrix (A: raw droplets)
    cells_loaded: int
    annotated_expected: int  # author-annotated cells for this specimen (A/B); NA for C
    annotated_matched: int  # of those, found in the matrix
    loaded_not_annotated: int  # loaded cells without an author annotation (B)
    n_genes: int
    integer_counts: bool
    min_count: float
    max_count: float
    total_counts: float


# --------------------------------------------------------------------------------------
# metadata


def sample_manifest() -> pd.DataFrame:
    m = pd.read_csv(config.SAMPLE_MANIFEST_PATH, sep="\t")
    m["donor_key"] = m["cohort"] + ":" + m["donor_id"].astype(str)
    return m


def population_lookup() -> tuple[dict[str, str], dict[tuple[str, str], str]]:
    """Suggested major populations: A keyed by broad label, B by (broad, detailed)."""
    lm = pd.read_csv(config.LABEL_MAPPING_PATH, sep="\t")
    a = lm[lm["cohort"] == "A"].set_index("broad_label")["population"].to_dict()
    b_rows = lm[lm["cohort"] == "B"]
    b = {(r.broad_label, r.detailed_label): r.population for r in b_rows.itertuples(index=False)}
    return a, b


def anchor_annotations() -> pd.DataFrame:
    return pd.read_csv(config.ANCHOR_ANNOTATIONS_PATH, sep="\t")


def cohort_b_annotations() -> pd.DataFrame:
    ann = pd.read_csv(config.COHORT_B_ANNOTATIONS_PATH, index_col=0)
    ann.index.name = "cell_id"
    return ann


# --------------------------------------------------------------------------------------
# helpers


def _count_checks(x: sp.spmatrix) -> tuple[bool, float, float, float]:
    data = x.data
    if data.size == 0:
        return True, 0.0, 0.0, 0.0
    is_int = bool(np.all(np.mod(data, 1) == 0))
    return is_int, float(data.min()), float(data.max()), float(data.sum(dtype=np.float64))


def _to_float32_csr(x: sp.spmatrix) -> sp.csr_matrix:
    x = sp.csr_matrix(x)
    if x.data.size and x.data.max() >= FLOAT32_EXACT:
        raise ValueError("counts too large to store exactly as float32")
    x.data = x.data.astype(np.float32)
    x.sum_duplicates()
    x.eliminate_zeros()
    return x


def _base_obs(row: pd.Series, barcodes: pd.Index) -> pd.DataFrame:
    obs = pd.DataFrame(
        {
            "barcode": barcodes.astype(str),
            "sample_id": row["sample_id"],
            "cohort": row["cohort"],
            "donor_id": str(row["donor_id"]),
            "donor_key": row["donor_key"],
            "condition": row["condition"],
            "assay": row["assay"],
        }
    )
    obs.index = (row["sample_id"] + "_" + obs["barcode"]).values
    obs.index.name = "cell_id"
    return obs


def specimen_path(cohort: str, sample_id: str) -> Path:
    return config.SPECIMEN_CKPT_DIR / cohort / f"{sample_id}.h5ad"


# --------------------------------------------------------------------------------------
# cohort A: 10x Matrix Market raw droplet matrices inside tar.gz archives


def load_a_specimen(row: pd.Series, annot: pd.DataFrame) -> tuple[ad.AnnData, SpecimenRecord]:
    """Restrict the raw droplet matrix to the author-annotated nuclei of this specimen."""
    sid = row["sample_id"]
    archive = config.PROJECT_ROOT / row["file"]
    config.TMP_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=config.TMP_DIR) as tmp:
        tmp = Path(tmp)
        with tarfile.open(archive, "r:gz") as tar:
            members = [
                m for m in tar.getmembers() if m.isfile() and "/raw_feature_bc_matrix/" in m.name
            ]
            names = sorted(Path(m.name).name for m in members)
            if names != ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
                raise ValueError(f"{sid}: unexpected archive contents {names}")
            for m in members:
                with tar.extractfile(m) as src, (tmp / Path(m.name).name).open("wb") as dst:
                    shutil.copyfileobj(src, dst)
        features = pd.read_csv(
            tmp / "features.tsv.gz",
            sep="\t",
            header=None,
            names=["gene_id", "symbol", "feature_type"],
        )
        barcodes = pd.read_csv(tmp / "barcodes.tsv.gz", header=None)[0]
        mat = scipy.io.mmread(tmp / "matrix.mtx.gz")  # genes x barcodes, COO

    if mat.shape != (len(features), len(barcodes)):
        raise ValueError(f"{sid}: matrix shape {mat.shape} != features/barcodes")
    mine = annot[annot["orig.ident"] == sid]
    idx = pd.Index(barcodes).get_indexer(mine["individual_barcode"])
    matched = idx >= 0
    mine, idx = mine[matched], idx[matched]

    x = sp.csc_matrix(mat)[:, idx].T  # annotated nuclei x genes
    del mat
    is_int, mn, mx, tot = _count_checks(x)
    x = _to_float32_csr(x)

    obs = _base_obs(row, pd.Index(mine["individual_barcode"]))
    pop_a, _ = population_lookup()
    obs["author_broad"] = mine["Cell_type.shorter"].values
    obs["author_detailed"] = mine["Cell_type.detailed"].values
    obs["author_group"] = mine["Cell_group_w_epithelialcelltypes"].values
    obs["author_cluster"] = mine["ident"].astype(str).values
    obs["author_umap_1"] = mine["UMAP_1"].values
    obs["author_umap_2"] = mine["UMAP_2"].values
    obs["population"] = obs["author_broad"].map(pop_a).fillna("Unmapped").values
    obs["in_author_annotation"] = True

    var = features.set_index("gene_id", drop=False)
    var.index.name = None
    adata = ad.AnnData(X=x, obs=obs, var=var)
    rec = SpecimenRecord(
        cohort="A",
        sample_id=sid,
        donor_key=row["donor_key"],
        condition=row["condition"],
        assay=row["assay"],
        manifest_cells=int(row["paper_annotated_cells"]),
        raw_barcodes=len(barcodes),
        cells_loaded=adata.n_obs,
        annotated_expected=int(matched.size),
        annotated_matched=int(matched.sum()),
        loaded_not_annotated=0,
        n_genes=adata.n_vars,
        integer_counts=is_int,
        min_count=mn,
        max_count=mx,
        total_counts=tot,
    )
    return adata, rec


# --------------------------------------------------------------------------------------
# cohort B: dense gene x cell count CSVs inside one tar


def _read_count_csv(fh, chunk_genes: int = 2000) -> tuple[sp.csr_matrix, pd.Index, pd.Index]:
    """Stream a gzip CSV (genes in rows, cells in columns) into a sparse genes x cells matrix."""
    blocks, genes, cells = [], [], None
    with gzip.open(fh, "rt") as text:
        for chunk in pd.read_csv(text, index_col=0, chunksize=chunk_genes):
            if cells is None:
                cells = chunk.columns
            blocks.append(sp.csr_matrix(chunk.to_numpy()))
            genes.extend(chunk.index.astype(str))
    return sp.vstack(blocks, format="csr"), pd.Index(genes), pd.Index(cells)


def load_b_specimen(
    row: pd.Series, tar: tarfile.TarFile, annot: pd.DataFrame
) -> tuple[ad.AnnData, SpecimenRecord]:
    sid = row["sample_id"]
    with tar.extractfile(row["archive_member"]) as fh:
        mat, genes, cells = _read_count_csv(fh)
    if genes.has_duplicates:
        raise ValueError(f"{sid}: duplicated gene symbols in CSV rows")

    # CSV headers are expected to be '<specimen>_<barcode>'; prefix if not.
    cells = pd.Index([c if c.startswith(f"{sid}_") else f"{sid}_{c}" for c in cells])
    x = mat.T.tocsr()
    is_int, mn, mx, tot = _count_checks(x)
    x = _to_float32_csr(x)

    barcodes = pd.Index([c[len(sid) + 1 :] for c in cells])
    obs = _base_obs(row, barcodes)
    mine = annot[annot.index.str.startswith(f"{sid}_")]
    joined = mine.reindex(obs.index)
    in_ann = obs.index.isin(mine.index)
    _, pop_b = population_lookup()
    obs["author_broad"] = joined["ano.l1"].fillna("NA").values
    obs["author_detailed"] = joined["ano.l2"].fillna("NA").values
    for src, dst in [
        ("TotalUMI", "author_total_umi"),
        ("gene.num", "author_n_genes"),
        ("MT.ratio", "author_mt_ratio"),
        ("doublest.score", "author_doublet_score"),
    ]:
        obs[dst] = joined[src].values
    obs["population"] = [
        pop_b.get((b, d), "Unmapped") if ok else "Not_in_author_annotation"
        for b, d, ok in zip(obs["author_broad"], obs["author_detailed"], in_ann, strict=True)
    ]
    obs["in_author_annotation"] = in_ann

    var = pd.DataFrame({"symbol": genes}, index=genes)  # Ensembl IDs assigned by genes.py
    adata = ad.AnnData(X=x, obs=obs, var=var)
    rec = SpecimenRecord(
        cohort="B",
        sample_id=sid,
        donor_key=row["donor_key"],
        condition=row["condition"],
        assay=row["assay"],
        manifest_cells=int(row["paper_annotated_cells"]),
        raw_barcodes=len(cells),
        cells_loaded=adata.n_obs,
        annotated_expected=len(mine),
        annotated_matched=int(in_ann.sum()),
        loaded_not_annotated=int((~in_ann).sum()),
        n_genes=adata.n_vars,
        integer_counts=is_int,
        min_count=mn,
        max_count=mx,
        total_counts=tot,
    )
    return adata, rec


# --------------------------------------------------------------------------------------
# cohort C: legacy 10x HDF5


def load_c_specimen(row: pd.Series) -> tuple[ad.AnnData, SpecimenRecord]:
    sid = row["sample_id"]
    with h5py.File(config.PROJECT_ROOT / row["file"], "r") as f:
        g = f["GRCh38"]
        shape = tuple(int(v) for v in g["shape"][()])
        mat = sp.csc_matrix((g["data"][()], g["indices"][()], g["indptr"][()]), shape=shape)
        gene_ids = [v.decode() for v in g["genes"][()]]
        symbols = [v.decode() for v in g["gene_names"][()]]
        barcodes = pd.Index([v.decode() for v in g["barcodes"][()]])
    x = mat.T.tocsr()  # cells x genes
    is_int, mn, mx, tot = _count_checks(x)
    x = _to_float32_csr(x)

    obs = _base_obs(row, barcodes)
    obs["population"] = "Unlabelled"
    obs["in_author_annotation"] = False
    var = pd.DataFrame({"gene_id": gene_ids, "symbol": symbols}, index=gene_ids)
    adata = ad.AnnData(X=x, obs=obs, var=var)
    rec = SpecimenRecord(
        cohort="C",
        sample_id=sid,
        donor_key=row["donor_key"],
        condition=row["condition"],
        assay=row["assay"],
        manifest_cells=int(row["paper_annotated_cells"]),
        raw_barcodes=len(barcodes),
        cells_loaded=adata.n_obs,
        annotated_expected=0,
        annotated_matched=0,
        loaded_not_annotated=0,
        n_genes=adata.n_vars,
        integer_counts=is_int,
        min_count=mn,
        max_count=mx,
        total_counts=tot,
    )
    return adata, rec


# --------------------------------------------------------------------------------------
# cohort driver


def _write_specimen(adata: ad.AnnData) -> Path:
    cohort, sid = adata.obs["cohort"].iloc[0], adata.obs["sample_id"].iloc[0]
    path = specimen_path(cohort, sid)
    path.parent.mkdir(parents=True, exist_ok=True)
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].astype(str)
    adata.write_h5ad(path, compression="gzip")
    return path


def load_cohort(cohort: str, progress=None) -> pd.DataFrame:
    """Load every specimen of a cohort, write per-specimen + concatenated checkpoints.

    Returns the reconciliation table (one row per specimen).
    """
    manifest = sample_manifest()
    rows = manifest[manifest["cohort"] == cohort].reset_index(drop=True)
    records, paths, var_ref = [], [], None

    def handle(adata: ad.AnnData, rec: SpecimenRecord) -> None:
        nonlocal var_ref
        if var_ref is None:
            var_ref = adata.var_names
        elif not adata.var_names.equals(var_ref):
            raise ValueError(f"{rec.sample_id}: gene set differs from first specimen in cohort")
        paths.append(_write_specimen(adata))
        records.append(asdict(rec))
        if progress:
            progress(rec)

    if cohort == "A":
        annot = anchor_annotations()
        for _, row in rows.iterrows():
            handle(*load_a_specimen(row, annot))
    elif cohort == "B":
        annot = cohort_b_annotations()
        with tarfile.open(config.PROJECT_ROOT / rows["file"].iloc[0]) as tar:
            for _, row in rows.iterrows():
                handle(*load_b_specimen(row, tar, annot))
    elif cohort == "C":
        for _, row in rows.iterrows():
            handle(*load_c_specimen(row))
    else:
        raise ValueError(cohort)

    out = config.cohort_raw_path(cohort)
    ad.experimental.concat_on_disk(paths, out, join="inner", merge="same")
    shutil.rmtree(config.TMP_DIR, ignore_errors=True)
    return pd.DataFrame(records)
