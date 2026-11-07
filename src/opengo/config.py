from pathlib import Path

# Base Paths
PROJECT_ROOT = Path(__file__).resolve().parents[2]  # src/opengo/config.py -> project root
DATA_DIR = PROJECT_ROOT / "data"
METADATA_DIR = PROJECT_ROOT / "metadata"
OUTPUTS_DIR = PROJECT_ROOT / "outputs"
FIGURES_DIR = OUTPUTS_DIR / "figures"
TABLES_DIR = OUTPUTS_DIR / "tables"
CHECKPOINTS_DIR = OUTPUTS_DIR / "checkpoints"

# Specific Data Paths
COHORT_A_DIR = DATA_DIR / "rna" / "cohort_A_GSE227898"
COHORT_B_TAR = DATA_DIR / "rna" / "cohort_B_GSE178481" / "GSE178481_RAW.tar"
COHORT_C_DIR = DATA_DIR / "rna" / "cohort_C_GSE159115"
BULK_DIR = DATA_DIR / "bulk"

# Metadata Paths
SAMPLE_MANIFEST_PATH = METADATA_DIR / "sample_manifest.tsv"
MULTIMODAL_CROSSWALK_PATH = METADATA_DIR / "multimodal_crosswalk.tsv"
LABEL_MAPPING_PATH = METADATA_DIR / "label_mapping.tsv"
ECM_PROGRAM_PATH = METADATA_DIR / "ecm_program.tsv"
BULK_COLUMN_MANIFEST_PATH = METADATA_DIR / "bulk_column_manifest.tsv"
BULK_SAMPLES_PATH = METADATA_DIR / "bulk_samples.tsv"
ORIGINAL_METADATA_DIR = METADATA_DIR / "original"
ANCHOR_ANNOTATIONS_PATH = ORIGINAL_METADATA_DIR / "anchor_annotations.tsv.gz"
ANCHOR_SAMPLES_PATH = ORIGINAL_METADATA_DIR / "anchor_samples.tsv"
COHORT_B_ANNOTATIONS_PATH = ORIGINAL_METADATA_DIR / "cohort_b_annotations.csv"
COHORT_C_SAMPLES_PATH = ORIGINAL_METADATA_DIR / "cohort_C_samples.tsv"

# ECM 10-gene Panel
ECM_GENES = [
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "COL5A1",
    "COL6A1",
    "COL6A2",
    "FN1",
    "SPARC",
    "DCN",
    "LUM",
]

# Random Seed for Reproducibility
RANDOM_SEED = 42
