# Generate more configuration keys
SCRIPT_FOLDER = os.path.join(config["pipeline_folder"], "workflow", "scripts")
WD = os.getcwd()

# Guaranteed-available scratch fallback for tools that need a writable tmp directory:
# nvme_storage/tmp_storage (config.yaml) may point at per-job cluster paths that aren't
# necessarily visible inside a container, so rules that use this probe those first and fall
# back to this project-local directory, which always exists and is on the same storage as
# results/.
PROJECT_TMP = Path("tmp/")

# reads
READS = Path("results/reads/")

# reference
REFERENCE = Path("results/reference/")
STAR_INDEX_GENOME = REFERENCE / "star_index/genome"
RSEM_INDEX = REFERENCE / "rsem_index"

# preprocess
PRE = Path("results/preprocess/")
CONCATENATED = PRE / "concatenate/"
TRIMMED = PRE / "fastp/"

# decontaminate
DECONTAM = Path("results/decontaminate/")
DECONTAM_INDEX = DECONTAM / "star_index/"
DECONTAM_STAR = DECONTAM / "star/"

# align
ALIGN = Path("results/align/")
GENOME_STAR = ALIGN / "star/genome/"

# assemble
ASSEMBLE = Path("results/assemble/")
STRINGTIE = ASSEMBLE / "stringtie/"

# quantify
QUANT = Path("results/quantify/")
QUANT_FEATURECOUNTS_REFERENCE = QUANT / "featurecounts/reference/"
QUANT_FEATURECOUNTS_STRINGTIE = QUANT / "featurecounts/stringtie/"
QUANT_RSEM = QUANT / "rsem/"

# reporting
PIPELINE_REPORT = Path("Rreports/")
