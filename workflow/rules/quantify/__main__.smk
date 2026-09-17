include: "featurecounts.smk"
include: "rsem.smk"

rule quantify:
    """Quantify every sample's aligned reads against the reference annotation and the merged
    StringTie assembly (featureCounts)

    RSEM's quantify__rsem__calculate_expression is kept available but not part of this
    default target, mirroring reference__rsem_index__build.
    """
    input:
        expand(str(QUANT_FEATURECOUNTS_REFERENCE / "{sample_id}_fc.txt"), sample_id=SAMPLES),
        expand(str(QUANT_FEATURECOUNTS_STRINGTIE / "{sample_id}_fc.txt"), sample_id=SAMPLES),
