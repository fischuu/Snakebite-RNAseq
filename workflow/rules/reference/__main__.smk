include: "star_index.smk"
include: "rsem_index.smk"

rule reference:
    """Build the STAR genome index used by the pipeline

    RSEM's index (reference__rsem_index__build) is kept available but is not part of the
    default target - build it on demand if you need quantify__rsem__calculate_expression.
    """
    input:
        rules.reference__star_index__genome.output,
