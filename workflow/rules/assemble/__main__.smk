include: "stringtie.smk"

rule assemble:
    """Assemble and merge transcripts for every sample (StringTie)"""
    input:
        rules.assemble__stringtie__merge.output,
