def _reads_sample_row(wildcards):
    """Look up the samples.tsv row matching a sample_id/library_id wildcard pair"""
    matches = samples.query(
        "sample_id == @wildcards.sample_id and library_id == @wildcards.library_id"
    )
    return matches.iloc[0]

def get_raw_read1(wildcards):
    """Path to the raw read1 FASTQ listed in samples.tsv for this sample_id/library_id"""
    return _reads_sample_row(wildcards)["read1"]

def get_raw_read2(wildcards):
    """Path to the raw read2 FASTQ listed in samples.tsv for this sample_id/library_id"""
    return _reads_sample_row(wildcards)["read2"]
