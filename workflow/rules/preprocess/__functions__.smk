def get_concatenation_inputs_read1(wildcards):
    """Every raw lane's read1 belonging to one sample_id, to be concatenated together"""
    rows = samples.query("sample_id == @wildcards.sample_id")
    return [row.read1 for row in rows.itertuples()]

def get_concatenation_inputs_read2(wildcards):
    """Every raw lane's read2 belonging to one sample_id, to be concatenated together"""
    rows = samples.query("sample_id == @wildcards.sample_id")
    return [row.read2 for row in rows.itertuples()]
