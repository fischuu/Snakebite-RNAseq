# Decontamination is a chain: each reference's input is the previous reference's unmapped
# ("still contaminated-free") output, and the very first reference reads straight from the
# trimmed (fastp) reads. CONTAM_NAMES (from features["decontamination"], order-preserving) is
# defined in the Snakefile.

def _decontam_chain_input(wildcards, end):
    """Reads feeding into the mapping step for one reference in the decontamination chain"""
    assert end in (1, 2)
    step = CONTAM_NAMES.index(wildcards.contam)
    if step == 0:
        return str(TRIMMED / f"{wildcards.sample_id}_R{end}.trimmed.fastq.gz")
    previous = CONTAM_NAMES[step - 1]
    return str(DECONTAM_STAR / previous / f"{wildcards.sample_id}_R{end}.fastq.gz")

def get_decontam_chain_r1(wildcards):
    """Read1 feeding into one decontamination-mapping step"""
    return _decontam_chain_input(wildcards, end=1)

def get_decontam_chain_r2(wildcards):
    """Read2 feeding into one decontamination-mapping step"""
    return _decontam_chain_input(wildcards, end=2)

# The final, fully decontaminated (or trimmed-only, if no references configured) reads that
# feed the align module.
def _final_preprocessed_read(wildcards, end):
    assert end in (1, 2)
    if not CONTAM_NAMES:
        return str(TRIMMED / f"{wildcards.sample_id}_R{end}.trimmed.fastq.gz")
    return str(DECONTAM_STAR / CONTAM_NAMES[-1] / f"{wildcards.sample_id}_R{end}.fastq.gz")

def get_final_read1(wildcards):
    """Read1 after the full preprocessing/decontamination chain"""
    return _final_preprocessed_read(wildcards, end=1)

def get_final_read2(wildcards):
    """Read2 after the full preprocessing/decontamination chain"""
    return _final_preprocessed_read(wildcards, end=2)
