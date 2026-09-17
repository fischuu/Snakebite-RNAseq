include: "__functions__.smk"
include: "star.smk"

rule decontaminate:
    """Run the full sequential decontamination chain (index + map, every reference, every sample)

    A no-op when features["decontamination"] is empty - align then reads straight from the
    trimmed (fastp) reads instead, see get_final_read1/get_final_read2.
    """
    input:
        [
            DECONTAM_STAR / CONTAM_NAMES[-1] / f"{sample_id}_R{end}.fastq.gz"
            for sample_id in SAMPLES
            for end in (1, 2)
        ] if CONTAM_NAMES else [],
