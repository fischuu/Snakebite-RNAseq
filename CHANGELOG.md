################################################################################
#####                                                                      #####
#####                         CHANGELOG                                    #####
#####                                                                      #####
################################################################################

#### PLANNED UPDATES FOR FUTURE VERSIONS
*  rule build_contamination_indices:
    - Output file is created as soon as one index is created. However, it would
        be good to report any index that was created successfully and then only
        create the output file in case that all indices were created

0.1.*: UPDATES
--------------------------------------------------------------------------------
* 15: Config options better explained
* 14: Output file dependencied from STAR alignments for metatranscriptomics further adjusted
* 13: Adjusted the output files from STAR alignments for metatranscriptomics
* 12: Added the 'mode' option to account for metagenomic studies
* 11: Option limitGenomeGenerateRAM added to the STAR genome builder
* 10: Added the genomeChrBinNbits option to STAR, to account for larger genomes
* 9: Bugfix, -t option in rule featureCounts_quantify_merged
* 8: Added to STAR index the parameter genomeSAindexNbases
* 7: FeatureCounts merged attribute types fixed
* 6: Bugfix, filenames
* 5: Attribute type and feature type added for featureCounts
* 4: Decontamination works now properly
* 3: Contamination statistics added to the final report
* 2: Added the decontamination step
* 1: Inital submission
