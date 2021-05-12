process bolt_lmm {
    label "bolt"
    publishDir "results/bolt/", mode: 'copy'
    tag "all"
    input:
        path(hardcalls)
        path(imputed)
        path(sample)
        path(fam)
        path(pheno)
    output:
        path("bolt_all_chr.*"), emit: bolt_results
    script:
        """
        for bgen in *bgen
            do
                prefix=\$(basename \$bgen .bgen)
                echo \$bgen ${sample} >> list_bgen.txt
            done
        bolt \
            --bed chr{1:22}.bed \
            --bim chr{1:22}.bim \
            --fam "${fam}" \
            --phenoFile "${pheno}" \
            --phenoCol standing_height \
            --covarFile "${pheno}" \
            --covarCol genetic_sex \
            --qCovarCol PC{1:15} \
            --LDscoresFile \$BOLT_LMM_ROOT/tables/LDSCORE.1000G_EUR.tab.gz \
            --geneticMapFile \$BOLT_LMM_ROOT/tables/genetic_map_hg19_withX.txt.gz \
            --lmmForceNonInf \
            --numThreads ${task.cpus} \
            --statsFile bolt_all_chr.hardcalls.txt \
            --bgenSampleFileList list_bgen.txt \
            --statsFileBgenSnps bolt_all_chr.bgen.txt \
            --verboseStats
        """
}

process filter_sample {
    label "R"
    input:
        path(fam)
        path(sample)
    output:
        path("ukb_all_chr.fam"), emit: fam
    script:
        """
        #!/usr/bin/env Rscript
        library(tidyverse)
        fam_file <- read_delim("${fam}",
                        delim = "\t",
                        col_names = FALSE)
        sample_file <- read_delim("${sample}",
                                delim = " ")

        sample_file_ids <- sample_file\$ID_1
        to_keep <- fam_file %>%
            filter(`X1` %in% sample_file_ids)

        write_delim(to_keep,
                    "ukb_all_chr.fam",
                    delim = "\t",
                    col_names = FALSE)
        """
}