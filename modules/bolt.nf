process bolt_lmm {
    label "bolt"
    tag "${prefix}"
    input:
        path(hardcalls)
        path(imputed)
        path(pheno)
    output:
        path("bolt_all_chr.*"), emit: bolt-results
    script:
        """
        for bgen in *bgen
            do
                prefix=\$(basename \$bgen .bgen)
                echo \$bgen \$prefix.sample >> list_bgen.txt
            done
        # --bgenFile chr{1:22}.bgen \
        # --sampleFile chr1.sample \
        bolt \
            --bed hardcalls_chr{1:22}.bed \
            --bim hardcalls_chr{1:22}.bim \
            --fam hardcalls_chr{1:22}.fam \
            --phenoFile "${pheno}" \
            --phenoCol standing_height \
            --covarFile "${pheno}" \
            --covarCol genetic_sex \
            --qCovarCol PC{1:15} \
            --LDscoresFile \$BOLTtables/LDSCORE.1000G_EUR.tab.gz \
            --geneticMapFile tables/genetic_map_hg19.txt.gz \
            --lmmForceNonInf \
            --numThreads ${task.cpus} \
            --statsFile bolt_all_chr.hardcalls.txt \
            --bgenSampleFileList list_bgen.txt
            --statsFileBgenSnps bolt_all_chr.bgen.txt \
            --verboseStats
        """
}