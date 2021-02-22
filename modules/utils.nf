nextflow.enable.dsl=2

process make_phenotypes {
    label "R"
    tag "phenotypes"
    publishDir "results/", mode: 'copy'

    input:
        path(phenotypes)
        path(to_include)
    output:
        path("phenotypes.txt"), emit: pheno
    script:
        """
        #!/usr/bin/env Rscript
        library(vroom)
        library(tidyverse)
        to_include <- read_delim("$to_include",
                                 col_names = FALSE,
                                 delim = " ",
                                 trim_ws = TRUE) %>% 
            mutate(ids = paste0("f.", X1, ".0.0"))
        all_phenotypes <- vroom("$phenotypes",
                                col_names = TRUE,
                                num_threads = $task.cpus,
                                col_select = c(f.eid,
                                               all_of(to_include\$ids))) %>% 
            rename(FID = f.eid) %>%
            mutate(IID = FID) %>%
            relocate(IID, .after = FID) %>%
            rename_at(vars(to_include\$ids), ~ to_include\$X2) %>%
            mutate(phenotype = case_when(is.na(age_astma_diagnosed) ~ 1,
                                         TRUE ~ 2)) %>%
            mutate(phenotype = as_factor(phenotype))
        
        write_tsv(all_phenotypes, path = "phenotypes.txt")
        """
}


process filter_cohort {
    label "plink2"
    publishDir "results/genotypes", mode: "copy"
    input:
        tuple val(prefix), path(genotypes), path(fam), path(to_include), path(phenotypes)
    output:
        tuple val(prefix), path("out/${prefix}.{bim,bed}"), path(fam)
    script:
        """
        mkdir out
        plink2 --threads "${task.cpus}" --bpfile "${prefix}" \
            --fam "${fam}" --keep "${to_include}" --maf 0.01 \
            --hwe 1e-20 --geno 0.05 --max-alleles 2 --make-bed \
            --out "out/${prefix}"
        """
}


process filter_hardcalls {
    label "plink2"
    publishDir "results/genotypes/hardcalls", mode: "copy"
    input:
        tuple val(prefix), path(genotypes), path(fam), path(to_include), path(phenotypes), path(hardcalls)
    output:
        path("out/${prefix}.{bim,bed}"), emit: genotypes_hardcalls_filtered
    script:
        """
        mkdir out
        plink2 --threads "${task.cpus}" --bfile "${prefix}" \
            --fam "${fam}" --keep "${to_include}" --maf 0.01 \
            --hwe 1e-20 --geno 0.05 --max-alleles 2 --make-bed \
            --extract "${hardcalls}" --out "out/${prefix}"
        """
}


process unpack_hard_calls {
    input:
        path(hardcalls)
    output:
        path("*.txt"), emit: hardcalls_list
    script:
        """
        tar xf "${hardcalls}"
        for bim in *bim
            do
                base=\$(basename "\${bim}" _v2.bim)
                prefix="\${base#ukb_snp_}"
                cut -f2 "\${bim}" > "\${prefix}".txt
            done
        """
}


process generate_exclude_non_biallelic {
    label "plink1"
    input:
        path(fam)
        path(genotypes_bim_bed)
    output:
        path("ukbb_all_chrs-merge.missnp"), emit: non_biallelic_to_exclude
    script:
        """
        for chr in {2..22}
            do
                echo "chr\${chr}.bed chr\${chr}.bim" "${fam}" >> list_beds.txt
            done
        trap 'if [[ \$? == 3 ]]; then echo OK; exit 0; fi' EXIT
        plink \
            --threads "${task.cpus}" \
            --bed chr1.bed \
            --bim chr1.bim \
            --fam "${fam}" \
            --merge-list list_beds.txt \
            --make-bed --out ukbb_all_chrs
        """
}


process merge_chromosomes {
    label "plink1"
    publishDir "results/genotypes/hardcalls", mode: "copy"
    input:
        path(fam)
        path(non_biallelic_to_exclude)
        path(genotypes_bim_bed)
    output:
        path("ukbb_all_chrs*.{bim,bed}"), emit: merged_chr
    script:
        """
        for chr in {2..22}
            do
                echo "chr\${chr}.bed chr\${chr}.bim" "${fam}" >> list_beds.txt
            done
        plink \
            --threads "${task.cpus}" \
            --bed chr1.bed \
            --bim chr1.bim \
            --fam "${fam}" \
            --merge-list list_beds.txt \
            --exclude ukbb_all_chrs-merge.missnp \
            --make-bed --out ukbb_all_chrs
        """
}