nextflow.enable.dsl=2

process make_phenotypes {
    label "R"

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
                                num_threads = $params.cpus,
                                col_select = c(f.eid,
                                               all_of(to_include\$ids))) %>% 
            rename(FID = f.eid) %>% 
            rename_at(vars(to_include\$ids), ~ to_include\$X2) %>%
            mutate(phenotype = case_when(is.na(age_astma_diagnosed) ~ "control",
                                         TRUE ~ "asthma")) %>%
            mutate(phenotype = as_factor(phenotype))
        
        write_tsv(all_phenotypes, path = "phenotypes.txt")
        """
}


process remove_not_biallelic {
    label "plink2"

    input:
        path(fam)
        path(genotypes_bim_bed)
    output:
        tuple val(prefix), path(fam), path("out/${prefix}*"), emit: genotypes_biallelic
    script:
        """
        mkdir out
        plink2 --threads "${task.cpus}" --bfile "${prefix}" \
            --fam "${fam}" --max-alleles 2 --make-bed --out "out/${prefix}"
        """
}


process merge_chromosomes {
    label "plink1"

    input:
        tuple val(prefix), path(fam), path(genotypes_bim_bed)
    output:
        tuple path("ukbb__all_chrs.bed"), path("ukbb__all_chrs.bim"), emit: geno_all_chrs
    script:
        """
        for chr in {2..22}
            do
                echo "\${chr}.bed \${chr}.bim" "${fam}" >> list_beds.txt
            done

        plink \
            --threads "${task.cpus}" \
            --bed chr1.bed \
            --bim chr1.bim \
            --fam "${fam}"
            --merge-list list_beds.txt \
            --make-bed --out ukbb_all_chrs
        """
}


process plink_qc {
    input:
    output:
    script:
        """
        plink2 \
            --bfile ukb_cal_allChrs \
            --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
            --mind 0.1 \
            --write-snplist --write-samples --no-id-header \
            --out qc_pass
        """
}