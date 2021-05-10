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
        
        # add PCs
        pcs <- paste0("f.22009.0.", seq(1, 15))
        pcs_name <- paste0("PC", seq(1, 15))
        pcs_tibble <- tibble(X1 = rep(22009, 15), X2 = pcs_name, ids = pcs)
        to_include <- to_include %>% bind_rows(pcs_tibble)

        all_phenotypes <- vroom("$phenotypes",
                                col_names = TRUE,
                                num_threads = $task.cpus,
                                col_select = c(f.eid,
                                               all_of(to_include\$ids))) %>% 
            rename(FID = f.eid) %>%
            mutate(IID = FID) %>%
            relocate(IID, .after = FID) %>%
            rename_at(vars(to_include\$ids), ~ to_include\$X2) %>%
            drop_na()
        write_tsv(all_phenotypes, path = "phenotypes.txt")
        """
}


process make_train_test {
    label "R"
    tag "ids"
    publishDir "results/", mode: "copy"
    input:
        path(phenotypes)
    output:
        path("train_ids_to_include.txt"), emit: train_ids_to_include
        path("test_ids_to_include.txt"), emit: test_ids_to_include
        path("phenotypes_split.rds"), emit: splits_rds
    script:
        """
        #!/usr/bin/env Rscript
        library(tidyverse)
        library(tidymodels)

        set.seed(1)
        phenotypes <- read_tsv("$phenotypes")
        splits <- initial_split(phenotypes,
                                strata = "genetic_sex")
        train_data <- training(splits)
        test_data <- testing(splits)
        train_ids_to_include <- train_data %>%
            select(FID, IID)
        test_ids_to_include <- test_data %>%
            select(FID, IID)
        
        write_tsv(train_ids_to_include, path = "train_ids_to_include.txt")
        write_tsv(test_ids_to_include, path = "test_ids_to_include.txt")
        saveRDS(splits, "phenotypes_split.rds")
        """
}


process filter_cohort {
    label "plink2"
    publishDir "results/", mode: "copy"
    input:
        tuple val(prefix), path(genotypes), path(fam), path(to_include), path(phenotypes)
    output:
        tuple val(prefix), path("genotypes/${prefix}.{bed,bim,fam}"), emit: genotypes_filtered
    script:
        """
        mkdir genotypes
        plink2 --threads "${task.cpus}" --bpfile "${prefix}" \
            --fam "${fam}" --keep "${to_include}" --maf 0.01 \
            --hwe 1e-20 --geno 0.05 --max-alleles 2 --make-bed \
            --out "genotypes/${prefix}"
        """
}


process filter_hardcalls {
    label "plink2"
    publishDir "results/genotypes/", mode: "copy"
    input:
        tuple val(prefix), path(genotypes), path(fam), path(to_include), path(phenotypes), path(hardcalls)
    output:
        path("hardcalls/${prefix}.{bed,bim,fam}"), emit: genotypes_hardcalls_filtered
    script:
        """
        mkdir hardcalls
        plink2 --threads "${task.cpus}" --bfile "${prefix}" \
            --fam "${fam}" --keep "${to_include}" --maf 0.01 \
            --hwe 1e-20 --geno 0.05 --max-alleles 2 --make-bed \
            --extract "${hardcalls}" --out "hardcalls/${prefix}"
        """
}


process make_bgen {
    label "plink2"
    publishDir "results/genotypes", mode: "copy"
    input:
        tuple val(prefix), path(genotypes)
    output:
        path("${prefix}.{bgen,sample}"), emit: bgen_full
    script:
        """
        plink2 --threads "${task.cpus}" --bfile "${prefix}" \
            --export bgen-1.2 bits=8 --out "${prefix}"
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


process merge_chromosomes {
    label "plink_dev"
    input:
        path(fam)
        path(genotypes)
    output:
        path("ukbb_merged.*"), emit: ukbb_merged
    script:
        """
        # this is a placeholder containing the correct commands, but this step 
        # is not implemented in the pipeline yet
        # TODO automate the merging. it will need splitting into several
        # processes
        # 1. generate the exclude missnp list. Plink errors at this step
        # the error code will need to be ignored (with a trap?)
        for chr in {2..22}
            do
                echo "chr\${chr}" >> temp_list_beds_to_exclude.txt
            done
        plink --threads "${task.cpus}" \
            --bfile chr1 \
            --merge-list temp_list_beds_to_exclude.txt \
            --make-bed --out to_exclude
        
        # 2. exclude the multiallelic snps from all chromosomes
        for chr in {1..22}
            do
                plink --threads "${task.cpus}" \
                --bfile chr\${chr} \
                --exclude to_exclude-merge.missnp \
                --make-bed \
                --out chr\${chr}_dedup
            done
        
        # 3. merge the files
        for chr in {2..22}
            do
                echo "chr\${chr}_dedup" >> list_beds.txt
            done
        plink --threads "${task.cpus}" \
            --bfile chr1_dedup \
            --merge-list list_beds.txt \
            --make-bed \
            --out hardcalls_merged
        """
}
