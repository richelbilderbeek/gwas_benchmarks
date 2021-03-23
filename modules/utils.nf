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
            drop_na()
        write_tsv(all_phenotypes, path = "phenotypes.txt")
        """
}


process filter_cohort {
    label "plink2"
    publishDir "results/", mode: "copy"
    input:
        tuple val(prefix), path(genotypes), path(fam), path(to_include), path(phenotypes)
    output:
        tuple val(prefix), path("genotypes/${prefix}.{bed,bim}"), path(fam)
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
        path("hardcalls/${prefix}.{bed,bim}"), emit: genotypes_hardcalls_filtered
    script:
        """
        mkdir hardcalls
        plink2 --threads "${task.cpus}" --bfile "${prefix}" \
            --fam "${fam}" --keep "${to_include}" --maf 0.01 \
            --hwe 1e-20 --geno 0.05 --max-alleles 2 --make-bed \
            --extract "${hardcalls}" --out "hardcalls/${prefix}"
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
        for chr in {2..22}
            do
                echo "chr\${chr}.pgen chr\${chr}.bim" "${fam}" >> list_beds.txt
            done
        plink2 --threads "${task.cpus}" \
            --bpfile chr1 \
            --fam "${fam}" \
            --pmerge-list list_beds.txt \
            --make-bed --out ukbb_all_chrs
        """
}
