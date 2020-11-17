nextflow.enable.dsl=2

process plink2 {
    tag "$prefix"
    label "plink2"

    input:
        tuple val(prefix), path(genotypes), path(fam), path(to_include)
    output:
        tuple val(prefix), path("${prefix}.raw.gz"), emit: genotypes_hard_calls_only
    script:
        """
        plink2 --bpfile "${prefix}" --fam "${fam}" --keep "${to_include}" \
            --extract --out "${prefix}" --threads "${task.cpus}"
            --ci 0.95 --glm hide-covar cols=+a1freq,+ax --maf 0.01 \
            --hwe 1e-20 --geno 0.05 --pheno phenotypes_to_include.txt \
            --covar phenotypes_to_include.txt --pheno-name phenotype \
            --covar-name birth_weight breastfed maternal_smoking ever_smoked bmi
        """
}