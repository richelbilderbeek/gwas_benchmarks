process plink2 {
    tag "$prefix"
    label "plink2"
    publishDir "results/plink2/imputed", mode: 'copy'
    input:
        tuple val(prefix), path(genotypes), path(fam), path(to_include), path(pheno)
    output:
        tuple val(prefix), path("${prefix}.phenotype.glm.logistic"), emit:  plink_results
    script:
        """
        plink2 --bpfile "${prefix}" --fam "${fam}" --keep "${to_include}" \
            --out "${prefix}" --threads "${task.cpus}" \
            --ci 0.95 --logistic hide-covar cols=+a1freq,+ax --maf 0.01 \
            --hwe 1e-20 --geno 0.05 --pheno ${pheno} \
            --covar ${pheno} --pheno-name phenotype \
            --covar-name birth_weight breastfed maternal_smoking ever_smoked bmi
        """
}


process plink2_hardcalls {
    tag "$prefix"
    label "plink2"
    publishDir "results/plink2/hardcalls", mode: 'copy'
    input:
        tuple val(prefix), path(genotypes), path(fam), path(to_include), path(pheno)
    output:
        tuple val(prefix), path("${prefix}.phenotype.glm.logistic"), emit:  plink_results_hardcalls
    script:
        """
        # --bfile allows to drop imputed
        plink2 --bfile "${prefix}" --fam "${fam}" --keep "${to_include}" \
            --out "${prefix}" --threads "${task.cpus}" \
            --ci 0.95 --logistic hide-covar cols=+a1freq,+ax --maf 0.01 \
            --hwe 1e-20 --geno 0.05 --pheno ${pheno} \
            --covar ${pheno} --pheno-name phenotype \
            --covar-name birth_weight breastfed maternal_smoking ever_smoked bmi
        """
}