nextflow.enable.dsl=2

process plink2 {
    tag "$prefix"

    input:
        tuple val(prefix), file(genotypes), file(fam), file(to_include)
    output:
        tuple val(prefix), file("${prefix}.raw.gz"), emit: genotypes_hard_calls_only
    script:
        """
        plink2 --bpfile "${prefix}" --fam "${fam}" --keep "${to_include}" \
            --extract --out "${prefix}" \
            --ci 0.95 --glm hide-covar cols=+a1freq,+ax --maf 0.01 \
            --hwe 1e-20 --geno 0.05 --pheno data/phenotypes_to_include.txt \
            --pheno-name 
        """
}