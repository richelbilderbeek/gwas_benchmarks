process saige_null_fitting {
    label "saige"
    input:
        tuple val(hardcalls), path(hardcalls_files), path(pheno)
    output:
        path("null_model.rda"), emit: null_model
        path("null_model.varianceRatio.txt"), emit: variance_ratio
    script:
        """
        step1_fitNULLGLMM.R \
            --plinkFile ${hardcalls} \
            --phenoFile "${pheno}" \
            --phenoCol standing_height \
            --sampleIDColinphenoFile=IID \  
            --covarColList genetic_sex \
            --traitType quantitative \
            --invNormalize TRUE \
            --nThreads "${task.cpus}" \
            --LOCO FALSE \
            --tauInit 1,0 \
            --outputPrefix "null_model"
        """
}


process saige {
    tag "$prefix"
    label "saige"
    publishDir "results/saige/", mode: 'copy'
    input:
        path(null_model)
        path(variance_ratio)
        tuple val(prefix), path(genotypes)
    output:
        tuple val(prefix), path("${prefix}.out.txt"), emit: saige_results
    script:
        """
        # saige wants a headerless .sample file
        sed -e 1,2d < ${prefix}.sample > ${prefix}_headerless.sample
        step2_SPAtests.R  \
            --bgenFile ${prefix}.bgen \
            --sampleFile ${prefix}_headerless.sample \
            --chrom ${prefix} \
            --GMMATmodelFile ${null_model} \
            --varianceRatioFile ${variance_ratio} \
            --SAIGEOutputFile ${prefix}.out.txt \
            --numLinesOutput 2 \
            --LOCO FALSE
        """
}