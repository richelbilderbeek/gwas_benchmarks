#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { filter_sample; bolt_lmm } from './modules/bolt.nf'
include { plink2 } from './modules/plink.nf'
include { saige_null_fitting; saige_assoc } from './modules/saige.nf'
include { regenie_step_1; regenie_step_2 } from './modules/regenie.nf'
include { make_phenotypes; make_train_test } from './modules/utils.nf'

Channel
    .fromPath(params.fam)
    .set{ fam }
Channel
    .fromPath(params.ids_to_include)
    .set{ ids_to_include }
Channel
    .fromPath(params.phenotypes_to_include)
    .set{ phenotypes_to_include }
Channel
    .fromPath(params.phenotypes_file)
    .set{ phenotypes_file }


workflow prep {
    make_phenotypes(phenotypes_file, phenotypes_to_include)
    make_train_test(make_phenotypes.out.pheno)
}


workflow plink {
    Channel
        .fromPath(params.ids_to_include_train)
        .set{ids_to_include_train}
    Channel
        .fromPath(params.phenotypes_filtered)
        .set{phenotypes}
    Channel
        .fromPath(params.fam)
        .set{fam}
    Channel
        .fromPath(params.genotypes)
        .map { file -> tuple(file.baseName, file) }
        .groupTuple(by:0)
        .combine(fam)
        .combine(ids_to_include_train)
        .combine(phenotypes)
        .dump()
        .set{ plink2_input }
    plink2(plink2_input)
}


workflow saige {
    Channel
        .fromPath(params.fam)
        .set{fam}
    Channel
        .fromPath(params.phenotypes_filtered)
        .set{phenotypes}
    Channel
        .fromPath(params.hardcalls_merged)
        .collect()
        .combine(fam)
        .combine(phenotypes)
        .dump()
        .set{ saige_input_hardcalls }

    saige_null_fitting(saige_input_hardcalls)

    Channel
        .fromPath(params.bgen_sample)
        .set{sample}
    Channel
        .fromPath(params.genotypes_bgen)
        .map { file -> tuple(file.baseName, file) }
        .groupTuple(by:0)
        .combine(sample)
        .combine(saige_null_fitting.out.null_model)
        .combine(saige_null_fitting.out.variance_ratio)
        .dump()
        .set{ saige_input_assoc }

    saige_assoc(saige_input_assoc)
}


workflow bolt {
    Channel
        .fromPath(params.fam)
        .set{fam}
    Channel
        .fromPath(params.phenotypes_filtered)
        .set{phenotypes}
    Channel
        .fromPath(params.hardcalls)
        .collect()
        .set{hardcalls}
    Channel
        .fromPath(params.genotypes_bgen)
        .collect()
        .set{imputed}
    Channel
        .fromPath(params.bgen_sample)
        .set{sample}

    filter_sample(fam, sample)
    bolt_lmm(hardcalls, imputed, sample, filter_sample.out.fam, phenotypes)
}


workflow regenie {
    make_phenotypes(phenotypes_file, phenotypes_to_include)
    Channel
        .fromPath(params.hardcalls_merged)
        .collect()
        .combine(make_phenotypes.out.pheno)
        .dump()
        .set{ regenie_input_step_1 }

    regenie_step_1(regenie_input_step_1)

    Channel
        .fromPath(params.genotypes_bgen)
        .map { file -> tuple(file.baseName, file) }
        .groupTuple(by:0)
        .combine(make_phenotypes.out.pheno)
        .combine(regenie_step_1.out.regenie_predictions)
        .dump()
        .set{ regenie_input_step_2 }

    regenie_step_2(regenie_input_step_2)
}


workflow {
    plink2()
    regenie()
}