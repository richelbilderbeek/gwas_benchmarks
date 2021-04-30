#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { bolt_lmm } from './modules/bolt.nf'
include { plink2; plink2_hardcalls } from './modules/plink.nf'
include { saige_null_fitting; saige_assoc } from './modules/saige.nf'
include { regenie_step_1; regenie_step_2 } from './modules/regenie.nf'
include { make_phenotypes; make_train_test; filter_cohort; filter_hardcalls; unpack_hard_calls; merge_chromosomes; make_bgen } from './modules/utils.nf'

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
    unpack_hard_calls(params.hardcalls)
    unpack_hard_calls.out.hardcalls_list
        .flatten()
        .map{file -> tuple(file.baseName, file)}
        .set{hardcalls_list_with_names}

    make_phenotypes(phenotypes_file, phenotypes_to_include)
    make_train_test(make_phenotypes.out.pheno)
    Channel
        .fromPath(params.genotypes)
        .map{file -> tuple(file.baseName, file)}
        .groupTuple(by:0)
        .combine(fam)
        .combine(ids_to_include)
        .combine(make_phenotypes.out.pheno)
        .dump()
        .set{filter_input}

    Channel
        .fromPath(params.genotypes_bim_bed)
        .map{file -> tuple(file.baseName, file)}
        .groupTuple(by:0)
        .combine(fam)
        .combine(ids_to_include)
        .combine(make_phenotypes.out.pheno)
        .combine(hardcalls_list_with_names, by:0)
        .dump()
        .set{filter_hardcalls_input}

    filter_cohort(filter_input) 
    filter_hardcalls(filter_hardcalls_input)

    filter_cohort.out.genotypes_filtered
        .set{genotypes_bim_bed_fam}
    make_bgen(genotypes_bim_bed_fam)
        
    // merge_chromosomes(fam, filter_hardcalls.out.genotypes_hardcalls_filtered.collect())

}


workflow plink {
    make_phenotypes(phenotypes_file, phenotypes_to_include)
    Channel
        .fromPath(params.genotypes_filtered)
        .map { file -> tuple(file.baseName, file) }
        .groupTuple(by:0)
        .combine(ids_to_include)
        .combine(make_phenotypes.out.pheno)
        .dump()
        .set{ plink2_input }
    Channel
        .fromPath(params.hardcalls_filtered)
        .map { file -> tuple(file.baseName, file) }
        .groupTuple(by:0)
        .combine(ids_to_include)
        .combine(make_phenotypes.out.pheno)
        .dump()
        .set{ plink2_input_hardcalls }

    plink2(plink2_input)
    plink2_hardcalls(plink2_input_hardcalls)
}


workflow saige {
    make_phenotypes(phenotypes_file, phenotypes_to_include)
    Channel
        .fromPath(params.hardcalls_merged)
        .collect()
        .combine(make_phenotypes.out.pheno)
        .dump()
        .set{ saige_input_hardcalls }

    saige_null_fitting(saige_input_hardcalls)

    Channel
        .fromPath(params.genotypes_bgen)
        .map { file -> tuple(file.baseName, file) }
        .groupTuple(by:0)
        .combine(saige_null_fitting.out.null_model)
        .combine(saige_null_fitting.out.variance_ratio)
        .dump()
        .set{ saige_input_assoc }

    saige_assoc(saige_input_assoc)
}


workflow bolt {
    make_phenotypes(phenotypes_file, phenotypes_to_include)
    Channel
        .fromPath(params.hardcalls_filtered)
        .collect()
        .set{hardcalls}
    Channel
        .fromPath(params.genotypes_bgen)
        .collect()
        .set{imputed}

    bolt_lmm(hardcalls, imputed, make_phenotypes.out.pheno)
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