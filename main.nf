#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { plink2; plink2_hardcalls } from './modules/plink.nf'
include { make_phenotypes; filter_cohort; filter_hardcalls; unpack_hard_calls; merge_chromosomes } from './modules/utils.nf'

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
    merge_chromosomes(filter_hardcalls.out.genotypes_hardcalls_filtered)
}


workflow plink {
    make_phenotypes(phenotypes_file, phenotypes_to_include)
    Channel
        .fromPath(params.genotypes)
        .map { file -> tuple(file.baseName, file) }
        .groupTuple(by:0)
        .combine(fam)
        .combine(ids_to_include)
        .combine(make_phenotypes.out.pheno)
        .dump()
        .set{ plink2_input }
    Channel
        .fromPath(params.genotypes_bim_bed)
        .map { file -> tuple(file.baseName, file) }
        .groupTuple(by:0)
        .combine(fam)
        .combine(ids_to_include)
        .combine(make_phenotypes.out.pheno)
        .dump()
        .set{ plink2_input_hardcalls }

    plink2(plink2_input)
    plink2_hardcalls(plink2_input_hardcalls)
}


workflow regenie {
    Channel
        .fromPath(params.genotypes_bim_bed)
        .map { file -> tuple(file.baseName, file) }
        .groupTuple(by:0)
        .combine(fam)
        .dump()
        .set{ regenie_input_unfiltered }
    remove_not_biallelic(fam, regenie_input_unfiltered)
    merge_chromosomes(remove_not_biallelic.out.genotypes_biallelic.collect())
    plink_qc
    regenie
}


workflow {
    plink2()
    regenie()
}