#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { plink2 } from './modules/plink.nf'
include { make_phenotypes } from './modules/utils.nf'

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

workflow {
    make_phenotypes(phenotypes_file, phenotypes_to_include)
    
    Channel
    .fromPath(params.genotypes)
    .map { file -> tuple(file.baseName, file) }
    .groupTuple(by:0)
    .combine(fam)
    .combine(make_phenotypes.out.pheno)
    .dump()
    .set{ plink2_input }

    plink2(plink2_input)
}