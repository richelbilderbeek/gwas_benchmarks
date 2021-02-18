# gwas_benchmarks



## Prepare phenotypes and henotypes

```bash
nextflow run main.nf -resume -entry prep
```

will result in the following file in the `results/` directory

```
results/
    phenotypes.txt
    genotypes/
        chr*.{bim,bed}
        hardcalls/
            chr*.{bim, bed}
            merged.{bim, bed}
```