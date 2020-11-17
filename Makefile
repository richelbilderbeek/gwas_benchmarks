run:
	nextflow run main.nf -resume -profile bianca

plink:
	nextflow run main.nf -resume -profile bianca -entry plink2

regenie:
	nextflow run main.nf -resume -profile bianca -entry regenie

debug:
	nextflow run main.nf -dump-channels -ansi-log false

clean:
	rm -rf work/
	rm -rf .nextflow.log*