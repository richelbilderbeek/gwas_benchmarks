run:
	nextflow run main.nf

debug:
	nextflow run main.nf -dump-channels -ansi-log false

clean:
	rm -rf work/
	rm -rf .nextflow.log*