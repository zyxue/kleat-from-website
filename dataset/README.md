This documents how the test dataset is generated

# Prerequisites

2. Assembly
	* transABySS (http://www.bcgsc.ca/platform/bioinfo/software/trans-abyss)
	* Samtools (http://www.htslib.org/)
	* GMAP (http://www.molecularevolution.org/software/genomics/gmap)
	* BWA (http://bio-bwa.sourceforge.net/)

# Results

The results in the `./assembly` and `./KLEAT` folder were
generated using

* transABySS 1.5.2
* GMAP version 2014-01-21
* bwa 0.7.4-r385
* samtools v1.1
* python 2.7.9
* pysam 0.8.1
* blat v. 34

# Assembly and Alignment generation

Use the included script (`TA.sh`), and run it inside a container to generate
the input necessary for KLEAT.

```
./TA.sh \
	-a <in1.fq.gz> \
	-b <in2.fq.gz> \
	-n <sample name> \
	-k <space-delimited k-values> \
	-o /output/path/for/assembly \
	-t <number threads> \
	-m <memory for sorting bam>
```

e.g.

```
./TA.sh -a sample_1.fq.gz -b sample_2.fq.gz -n DHX30 -k "32 52 72" -o test/assembly -t 6 -m 15G
```

This should generate the appropriate files for KLEAT input.

