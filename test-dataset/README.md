This documents how the test dataset is generated

# Prerequisites (versions used at time of writing)

* Trans-ABySS (http://www.bcgsc.ca/platform/bioinfo/software/trans-abyss)
* Samtools (http://www.htslib.org/)
* BWA (http://bio-bwa.sourceforge.net/)
				
The prerequisites has been installed in the Dockerfile, as well.

# Results

The results in the `../assembly` folder were generated using

* Trans-ABySS 1.5.3
* bwa 0.7.13-r1126
* samtools-1.3

# Assembly and Alignment generation

Use the included script (`TA.sh`), and run it inside a container to generate
the input necessary for KLEAT.

```
./build_input.sh \
	-a <in1.fq.gz> \
	-b <in2.fq.gz> \
	-n <sample name> \
	-k <space-delimited k-values> \
	-o /output/path/for/assembly \
	-t <number threads> \
	-m <memory for sorting bam>
```

e.g. the command used to generate the KLEAT input in the `../assembly`
directory:

```
./build_input.sh \
    -a sample_input_reads/sample_1.fq.gz \
    -b sample_input_reads/sample_2.fq.gz \
    -n DHX30 \
    -k "25 35" \
    -o assembly \
    -t 8 \
    -m '10G'
```

This should generate the appropriate files for KLEAT input.

