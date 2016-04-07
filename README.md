# Prerequisites

* BLAT (http://www.kentinformatics.com/)
* pysam (https://github.com/pysam-developers/pysam)

# Run KLEAT

To run KLEAT, simply have python installed along with the python module
pysam. You may then run KLEAT as below:

```
python KLEAT.py \
	<assembly_dir>/c2g.bam <assembly_dir>/sample-merged.fa \
	<reference_genome.fa> \
	<gtf> \
	<assembly_dir>/r2c_sorted.bam \
	/output/path/prefix \
	-k <track_name> \
	<track_description> \
	-ss
```

e.g.

```
python KLEAT.py \
	test/assembly/c2g.bam \
	test/assembly/merged/DHX30-merged.fa \
	hg19.fa \
	ensembl.fixed.sorted.gz \
	test/assembly/r2c_sorted.bam \
	test/KLEAT/DHX30 \
	-k KLEAT_test "KLEAT cleavage sites" \
	-ss
```

You may need to alter some of the options, for example the `-ss` option is
indicating that your data is strand-specific.

# KLEAT output

If your output is specified as `/path/to/output/prefix` then all files will be
output in the `/path/to/output/` directory with various filenames beginning
with `'prefix'`. Below are the files output by KLEAT using the above options.

```
prefix.KLEAT
    This is the main output file that KLEAT produces.
    Inside you can find rows of data corresponding to
    various cleavage site locations as well as other
    useful information.

prefix.-.bg
    This is a bedgraph track of the cleavage sites
    associated with the negatively stranded contigs.

prefix.+.bg
    Like above, but for positive strands.

prefix.HEXAMERS.bed
    A bed track consisting of all the hexamer binding
    sites identified by KLEAT which provide supporting
    evidence for cleavage sites. These hexamers are 6
    nucleotides in length and function as a binding site
    for CPSF, which is involved in cleavage of the 3'
    UTR. There are 16 supported hexamers each with a 
    colour corresponding to their strength level. Below 
    are the hexamers in order of strongest to weakest, 
    with their respective colours 
    (visible in UCSC genome browser).
    AATAAA  =   255,0,0
    ATTAAA  =   255,100,100
    AGTAAA  =   255,150,150
    TATAAA  =   255,200,200
    CATAAA  =   0,255,0
    GATAAA  =   100,255,100
    AATATA  =   150,255,150
    AATACA  =   200,255,200
    AATAGA  =   0,0,255
    AAAAAG  =   100,100,255
    ACTAAA  =   150,150,255
    AAGAAA  =   200,200,255
    AATGAA  =   255,0,255
    TTTAAA  =   255,100,255
    AAAACA  =   255,150,255
    GGGGCT  =   255,200,255

prefix.stats
    A file containing summary statistics
```
