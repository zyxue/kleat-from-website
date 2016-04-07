#!/bin/bash

set -x

# need samtools-1.x NOT samtools-0.1

# need gmap 2014-03-28 (http://research-pub.gene.com/gmap/src/gmap-gsnap-2014-03-28.v2.tar.gz)

inPath () {
    if hash $1 2>/dev/null; then
        return 1
    else
        echo "${1} must be in your path before running this script!"
        exit 1
    fi
}

inPath samtools
inPath transabyss
inPath bwa
inPath gmap

red='\033[0;31m'
yellow='\033[1;33m'
green='\033[1;32m'
nc='\033[0m'

usage()
{
cat << EOF
usage: $0 options

This script will run trans-ABySS on an Intel demo computer (centos7).

OPTIONS:
    -a  REQUIRED    Read 1
    -b  REQUIRED    Read 2
    -n  REQUIRED    Name for job
    -k  REQUIRED    k-values to run
    -o  REQUIRED    Output directory
    -t  OPTIONAL    Threads to use (default=1)
    -m  OPTIONAL    Maximum memory to use for samtools sort (default=10G)
    -h  OPTIONAL    Show this message and exit
EOF
}

read1=
read2=
name=
ks=
outdir=
threads=1
maxmem=10
while getopts "a:b:n:k:o:t:m:h" OPTION
do
    case $OPTION in
        a) read1="$OPTARG" ;;
        b) read2="$OPTARG" ;;
        n) name="$OPTARG" ;;
        k) ks="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        m) maxmem="$OPTARG" ;;
        h) usage && exit 1 ;;
    esac
done

if [ -n "$read1" -a -n "$read2" -a -n "$name" -a -n "$ks" -a -n "$outdir" ]
then
    readexample=$(zcat $read1 | head | tail -1)
    readlen=${#readexample}
    ks=(${ks// / })
    mink=${ks[0]}
    maxk=${ks[${#ks[@]}-1]}
    if [ "$mink" -ge "$maxk" ]
    then
        echo "k values must be in ascending order!"
        exit 1
    fi

    merge=1
    # Run the assemblies
    for k in ${ks[@]}
    do
        stdo="$outdir"/k"$k"/"$k".ta.std.o
        stde="$outdir"/k"$k"/"$k".ta.std.e
        # Check if output directories exist and if not, create them
        if [ ! -d "$outdir"/k"$k" ]
        then
            echo "$outdir"/k"$k" does not exist! Creating.
            mkdir -p "$outdir"/k"$k"
            echo "Running transABySS with k="$k". Output logs in '"$stdo"' and '"$stde"'"
            time transabyss --SS --kmer "$k" --pe "$read1" "$read2" --outdir "$outdir"/k"$k" --name "$name" --threads "$threads" > "$stdo" 2> "$stde"
            error=$(grep -i error "$stde")
            if [ -n "$error" ]
            then
                merge=0
            fi
        else
            echo "==================================="
            echo "Assembly for k$k complete. Skipping"
            echo "==================================="
        fi
    done

    if [ "$merge" -eq 1 -a ! -d "$outdir/merged" ]
    then
        echo Merging...
        mkdir "$outdir"/merged
        # Put statement here to continue if merge has already been completed for this k
        contigs=
        prefixes=
        for k in ${ks[@]}
        do
            contigs="$contigs $outdir/k$k/$name-final.fa"
            prefixes="$prefixes k$k"
        done
        contigs=${contigs:1}
        prefixes=${prefixes:1}
        contigs=(${contigs// / })
        prefixes=(${prefixes// / })
        echo "Merging ${prefixes[@]}..."
        time transabyss-merge --SS --threads "$threads" --mink "$mink" --maxk "$maxk" --prefixes "${prefixes[@]}" --length "$readlen" "${contigs[@]}" --out "$outdir"/merged/"$name"-merged.fa --force > "$outdir"/merged/merge.o 2> "$outdir"/merged/merge.e
    elif [ "$merge" -ne 1 ]
    then
        echo "There was a problem with one of the assemblies, unable to merge. Please check log files for details."
        exit 0
    else
        echo "==================================="
        echo "         Merging complete."
        echo "==================================="
    fi

    # bwa index
    echo "Indexing merged fasta..."
    if [ ! -f "$outdir"/merged/"$name"-merged.fa.bwt ]
    then
        bwa index "$outdir"/merged/"$name"-merged.fa
    fi

    echo -e "${green}Indexing complete.${nc}"
    echo "Generating unsorted reads-to-contigs alignment (r2c.bam)..."
    if [ ! -f "$outdir"/r2c.bam -a ! -f ${outdir}/r2c_sorted.bam ]
    then
        bwa mem -t "$threads" "$outdir"/merged/"$name"-merged.fa <(zcat "$read1" "$read2") | samtools view -bhS - -o "$outdir"/r2c.bam
    fi

    echo -e "${green}r2c.bam generated.${nc}"
    echo "Sorting r2c.bam to r2c_ns.bam...(this may take a while)"
    if [ ! -f ${outdir}/r2c_ns.bam -a ! -f ${outdir}/r2c_sorted.bam ]
    then
        samtools sort -m ${maxmem} -n -o ${outdir}/r2c_ns.bam ${outdir}/r2c.bam
    fi
    echo -e "${green}r2c_ns.bam generated!${nc}"

    echo "Fixing mates in r2c_ns.bam to r2c_fm.bam..."
    if [ ! -f ${outdir}/r2c_fm.bam -a ! -f ${outdir}/r2c_sorted.bam ]
    then
        samtools fixmate "$outdir"/r2c_ns.bam "$outdir"/r2c_fm.bam
    fi
    echo -e "${green}r2c_fm.bam generated!${nc}"

    echo "Sorting r2c_fm.bam to r2c_sorted.bam..."
    if [ ! -f ${outdir}/r2c_sorted.bam -a ! -f ${outdir}/r2c_sorted.bam ]
    then
        samtools sort -m ${maxmem} -o ${outdir}/r2c_sorted.bam ${outdir}/r2c_fm.bam
    fi
    echo -e "${green}r2c_sorted.bam generated!${nc}"

    echo "Indexing r2c_sorted.bam..."
    if [ ! -f ${outdir}/r2c_sorted.bam.bai ]
    then
        samtools index "$outdir"/r2c_sorted.bam
    fi
    echo -e "${green}Index for r2c_sorted.bam generated!${nc}"
    rm "${outdir}"/r2c_ns*

#    echo "Removing extra files..."
#    files=("$outdir"/merged/"$name"-merged.fa.sa "$outdir"/merged/"$name"-merged.fa.amb "$outdir"/merged/"$name"-merged.fa.ann "$outdir"/merged/"$name"-merged.fa.bwt "$outdir"/merged/"$name"-merged.fa.pac "$outdir"/r2c.bam "$outdir"/r2c_ns.bam "$outdir"/r2c_fm.bam)
#    for f in ${files[@]}
#    do
#        if [ -f $f ]
#        then
#            rm $f
#        fi
#    done

    #c2g generation
    echo "Generating contig-to-genome alignment..."
    echo -e "${yellow}Executing:${nc} time gmap -d hg19 -D /projects/btl/arch/gmapdb_sarray/hg19 "$outdir"/merged/"$name"-merged.fa -t "$threads" -f samse -n 0 -x 10 | grep -v chrM | egrep '^[@k]' > "$outdir"/c2g.sam"

    time gmap \
	 -d hg19 \
	 -D experiment/lele/KLEAT-2.0/gmapdb_sarray/hg19 \
	 "$outdir"/merged/"$name"-merged.fa \
	 -t "$threads" \
	 -f samse \
	 -n 0 \
	 -x 10 \
	| grep -v chrM \
	| egrep '^[@k]' > "$outdir"/c2g.sam

    echo "Compressing sam to bam..."
    time samtools view -bhS "$outdir"/c2g.sam -o "$outdir"/c2g.bam
    rm "$outdir"/c2g.sam

    # bwa mem experiment/hg19/bwa-index/hg19.fa ${outdir}/merged/${name}-merged.fa \
    # 	| samtools view -h -F 2052 -S - \
    # 	| samtools sort -o ${outdir}/c2g.bam

    echo "Success!"
else
    usage
    exit 1
fi
