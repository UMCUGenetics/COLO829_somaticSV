#! /bin/sh

SAMBAMBA=/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba

#1,000,000 random positions
BED=human_hg19.bed


for TECH in nanopore pacbio
do
        for TYPE in tumor normal
                BAM=${TECH}.${TYPE}.bam
                $SAMBAMBA depth base -t 4 --min-coverage=0 -L $BED > ${BAM}.depth
        done
done

Rscript depth_plot.R