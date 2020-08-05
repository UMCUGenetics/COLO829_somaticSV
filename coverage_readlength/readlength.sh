#! /bin/sh

###1###
#Use python script to get readlengths per chromosome for each technology (tumor and normal) (This is designed to submit jobs to an HPC). 
for TECH in nanopore pacbio
do
	for TYPE in tumor normal
		BAM=${TECH}.${TYPE}.bam
		for CHROM in {1..22} X
		do
			python readlength_from_bam.py $BAM $CHROM > ${BAM/.bam/.$CHROM.readlength}
		done
	done
done

#10X is a bit different, need to compute the MI tags 
for TYPE in tumor normal
do
	BAM=tenx.${TYPE}.bam
        for CHROM in {1..22} X
        do
		python tenX_readlength_MItags.py $BAM $CHROM > ${BAM/.bam/.$CHROM.readlength}
	done
done

#Bionano is also different
for TECH in bionano
do
  	for TYPE in tumor normal
        do
        BNX=${TECH}.${TYPE}.bnx
        echo $BNX
        python bionano_generate_hist.py $BNX > ${TECH}.${TYPE}.merged.readlength
        done
done

#Get insert for illumina
for TECH in illumina
do
  	for TYPE in tumor normal
        do
        java -jar /hpc/local/CentOS7/cog_bioinf/picard-tools-1.141/picard.jar CollectInsertSizeMetrics I=${TECH}.${TYPE}.bam H=${TECH}.${TYPE}.insert.picard.pdf O=${TECH}.${TYPE}.insert.picard.txt
        awk '/## HISTOGRAM/,0' ${TECH}.${TYPE}.insert.picard.txt  | tail -n+3 | grep -v "^$" | tr '\t' ',' | awk -v TYPE=$TYPE '{print $0",illumina,"TYPE}'> ${TECH}.${TYPE}.merged.readlength.hist
        done
done



###2###
#Merge all the readlength files per technology (remove the header) bionano is already merged, also illumina
for TECH in nanopore pacbio tenx
do
	for TYPE in tumor normal
	do
		cat ${TECH}.${TYPE}.*.readlength | grep -v "tech" > ${TECH}.${TYPE}.merged.readlength
		
	done
done


###3###
#Make histogram file
#1kb bins for techs, Illumina has it own from picard
for TECH in nanopore pacbio tenx bionano
do
        for TYPE in tumor normal
        do
                INPUTFILE=${TECH}.${TYPE}.merged.readlength
                echo $INPUTFILE
                awk -F "," -v TECH=$TECH -v TYPE=$TYPE 'BEGIN { MIN=9999999999999; MAX=-MIN; OFS=","; BINSIZE=1000;}
                                {        A=sprintf("%d", $4/BINSIZE);
                                        BIN[A*BINSIZE]++;
                                        }
                                END { for(X in BIN) print X, BIN[X], TECH, TYPE;}' $INPUTFILE > ${INPUTFILE}.hist
        done
done

###Plot
Rscript readlength_plot.R