all: S7_rep1_1.txt S7_rep1_2.txt S7_rep2_1.txt S7_rep2_2.txt

rep1:
	../preprocess.py S7specs.txt \
		/data/GF7_1_15019_CCTTCA_read1.fastq.gz \
		/data/GF7_1_15019_CCTTCA_read2.fastq.gz rep1
#		/data/S07_rep1_read1.fastq.gz \
#		/data/S07_rep1_read2.fastq.gz rep1

rep2:
	../preprocess.py S7specs.txt \
		/data/S07_rep2_read1.fastq.gz \
		/data/S07_rep2_read2.fastq.gz rep2
#		/data/GF7_2_15020_TTCAGC_read1.fastq.gz \
#		/data/GF7_2_15020_TTCAGC_read2.fastq.gz rep2

S7_rep1_1.txt S7_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S7_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S7_rep1_2.txt

S7_rep2_1.txt S7_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S7_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S7_rep2_2.txt

cleantmp:
	rm -r rep1 rep2
