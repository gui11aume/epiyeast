all: S5_rep1_1.txt S5_rep1_2.txt S5_rep2_1.txt S5_rep2_2.txt

rep1:
	../preprocess.py S5specs.txt \
		/data/GF5_1_12302_AAGACG_read1.fastq.gz \
		/data/GF5_1_12302_AAGACG_read2.fastq.gz rep1
rep2:
	../preprocess.py S5specs.txt \
		/data/GF5_2_12303_GGATGT_read1.fastq.gz \
		/data/GF5_2_12303_GGATGT_read2.fastq.gz rep2

S5_rep1_1.txt S5_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S5_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S5_rep1_2.txt

S5_rep2_1.txt S5_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S5_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S5_rep2_2.txt

cleantmp:
	rm -r rep1 rep2
