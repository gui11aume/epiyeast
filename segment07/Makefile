all: S7_rep1_1.txt S7_rep1_2.txt S7_rep2_1.txt S7_rep2_2.txt

SRR5680240_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/000/SRR5680240/SRR5680240_1.fastq.gz
	echo "1cd3905816a9a62dcc8fb546bb9cd28a  SRR5680240_1.fastq.gz" | md5sum -c

SRR5680240_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/000/SRR5680240/SRR5680240_2.fastq.gz
	echo "8f2d9b2b92e3afa0570e89dd4d128b80  SRR5680240_2.fastq.gz" | md5sum -c

SRR5680241_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/001/SRR5680241/SRR5680241_1.fastq.gz
	echo "3ddf16213944a587ae341eda32ff666c  SRR5680241_1.fastq.gz" | md5sum -c

SRR5680241_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/001/SRR5680241/SRR5680241_2.fastq.gz
	echo "ac0477891ff0e718be2562a21441560d  SRR5680241_2.fastq.gz" | md5sum -c

rep1: SRR5680240_1.fastq.gz SRR5680240_2.fastq.gz
	../preprocess.py S7specs.txt \
		SRR5680240_1.fastq.gz \
		SRR5680240_2.fastq.gz rep1

rep2: SRR5680241_1.fastq.gz SRR5680241_2.fastq.gz
	../preprocess.py S7specs.txt \
		SRR5680241_1.fastq.gz \
		SRR5680241_2.fastq.gz rep2

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
