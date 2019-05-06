all: S5_rep1_1.txt S5_rep1_2.txt S5_rep2_1.txt S5_rep2_2.txt

SRR5680236_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/006/SRR5680236/SRR5680236_1.fastq.gz
	echo "6bc33c10b60c775f9ca5eb800145d401  SRR5680236_1.fastq.gz" | md5sum -c

SRR5680236_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/006/SRR5680236/SRR5680236_2.fastq.gz
	echo "51e9b3e420d77ee878188ff05acb1c62  SRR5680236_2.fastq.gz" | md5sum -c

SRR5680237_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/007/SRR5680237/SRR5680237_1.fastq.gz
	echo "96e9ed4280b61baef6647fb307d94c92  SRR5680237_1.fastq.gz" | md5sum -c

SRR5680237_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/007/SRR5680237/SRR5680237_2.fastq.gz
	echo "aa5c2704fe9b8e6c809c1a6308d72852  SRR5680237_2.fastq.gz" | md5sum -c

rep1: SRR5680236_1.fastq.gz SRR5680236_2.fastq.gz
	../preprocess.py S5specs.txt \
		SRR5680236_1.fastq.gz \
		SRR5680236_2.fastq.gz rep1

rep2: SRR5680237_1.fastq.gz SRR5680237_2.fastq.gz
	../preprocess.py S5specs.txt \
		SRR5680237_1.fastq.gz \
		SRR5680237_2.fastq.gz rep2

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
