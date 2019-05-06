all: S3_rep1_1.txt S3_rep1_2.txt S3_rep2_1.txt S3_rep2_2.txt

SRR5680232_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/002/SRR5680232/SRR5680232_1.fastq.gz
	echo "bca113392c943975561865a81873eb2a  SRR5680232_1.fastq.gz" | md5sum -c

SRR5680232_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/002/SRR5680232/SRR5680232_2.fastq.gz
	echo "46ea72875a31e0d04d98f0e831ee4afd  SRR5680232_2.fastq.gz" | md5sum -c

SRR5680233_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/003/SRR5680233/SRR5680233_1.fastq.gz
	echo "d9abef6753592724eaeb5cc20930c648  SRR5680233_1.fastq.gz" | md5sum -c

SRR5680233_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/003/SRR5680233/SRR5680233_2.fastq.gz
	echo "c7e607dfd479d0d432900d955b425ad7  SRR5680233_2.fastq.gz" | md5sum -c

rep1: SRR5680232_1.fastq.gz SRR5680232_2.fastq.gz
	../preprocess.py S3specs.txt \
		SRR5680232_1.fastq.gz \
		SRR5680232_2.fastq.gz rep1

rep2: SRR5680233_1.fastq.gz SRR5680233_2.fastq.gz
	../preprocess.py S3specs.txt \
		SRR5680233_1.fastq.gz \
		SRR5680233_2.fastq.gz rep2

S3_rep1_1.txt S3_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S3_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S3_rep1_2.txt

S3_rep2_1.txt S3_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S3_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S3_rep2_2.txt

cleantmp:
	rm -r rep1 rep2

