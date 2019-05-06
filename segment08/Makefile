all: S8_rep1_1.txt S8_rep1_2.txt S8_rep2_1.txt S8_rep2_2.txt

SRR8957086_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR895/006/SRR8957086/SRR8957086_1.fastq.gz
	echo "bf69e116511ee22a69617bffb7c802ff  SRR8957086_1.fastq.gz" | md5sum -c

SRR8957086_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR895/006/SRR8957086/SRR8957086_2.fastq.gz
	echo "5538172cc7f4a63bfdd46a2cf89b2e56  SRR8957086_2.fastq.gz" | md5sum -c

SRR8957087_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR895/007/SRR8957087/SRR8957087_1.fastq.gz
	echo "1752f1581180349d78f4143b44b3a64b  SRR8957087_1.fastq.gz" | md5sum -c

SRR8957087_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR895/007/SRR8957087/SRR8957087_2.fastq.gz
	echo "bf396368241adf5a9a9dfbe78409c89e  SRR8957087_2.fastq.gz" | md5sum -c

rep1: SRR8957086_1.fastq.gz SRR8957086_2.fastq.gz
	../preprocess.py S8specs.txt \
		SRR8957086_1.fastq.gz \
		SRR8957086_2.fastq.gz rep1

rep2: SRR8957087_1.fastq.gz SRR8957087_2.fastq.gz
	../preprocess.py S8specs.txt \
		RR8957087_1.fastq.gz
		SRR8957087_2.fastq.gz rep2

S8_rep1_1.txt S8_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S8_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S8_rep1_2.txt

S8_rep2_1.txt S8_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S8_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S8_rep2_2.txt

cleantmp:
	rm -r rep1 rep2
