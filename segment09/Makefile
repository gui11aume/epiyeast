all: S9_rep1_1.txt S9_rep1_2.txt S9_rep2_1.txt S9_rep2_2.txt

SRR5680244_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/004/SRR5680244/SRR5680244_1.fastq.gz
	echo "dc35b2767b1a511c35225c7a6aa75770  SRR5680244_1.fastq.gz" | md5sum -c

SRR5680244_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/004/SRR5680244/SRR5680244_2.fastq.gz
	echo "b9fcd7123839b4d8f9c1ee610032e485  SRR5680244_2.fastq.gz" | md5sum -c

SRR5680245_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/005/SRR5680245/SRR5680245_1.fastq.gz
	echo "7fe135333c3619545837a05690adcd2b  SRR5680245_1.fastq.gz" | md5sum -c

SRR5680245_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/005/SRR5680245/SRR5680245_2.fastq.gz
	echo "9cdc2fa6224424743d90f9017976871e  SRR5680245_2.fastq.gz" | md5sum -c

rep1: SRR5680244_1.fastq.gz SRR5680244_2.fastq.gz
	../preprocess.py S9specs.txt \
		SRR5680244_1.fastq.gz \
		SRR5680244_2.fastq.gz rep1

rep2: SRR5680245_1.fastq.gz SRR5680245_2.fastq.gz:
	../preprocess.py S9specs.txt \
		SRR5680245_1.fastq.gz \
		SRR5680245_2.fastq.gz rep2

S9_rep1_1.txt S9_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S9_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S9_rep1_2.txt

S9_rep2_1.txt S9_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S9_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S9_rep2_2.txt

cleantmp:
	rm -r rep1 rep2
