all: S4_rep1_1.txt S4_rep1_2.txt S4_rep2_1.txt S4_rep2_2.txt

SRR5680234_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/004/SRR5680234/SRR5680234_1.fastq.gz
	echo "f6b7381ff8450f964fec7dbbb7ab6c08  SRR5680234_1.fastq.gz" | md5sum -c

SRR5680234_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/004/SRR5680234/SRR5680234_2.fastq.gz
	echo "d1186eff8e3fbdb915fb9ea2a1a52124  SRR5680234_2.fastq.gz" | md5sum	-c

SRR5680235_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/005/SRR5680235/SRR5680235_1.fastq.gz
	echo "7715c9c52a7337e92704479c295470e3  SRR5680235_1.fastq.gz" | md5sum -c

SRR5680235_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/005/SRR5680235/SRR5680235_2.fastq.gz
	echo "50892c639708f097219cdc4144c7c1d4  SRR5680235_2.fastq.gz" | md5sum -c

rep1: SRR5680234_1.fastq.gz SRR5680234_2.fastq.gz
	../preprocess.py S4specs.txt \
		SRR5680234_1.fastq.gz \
		SRR5680234_2.fastq.gz rep1

rep2: SRR5680235_1.fastq.gz SRR5680235_2.fastq.gz
	../preprocess.py S4specs.txt \
		SRR5680235_1.fastq.gz \
		SRR5680235_2.fastq.gz rep2

S4_rep1_1.txt S4_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S4_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S4_rep1_2.txt

S4_rep2_1.txt S4_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S4_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S4_rep2_2.txt

cleantmp:
	rm -r rep1 rep2

