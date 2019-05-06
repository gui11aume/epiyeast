all: S12_rep1_1.txt S12_rep1_2.txt S12_rep2_1.txt S12_rep2_2.txt

SRR5680250_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/000/SRR5680250/SRR5680250_1.fastq.gz
	echo "36ceb2350d337dcd443d34f813e2fa7c  SRR5680250_1.fastq.gz" | md5sum -c

SRR5680250_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/000/SRR5680250/SRR5680250_2.fastq.gz
	echo "cb8aa4f415ba0c31b575d277d1cc0099  SRR5680250_2.fastq.gz" | md5sum -c

SRR5680251_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/001/SRR5680251/SRR5680251_1.fastq.gz
	echo "464a04be63605e3afeb2dca0c3dd63df  SRR5680251_1.fastq.gz" | md5sum -c

SRR5680251_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/001/SRR5680251/SRR5680251_2.fastq.gz
	echo "2640b755593fe7f12deea54450b9330e  SRR5680251_2.fastq.gz" | md5sum -c

rep1: SRR5680250_1.fastq.gz SRR5680250_2.fastq.gz
	../preprocess.py S12specs.txt \
		SRR5680250_1.fastq.gz \
		SRR5680250_2.fastq.gz rep1

rep2: SRR5680251_1.fastq.gz SRR5680251_2.fastq.gz
	../preprocess.py S12specs.txt \
		SRR5680251_1.fastq.gz \
		SRR5680251_2.fastq.gz rep2

S12_rep1_1.txt S12_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S12_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S12_rep1_2.txt

S12_rep2_1.txt S12_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S12_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S12_rep2_2.txt

cleantmp:
	rm -r rep1 rep2

