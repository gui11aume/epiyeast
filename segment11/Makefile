all: S11_rep1_1.txt S11_rep1_2.txt S11_rep2_1.txt S11_rep2_2.txt

SRR5680248_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/008/SRR5680248/SRR5680248_1.fastq.gz
	echo "c816bb5f552e402143dca086eefd4fbf  SRR5680248_1.fastq.gz" | md5sum -c

SRR5680248_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/008/SRR5680248/SRR5680248_2.fastq.gz
	echo "8a6fd14bc3ae46b8e5ae899df88db030  SRR5680248_2.fastq.gz" | md5sum -c

SRR5680249_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/009/SRR5680249/SRR5680249_1.fastq.gz
	echo "d8397122c74003757aa776edc381e92f  SRR5680249_1.fastq.gz" | md5sum -c

SRR5680249_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/009/SRR5680249/SRR5680249_2.fastq.gz
	echo "1fab621ce54629eea68002a0dbf91b46  SRR5680249_2.fastq.gz" | md5sum -c
	

rep1: SRR5680248_1.fastq.gz SRR5680248_2.fastq.gz
	../preprocess.py S11specs.txt \
		SRR5680248_1.fastq.gz \
		SRR5680248_2.fastq.gz rep1

rep2: SRR5680249_1.fastq.gz SRR5680249_2.fastq.gz
	../preprocess.py S11specs.txt \
		SRR5680249_1.fastq.gz \
		SRR5680249_2.fastq.gz rep2

S11_rep1_1.txt S11_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S11_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S11_rep1_2.txt

S11_rep2_1.txt S11_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S11_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S11_rep2_2.txt

cleantmp:
	rm -r rep1 rep2
