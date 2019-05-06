all: S10_rep1_1.txt S10_rep1_2.txt S10_rep2_1.txt S10_rep2_2.txt

SRR5680246_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/006/SRR5680246/SRR5680246_1.fastq.gz
	echo "7c5c6bb3d2803ca154b6d833e72c2f8e  SRR5680246_1.fastq.gz" | md5sum -c

SRR5680246_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/006/SRR5680246/SRR5680246_2.fastq.gz
	echo "fee52058be176c272808fe1f60617412  SRR5680246_2.fastq.gz" | md5sum -c

SRR5680247_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/007/SRR5680247/SRR5680247_1.fastq.gz
	echo "a08ee5b0b24b8d4407891063b14f1a4a  SRR5680247_1.fastq.gz" | md5sum -c

SRR5680247_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/007/SRR5680247/SRR5680247_2.fastq.gz
	echo "2d338fdea2f04e430889197054c46e16  SRR5680247_2.fastq.gz" | md5sum -c

rep1: SRR5680246_1.fastq.gz SRR5680246_2.fastq.gz
	../preprocess.py S10specs.txt \
		SRR5680246_1.fastq.gz \
		SRR5680246_2.fastq.gz rep1

rep2: SRR5680247_1.fastq.gz SRR5680247_2.fastq.gz
	../preprocess.py S10specs.txt \
		SRR5680247_1.fastq.gz \
		SRR5680247_2.fastq.gz rep2

S10_rep1_1.txt S10_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S10_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S10_rep1_2.txt

S10_rep2_1.txt S10_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S10_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S10_rep2_2.txt

cleantmp:
	rm -r rep1 rep2
