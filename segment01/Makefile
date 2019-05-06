all: S1_rep1_1.txt S1_rep1_2.txt S1_rep2_1.txt S1_rep2_2.txt

SRR5680228_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/008/SRR5680228/SRR5680228_1.fastq.gz
	echo "c145181de18a466ffe684223dd82822a  SRR5680228_1.fastq.gz" | md5sum -c

SRR5680228_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/008/SRR5680228/SRR5680228_2.fastq.gz
	echo "6d43b781ada3d4daa265ab936afdb71f  SRR5680228_2.fastq.gz" | md5sum -c

SRR5680229_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/009/SRR5680229/SRR5680229_1.fastq.gz
	echo "7541ee58936c43aadbfe67982cf825a8  SRR5680229_1.fastq.gz" | md5sum -c

SRR5680229_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/009/SRR5680229/SRR5680229_2.fastq.gz
	echo "867402d54cd36fabbbc486ebcb615584  SRR5680229_2.fastq.gz" | md5sum -c

rep1: SRR5680228_1.fastq.gz SRR5680228_2.fastq.gz
	../preprocess.py S1specs.txt \
		SRR5680228_1.fastq.gz \
		SRR5680228_2.fastq.gz rep1

rep2: SRR5680229_1.fastq.gz SRR5680229_2.fastq.gz
	../preprocess.py S1specs.txt \
		SRR5680229_1.fastq.gz \
		SRR5680229_2.fastq.gz rep2

S1_rep1_1.txt S1_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S1_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S1_rep1_2.txt

S1_rep2_1.txt S1_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S1_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S1_rep2_2.txt

cleantmp:
	rm -r rep1 rep2
