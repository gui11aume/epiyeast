all: S2_rep1_1.txt S2_rep1_2.txt S2_rep2_1.txt S2_rep2_2.txt

SRR5680230_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/000/SRR5680230/SRR5680230_1.fastq.gz
	echo "c4160de5393a82354722c0e327a0558c  SRR5680230_1.fastq.gz" | md5sum -c

SRR5680230_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/000/SRR5680230/SRR5680230_2.fastq.gz
	echo "6b1e84d9b1c1ee5718abcf3f2ad0f9c0  SRR5680230_2.fastq.gz" | md5sum -c

SRR5680231_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/001/SRR5680231/SRR5680231_1.fastq.gz
	echo "2e82589a71ebf2d18d4a086c4431a926  SRR5680231_1.fastq.gz" | md5sum -c

SRR5680231_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/001/SRR5680231/SRR5680231_2.fastq.gz
	echo "2d53ddf639719601e5e87ba0e1ce1fd5  SRR5680231_2.fastq.gz" | md5sum -c

rep1: SRR5680230_1.fastq.gz SRR5680230_2.fastq.gz
	../preprocess.py S2specs.txt \
		SRR5680230_1.fastq.gz \
		SRR5680230_2.fastq.gz rep1

rep2: SRR5680231_1.fastq.gz SRR5680231_2.fastq.gz
	../preprocess.py S2specs.txt \
		SRR5680231_1.fastq.gz \
		SRR5680231_2.fastq.gz rep2

S2_rep1_1.txt S2_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S2_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S2_rep1_2.txt

S2_rep2_1.txt S2_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S2_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S2_rep2_2.txt

cleantmp:
	rm -r rep1 rep2
