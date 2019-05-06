all: S6_rep1_1.txt S6_rep1_2.txt S6_rep2_1.txt S6_rep2_2.txt

SRR5680238_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/008/SRR5680238/SRR5680238_1.fastq.gz
	echo "2aba3affdafa4ad43a804b1261717c58  SRR5680238_1.fastq.gz" | md5sum -c

SRR5680238_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/008/SRR5680238/SRR5680238_2.fastq.gz
	echo "1e7286dc5e287155ab4bedb3d1e2634d  SRR5680238_2.fastq.gz" | md5sum -c

SRR5680239_1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/009/SRR5680239/SRR5680239_1.fastq.gz
	echo "e590f9dabf3f4cc49f317fec7c3bb548  SRR5680239_1.fastq.gz" | md5sum -c

SRR5680239_2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/009/SRR5680239/SRR5680239_2.fastq.gz
	echo "321113e967948d43a4b17c699497f723  SRR5680239_2.fastq.gz" | md5sum -c

rep1: SRR5680238_1.fastq.gz SRR5680238_2.fastq.gz
	../preprocess.py S6specs.txt \
		SRR5680238_1.fastq.gz \
		SRR5680238_2.fastq.gz rep1

rep2: SRR5680239_1.fastq.gz SRR5680239_2.fastq.gz
	../preprocess.py S6specs.txt \
		SRR5680239_1.fastq.gz \
		SRR5680239_2.fastq.gz rep2

S6_rep1_1.txt S6_rep1_2.txt: rep1
	../do_starcode.sh rep1
	../merge.py rep1/1t0.stc rep1/1t1.stc rep1/1t14.stc > S6_rep1_1.txt
	../merge.py rep1/2t0.stc rep1/2t1.stc rep1/2t14.stc > S6_rep1_2.txt

S6_rep2_1.txt S6_rep2_2.txt: rep2
	../do_starcode.sh rep2
	../merge.py rep2/1t0.stc rep2/1t1.stc rep2/1t14.stc > S6_rep2_1.txt
	../merge.py rep2/2t0.stc rep2/2t1.stc rep2/2t14.stc > S6_rep2_2.txt

cleantmp:
	rm -r rep1 rep2
