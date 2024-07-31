# silverfish
Local de novo assembly project, De bruijn graph asssembler using the networkx module. Should work on regions < 2kbp.

# Run
Silverfish supports fasta or bam input. Contigs are written to stdout 

	python silverfish.py --bam <bam>

	python silverfish.py --fasta <fasta>

Optionally change the kmer lenght or minimum read support.

	python silverfish.py --bam <bam> -k 1337

	python silverfish.py --bam <bam> -s 3

Example:

python silverfish.py test_data/1_1598414_1598580_0.bam

note: when reading fasta, I suggest reverse complementing reads on the negative strand before running.
When providing bam file, SIlverfish will skip duplicates, secondary alignments, and supplementary alignments. 

# Dependencies

python3, networkx, matplotlib, pysam
