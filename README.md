# silverfish
Local de novo assembly project, De bruijn graph asssembler using the networkx module

# Run

python silverfish.py <input_fasta> > <output_fasta>

Example:

python silverfish.py test_data/1_1598414_1598580_0.fasta

note the assmebler does not handle reads from different stands, I suggest reverse complementing reads on the negative strand before running.

# Dependencies

python3, networkx, matplotlib
