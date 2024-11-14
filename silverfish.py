import time
import sys
import pysam
import argparse
import graphlib
import copy

def read_fasta(fasta_path):
	reads={}

	readname=0
	for line in open(fasta_path):
		if line[0] == ">":
			readname+=1		
		else:
			reads[str(readname)]=line.strip()

	return(reads)

def build_kmer_hist(seq,kmer_hist,k):
	seq_len=len(seq)
	kmers=[]
	for i in range(0,seq_len-k+1):

		kmer=seq[i:i+k]
		if len(kmer) < k:
			break
		kmers.append(kmer)

		if not kmer in kmer_hist:
			kmer_hist[kmer]=0

		kmer_hist[kmer]+=1

	return(kmers,kmer_hist)	

def read_bam(bam_path):
	reads={}
	readname=0
	samfile = pysam.AlignmentFile(bam_path, "rb")
	for read in samfile:
		if read.is_duplicate or read.is_supplementary or read.is_secondary:
			continue

		reads[readname]=read.query_sequence
		readname+=1

	return(reads)

def find_chain(graph,start,ends):
	chain=[start]
	current_node=start
	if start in ends:
		return(chain)
		
	while True:
		current_node=graph.sucessors[current_node]

		if not current_node or len(current_node) > 1 or current_node == start or current_node in ends:
			return(chain)
		current_node=list(current_node)[0]
		chain.append(current_node)
		if current_node in ends:
			return(chain)

def drop_kmers(graph,min_support):
	kmers=list(graph.kmers.keys())
	for kmer in kmers:
		if len(graph.kmers[kmer]) < min_support:
			graph.delete_kmer(kmer)
	return(graph)

def trim_edges(graph,min_weight):
	edge_list=list(graph.vertice_set)
	for edge in edge_list:
		if len(graph.vertices[edge[0]][edge[1]]) < min_weight:
			graph.delete_vertice(edge[0],edge[1])
	return(graph)

def remove_tips(graph,min_tip_length):
	branch_start=graph.in_branch_points
	branch_end=graph.out_branch_points
	starting_point=graph.starting_points

	switches=branch_end.union(branch_start)
	for start in starting_point.union(branch_start,branch_end):	
		chains=[]
		for node in graph.sucessors[start]:
			chains.append([start]+find_chain(graph,node,switches))

		for chain in chains:
			if len(chain) < 20 and chain[-1] in graph.end_points:
				for node in chain:
					graph.delete_kmer(node)

	return(graph)


def chain_typer(chain,graph):

	if chain[0] in graph.starting_points:
		return("starting_point")
	elif chain[-1] in graph.end_points:
		return("end_point")

	elif chain[0] in graph.in_branch_points:
		if chain[-1] in graph.out_branch_points:
			return("in_out")
		elif chain[-1] in graph.in_branch_points:
			return("in_in")

	elif chain[0] in graph.out_branch_points:
		if chain[-1] in graph.out_branch_points:
			return("out_out")

		elif chain[-1] in graph.in_branch_points:
			return("out_in")

	return("unknown")

def forward_scaffold(scaffold,chains,graph,chain_numbers):
	results=[]
	
	for i in range(0,len(chains)):
		if i in chain_numbers:
			continue

		if chains[i][0][0] == chains[scaffold][0][-1]:
			r=forward_scaffold(i,chains,graph,set([i]) | chain_numbers )

			for j in range(0,len(r)):
				results.append( [ chains[scaffold][0]+r[j][0][1:],r[j][1] | set([scaffold]) ] )
			
	if not results:
		results=[ [chains[scaffold][0], set([scaffold]) ] ]

	return(results)

def backward_scaffold(scaffold,chains,graph,chain_numbers):
	results=[]
	
	for i in range(0,len(chains)):
		if i in chain_numbers:
			continue

		if chains[i][0][-1] == chains[scaffold][0][0]:
			r=backward_scaffold(i,chains,graph,set([i]) | chain_numbers )

			for j in range(0,len(r)):
				results.append( [ r[j][0]+chains[scaffold][0][1:],r[j][1] | set([scaffold]) ] )
			
	if not results:
		results=[ [chains[scaffold][0], set([scaffold]) ] ]

	return(results)

def main(reads,k,min_support):

	time_all=time.time()

	kmers={}
	time_kmerize=time.time()
	graph = graphlib.graph()

	kmer_hist={}
	for read in reads:
		if len(reads[read]) < k:
			continue

		read_kmers,kmer_hist=build_kmer_hist(reads[read],kmer_hist,k)
		kmers[read]=read_kmers

	for read in kmers:
		if len(reads[read]) < k+1:
			continue


		for i in range(1,len(kmers[read])):

			if kmer_hist[kmers[read][i-1]] < min_support and kmer_hist[kmers[read][i]] < min_support:
				continue

			if kmer_hist[kmers[read][i]] < min_support and kmer_hist[kmers[read][i-1]] >= min_support:
				graph.add_kmer(kmers[read][i-1],read)

			elif kmer_hist[kmers[read][i]] >= min_support and kmer_hist[kmers[read][i-1]] < min_support:
				graph.add_kmer(kmers[read][i],read)

			if kmer_hist[kmers[read][i]] >= min_support and kmer_hist[kmers[read][i]] >= min_support:
				graph.add_vertice(kmers[read][i-1],kmers[read][i],read)		

	print("kmerized:",time.time()-time_kmerize)

	if not reads:
		print("no reads found")
		print("make sure k is shorter than read lenght, and that the input contains reads")
		quit()
	

	time_graph=time.time()

	print("graph construction:",time.time()-time_graph)

	time_clean=time.time()
	graph=drop_kmers(graph,min_support)
	graph=trim_edges(graph,min_support)
	graph=remove_tips(graph,10)
	print("cleaned graph:",time.time()-time_clean)

	branch_start=graph.in_branch_points
	branch_end=graph.out_branch_points
	starting_point=graph.starting_points

	chains=[]
	time_define_chains=time.time()

	switches=branch_end.union(branch_start)
	for start in starting_point.union(branch_start,branch_end):	
		for node in graph.sucessors[start]:

			chain=[start]+find_chain(graph,node,switches)
			chain_type=chain_typer(chain,graph)
			chains.append([chain,chain_type])

	print("defined chains:", time.time()-time_define_chains)

	scaffolds=[]
	build_scaffolds=time.time()

	for i in range(0,len(chains)):
		chain=chains[i][0]
		start=chain[0]
		end=chain[-1]

		scaffold=[]	

		if chains[i][1] == "end_point":
			results=backward_scaffold(i,chains,graph,set([i]) )
			#print(results[0][1],i)
			#print(results[0][0][-1] == chains[i][0][-1])
			scaffolds+=results

		elif chains[i] == "start_point":
			results=forward_scaffold(i,chains,graph,set([i]) )
			scaffolds+=results
		else:
			forward=forward_scaffold(i,chains,graph,set([i]) )
			for forward_result in forward:
				backward_result=backward_scaffold(i,chains,graph,forward_result[1] )
				for result in backward_result:
					scaffolds.append([  result[0]+forward_result[0][len(chains[i][0])-1:],forward_result[1] | result[1]])
	#print(scaffolds)
	print("build scaffolds",time.time()-build_scaffolds)

	for i in range(0,len(scaffolds)):
	
		skip=False
		for j in range(0,len(scaffolds)):
			if j ==i or j < i:
				continue
			if not len(scaffolds[i][-1]-scaffolds[j][-1]):
				skip=True

		if skip:
			continue

		scaffolded_chains=list(map(str,scaffolds[i][-1]))

		print(f">scaffold_{i} {','.join(scaffolded_chains)}")
		out=[]

		for j in range(1,len(scaffolds[i][0])):
			out.append(scaffolds[i][0][j][-1])

		print(scaffolds[i][0][0]+"".join(out))

min_branch_length=2
min_overlap=0.2
max_overlap=0.8


parser = argparse.ArgumentParser(prog='Silverfish',description='Local de novo assembler')
parser.add_argument('-b', '--bam')
parser.add_argument('-f', '--fasta')
parser.add_argument('-k','--kmer-length',default=91,type=int)
parser.add_argument('-s','--min-support',default=4,type=int)

args = parser.parse_args()

if args.bam and args.fasta:
	print("input either bam or fasta, not both")
	quit()

if not args.bam and not args.fasta:
	print("bam or fasta is required")
	quit()

if args.bam:
	reads=read_bam(args.bam)
else:
	reads=read_fasta(args.fasta)

main(reads,args.kmer_length,args.min_support)

