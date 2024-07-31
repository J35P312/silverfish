import networkx as nx
import matplotlib.pyplot as plt
import time
import sys
import pysam
import argparse

def read_fasta(fasta_path):
	reads={}

	readname=0
	for line in open(fasta_path):
		if line[0] == ">":
			readname+=1		
		else:
			reads[str(readname)]=line.strip()

	return(reads)

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

	successor_node=list(graph.successors(current_node))
	if not successor_node or len(successor_node) > 1:
		return(chain)	
	successor_node=successor_node[0]
		
	while successor_node:
		current_node=graph.successors(current_node)
		current_node=list(current_node)

		if not current_node:
			break
		current_node=current_node[0]
		chain.append(current_node)

		successor_node=list(graph.successors(current_node))
		if not successor_node or len(successor_node) > 1:
			return(chain)

		if successor_node[0] in ends:
			chain.append(successor_node[0])
			return(chain)

		successor_node=successor_node[0]
	return(chain)

def kmerizer(seq,k):
	read_graph = nx.DiGraph()

	seq_len=len(seq)
	kmers=[]
	for i in range(0,seq_len-k+1):
		read_graph.add_node( seq[i:i+k] )
		kmers.append(seq[i:i+k])

	for i in range(0,len(kmers)-1):
		read_graph.add_edge(kmers[i], kmers[i+1])
		
	return([read_graph,kmers[0],kmers[-1]])

def graph_to_seq(graph):
	seq=[]

	nodes=graph.nodes()
	starting_point=False

	for node in nodes:
		successors=list(graph.successors(node))
		predecessors=list(graph.predecessors(node))

		if not predecessors and starting_point:
			print("error, multiple starting points")
			quit()

		if not predecessors:
			starting_point=node

	node=starting_point
	seq=[node]
	while True:
		node=list(graph.successors(node))
		if not node:
			break
		seq.append(node[0][-1])
		node=node[0]

	return("".join(seq))

def add_edge_weight(reads,graph):
	for edge in graph.edges:
			graph[edge[0]][edge[1]]["weight"]=0

	for read in reads:
		for edge in reads[read][0].edges:
			graph[edge[0]][edge[1]]["weight"]+=1

	return(graph)

def trim_edges(graph,min_weight):
	edge_list=list(graph.edges)
	for edge in edge_list:
		if graph[edge[0]][edge[1]]["weight"] < min_weight:
			graph.remove_edge(edge[0],edge[1])
	return(graph)

def main(reads,k,min_support):

	kmers={}
	for read in reads:
		if len(reads[read]) < k+1:
			continue

		kmers[read]=kmerizer(reads[read],k)

	if not reads:
		print("no reads found")
		print("make sure k is shorter than read lenght, and that the input contains reads")
		quit()
	
	graph = nx.DiGraph()
	#add all nodes
	for read in kmers:
		graph.update(kmers[read][0])
	nodes=set(graph.nodes)


	graph=add_edge_weight(kmers,graph)
	graph=trim_edges(graph,min_support)

	starting_point=set([])
	end_point=set([])
	branch_start=set([])
	branch_end=set([])

	for node in nodes:
		successors=list(graph.successors(node))
		predecessors=list(graph.predecessors(node))

		if not successors:
			end_point.add(node)
		if not predecessors:
			starting_point.add(node)

		if len(successors) > 1:
			branch_start.add(node)

		if len(predecessors) > 1:
			branch_end.add(node)

	stretches=[]
	branch_start_chains={}
	branch_end_chains={}

	for start in starting_point.union(branch_start,branch_end):	
		for node in graph.successors(start):
			node=str(node)
			chain=[start]+find_chain(graph,node,branch_end.union(branch_start))
			if len(chain) > min_branch_length:
				stretches.append(chain)

	for i in range(0,len(stretches)):
		if stretches[i][0] in branch_end:
			branch_end_chains[stretches[i][0]]=i
		if stretches[i][-1] in branch_start:
			branch_start_chains[stretches[i][-1]]=i

	stretch_graph=[]
	bubble_arcs=[]
	
	bridging_reads={}
	scaffold_graph=nx.DiGraph()
	for i in range(0,len(stretches)):
		#stretch_graph.add_node( i )
		if stretches[i][0] in branch_start and stretches[i][-1] in branch_end:
			bubble_arcs.append(i)

		bridging_reads[i]=set([])
		stretch_graph.append( nx.subgraph(graph,stretches[i]) )

	scaffolds={}
	scaffolded_arcs=set([])

	for i in range(0,len(stretches)):
		chain=stretches[i]

		scaffold=nx.DiGraph(stretch_graph[i])
		scaffolded=set([])

		if chain[0] in starting_point:
			continue

		if chain[0] in branch_start and chain[0] in branch_start_chains:		
			scaffold.update(stretch_graph[branch_start_chains[chain[0]]])
			scaffolded_arcs.add(branch_start_chains[chain[0]])
		if chain[-1] in branch_end and chain[-1] in branch_end_chains:		
			scaffold.update(stretch_graph[branch_end_chains[chain[-1]]])
			scaffolded_arcs.add(branch_end_chains[chain[-1]])

		if scaffold:	
			scaffolds[i]=scaffold
		
	final_scaffolds=[]	
	for i in range(0,len(stretches)):
		if i in scaffolds:
			final_scaffolds.append(scaffolds[i])
		elif i in scaffolded_arcs:
			continue
		else:
			final_scaffolds.append(stretch_graph[i])

	scaffold_id=1
	for scaffold in final_scaffolds:
		print(f">scaffold_{scaffold_id}")
		print(graph_to_seq(scaffold))
		scaffold_id+=1

	#subax1 = plt.subplot(131)
	nx.draw(final_scaffolds[0], with_labels=True, font_weight='bold')
	#subax2 = plt.subplot(132)
	#nx.draw(final_scaffolds[1], with_labels=True, font_weight='bold')
	#subax4 = plt.subplot(133)
	#nx.draw(graph, with_labels=True, font_weight='bold')

	#plt.show()

min_branch_length=2
min_overlap=0.2
max_overlap=0.8

parser = argparse.ArgumentParser(prog='Silverfish',description='Local de novo assembler')
parser.add_argument('-b', '--bam')
parser.add_argument('-f', '--fasta')
parser.add_argument('-k','--kmer-length',default=101)
parser.add_argument('-s','--min-support',default=3)

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

