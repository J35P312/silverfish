import networkx as nx
import matplotlib.pyplot as plt
import time

k=5
min_branch_length=2
min_coverage=2

min_overlap=0.2
max_overlap=0.8


read="ATTGCCGATATTA"
read3="ATTGCCGATATTA"
read2="ATATTACCCCAAATTG"
read5="ATTGCCGATCCCCAAATTG"

t=time.time()

reads={"1":read,"2":read2,"3":read3,"4":read2[0:10],"5":read5,"6":read5,"7":read5,"8":read,"9":read2}


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

		#print(current_node)
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

kmers={}
for read in reads:
	print(reads[read])
	kmers[read]=kmerizer(reads[read],k)
	
graph = nx.DiGraph()
#add all nodes
for read in kmers:
	graph.update(kmers[read][0])
nodes=set(graph.nodes)

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

#print(starting_point)
#print(end_point)
#print(branch_start)
#print(branch_end)


stretches=[]
print(starting_point.union(branch_start))

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



#stretch_graph= nx.DiGraph()
stretch_graph=[]
bubble_arcs=[]

	
bridging_reads={}
scaffold_graph=nx.DiGraph()
for i in range(0,len(stretches)):
	#stretch_graph.add_node( i )
	if stretches[i][0] in branch_start and stretches[i][-1] in branch_end:
		bubble_arcs.append(i)
		#print(i)
		#print(stretches[i])

	bridging_reads[i]=set([])
	stretch_graph.append( nx.subgraph(graph,stretches[i]) )



print(branch_start_chains)
print(branch_end_chains)

scaffolds={}
scaffolded_arcs=set([])

for i in range(0,len(stretches)):
	chain=stretches[i]

	scaffold=nx.DiGraph(stretch_graph[i])
	scaffolded=set([])

	if chain[0] in starting_point:
		continue

	if chain[0] in branch_start:		
		scaffold.update(stretch_graph[branch_start_chains[chain[0]]])
		scaffolded_arcs.add(branch_start_chains[chain[0]])
	if chain[-1] in branch_end:		
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



print(time.time()-t)

#subax1 = plt.subplot(131)
#nx.draw(final_scaffolds[0], with_labels=True, font_weight='bold')
#subax2 = plt.subplot(132)
#nx.draw(final_scaffolds[1], with_labels=True, font_weight='bold')
#subax4 = plt.subplot(133)
#nx.draw(graph, with_labels=True, font_weight='bold')

#plt.show()

