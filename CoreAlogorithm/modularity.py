import networkx as nx

def create_network(data):
    G = nx.Graph()
    f = open(data, 'r')
    line = f.readline()
    while line[0] == '#':
        print(line)
        line = f.readline()

    while line:
        edge = line.split()
        # print(int(edge[0]), (edge[1]))
        # print(line)
        # print(edge)
        # print(edge[0], edge[1])
        G.add_edge(edge[0], edge[1])
        line = f.readline()
    f.close()
    print('number of nodes:',G.number_of_nodes())
    return G

def get_community_node(nodefile, n):
	f = open(nodefile, 'r')
	node = list()

	line = f.readline()
	while line[0] == '#':
		print(line)
		line = f.readline()

	while line:
		edge = line.split()
		if int(edge[1]) == len(node):
			node.append([])
		node[int(edge[1])].append(edge[0])
		line = f.readline()

	f.close()
	return node

def get_degree(G, node):
	k = 0
	for i in node:
		k = k + G.degree(i)
	# print(k)
	return k

def compute_modularity(data, nodefile, n):
	G = create_network(data)
	node = get_community_node(nodefile, n);
	m = G.number_of_edges()
	# print(m)
	Q = 0
	for i in range(len(node)):
		h = G.subgraph(node[i])
		ei = h.number_of_edges()
		k = get_degree(G, node[i])
		# print(ei, k, m)
		# print('ei/m: ', ei/m)
		# print('k/(2*m):', k/(2*m))
		Q = Q + float(ei)/m - pow(float(k)/(2*m),2)
	return Q

def main():
	data = '/Users/dencesun/Desktop/snap/examples/as20graph.txt'
	nodefile = 'communities.txt'
	Q = compute_modularity(data, nodefile, 2)
	print(Q)
	# h = G.subgraph(node[0])
	# print(G.nodes())
	# print(h.nodes())

if __name__ == '__main__':
	main()
