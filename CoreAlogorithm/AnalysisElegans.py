# -*- coding: utf-8

"""
    analysis C.elegans network's property

    author: Xiaowei Sun
    date: 8.24.2017

"""

import networkx as nx
import matplotlib.pyplot as plt
# from collections import Counter


def create_network(data):
    DG = nx.DiGraph()
    with open(data, 'r') as f:
        line = f.readline()
        while line[0] == '#':
            print(line)
            line = f.readline()

        while line:
            edge = line.split()
            # print(int(edge[0]), (edge[1]))
            DG.add_edge(int(edge[0]), int(edge[1]))
            line = f.readline()
        f.close()
    return DG


def draw_network(DG):
    nx.draw_networkx(DG, pos=nx.spring_layout(DG), node_color = 'red')
    plt.axis('off') # turn off axis
    plt.savefig('c-elegans.svg')
    # plt.show()


def degree_distribution(DG):
    in_degree = DG.in_degree()
    out_degree = DG.out_degree()

    # print(in_degree)
    in_histogram = list([0]) * (in_degree[max(in_degree, key=in_degree.get)] + 1)
    out_histogram = list([0]) * (out_degree[max(out_degree, key=out_degree.get)] + 1)

    for k, v in in_degree.items():
        in_histogram[v] += 1

    # in_histogram = list([0])
    # in_degree_count = Counter(in_degree_sequence)
    # print(sum(in_histogram))

    in_pk = [z/sum(in_histogram) for z in in_histogram]
    # print(in_histogram)
    # print(in_pk)

    plt.title('in_degree')
    plt.xlabel('K+1')
    plt.ylabel('Pk')
    plt.loglog(range(len(in_histogram)), in_pk, 'bo')
    plt.savefig('c-elegans_in_degree_distribution.svg')
    # plt.show()
    for k, v in out_degree.items():
        out_histogram[v] += 1
    # print(len(out_degree_sequence))
    out_pk = [z/sum(out_histogram) for z in out_histogram]
    plt.title('out_degree_distribution')
    plt.xlabel('K+1')
    plt.ylabel('Pk')
    plt.loglog(range(len(out_histogram)), out_pk, 'bo')
    plt.savefig('c-elegans_out_degree_distribution.svg')


def clustering_coefficient(DG):
    pass


def centrality(DG):
    in_degree_centrality = nx.in_degree_centrality(DG)
    out_degree_centrality = nx.out_degree_centrality(DG)
    with open('/home/sun/PycharmProjects/Network/in_degree_centrality.csv', 'w') as f:
        for k, v in in_degree_centrality.items():
            f.write(str(k) + ': ' + str(v) + '\n')
        f.close()

    with open('/home/sun/PycharmProjects/Network/out_degree_centrality.csv', 'w') as f:
        for k, v in out_degree_centrality.items():
            f.write(str(k) + ': ' + str(v) + '\n')
        f.close()


# def main():
#     data = '/home/sun/PycharmProjects/Network/C-elegans-frontal.txt'
#     # data = 'www.adj'
#     DG = create_network(data)
#
#     # draw_network(DG)
#     # clustering_coefficient(DG)
#     # centrality(DG)
#     degree_distribution(DG)
#
# if __name__ == '__main__':
#     main()
#
#     # DG = nx.DiGraph()
#     # DG.add_edge(1,2)
#     # print(DG.edges())
#     # # pos = nx.nx_agraph.graphviz_layout(DG)
#     # nx.draw_networkx(DG, pos = nx.spring_layout(DG))
#     # plt.show()
#     # plt.ishold()
#     # plt.draw(DG)