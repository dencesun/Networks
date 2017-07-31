# -*- coding: utf-8 -*-

"""
    reference: M. E. J. Newman(2006) Modularity and community structure in networks
    reference URL: http://www.pnas.org/content/103/23/8577.full.pdf

    author: Xiaowei Sun
    e-mail: dencesun@gmail.com
    date: 7.31.2017
"""

import networkx as nx
import numpy as np
from numpy import linalg as LA


class Modularity(object):

    def __init__(self, graph):
        self.graph = graph

    def modularity(self):
        # convert to numpy adjacency matrix
        A = nx.to_numpy_matrix(self.graph)

        # compute adjacency matrix A's degree centrality
        degree_centrality = np.sum(A, axis=0, dtype=int)
        m = np.sum(degree_centrality, dtype=int)/2

        # compute matrix B


        B = np.zeros(A.shape, dtype=float)

        for i in range(len(A)):
            for j in range(len(A)):
                B[i, j] = A[i, j] - (degree_centrality[0, i]*degree_centrality[0, j])/float(2*m)

        # compute A's eigenvector
        w, v = LA.eig(B)
        wmax = np.argmax(w)
        s = np.zeros((len(A), 1), dtype=float)

        for i in range(len(A)):
            if v[i, wmax] < 0:
                s[i, 0] = -1
            else:
                s[i, 0] = 1

        Q = s.T.dot(B.dot(s))/float(4*m)
        return Q[0, 0], s


def main():
    g = nx.karate_club_graph()

    Q = Modularity(g)
    q, s = Q.modularity()

    print(q)

if __name__ == '__main__':
    main()