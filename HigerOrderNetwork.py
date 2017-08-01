# -*- coding: utf-8 -*-

"""
    reference: Jure Leskovec(2016) Higher-order organization of complex networks
    reference URL: http://science.sciencemag.org/content/353/6295/163.full.pdf

    author: Xiaowei Sun
    e-mail: dencesun@gmail.com
    date: 8.1.2017

"""

import networkx as nx
import numpy as np
from numpy import linalg as LA


class HigerOrderNetwork(object):
    def __init__(self, graph, motif=None):
        self.graph = nx.to_numpy_matrix(graph, dtype=int)
        self.motif = motif

    def MotifAdjacency(self):
        A = self.graph
        print(A)
        # ignore diagonals and weights
        A = A - np.diag(np.diag(A))
        A[A > 1] = 1

        # compute get matrix W
        try:
            lmotif = str.lower(self.motif)
        except TypeError:
            print("please input a right motif module, for example motif = 'm1'....")
            return None

        if lmotif == 'm1':
            W = self.M1(A)
        elif lmotif == 'm2':
            W = self.M2(A)
        elif lmotif == 'm3':
            W = self.M3(A)
        elif lmotif == 'm4':
            W = self.M4(A)
        elif lmotif == 'm5':
            W = self.M5(A)
        elif lmotif == 'm6':
            W = self.M6(A)
        elif lmotif == 'm7':
            W = self.M7(A)
        elif lmotif == 'm8':
            W = self.M8(A)
        elif lmotif == 'm9':
            W = self.M9(A)
        elif lmotif == 'm10':
            W = self.M10(A)
        elif lmotif == 'm11':
            W = self.M11(A)
        elif lmotif == 'm12':
            W = self.M12(A)
        elif lmotif == 'm13':
            W = self.M13(A)
        elif lmotif == 'bifan':
            W = self.Bifan(A)
        elif lmotif == 'edge':
            x = self.DirectionalBreakup(A)
            W = x[2]
        else:
            print('Error.\nUnknown motif: ', self.motif)
            return None

        return W
        # print(A)
        # x = np.random.randn(3, 3)
        # x[1, 1] = 0
        # y = x.copy()
        # print(y)
        # x[x != 0] = 1
        # print(x)
        # print(x-y)

    def M1(self, A):
        X = self.DirectionalBreakup(A)
        W = (X[1]*X[1]).dot(X[1])
        return W

    def M2(self, A):
        pass

    def M3(self, A):
        pass

    def M4(self, A):
        pass

    def M5(self, A):
        pass

    def M6(self, A):
        pass

    def M7(self, A):
        pass

    def M7(self, A):
        B, U, G = self.DirectionalBreakup(A)
        print('B\n', B)
        print('U\n', U)
        print('G\n', G)
        C1 = (U.T * B).dot(U.T)
        C1 = C1 + C1.T
        C2 = (U*U.T).dot(B)
        W = C1+C2
        return W

    def M8(self, A):
        pass

    def M9(self, A):
        pass

    def M10(self, A):
        pass

    def M11(self, A):
        pass

    def M12(self, A):
        pass

    def M13(self, A):
        pass

    def Bifan(self, A):
        pass

    def DirectionalBreakup(self, A):
        A[A != 0] = 1
        B = np.logical_and(A, A.T)*1
        U = A - B
        G = np.logical_or(A, A.T)*1
        return [B, U, G]


def main():
    g = nx.DiGraph()
    g.add_edge(1, 2)
    g.add_edge(1, 3)
    g.add_edge(1, 4)
    g.add_edge(1, 5)
    g.add_edge(1, 6)
    g.add_edge(2, 1)
    g.add_edge(2, 3)
    g.add_edge(2, 4)
    g.add_edge(2, 5)
    g.add_edge(6, 2)
    g.add_edge(6, 7)
    g.add_edge(7, 2)
    g.add_edge(7, 6)
    g.add_edge(8, 2)
    g.add_edge(8, 6)
    g.add_edge(8, 9)
    g.add_edge(8, 10)
    g.add_edge(9, 6)
    g.add_edge(9, 7)
    g.add_edge(9, 8)
    g.add_edge(9, 10)

    test = HigerOrderNetwork(g, 'm7')
    W = test.MotifAdjacency()
    print('W\n', W)

if __name__ == '__main__':
    main()