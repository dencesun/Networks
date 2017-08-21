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
from scipy import linalg as SLA
from numpy import linalg as NLA

class HigerOrderNetwork(object):
    def __init__(self, graph):
        self.graph = nx.to_numpy_matrix(graph, dtype=int)

    def MotifAdjacency(self, motif = None):
        A = self.graph

        # ignore diagonals and weights
        A = A - np.diag(np.diag(A))
        A[A > 1] = 1

        # compute get matrix W
        try:
            lmotif = str.lower(motif)
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

    def DirectionalBreakup(self, A):
        A[A != 0] = 1
        B = np.logical_and(A, A.T)*1
        U = A - B
        G = np.logical_or(A, A.T)*1
        return [B, U, G]

    def M1(self, A):
        B, U, G = self.DirectionalBreakup(A)
        W = np.multiply(U*U, U)
        return W

    def M2(self, A):
        B, U, G = self.DirectionalBreakup(A)
        C = np.multiply(B*U, U.T) + np.multiply(U*B, U.T) + np.multiply(U*U, B)
        W = C + C.T
        return W

    def M3(self, A):
        B, U, G = self.DirectionalBreakup(A)
        C = np.multiply(B*B, U) + np.multiply(B*U, B) + np.multiply(U*B, B)
        W = C + C.T
        return W

    def M4(self, A):
        B, U, G = self.DirectionalBreakup(A)
        W = np.multiply(B*B, B)
        return W

    def M5(self, A):
        B, U, G = self.DirectionalBreakup(A)
        T1 = np.multiply(U*U, U)
        T2 = np.multiply((U.T)*U, U)
        T3 = np.multiply(U*(U.T), U)
        C = T1 + T2 + T3
        W = C + C.T
        return W

    def M6(self, A):
        B, U, G = self.DirectionalBreakup(A)
        C1 = np.multiply(U*B, U)
        C1 = C1 + C1.T
        C2 = np.multiply((U.T)*U, B)
        W = C1 + C2
        return W

    def M7(self, A):
        B, U, G = self.DirectionalBreakup(A)
        C1 = np.multiply((U.T)*B, U.T)
        C1 = C1 + C1.T
        C2 = np.multiply(U*(U.T), B)
        W = C1 + C2
        return W

    def M8(self, A):
        B, U, G = self.DirectionalBreakup(A)
        W = np.zeros(G.shape)
        W = np.asmatrix(W)
        N = G.shape[0]
        # print('U\n', U)
        for i in range(N):
            J = np.nonzero(U[i, :])[1]
            # print(J)
            for j1 in range(len(J)):
                for j2 in range(j1+1, len(J)):
                    k1 = J[j1]
                    k2 = J[j2]
                    # print(k1, k2)
                    if A[k1, k2] == 0 and A[k2, k1] == 0:
                        W[i, k1] = W[i, k1] + 1
                        W[i, k2] = W[i, k2] + 1
                        W[k1, k2] = W[k1, k2] + 1

        W = W + W.T
        # matlab use W = sparse(W + W')
        # I think it is properly use W = W+W'T
        return W

    def M9(self, A):
        B, U, G = self.DirectionalBreakup(A)
        W = np.zeros(G.shape)
        N = G.shape[0]
        # print(np.nonzero(U))
        for i in range(N):
            J1 = np.nonzero(U[i, :])[1]
            # np.nonzero(x)中x不能是列向量,只能是行向量或者矩阵
            J2 = np.nonzero(U[:, i].T)[1]
            for j1 in range(len(J1)):
                for j2 in range(len(J2)):
                    k1 = J1[j1]
                    k2 = J2[j2]
                    if A[k1, k2] == 0 and A[k2, k1] == 0:
                        W[i, k1] = W[i, k1] + 1
                        W[i, k2] = W[i, k2] + 1
                        W[k1, k2] = W[k1, k2] + 1
        W = W + W.T
        return W

    def M10(self, A):
        W = self.M8(A.T)
        return W

    def M11(self, A):
        B, U, G = self.DirectionalBreakup(A)
        W = np.zeros(G.shape)
        N = G.shape[0]
        for i in range(N):
            J1 = np.nonzero(B[i, :])[1]
            J2 = np.nonzero(U[i, :])[1]
            for j1 in range(len(J1)):
                for j2 in range(len(J2)):
                    k1 = J1[j1]
                    k2 = J2[j2]
                    if A[k1, k2] == 0 and A[k2, k1] == 0:
                        W[i, k1] = W[i, k1] + 1
                        W[i, k2] = W[i, k2] + 1
                        W[k1, k2] = W[k1, k2] + 1
        W = W + W.T
        return W

    def M12(self, A):
        W = self.M11(A.T)
        return W

    def M13(self, A):
        B, U, G = self.DirectionalBreakup(A)
        print('B\n', B)
        W = np.zeros(G.shape)
        N = G.shape[0]
        for i in range(N):
            J = np.nonzero(B[i, :])[1]
            print(J)
            for j1 in range(len(J)):
                for j2 in range(j1+1, len(J)):
                    print(j1, j2)
                    k1 = J[j1]
                    k2 = J[j2]
                    if A[k1, k2] == 0 and A[k2, k1] == 0:
                        W[i, k1] = W[i, k1] + 1
                        W[i, k2] = W[i, k2] + 1
                        W[k1, k2] = W[k1, k2] + 1
        W = W + W.T
        return W

    # 未完成
    def Bifan(self, A):
        B, U, G = self.DirectionalBreakup(A)
        NA = np.logical_and((A==0)*1, (A.T == 0)*1)*1
        # print(NA)
        W = np.zeros(G.shape)

        return None

    def SpectralPartitioning(self, A):
        part_vec, x= self.nfiedler(A)
        # print(part_vec)
        # print(np.argsort(part_vec))

    def nfiedler(self, A = None, tol = 1e-12):
        # if A == None:
        #     print('Error! matrix A is None..')
        #     return None
        L = self.nlaplacian(A)
        n = A.shape[0]
        # print(n)
        eigvalue, eigvector = SLA.eigh(L + np.eye(n))
        # print(L + np.eye(n))
        print('eigvalue\n', eigvalue)
        # print('eigvalue\n', eigvalue[:2])
        # print('eigvector\n', eigvector[0:1])
        # print(eigvalue)
        x = eigvector[1]
        # print(x)
        x = x/(np.sqrt(np.sum(A, 1)).T)
        # print(x)
        eigvalue  = eigvalue[1] - 1
        print(eigvalue)
        print(np.argsort(eigvector[1]))
        print(eigvector[1])

        return x, eigvalue

    def nlaplacian(self, A):
        d = A.sum(axis = 1) # axis = 0 求每一列的和，axis = 1 求每一行的和
        # print(A)
        # print('np.sqrt(d[np.nonzero(d)])\n', np.sqrt(d[np.nonzero(d)]))
        # print('1/np.sqrt(d[np.nonzero(d)])\n', 1/(np.sqrt(d[np.nonzero(d)])).T)
        d = 1/np.sqrt(d[np.nonzero(d)]).T
        # print('d \n', d)
        # print('type(d)\n', type(d))
        # print(np.nonzero(A)[0].T)
        # print(A[np.nonzero(A)])
        x = np.nonzero(A)
        i = x[0]
        j = x[1]
        v = A[x]
        # print(i.T, '\n', j.T, '\n', v.sum(), '\n')
        m, n = A.shape
        # print('m, n\n', m, n)
        # test = sparse.csc_matrix((v, (i, j)), shape=(m,n))
        v = np.multiply(-v.T, np.multiply(d[i], d[j]))
        # print(v)
        # print(v.shape)
        # print(v.sum())
        L = np.zeros((m, n))
        # print(L.shape)
        ziped = zip(i.T, j.T)
        ziped = zip(ziped, list(v))
        for z in ziped:
            L[z[0][0], z[0][1]] = z[1][0]

        L = L + np.mat(np.eye(n))
        return L


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

    test = HigerOrderNetwork(g)
    W = test.MotifAdjacency('m7')
    # print('W\n', W)
    # test.nfiedler(W)
    test.SpectralPartitioning(W)

if __name__ == '__main__':
    main()