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
import math
import matplotlib.pyplot as plt
from scipy.io import savemat
from scipy.sparse.csgraph import connected_components

np.seterr(divide='ignore', invalid='ignore')

class HigerOrderNetwork(object):
    def __init__(self, graph):
        self.graph = nx.to_numpy_matrix(graph, dtype=int)

    def MotifAdjacency(self, motif = None):
        """
        MotifAdjacency forms the motif adjacency matrix for the adjacency matrix A and the specified motif

        :param motif: motif is one of
        m1, m2, m3, m4, m5, m6, m7, m8. m9, m10, m11, m12, m13, bifan, etc
        :return: W = MotifAdjacency(motif)
                 W: the motif adjacency matrix

        """
        A = self.graph
        # print('A\n', A)
        # print(type(A))
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

        # print(A)
        # x = np.random.randn(3, 3)
        # x[1, 1] = 0
        # y = x.copy()
        # print(y)
        # x[x != 0] = 1
        # print(x)
        # print(x-y)

        return W

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
        tmp = A.T
        NA = np.logical_and((A==0)*1, (tmp == 0)*1)*1

        # print('bifan U\n', U)
        # print('bifan G\n', G)
        # print('NA\n', NA)
        W = np.zeros(G.shape)
        W = np.mat(W)
        # print('bifan W\n', W)
        # print('bifan triu(NA, 1)\n', np.triu(NA, 1))
        nzero_ind = np.nonzero(np.triu(NA, 1))
        # print('bifan nzero_ind\n', nzero_ind[1].shape)

        for ind in range(nzero_ind[0].shape[0]):
            x = nzero_ind[0][ind]
            y = nzero_ind[1][ind]
            # print(x, y)
            xout = np.nonzero(U[x, :])
            yout = np.nonzero(U[y, :])
            # print('bifan np.nonzero(U[x, :])\n', np.nonzero(U[x, :]))
            # print('bifan np.nonzero(U[y, :])\n', np.nonzero(U[y, :]))
            # print('bifan xout\n', xout)
            # print('bifan yout\n', yout)
            common = np.intersect1d(xout[1], yout[1])
            # print('conmon\n',common)
            nc = len(common)
            # print('bifan nc', nc)
            # print(common[0, 0])
            for i in range(nc):
                for j in range(i+1, nc):
                    w = common[i]
                    v = common[j]
                    if NA[w, v] == 1:
                        W[x, y] = W[x, y] + 1
                        W[x, w] = W[x, w] + 1
                        W[x, v] = W[x, v] + 1
                        W[y, w] = W[y, w] + 1
                        W[y, v] = W[y, v] + 1
                        W[w, v] = W[w, v] + 1

        # print(nzero_ind[0])
        # print(nzero_ind[1])
        # print(type(nzero_ind[0]))
        # print('U\n', U)
        # print('G\n', G)
        # print('NA\n', NA)
        # print('W\n', W)
        # print('nzero_ind\n', nzero_ind)
        # print(nzero_ind[0].shape[0])

        W = W + W.T

        # print(type(W))
        # print(W)
        
        return W

    def SpectralPartitioning(self, A):
        """
        SpectralPartitioning performs a spectral partitioning  of A and returns several relevant quantities. It assumes
        that A is undirected and connected

        :param A: motif adjacency matrix
        :return:
                cluster, condv, condc, order = SpectralPartitioning(A)
                cluster: vector of nodes in the smaller side of partition
                condv: the sweep conductance vector
                condc: the conductance of cluster
                order: the spectral ordering

        """
        part_vec, x= self.nfiedler(A)
        # print(part_vec)
        order = np.argsort(part_vec) # return an index matrix order

        # order = np.array([9,7,8, 5,6, 1, 4, 0, 3, 2])
        # order = np.mat(order)
        # print('\norder.shape\n', order.shape)
        # print(A)

        n = order.shape[1]
        crtesianProduct = [(order[0, i], order[0, j]) for i in range(n) for j in range(n)]
        # compute the conductance values (vectorized)
        B = np.zeros([n, n])
        for i in range(len(crtesianProduct)):
            # print(int(i/10), i%10)
            # print(crtesianProduct[i][0], crtesianProduct[i][1])
            # B[int(i/10), i%10] = A[crtesianProduct[i][0], crtesianProduct[i][1]]
            B[int(i/n), i%n] = A[crtesianProduct[i][0], crtesianProduct[i][1]]
        # print('B\n', B)

        B_lower = np.tril(B)
        # print('B_lower\n', B_lower)
        B_sums = np.mat(B.sum(axis = 1)).T
        # print('B_sums\n', B_sums)
        B_lower_sums = np.mat(B_lower.sum(axis = 1)).T
        # print('B_lower_sums\n', B_lower_sums)
        volumes = np.cumsum(B_sums, axis = 0)

        # print('\nB_sums\n', B_sums)
        # print('\nvolumes\n', volumes)
        # print('B_sums-2*B_lower_sums\n', B_sums-2*B_lower_sums)

        num_cut = np.cumsum(B_sums-2*B_lower_sums, axis = 0)
        # print('num_cut\n', num_cut)
        total_vol = np.sum(A)
        # print('total_vol\n', total_vol)
        volumes_other = total_vol * np.ones((order.shape[1], 1)) - volumes
        # print('\nvolumes_other\n', volumes_other)
        vols = np.fmin(volumes, volumes_other)
        # print('vols\n', vols)
        scores = num_cut / vols
        # print('\nscores\n', scores)
        scores = scores[0:-1, :]
        # print('scores\n', scores)
        min_ind = np.argmin(scores)
        condc = scores[min_ind, 0]
        # print(min_ind, condc)
        n = A.shape[0]
        # the output cluster is the smaller of the two sides of the partition
        if min_ind <= math.floor(n/2):
            cluster = order[0, 0:min_ind+1]
        else:
            cluster = order[0, min_ind+1:n]

        # print('scores\n', scores)
        # print('scores[::-1, 0]\n', scores[::-1, 0])

        condv = np.fmin(scores, scores[::-1, :])
        # print('condv\n', condv)
        condv = condv[0:math.ceil(scores.shape[0]/2), :]

        # print(condv)
        # 画图验证正确性
        # x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        # plt.plot(x, scores[:, 0])
        # # plt.show()
        # plt.savefig('demo.svg')

        return cluster, condv, condc, order

    def nfiedler(self, A = None, tol = 1e-12):
        """
        :param A: motif adjacency matrix A
        :param tol: no uesd in this program
        :return: the fiedler vector of the normalized laplacian of A
                 (fiedler vector: eigenvector corresponding to the second smallest eigenvalue)

        """
        if A is None:
            print('Error! matrix A is None..')
            return None

        L = self.nlaplacian(A)
        # print('L\n', L, '\n')
        n = A.shape[0]
        # print(n)
        eigvalue, eigvector = SLA.eigh(L + np.eye(n))

        # print('eigvalue\n', eigvalue)
        # print('eigvector\n', eigvector[:, 0], '\n\n', eigvector[:, 1])
        # print(L + np.eye(n))
        # print('eigvalues\n', eigvalue)
        # print('eigvalue\n', eigvalue[:2])
        # print('eigvector\n', eigvector[0:1])
        # print(eigvalue)

        x = eigvector[:, 1]

        # print('\n\nx\n', x, '\n')
        # print(x)

        x = x/(np.sqrt(np.sum(A, 1)).T)

        # x = x.T
        # print('\nx\n', x.T)

        eigvalue  = eigvalue[1] - 1

        # print('\neigvalue\n', eigvalue)
        # print('\nnp.argsort(e+igvector[1])\n', np.argsort(eigvector[1]), '\n')
        # print('eigvector[1]\n', eigvector[1])
        # print('\nx\n', x)
        # print('\neigvalue\n', eigvalue)

        return x, eigvalue

    def nlaplacian(self, A):
        """
        :param A: matrix A
        :return: the normalized laplacian of A

        """
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

    def DirectionalBreakup(self, A):
        """
        :param A:
        :return: B, U, G = DirectionalBreakup(A)
                B: the bidirectional subgraph
                U: the undirectional subgraph
                G: the undirected graph

        """

        A = (A!=0)*1 # replace nonzero elements with one
        B = np.logical_and(A, A.T)*1
        U = A - B
        G = np.logical_or(A, A.T)*1

        # print('A\n', A)
        # print('B\n', B)
        # print('U\n', U)
        # print('G\n', G)

        return B, U, G

    def LargestConnectComponent(self, A):
        """
        LargestConnectComponent gets the largest connnected component of A
        A is assumed to be undirected

        :param A: motif adjacency matrix
        :return: LCC, lcc_inds, ci, sizes = LargestConnectComponent(A)
                LCC: the largest connected component
                lcc_inds: the indices in A corresponding to the connected component
                ci: the component indices of the nodes in A
                sizes: the sizes of the connected components

        """
        x, y = connected_components(A)
        mask = list()

        for i in range(y.shape[0]):
            if y[i] != 0:
                A[i, :] = 0
                A[:, i] = 0
                mask.append(i)

        A = np.delete(A, mask, 0)
        A = np.delete(A, mask, 1)

        return A


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


def main():
    """
    graph data g come from Jure Leskovec(2016) Higher-order organization of complex networks Fig 1

    :return: None

    """

    # g = nx.DiGraph()
    # g.add_edge(1, 2)
    # g.add_edge(1, 3)
    # g.add_edge(1, 4)
    # g.add_edge(1, 5)
    # g.add_edge(1, 6)
    # g.add_edge(2, 1)
    # g.add_edge(2, 3)
    # g.add_edge(2, 4)
    # g.add_edge(2, 5)
    # g.add_edge(6, 2)
    # g.add_edge(6, 7)
    # g.add_edge(7, 2)
    # g.add_edge(7, 6)
    # g.add_edge(8, 2)
    # g.add_edge(8, 6)
    # g.add_edge(8, 9)
    # g.add_edge(8, 10)
    # g.add_edge(9, 6)
    # g.add_edge(9, 7)
    # g.add_edge(9, 8)
    # g.add_edge(9, 10)
    # #
    # test = HigerOrderNetwork(g)
    # W = test.MotifAdjacency('bifan')


    data = '/home/sun/PycharmProjects/Network/C-elegans-frontal.txt'
    DG = create_network(data)
    test = HigerOrderNetwork(DG)
    W = test.MotifAdjacency('bifan')
    W = test.LargestConnectComponent(W)
    print('W.shape', W.shape)
    # print('y.shape', y.shape)
    # print('test x\n', x)
    # print('test y\n', y)
    # save as matlab file format
    # savemat('W.mat', {'W': W})
    cluster, condv, condc, order = test.SpectralPartitioning(W)
    print('condc: ', condc)
    print('cluster:',  cluster)

if __name__ == '__main__':
    main()
