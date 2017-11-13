# -*- coding: utf-8 -*-
from tests import *
from time import time


if __name__ == '__main__':
    # start1 = time()
    # test_modularity.test_modularity()
    # end1 = time()
    # print('modularity runtime = ', str(end1 - start1) + 's')

    start2 = time()
    test_HigherOrderNetwork.higher_order_network()
    end2 = time()
    print('HigerOrderNetwork runtime: ', str(end2 - start2) + 's')
    # print('all module runtime: ', str(end2 - start1) + 's')