import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import argparse

G = 6.67e-8
mu = 2.31
mp = 1.67e-24
gamma = 5/3
kb = 1.3806e-16
lambda_0 = 2.671e-33
m = 1.3
n = 2.3

def func(T, *args):
    nH, m, n, lambda_0 = args
    return np.sqrt(3*np.pi / 32.0 / G/ mu / mp / nH) - nH * kb * T / (gamma - 1) / lambda_0 / nH ** m / T ** n

def main():
    nH_start = -2
    nH_end = 10
    nH_count = 1000
    nH_array = np.logspace(nH_start, nH_end, nH_count)
    root = []
    lambda_0 = 2.671e-33
    m = 1.3
    n = 2.3
    for nH in nH_array:
        root.append(fsolve(func, 1.0, args=(nH, m, n, lambda_0))[0])
    print(root)
    # plt.plot(np.log10(nH_array), np.log10(root))
    # plt.show()
    return


if __name__ == "__main__":

    # parser = argparse.ArgumentParser()
    # parser.add_argument("-m", "--m", help="power index of density")
    # parser.add_argument("-n", "--n", help="power index of temperature")
    # parser.add_argument("-l", "--lambda_0", help="initial cooling rate")
    # args = parser.parse_args()

    # if arg.m != ""
    #     m = float(m)
    # if arg.n != ""
    #
    # if args.row != "":
    #     row = int(args.row)
    # else:
    #     raise Exception("row must > 0")
    #
    # if args.column != "":
    #     column = int(args.column)
    # else:
    #     raise Exception("column must > 0")
    #
    # if args.input != "":
    #     index_list = args.input.split(",")
    #
    #     if len(index_list) != (row * column):
    #         raise Exception("len of index list is %d. Not compatiable in %d x %d grid" % (len(index_list), row, column))
    # else:
    #     raise Exception("must be specific the index list")
    #
    # if args.folder != "":
    #     folder_path = args.folder
    # else:
    #     raise Exception("must be specific the root folder")
    main()
