#!/usr/bin/env python3

from kayadata import *
import os

REGIONS = ["World", "EU", "USA", "China", "MENA", "AWC", "AES"]


def main():

    os.system("rm -rf results")
    os.system("mkdir results")

    for r in REGIONS:

        os.system("mkdir results/" + r)
        os.system("mkdir results/" + r + "/sda")
        os.system("mkdir results/" + r + "/ida")

        (Years, Names, NDataY) = import_data("Data/" + r + "Data.csv")
        for i in range(len(Years) - 1):
            Kaya = kayaData(NDataY[Years[i]], NDataY[Years[i+1]], (Years[i], Years[i+1]), Names)
            Kaya.sdaGlobal()
            Kaya.save_sda_weights_rankings("results/" + r + "/sda/")
        



if __name__ == "__main__":
    main()
