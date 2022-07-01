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
        print()
        print(r)

        (Years, Names, DataY, NDataY) = import_data("Data/" + r + "Data.csv")
        for i in range(len(Years) - 1):
            print(str(Years[i]) + " - " + str(Years[i+1]))
            Kaya = kayaData(DataY[Years[i]], DataY[Years[i+1]], (Years[i], Years[i+1]), Names)

            Kaya.sdaGlobal()
            Kaya.save_sda_weights_rankings("results/" + r + "/sda/")
            Kaya.save_sda_coefficients_rankings("results/" + r + "/sda/")

            Kaya.idaGlobal()
            Kaya.save_ida_rankings("results/" + r +"/ida/")
        



if __name__ == "__main__":
    main()
