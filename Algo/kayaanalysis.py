#!/usr/bin/env python3

from tkinter import W
from kayadata import *
import os

REGIONS = ["World"] #, "EU", "USA", "China", "MENA", "AWC", "AES"]
# METHODS = [("sda weights", 4), ("sda coefficients", 4), ("ida_mult_param_two", 5), ("ida_add_param_one",5) ("ida_add_non_param_one", 5), ("ida_add_param_two", 5)]

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

        region_results = {}
        region_evolution = {}
        com = 0

        for i in range(len(Years) - 20):
            com += 1
            print(str(Years[i]) + " - " + str(Years[i+1]))
            Kaya = kayaData(DataY[Years[i]], DataY[Years[i+1]], (Years[i], Years[i+1]), Names)

            Kaya.sdaGlobal()
            Kaya.save_sda_weights_rankings("results/" + r + "/sda/")
            Kaya.save_sda_coefficients_rankings("results/" + r + "/sda/")
            w_rankings = Kaya.sda_send_weights_rankings()
            c_rankings = Kaya.sda_send_coefficients_rankings()
            if w_rankings[0] not in region_results.keys():
                region_results[w_rankings[0]] = w_rankings[1].copy()
                region_evolution[w_rankings[0]] = [[]]
                for w in w_rankings[1]:
                    region_evolution[w_rankings[0]][0].append(w.copy())
            else:
                region_evolution[w_rankings[0]].append(w_rankings[1].copy())
                for i in range(len(w_rankings[1])):
                    for j in range(len(w_rankings[1][i])):
                        region_results[w_rankings[0]][i][j] += w_rankings[1][i][j]
            print(region_evolution["sda weights"])
            if c_rankings[0] not in region_results.keys():
                region_results[c_rankings[0]] = c_rankings[1].copy()
                region_evolution[c_rankings[0]] = [[]]
                for c in c_rankings[1]:
                    region_evolution[c_rankings[0]][0].append(c.copy())
            else:
                region_evolution[c_rankings[0]].append(c_rankings[1].copy())
                for i in range(len(c_rankings[1])):
                    for j in range(len(c_rankings[1][i])):
                        region_results[c_rankings[0]][i][j] += c_rankings[1][i][j]

            Kaya.idaGlobal()
            Kaya.save_ida_rankings("results/" + r +"/ida/")
            ida_rankings = Kaya.ida_send_rankings()
            for method in ida_rankings.keys():
                if method not in region_results.keys():
                    region_results[method] = ida_rankings[method].copy()
                    region_evolution[method] = [[]]
                    for rk in ida_rankings[method]:
                        region_evolution[method][0].append(rk.copy())
                else:
                    region_evolution[method].append(ida_rankings[method].copy())
                    for i in range(len(ida_rankings[method])):
                        for j in range(len(ida_rankings[method][i])):
                            region_results[method][i][j] += ida_rankings[method][i][j]


        print()
        print()
        pprint.pprint(region_results)

        # for method in region_results.keys():
        #     print("\n" + method)
        #     for c in region_results[method]:
        #         print(sum(c))

        pprint.pprint(region_evolution)
        print(com, len(region_evolution["sda weights"]))
        



if __name__ == "__main__":
    main()
