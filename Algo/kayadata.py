#!/usr/bin/env python3

from analysisfunctions import *
from math import factorial
import numpy as np
import matplotlib.pyplot as plt
import csv

COLORS = ["r", "b", "g", "m", "y", "c"]


def import_data(namefile):
    with open(namefile) as file:
        years = []
        coefficient_names = []
        data_values = []
        first = True
        for row in csv.reader(file):
            if first:
                for cell in row:
                    try:
                        years.append(int(cell))
                    except ValueError:
                        continue
                first = False
            else:
                data_values.append([])
                for cell in row:
                    try:
                        data_values[-1].append(float(cell))
                    except ValueError:
                        coefficient_names.append(cell)

    normalized_data = []
    for x in range(len(data_values)):
        normalized_data.append([])
        for i in range(len(data_values[x])):
            normalized_data[-1].append(data_values[x][i] / data_values[x][0])

    normalized_data_per_year = {}
    for y in years:
        normalized_data_per_year[y] = []
        i = years.index(y)
        for nd in normalized_data:
            normalized_data_per_year[y].append(nd[i])
    
    return(years, coefficient_names, data_values, normalized_data, normalized_data_per_year)




class kayaData:

    def __init__(self, dataT1, dataT2, years = (0,0), coefficient_names = []):
        self.dataT1 = dataT1
        self.dataT2 = dataT2
        self.years = years
        self.dataDelta = []
        self.n = len(dataT1) #The number of coefficients
        self.coefficients_names = []
        self.reduced_c_names = []
        if coefficient_names == []:
            for i in range(self.n):
                self.coefficient_names.append("C" + str(i))
                self.reduced_c_names.append("C" + str(i))
        else:
            for i in range(len(coefficient_names)):
                if not i % 2:
                    self.coefficients_names.append(coefficient_names[i])
                else:
                    self.reduced_c_names.append(coefficient_names[i])
        self.sdaMin = [0] * self.n
        self.sdaMax = [0] * self.n
        self.sda_range = [0] * self.n
        self.sdaStd = []
        self.coefficients = []
        self.mean_coefficients = []
        self.rank_coefficients = []
        for i in range(self.n):
            self.coefficients.append([])
            self.rank_coefficients.append([0] * self.n)
        for T in zip(dataT1, dataT2):
            self.dataDelta.append(T[1] - T[0])

    def deltaYi(self, i):
        Y = self.dataT2[i]
        for j in range(self.n):
            if j != i:
                Y *= self.dataT1[j]
        return(Y)

    # SDA

    def sda(self, indexes):
        return sda(self.dataT1, self.dataT2, self.dataDelta, indexes)

    def sdaMean(self):
        if self.coefficients[0] == []:
            return -1
        meanSda = []
        for cl in self.coefficients:
            Mean = 0
            for c in cl:
                Mean += c
            meanSda.append(Mean/len(cl))
        self.mean_coefficients = meanSda

    def sdaGlobal(self):
        allIndexes = genOrders(self.n)
        n = 0
        for indexes in allIndexes:
            n+=1
            coefficients = self.sda(indexes)
            ranking = {}
            for i in range(self.n):
                self.coefficients[i].append(coefficients[i])
                ranking[i] = coefficients[i]
            ranks = rank_dict(ranking, self.n)
            for r in zip(ranks, range(self.n)):
                self.rank_coefficients[r[0]][r[1]] += 1
        for i in range(self.n):
            self.sdaMax[i] = max(self.coefficients[i])
            self.sdaMin[i] = min(self.coefficients[i])
        for cl in self.coefficients:
            self.sdaStd.append(np.std(cl))
        self.sdaMean()
    
    def print_sda(self):
        if self.coefficients[0] == []:
            print("SDA non effectuée")
        else:
            weights = []
            results = []

            for i in range(self.n):
                weights.append(100 * self.mean_coefficients[i] / sum(self.mean_coefficients))
                self.sda_range[i] = self.sdaMax[i] - self.sdaMin[i]
                results.append(2.932960 / self.mean_coefficients[i])
            print("Mean values : ", self.mean_coefficients)
            print("Weights : ", weights, "  /  Sum of weights : ", sum(weights))
            print("Min values : ", self.sdaMin)
            print("Max values : ", self.sdaMax)
            print("ranges : ", self.sda_range)
            print("results : ", results)
            print("Classements de C : ", self.rank_coefficients[0])
            print("Classements de E : ", self.rank_coefficients[1])
            print("Classements de A : ", self.rank_coefficients[2])
            print("Classements de P : ", self.rank_coefficients[3])

    def show_sda(self):
        fig_sda = plt.subplot()
        legend_text = []
        x_labels = []
        for x in range(self.n):
            fig_sda.plot(x + 1, self.sdaMin[x], "vk", label = '_nolegend_')
            fig_sda.plot(x + 1, self.sdaMax[x], "^k", label = '_nolegend_')
            bars = fig_sda.bar(x + 1, self.mean_coefficients[x], yerr = self.sdaStd[x], color = COLORS[x], width = .95, linewidth = .7, edgecolor = "black", capsize = 10)
            fig_sda.bar_label(bars, label_type = "center")
            legend_text.append(self.coefficients_names[x] + " | range : " + str(round(self.sdaMin[x], 3)) + " -> " + str(round(self.sdaMax[x], 3)) + " | std : " + str(round(self.sdaStd[x], 3)) + " (" + str(round(100 * self.sdaStd[x] / self.mean_coefficients[x], 1)) + "%)")
            x_labels.append(self.coefficients_names[x])
        fig_sda.set_xticks(np.arange(1, self.n + 1), x_labels)
        fig_sda.set_title("Valeur moyenne des coefficients sur les {}! décompositions classiques pour la variation de {} à {}".format(self.n, self.years[0], self.years[1]), fontsize = 15)
        fig_sda.legend(legend_text)
        fig_sda.set_xlabel("Coefficients à décomposer", fontsize = 12)
        fig_sda.set_ylabel("Valeurs des coefficients", fontsize = 12)
        
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()
        plt.show()


    # IDA

    def ida_mult_param_two(self, alpha = None):
        return mult_parametric_method_two(self.dataT1, self.dataT2, alpha)

    def ida_add_param_one(self, alpha = 0.5):
        return add_parametric_method_one(self.dataT1, self.dataT2, alpha)
    
    def ida_add_non_param_one(self):
        return add_non_parametric_method_one(self.dataT1, self.dataT2)
    
    def ida_add_param_two(self, alpha = None):
        return add_parametric_method_two(self.dataT1, self.dataT2, alpha)


    # Results
    
    def show_rankings(self):
        fig_ranks = plt.subplot()
        d = 1
        max_rank = 0
        x_labels = []
        for i in range(self.n):
            x = []
            for r in range(self.n):
                if max_rank < self.rank_coefficients[i][r]:
                    max_rank = self.rank_coefficients[i][r]
                x.append( d + r)
                x_labels.append(self.reduced_c_names[i] + "\nR" + str(r+1))
            d += self.n 
            bars = fig_ranks.bar(x, self.rank_coefficients[i], width = .95, color = COLORS[i], linewidth = .7, edgecolor = "black")
            fig_ranks.bar_label(bars)
        fig_ranks.set(xlim = (0, self. n * self.n + 1), ylim = (0, max_rank + 1), yticks = np.arange(1, max_rank + 2))
        fig_ranks.set_xticks(np.arange(1, self.n * self.n + 1), x_labels)
        fig_ranks.set_title("Nombre de fois où un coefficient se trouve au rang \"R\" en terme d'importance pour la variation de {} à {}".format(self.years[0], self.years[1]))
        fig_ranks.set_xlabel("Nom et rang \"R\" des coefficients")
        fig_ranks.set_ylabel("Nombre de décomposition donnant le coefficient au rang \"R\"")
        
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()
        plt.show()
            





def main():
    # T1 = [2.39963380951199, 1/3.49203330338253, 5558.11979126848, 5280062644]
    # T2 = [2.39146347617893, 1/5.04044509019937, 8016.70518839688, 6114324044]
    # T3 = []
    # for k in range(4):
    #     T3.append(T2[k] / T1[k])

    
    
    (Years, Names, Data, NData, NDataY) = import_data("kayaData.csv")

    Kaya9000 = kayaData(NDataY[1990], NDataY[2000], (1990, 2000), Names)
    Kaya0005 = kayaData(NDataY[2000], NDataY[2005], (2000, 2005), Names)
    Kaya0510 = kayaData(NDataY[2005], NDataY[2010], (2005, 2010), Names)
    Kaya1014 = kayaData(NDataY[2010], NDataY[2014], (2010, 2014), Names)
    KayaGlobal = kayaData(NDataY[1990], NDataY[2014], (1990, 2014), Names)
    
    Kaya9000.sdaGlobal()
    # Kaya9000.show_sda()
    # Kaya9000.show_rankings()
    P = np.prod(Kaya9000.dataT2)
    print(P)
    print("mult_param_two")
    print(Kaya9000.ida_mult_param_two(), np.prod(Kaya9000.ida_mult_param_two()), "res = ", P / np.prod(Kaya9000.ida_mult_param_two()))
    print(Kaya9000.ida_mult_param_two(0), np.prod(Kaya9000.ida_mult_param_two(0)), "res = ", P /  np.prod(Kaya9000.ida_mult_param_two(0)))
    print(Kaya9000.ida_mult_param_two(0.5), np.prod(Kaya9000.ida_mult_param_two(.5)), "res = ", P /  np.prod(Kaya9000.ida_mult_param_two(.5)))
    print(Kaya9000.ida_mult_param_two(1), np.prod(Kaya9000.ida_mult_param_two(1)), "res = ", P /  np.prod(Kaya9000.ida_mult_param_two(1)))
    
    S = P - 1
    print("\n", S, sep ="")
    print("add_param_one")
    print(Kaya9000.ida_add_param_one(0), sum(Kaya9000.ida_add_param_one(0)), "res = ", S - sum(Kaya9000.ida_add_param_one(0)))
    print(Kaya9000.ida_add_param_one(), sum(Kaya9000.ida_add_param_one()), "res = ", S - sum(Kaya9000.ida_add_param_one()))
    print(Kaya9000.ida_add_param_one(1), sum(Kaya9000.ida_add_param_one(1)), "res = ", S - sum(Kaya9000.ida_add_param_one(1)))
    print("add_non_param_one")
    print(Kaya9000.ida_add_non_param_one(), sum(Kaya9000.ida_add_non_param_one()), "res = ", S - sum(Kaya9000.ida_add_non_param_one()))
    print("add_param_two")
    print(Kaya9000.ida_add_param_two(), sum(Kaya9000.ida_add_param_two()), "res = ", S - sum(Kaya9000.ida_add_param_two()))
    print(Kaya9000.ida_add_param_two(0), sum(Kaya9000.ida_add_param_two(0)), "res = ", S - sum(Kaya9000.ida_add_param_two(0)))
    print(Kaya9000.ida_add_param_two(0.5), sum(Kaya9000.ida_add_param_two(.5)), "res = ", S - sum(Kaya9000.ida_add_param_two(.5)))
    print(Kaya9000.ida_add_param_two(1), sum(Kaya9000.ida_add_param_two(1)), "res = ", S - sum(Kaya9000.ida_add_param_two(1)))
    print("\nsda")
    print(Kaya9000.mean_coefficients)

    # Kaya0005.sdaGlobal()
    # Kaya0005.show_sda()
    # Kaya0005.show_rankings()

    # Kaya0510.sdaGlobal()
    # Kaya0510.show_sda()
    # Kaya0510.show_rankings()

    # Kaya1014.sdaGlobal()
    # Kaya1014.show_sda()
    # Kaya1014.show_rankings()

    # KayaGlobal.sdaGlobal()
    # KayaGlobal.show_sda()
    # KayaGlobal.show_rankings()


if __name__ == "__main__":
    main()