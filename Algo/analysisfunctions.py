#!/usr/bin/env python3

from itertools import permutations, combinations
from xmlrpc.client import MAXINT #fonction combinations à utiliser pour l'optimisation de Rørmose et Olsen
import numpy as np


## General functions


def rank_dict(d, n):
    excl = []
    ranking = []
    for i in range(n):
        Max = None
        max_k = None
        for key in d:
            if key not in excl:
                if Max is None or abs(d[key]) > Max:
                    Max = abs(d[key])
                    max_k = key
        ranking.append(max_k)
        excl.append(max_k)
    return ranking

def min_dict_key(d):
    Min = float('inf')
    for k in d.keys():
        if d[k] < Min:
            min_k = k
            Min = d[k]
    return min_k




## Functions for SDA 


def sda(dataT1, dataT2, dataDelta, indexes):
    """
    Functions that perform a structural decomposition analysis given an indexes list that gives one order for the n coefficients among the n! possibilities for that order.
    """
    n = len(indexes)
    coefficientWeights = [0] * n
    coefficientFactors = [0] * n
    for i in zip(indexes, range(n)):
        c = 1
        for k in range(n):
            if k < i[1]:
                c *= dataT1[indexes[k]]
            elif k > i[1]:
                c *= dataT2[indexes[k]]
        coefficientFactors[i[0]] = c * dataDelta[i[0]]
        coefficientWeights[i[0]] = c
    return coefficientWeights, coefficientFactors

def genOrders(n):
    l = list(range(n))
    return list(permutations(l))

def sda_test(indexes):
    n = len(indexes)
    chars = ["a", "b", "c", "d"]
    print(chars[indexes[0]] + " " + chars[indexes[1]] + " " + chars[indexes[2]] + " " + chars[indexes[3]])
    for i in zip(indexes, range(n)):
        c = chars[i[0]] + " := 1"
        for k in range(n):
            if k < i[1]:
                c += " * " + chars[indexes[k]] + "0"
            elif k > i[1]:
                c += " * " + chars[indexes[k]] + "1"
            else:
                c += " * Δ" + chars[indexes[k]]
        print(c)


## Functions for IDA

# Multiplicative (Y2 / Y1) = prod

def mult_parametric_method_two(dataT1, dataT2, result, alpha = None):
    n = len(dataT1)
    D_factors = []
    Y1 = np.prod(dataT1)
    Y2 = np.prod(dataT2)
    for i in range(n):
        if alpha is None:
            alpha = (dataT2[i] - dataT1[i]) / dataT1[i] - np.log(dataT2[i] / dataT1[i]) * Y1
            alpha /= 0 - (1 / dataT2[i] - 1 / dataT1[i]) * (dataT2[i] - dataT1[i])
        D_factors.append(np.exp((dataT2[i] - dataT1[i]) * (1 / dataT1[i] + alpha * (1 / dataT2[i] - 1 / dataT1[i]))))
    return D_factors.append(result / np.prod(D_factors))

# Additive (Y2 - Y1) = sum

def add_parametric_method_one(dataT1, dataT2, result, alpha = 0.5):
    if alpha is None:
        return None
    n = len(dataT1)
    D_factors = []
    Y1 = np.prod(dataT1)
    Y2 = np.prod(dataT2)
    for i in range(n):
        D_factors.append(np.log(dataT2[i] / dataT1[i]) * (Y1 + alpha * (Y2 - Y1)))
    return D_factors.append(result - sum(D_factors))

def add_non_parametric_method_one(dataT1, dataT2, result):
    n = len(dataT1)
    D_factors = []
    Y1 = np.prod(dataT1)
    Y2 = np.prod(dataT2)
    for i in range(n):
        D_factors.append(np.log(dataT2[i] / dataT1[i]) * (Y2 - Y1) / np.log(Y2 / Y1))
    return D_factors.append(result - sum(D_factors))

def add_parametric_method_two(dataT1, dataT2, result, alpha = None):
    n = len(dataT1)
    D_factors = []
    Y1 = np.prod(dataT1)
    Y2 = np.prod(dataT2)
    for i in range(n):
        if alpha is None:
            alpha = ((dataT2[i] - dataT1[i]) * (Y1 / dataT1[i]) - np.log(dataT2[i] / dataT1[i]) * Y1)
            alpha /= ((Y2 - Y1) * np.log(dataT2[i] / dataT1[i]) - (Y2 / dataT2[i] - Y1 * dataT1[i]) * (dataT2[i] - dataT1[i]))
        D_factors.append((dataT2[i] - dataT1[i]) * ((Y1 / dataT1[i]) + alpha * (Y2 / dataT2[i] - Y1 / dataT1[i])))
    return D_factors.append(result - sum(D_factors))
        