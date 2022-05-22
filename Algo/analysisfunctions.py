#!/usr/bin/env python3

from itertools import permutations, combinations #fonction combinations à utiliser pour l'optimisation de Rørmose et Olsen

def sda(dataT1, dataT2, dataDelta, indexes):
    """
    Functions that perform a structural decomposition analysis given an indexes list that gives one order for the n coefficients among the n! possibilities for that order.
    """
    n = len(indexes)
    coefficientsFactors = [0] * n
    for i in zip(indexes, range(n)):
        c = 1
        for k in range(n):
            if k < i[1]:
                c *= dataT1[indexes[k]]
            elif k > i[1]:
                c *= dataT2[indexes[k]]
        # coefficientsFactors[i[0]] = c*dataDelta[i[0]]
        coefficientsFactors[i[0]] = c
    return coefficientsFactors

def genOrders(n):
    l = list(range(n))
    return list(permutations(l))

def rank_dict(d, n):
    excl = []
    ranking = []
    for i in range(n):
        max = None
        max_k = None
        for key in d:
            if key not in excl:
                if max is None or d[key] > max:
                    max = d[key]
                    max_k = key
        ranking.append(max_k)
        excl.append(max_k)
    return ranking



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
        print(c)

