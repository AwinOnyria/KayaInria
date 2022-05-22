#!/usr/bin/env python3

from analysisfunctions import sda_test, genOrders

def main():
    allIndexes = genOrders(4)
    for indexes in allIndexes:
        print()
        sda_test(indexes)

if __name__ == "__main__":
    main()