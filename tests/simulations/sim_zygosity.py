#!/usr/bin/env python


# Script: sim_zygosity.py
# Author: David Twesigomwe
# Date modified: 15-Dec-2020
# Purpose: This script simulates the zygosity of variants after concatenating 2 haplotype VCF files.
# Usage: python3 sim_zygosity.py in_database.txt > output.txt

import os
import sys

in_file = sys.argv[1]

f = open (in_file, "r")


for line in f:
    list1 = line.strip().split(";")
    new_list = []
    for i in list1:
        x = i + '~0/1'
        y = i + '~1/1'
        if x not in new_list:
            new_list.append(x)
            
        else:
            new_list = list(map(lambda st: str.replace(st, x, y), new_list))

    print(";".join(new_list))

f.close()
