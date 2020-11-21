#!/usr/bin/env python3


import os
import sys


def get_backgroud_alleles(database, core_vars):

    dbs = []
    dbs_temp = []

    core_vars_list = core_vars.split(";")
    core_temp1 = core_vars_list[-1][:-4]
    core_temp2 = core_vars_list[0][:-4]

    for line in open(database, "r"):
        line = line.strip().split("\t")
        dbs.append(line)

    for record in dbs:
        temp_rec = record[1]
        # record_core_var = record[1].split(";")
        
        if core_temp1 and core_temp2 in temp_rec:
            dbs_temp.append(record)

            
    scores = []

    for elem in dbs_temp:
        record_core_var = elem[1].split(";")
        
        counter = 0

        for i in record_core_var:
            if i in core_vars_list:
                counter += 3
            elif i[:-4] in core_vars:
                counter += 1
            else:
                counter += -1

        scores.append(counter)

    max_score = max(scores)
    max_index = scores.index(max_score)

    diplo1 = dbs_temp[max_index][0]

    res1 = [i for i in range(len(diplo1)) if diplo1.startswith("_", i)]
    res2 = [i for i in range(len(diplo1)) if diplo1.startswith(".", i)]
    hap1 = "*" + str (diplo1[:res2[0]])
    hap2 = "*" + str (diplo1[res1[0]+1:res2[1]])
    allele_res =  hap1 + "/" + hap2 
    return allele_res
