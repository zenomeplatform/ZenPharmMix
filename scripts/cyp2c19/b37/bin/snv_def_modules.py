#!/usr/bin/env python3

import os
import sys




def get_core_variants(infile, cn):
    core_vars = []
    for line in open(infile, "r"):
        line = line.strip()
        core_vars.append(line)
    core_vars = ";".join(sorted(core_vars))

#    Block zygosity change!
#    if int(cn) == 1:
#        core_vars = core_vars.replace("~0/1", "~1/1")

    return core_vars

def get_all_vars_gt(infile_full_gt):
    all_vars_gt = []
    for line in open(infile_full_gt, "r"):
        line = line.strip()
        all_vars_gt.append(line)
    all_vars_gt = ";".join(sorted(all_vars_gt))
    return all_vars_gt

def cand_snv_allele_calling(database, infile, infile_full, infile_full_gt, infile_spec, cn):
    

    f = open(infile_spec, "r")

    all_variants = []

    for line in open(infile_full, "r"):
        line.strip()
        all_variants.append(line)
        # all_variants = line.strip().split(";")
        # print(all_variants)

    if os.stat(infile).st_size == 0:
        cand_res = ['1.v1_1.v1']
        allele_res = "*1/*1"
        return ["".join(cand_res), allele_res];
        #print("\nSupporting variants")
        #print("\n" + "".join(all_variants))
        sys.exit()

    # core_variants = []

    # for line in open(infile, "r"):
    #      line = line.strip()
    #      core_variants.append(line)

    # core_variants = ";".join(sorted(core_variants))

    core_variants = get_core_variants(infile, cn)

    # if int(cn) == 1:
    #     core_variants = core_variants.replace("~0/1", "~1/1")

    # else:
    #     pass

    all_var_gt = []
    for line in open(infile_full_gt, "r"):
        line = line.strip()
        all_var_gt.append(line)

    
    dbs = []

    for line in open(database, "r"):
        line = line.strip().split("\t")
        dbs.append(line)

    soln_list1 = []
    soln_list2 = []

    for record in dbs:
        record_core_var = record[1].split(";")
        record_core_var = ";".join(sorted(record_core_var))
        if record_core_var == core_variants:
            diplo = record[0]
            full_dip = record[2]
            soln_list1.append(record[0])
            soln_list2.append(record[2])
        else:
            pass

    #return soln_list1

    #print("\nResult:")

    diff_alleles_check = False

    def chkList(lst):
        if len(lst) < 0 :
            diff_alleles_check = True
        diff_alleles_check = all(ele == lst[0] for ele in lst)

        if(diff_alleles_check):
            return("Equal")
        else:
            return("Not equal")


    if len(soln_list1) == 1:
        diplo = "".join(soln_list1)
        res1 = [i for i in range(len(diplo)) if diplo.startswith("_", i)]
        res2 = [i for i in range(len(diplo)) if diplo.startswith(".", i)]
        hap1 = "*" + str (diplo[:res2[0]])
        hap2 = "*" + str (diplo[res1[0]+1:res2[1]])
        allele_res = hap1 + "/" + hap2
        return [soln_list1, diplo, allele_res];
        #print ("\nSupporting variants:")
        #print ("\n" + core_variants + "\n")


    elif len(soln_list1) == 2:
        print(soln_list1)
        diplo1 = soln_list1[0]
        diplo2 = soln_list1[1]
        diplo1_supp_var = soln_list2[0].split(";")
        diplo2_supp_var = soln_list2[1].split(";")
        uniq_diplo1 = []
        uniq_diplo2 = []
        for i in all_variants:
            if i not in diplo1_supp_var:
                uniq_diplo1.append(i)
        
            if i not in diplo2_supp_var:
                uniq_diplo2.append(i)

        #print("\nUnique variants in soln 1: {}".format(len(uniq_diplo1)))
        #print("\nUnique variants in soln 2: {}".format(len(uniq_diplo2)))
            
        if len(uniq_diplo1) < len(uniq_diplo2):
            res1 = [i for i in range(len(diplo1)) if diplo1.startswith("_", i)]
            res2 = [i for i in range(len(diplo1)) if diplo1.startswith(".", i)]
            hap1 = "*" + str (diplo1[:res2[0]])
            hap2 = "*" + str (diplo1[res1[0]+1:res2[1]])
            allele_res =  hap1 + "/" + hap2 
            return [soln_list1, diplo1, allele_res];
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")

        elif len(uniq_diplo1) > len(uniq_diplo2):
            res1 = [i for i in range(len(diplo2)) if diplo2.startswith("_", i)]
            res2 = [i for i in range(len(diplo2)) if diplo2.startswith(".", i)]
            hap1 = "*" + str (diplo2[:res2[0]])
            hap2 = "*" + str (diplo2[res1[0]+1:res2[1]])
            allele_res =  hap1 + "/" + hap2 
            return [soln_list1, diplo2, allele_res];
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")

        # elif len(uniq_diplo1) == len(uniq_diplo2) and (diplo1 == "4.v11_74.v1" and diplo2 == "4.v12_1.v1"):
        #     res1 = [i for i in range(len(diplo2)) if diplo2.startswith("_", i)]
        #     res2 = [i for i in range(len(diplo2)) if diplo2.startswith(".", i)]
        #     hap1 = "*" + str (diplo2[:res2[0]])
        #     hap2 = "*" + str (diplo2[res1[0]+1:res2[1]])
        #     allele_res =  hap1 + "/" + hap2
        #     return [soln_list1, diplo2, allele_res];
    
        # elif len(uniq_diplo1) == len(uniq_diplo2) and diplo2 == "41.v1_65.v1":
        #     res1 = [i for i in range(len(diplo2)) if diplo2.startswith("_", i)]
        #     res2 = [i for i in range(len(diplo2)) if diplo2.startswith(".", i)]
        #     hap1 = "*" + str (diplo2[:res2[0]])
        #     hap2 = "*" + str (diplo2[res1[0]+1:res2[1]])
        #     allele_res =  hap1 + "/" + hap2 
        #     return [soln_list1, diplo2, allele_res];
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")

        # elif len(uniq_diplo1) == len(uniq_diplo2) and (diplo1 == "4.v1_6.v1" and diplo2 == "4.v4_6.v2") :
        #     res1 = [i for i in range(len(diplo1)) if diplo1.startswith("_", i)]
        #     res2 = [i for i in range(len(diplo1)) if diplo1.startswith(".", i)]
        #     hap1 = "*" + str (diplo1[:res2[0]])
        #     hap2 = "*" + str (diplo1[res1[0]+1:res2[1]])
        #     allele_res =  hap1 + "/" + hap2 
        #     return [soln_list1, diplo1, allele_res];
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")
    
    
        else:
            tiebreak1 = []
            tiebreak2 = []
            tiebreak3 = []
            score = []
            for line in f:
                line = line.strip().split()
                #print(line)
                if line[2] == core_variants:
                    tiebreak1.append(line[1])
                    tiebreak2.append(line[3])
                    tiebreak3.append(line[0])
            for full_dip in tiebreak2:
                diplo_supp_gt = full_dip.split(";")
                uniq_gt = []
                for i in all_var_gt:
                    if i not in diplo_supp_gt:
                        uniq_gt.append(i)
                score_dip = len(uniq_gt)
                score.append(score_dip)

            min_score = min(score)    
            #print(score)


            if chkList(score) == "Equal" and soln_list1[0] == "17.v1_4.v1":
                elem = "17.v1_4.v1"
                res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                hap1 = "*" + str (elem[:res2[0]])
                hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                result_dip = hap1 + "/" + hap2
                return [soln_list1, elem, result_dip];


            elif chkList(score) == "Equal":
                amb_soln_set = []
                for elem in soln_list1:
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    amb_soln_set.append(result_dip)
                    #elem_pos = tiebreak1.index(elem)
                    #print ("Solution " + str(elem_pos) + ": " + result_dip)
                allele_res =  " or ".join(amb_soln_set) 
                return [soln_list1, allele_res];
        

            # elif score.count(min_score) > 1 and soln_list1[0] == "11.v1_2.v2":
            #     elem = "11.v1_2.v2"
            #     res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
            #     res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
            #     hap1 = "*" + str (elem[:res2[0]])
            #     hap2 = "*" + str (elem[res1[0]+1:res2[1]])
            #     result_dip = hap1 + "/" + hap2
            #     return [soln_list1, elem, result_dip];



            elif score.count(min_score) > 1:
                index_scores = []
                amb_soln_set = []

                for i in score:
                    if i == min_score:
                        index_scores.append(score.index(i))
                        
                alt_solns = []
                for j in index_scores:
                    elem = tiebreak1[j]
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    alt_solns.append(result_dip)
                
                if chkList(alt_solns) == "Equal":
                    return[soln_list1, alt_solns[0]];
                else:
                    alt_solns = sorted(alt_solns)
                    amb_soln_set.append(alt_solns[0])
                    amb_soln_set.append(alt_solns[-1])
                    allele_res = " or ".join(amb_soln_set)
                    return [soln_list1, allele_res];

                # amb_soln_set = []
                # temp_set = []
                # temp_set.append(tiebreak1[0])
                # temp_set.append(tiebreak1[-1])

                # for elem in temp_set:
                #     res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                #     res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                #     hap1 = "*" + str (elem[:res2[0]])
                #     hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                #     result_dip = hap1 + "/" + hap2
                #     amb_soln_set.append(result_dip)
                    #elem_pos = tiebreak1.index(elem)                                                                                                             
                #print ("Solution " + str(elem_pos) + ": " + result_dip)                                                                                      
                # allele_res = " or ".join(amb_soln_set)
                # return [soln_list1, allele_res];
                #print ("\nSupporting core variants:")
                #print ("\n" + core_variants + "\n")


            else:
                minpos = score.index(min_score)
                best_diplo = tiebreak1[minpos]
                best_cand_haps = tiebreak3[minpos] 
                res1 = [i for i in range(len(best_diplo)) if best_diplo.startswith("_", i)]
                res2 = [i for i in range(len(best_diplo)) if best_diplo.startswith(".", i)]
                hap1 = "*" + str (best_diplo[:res2[0]])
                hap2 = "*" + str (best_diplo[res1[0]+1:res2[1]])
                allele_res =  hap1 + "/" + hap2 
                return [soln_list1, best_cand_haps, allele_res];
                #print ("Supporting core variants:")
                #print ("\n" + core_variants + "\n")


    elif len(soln_list1) == 3:
        diplo1 = soln_list1[0]
        diplo2 = soln_list1[1]
        diplo3 = soln_list1[2]
        diplo1_supp_var = soln_list2[0].split(";") 
        diplo2_supp_var = soln_list2[1].split(";")
        diplo3_supp_var = soln_list2[2].split(";")
        uniq_diplo1 = []
        uniq_diplo2 = []
        uniq_diplo3 = []

        for i in all_variants:
            if i not in diplo1_supp_var:
                uniq_diplo1.append(i)

            if i not in diplo2_supp_var:
                uniq_diplo2.append(i)

            if i not in diplo3_supp_var:
                uniq_diplo3.append(i)


        if len(uniq_diplo1) < len(uniq_diplo2) and len(uniq_diplo1) < len(uniq_diplo3):
            res1 = [i for i in range(len(diplo1)) if diplo1.startswith("_", i)]
            res2 = [i for i in range(len(diplo1)) if diplo1.startswith(".", i)]
            hap1 = "*" + str (diplo1[:res2[0]])
            hap2 = "*" + str (diplo1[res1[0]+1:res2[1]])
            allele_res = hap1 + "/" + hap2
            return [soln_list1, diplo1, allele_res];
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")

        elif len(uniq_diplo1) > len(uniq_diplo2) and len(uniq_diplo2) < len(uniq_diplo3):
            res1 = [i for i in range(len(diplo2)) if diplo2.startswith("_", i)]
            res2 = [i for i in range(len(diplo2)) if diplo2.startswith(".", i)]
            hap1 = "*" + str (diplo2[:res2[0]])
            hap2 = "*" + str (diplo2[res1[0]+1:res2[1]])
            allele_res = hap1 + "/" + hap2
            return [soln_list1, diplo2, allele_res]
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")

        elif len(uniq_diplo1) > len(uniq_diplo2) and len(uniq_diplo2) > len(uniq_diplo3):
            res1 = [i for i in range(len(diplo3)) if diplo3.startswith("_", i)]
            res2 = [i for i in range(len(diplo3)) if diplo3.startswith(".", i)]
            hap1 = "*" + str (diplo3[:res2[0]])
            hap2 = "*" + str (diplo3[res1[0]+1:res2[1]])
            allele_res = hap1 + "/" + hap2
            return [soln_list1, diplo3, allele_res]
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")


        # elif len(uniq_diplo1) == len(uniq_diplo2) == len(uniq_diplo3) and diplo3 == "39.v1_4.v4":
        #     res1 = [i for i in range(len(diplo3)) if diplo3.startswith("_", i)]
        #     res2 = [i for i in range(len(diplo3)) if diplo3.startswith(".", i)]
        #     hap1 = "*" + str (diplo3[:res2[0]])
        #     hap2 = "*" + str (diplo3[res1[0]+1:res2[1]])
        #     allele_res = hap1 + "/" + hap2
        #     return [soln_list1, diplo3, allele_res]
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")


        elif len(uniq_diplo1) == len(uniq_diplo2) == len(uniq_diplo3) or (len(uniq_diplo1) != len(uniq_diplo2) == len(uniq_diplo3)) or (len(uniq_diplo1) == len(uniq_diplo2) != len(uniq_diplo3)):

            tiebreak1 = []
            tiebreak2 = []
            tiebreak3 = []
            score = []
            for line in f:
                line = line.strip().split()
                #print(line)                                                                                                                  
                if line[2] == core_variants:
                    tiebreak1.append(line[1])
                    tiebreak2.append(line[3])
                    tiebreak3.append(line[0])
            for full_dip in tiebreak2:
                diplo_supp_gt = full_dip.split(";")
                uniq_gt = []
                for i in all_var_gt:
                    if i not in diplo_supp_gt:
                        uniq_gt.append(i)
                score_dip = len(uniq_gt)
                score.append(score_dip)

            min_score = min(score)
            # print(score)
        
            if chkList(score) == "Equal":
                amb_soln_set = []
                for elem in tiebreak1:
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    amb_soln_set.append(result_dip)
                    #elem_pos = tiebreak1.index(elem)                                                                                                             
                    #print ("Solution " + str(elem_pos) + ": " + result_dip)                                                                                      
                allele_res = " or ".join(amb_soln_set)
                return [soln_list1, tiebreak1, allele_res];
                #print ("\nSupporting core variants:")
                #print ("\n" + core_variants + "\n")


            else:
                minpos = score.index(min_score)
                best_diplo = tiebreak1[minpos]
                best_cand_haps = tiebreak3[minpos]
                res1 = [i for i in range(len(best_diplo)) if best_diplo.startswith("_", i)]
                res2 = [i for i in range(len(best_diplo)) if best_diplo.startswith(".", i)]
                hap1 = "*" + str (best_diplo[:res2[0]])
                hap2 = "*" + str (best_diplo[res1[0]+1:res2[1]])
                allele_res = hap1 + "/" + hap2
                return [soln_list1, best_cand_haps, allele_res];
                #print ("Supporting core variants:")
                #print ("\n" + core_variants + "\n")



    #print("\nFull diplotype variants:")
    #print("\n" + ";".join(all_var_gt))
