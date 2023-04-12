import os
import sys
import subprocess
from snv_def_modules import *
from bkg_modules import *
from sv_modules import *


print("--------------------------------------------\n")

print("CYP2D6 Star Allele Calling with StellarPGx\n")

print("--------------------------------------------\n")



database = sys.argv[1]
infile = sys.argv[2]
infile_full = sys.argv[3]
infile_full_gt = sys.argv[4]
infile_spec = sys.argv[5]
sv_del = sys.argv[6]
sv_dup = sys.argv[7]
cov_file = sys.argv[8]
hap_dbs = sys.argv[9]
act_score = sys.argv[10]


cn = get_total_CN(cov_file)[0]

print("Initially computed CN = {}".format(cn))

supp_core_vars = get_core_variants(infile, cn)

print("\nSample core variants:")
print(supp_core_vars)


snv_def_calls = cand_snv_allele_calling(database, infile, infile_full, infile_full_gt, infile_spec, cn)
bac_alleles = get_backgroud_alleles(database, supp_core_vars)

if snv_def_calls == None and "*139/*4" in bac_alleles:
    snv_def_calls = [['1.v1_4.v11'], '1.v1_4.v11', '*1/*4']
elif snv_def_calls == None and "*139/*3" in bac_alleles:
    snv_def_calls = [['1.v1_3.v1'], '1.v1_3.v1', '*1/*3']

if snv_def_calls == None:

    bac_alleles = get_backgroud_alleles(database, supp_core_vars)

    if int(cn) == 0:
        print("\nResult:")
        print("*5/*5")

    elif bac_alleles == None:
        print("\nResult:")
        print("Possible novel allele or suballele present: interpret with caution")

    elif bac_alleles != None and int(cn) < 2:
        bac_alleles = bac_alleles[0].split("/")
        bac_alleles1 = bac_alleles[0] + "/" + "*5"
        print("\nResult:")
        print("Possible novel allele or suballele present: interpret with caution; experimental validation and expert review through PharmVar is recommended")
        print("\nLikely background alleles:")
        print("[" + bac_alleles1 + "]")

    else:
        print("\nCandidate alleles:")
        print("[" + bac_alleles[-1] + "]")

        print("\nResult:")
        print("Possible novel allele or suballele present: interpret with caution; experimental validation and expert review through PharmVar is recommended")
        print("\nLikely background alleles:")
        print("[" + bac_alleles[0] + "]")

    print("\nActivity score:")
    print("Indeterminate")

    print("\nMetaboliser status:")
    print("Indeterminate")


    sys.exit()

best_diplos = snv_def_calls[0]

print("\nCandidate alleles:")
print(best_diplos)


snv_def_alleles = snv_def_calls[-1]

if "or" in snv_def_alleles:
    pass
else:
    snv_cand_alleles = snv_def_calls[1]

dip_variants = get_all_vars_gt(infile_full_gt)


print("\nResult:")


# cn = get_total_CN(cov_file)[0]
av_cov = get_total_CN(cov_file)[3]
cn_in1_3pr = get_total_CN(cov_file)[2]
cn_ex9_3pr = get_total_CN(cov_file)[4]
in1_3pr_float = get_total_CN(cov_file)[5]

cov_in4_3pr = get_total_CN(cov_file)[6]
cov_5pr_in4 = get_total_CN(cov_file)[7]
cn_2d7_ex9 = get_total_CN(cov_file)[8]
cn_2d7_in4_in8 = get_total_CN(cov_file)[9]
cov_2d7_ex2_in8 = get_total_CN(cov_file)[10]
cov_2d7_5pr_in1 = get_total_CN(cov_file)[11]
cov_2d7_hyb = get_total_CN(cov_file)[12]
cn_hyb_68 = get_total_CN(cov_file)[13]
cn_hyb_36 = get_total_CN(cov_file)[14]
cn_hyb_61 = get_total_CN(cov_file)[15]
cn_hyb_63 = get_total_CN(cov_file)[16]
cn_hyb_13 = get_total_CN(cov_file)[17]
cn_hyb_83 = get_total_CN(cov_file)[18]

# print(float(cn_ex9_3pr))

gene_alleles = ""


if snv_def_alleles != '*2/*2':
    in_list = dup_test_init(sv_dup, av_cov)



if cn == '2' and snv_def_alleles == '*4/*4':

    test_68 = hyb_test_5_68_4(sv_del, in1_3pr_float, av_cov)

    if test_68 == 'norm_art':
        pass
    elif test_68 == 'del_hyb':
        snv_def_alleles = (snv_def_alleles.replace('*4', '*5', 1)).replace('*4', '*68+*4')

    gene_alleles = snv_def_alleles
    print(snv_def_alleles)
    cn_recalib = cn_hyb_68
    if cn_recalib == int(0):
        print("\nRecalibrated CN = ", cn)
    else:
        print("\nRecalibrated CN = ", cn_recalib)



elif cn == '2':
    if 'or' in snv_def_alleles:

        print (snv_def_alleles)
    else:
        snv_def_alleles = snv_def_alleles.split("/")


        if snv_def_alleles[0] == '*2' or snv_def_alleles[1] == '*2':
            ind_star2 = snv_def_alleles.index('*2')
            ind_other = 1 - ind_star2
            def_1 = snv_def_alleles[0]
            def_2 = snv_def_alleles[1]
            lst = ['*2x2', '*4x2']
            lst2 = ['*13', '*13+*2']

            test_13_2_v1 = hybrid_13_2_v1(cov_in4_3pr, cov_5pr_in4)

            test_13_2_v2 = hybrid_13_2_v2(cov_2d7_ex2_in8, cov_2d7_5pr_in1)

            if test_13_2_v1 == 'norm_var':
                pass

            elif test_13_2_v1 == 'hyb_13_2' and cn_hyb_13 == '2':
                gene_alleles = snv_def_alleles[ind_other] + "/" + "*13+*2"
                print(gene_alleles)

            elif test_13_2_v2 == 'hyb_13_2_v2' and cn_hyb_13 == '2':
                gene_alleles = snv_def_alleles[ind_other] + "/" + "*13"
                print(gene_alleles)
            
            elif test_13_2_v1 == 'hyb_13_2' and cn_hyb_13 == '3':
                gene_alleles = snv_def_alleles[ind_other] + "/" + "*13+*2"
                print(gene_alleles)
                cn_recalib = cn_hyb_13
                if cn_recalib == int(0):
                    print("\nRecalibrated CN = ", cn)
                else:
                    print("\nRecalibrated CN = ", cn_recalib)

            elif test_13_2_v2 == 'hyb_13_2_v2' and cn_hyb_13 == '3':
                gene_alleles = snv_def_alleles[ind_other] + "/" + "*13+*2"
                print(gene_alleles)
                cn_recalib = cn_hyb_13
                if cn_recalib == int(0):
                    print("\nRecalibrated CN = ", cn)
                else:
                    print("\nRecalibrated CN = ", cn_recalib)

            test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

            if test_68 == 'norm_dup':
                if any([x in gene_alleles for x in lst2]):
                    pass
                else:
                    gene_alleles = "/".join(snv_def_alleles)
                    print(gene_alleles)
            
            elif test_68 == 'hyb_68':
                gene_alleles = snv_def_alleles[ind_star2] + "/" + "*68+*4"
                print(gene_alleles)
                cn_recalib = cn_hyb_68
                if cn_recalib == int(0):
                    print("\nRecalibrated CN = ", cn)
                else:
                    print("\nRecalibrated CN = ", cn_recalib)


        elif snv_def_alleles[0] == '*39' or snv_def_alleles[1] == '*39':
            ind_star2 = snv_def_alleles.index('*39')
            ind_other = 1 - ind_star2

            test_83_single = hybrid_test_83_single(sv_dup, cn, av_cov, cn_ex9_3pr)

            if test_83_single == 'norm_star39':
                gene_alleles = "/".join(snv_def_alleles)
                print(gene_alleles)

            elif test_83_single == 'hyb_83_single' and cn_hyb_83 == '2':
                gene_alleles = snv_def_alleles[ind_other] + "/" + "*83"
                print(gene_alleles)
            
            test_83 = hybrid_test_83(sv_dup, cn, av_cov, cn_ex9_3pr)
            
            if test_83_single == 'norm_star39':
                gene_alleles = "/".join(snv_def_alleles)
                print(gene_alleles)
                
            elif test_83 == 'hyb_83' and cn_hyb_83 == '3':
                gene_alleles = snv_def_alleles[ind_other] + "x2" + "/" + "*83"
                print(gene_alleles)
                cn_recalib = cn_hyb_83
                if cn_recalib == int(0):
                    print("\nRecalibrated CN = ", cn)
                else:
                    print("\nRecalibrated CN = ", cn_recalib)


        elif snv_def_alleles[0] == '*10' or snv_def_alleles[1] == '*10':
            ind_star2 = snv_def_alleles.index('*10')
            ind_other = 1 - ind_star2

            test_36_single = hybrid_test_36_single(sv_dup, cn, av_cov, cn_ex9_3pr)

            if test_36_single == 'norm_star10':
                gene_alleles = "/".join(snv_def_alleles)
                print(gene_alleles)

            elif test_36_single == 'hyb_36_single' and cn_hyb_36 == '2':
                gene_alleles = snv_def_alleles[ind_other] + "/" + "*36"
                print(gene_alleles)
            
            test_36 = hybrid_test_36(sv_dup, cn, av_cov, cn_ex9_3pr, cn_2d7_ex9, cn_2d7_in4_in8)
            
            if test_36 == 'norm_dup':
                gene_alleles = "/".join(snv_def_alleles)
                print(gene_alleles)
                
            elif test_36 == 'hyb_36_10' and cn_hyb_36 == '3':
                gene_alleles = snv_def_alleles[ind_other] + "/" + "*36+*10"
                print(gene_alleles)
                cn_recalib = cn_hyb_36
                if cn_recalib == int(0):
                    print("\nRecalibrated CN = ", cn)
                else:
                    print("\nRecalibrated CN = ", cn_recalib)
            
            elif test_36 == 'hyb_36_36' and cn_hyb_36 == '3':
                gene_alleles = snv_def_alleles[ind_other] + "/" + "*36x2"
                print(gene_alleles)
                cn_recalib = cn_hyb_36
                if cn_recalib == int(0):
                    print("\nRecalibrated CN = ", cn)
                else:
                    print("\nRecalibrated CN = ", cn_recalib)


        elif snv_def_alleles[0] == '*4' or snv_def_alleles[1] == '*4':
            def_1 = snv_def_alleles[0]
            def_2 = snv_def_alleles[1]
            ind_star2 = snv_def_alleles.index('*4')
            ind_other = 1 - ind_star2

            test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

            if test_68 == 'norm_dup':
                gene_alleles = "/".join(snv_def_alleles)
                print(gene_alleles)
            
            elif test_68 == 'hyb_68':
                gene_alleles = snv_def_alleles[ind_other] + "/" + "*68+*4"
                print(gene_alleles)
                cn_recalib = cn_hyb_68
                if cn_recalib == int(0):
                    print("\nRecalibrated CN = ", cn)
                else:
                    print("\nRecalibrated CN = ", cn_recalib)

        else:
            gene_alleles = "/".join(snv_def_alleles)
            print(gene_alleles)


elif cn == '0':
    del_confirm = del_test(sv_del)
    if del_confirm == '*5/*5':
        gene_alleles = del_confirm
        print (gene_alleles)

    elif del_confirm == '*5':
        gene_alleles = del_confirm + "/" + "*other"
        print(gene_alleles)

    else:
        gene_alleles = "*5/*5"
        print(gene_alleles)

elif cn == '1':
    del_confirm = del_test(sv_del)

    if "or" in snv_def_alleles and del_confirm == 'None':
        print (snv_def_alleles + "\t" + "Possible CYP2D6 gene deletion (*5) present")

    elif "or" not in snv_def_alleles and del_confirm == 'None':
        snv_def_alleles = snv_def_alleles.split("/")
        snv_cand_alleles = "".join(snv_cand_alleles)
        snv_cand_alleles = snv_cand_alleles.split("_")

        if snv_def_alleles[0] == snv_def_alleles[1]:
            gene_alleles = snv_def_alleles[0] + "/" + "*5"
            print(gene_alleles)

        elif snv_def_alleles[0] != snv_def_alleles[1]:
            samp_allele1 = del_adv_test(hap_dbs, snv_cand_alleles[0], snv_cand_alleles[1], snv_def_alleles[0], snv_def_alleles[1], supp_core_vars)

            gene_alleles = samp_allele1 + "/" + "*5"
            print(gene_alleles)

    else:
        snv_def_alleles = snv_def_alleles.split("/")
        snv_cand_alleles = "".join(snv_cand_alleles)
        snv_cand_alleles = snv_cand_alleles.split("_")

        if snv_def_alleles[0] == snv_def_alleles[1]:

            if del_confirm == "*5/*5":
                del_confirm = "*5"
            else:
                del_confirm = "*5"

            gene_alleles = del_confirm + "/" + snv_def_alleles[0]
            print(gene_alleles)

        elif snv_def_alleles[0] != snv_def_alleles[1]:
            samp_allele1 = del_adv_test(hap_dbs, snv_cand_alleles[0], snv_cand_alleles[1], snv_def_alleles[0], snv_def_alleles[1], supp_core_vars)

            if del_confirm == "*5/*5":
                del_confirm = "*5"
            else:
                del_confirm = "*5"

            gene_alleles = del_confirm + "/" + samp_allele1
            print(gene_alleles)

elif (int(cn) == 3 or int(cn) == 4) and snv_def_alleles != None:

    orig = snv_def_alleles

    if "or" in snv_def_alleles:
        print (snv_def_alleles + "\t" + "Duplication present")

    else:
        snv_def_alleles = snv_def_alleles.split("/")
        snv_cand_alleles = "".join(snv_cand_alleles)
        snv_cand_alleles = snv_cand_alleles.split("_")


        if snv_def_alleles[0] == '*90' or snv_def_alleles[1] == '*90':

            alt_allele_ind = 1 - snv_def_alleles.index('*90')
            alt_allele = snv_def_alleles[alt_allele_ind]
            sp_allele = tandem_90_1(in_list, alt_allele, cn)


            sp_allele1 = sp_allele.split("/")


            if "*10x2" in sp_allele1:

                test_36 = hybrid_test_36(sv_dup, cn, av_cov, cn_ex9_3pr, cn_2d7_ex9, cn_2d7_in4_in8)

                if test_36 == 'norm_dup':
                    pass

                elif test_36 == 'hyb_36_10':
                    sp_allele = sp_allele.replace('*10x2', '*36+*10')

                elif test_36 == 'hyb_36_36':
                    sp_allele = sp_allele.replace('*10x2', '*36x2')

            gene_alleles = sp_allele
            print(sp_allele)


        elif snv_def_alleles[0] == '*57' or snv_def_alleles[1] == '*57':

            alt_allele_ind = 1 - snv_def_alleles.index('*57')
            alt_allele = snv_def_alleles[alt_allele_ind]
            sp_allele = tandem_57_10(in_list, alt_allele, cn)

            print(sp_allele)


        elif snv_def_alleles[0] != snv_def_alleles[1]:


            phased_dup = dup_test_cn_3_4(sv_dup, hap_dbs, snv_cand_alleles[0], snv_cand_alleles[1], snv_def_alleles[0], snv_def_alleles[1], cn, av_cov, in_list)

            if phased_dup == 'check':
                phased_dup == 'No_call'

            else:
                pass

            phased_dup1 = phased_dup.split("/")


            if '*4' in phased_dup1:
                count1 = phased_dup1.count('*4')
                a_ind1 = phased_dup1.index('*4')
                a_ind2 = 1 - a_ind1
                other_hap = phased_dup1[a_ind2]
                lst = ['*4x2', '*4x3']

                if count1 == 1:

                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass
                    elif test_68 == 'hyb_68':
                        if any([x in phased_dup1 for x in lst]):
                            pass
                            
                        elif int(cn_hyb_68) >= 4 and 'x' not in other_hap:
                            phased_dup = phased_dup.replace('*4', '*68+*4')
                            phased_dup = phased_dup.replace(other_hap, (other_hap + 'x2'))
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                        else:
                            phased_dup = phased_dup.replace('*4', '*68+*4')
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                elif count1 == 2:
                    pass

            if '*4x2' in phased_dup1:
                count1 = phased_dup1.count('*4x2')
                a_ind1 = phased_dup1.index('*4x2')
                a_ind2 = 1 - a_ind1
                other_hap = phased_dup1[a_ind2]

                if count1 == 1:

                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass
                    
                    elif test_68 == 'hyb_68':
                            
                        if int(cn_hyb_68) >= 4 and 'x' not in other_hap:
                            phased_dup = phased_dup.replace('*4x2', '*68+*4')
                            phased_dup = phased_dup.replace(other_hap, (other_hap + 'x2'))
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                        else:
                            phased_dup = phased_dup.replace('*4x2', '*68+*4')
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                elif count1 == 2:
                    pass

            if '*4x3' in phased_dup1:
                count1 = phased_dup1.count('*4x3')
                a_ind1 = phased_dup1.index('*4x3')
                a_ind2 = 1 - a_ind1
                other_hap = phased_dup1[a_ind2]

                if count1 == 1:

                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass
                    elif test_68 == 'hyb_68':

                        if int(cn_hyb_68) >= 4 and 'x' not in other_hap:
                            phased_dup = phased_dup.replace('*4x3', '*68+*4')
                            phased_dup = phased_dup.replace(other_hap, (other_hap + 'x2'))
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                        else:
                            phased_dup = phased_dup.replace('*4x3', '*68+*4')
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                elif count1 == 2:
                    pass

            if '*10' in phased_dup1:
                count2 = phased_dup1.count('*10')
                b_ind1 = phased_dup1.index('*10')
                b_ind2 = 1 - b_ind1


                if count2 == 1:
                    test_36 = hybrid_test_36(sv_dup, cn, av_cov, cn_ex9_3pr)


                    if test_36 == 'norm_dup':
                        pass

                    elif test_36 == 'hyb_36_10' and cn == '3' and cn_hyb_36 == '4':
                        phased_dup = phased_dup.replace('*10', '*36+*10')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

                    elif test_36 == 'hyb_36_36' and cn == '3' and cn_hyb_36 == '4':
                        phased_dup = phased_dup.replace('*10', '*36x2')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

            if '*10x2' in phased_dup1:
                count2 = phased_dup1.count('*10x2')
                b_ind1 = phased_dup1.index('*10x2')
                b_ind2 = 1 - b_ind1


                if count2 == 1:
                    test_36 = hybrid_test_36(sv_dup, cn, av_cov, cn_ex9_3pr)


                    if test_36 == 'norm_dup':
                        pass

                    elif test_36 == 'hyb_36_10':
                        phased_dup = phased_dup.replace('*10x2', '*36+*10')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

                    elif test_36 == 'hyb_36_36':
                        phased_dup = phased_dup.replace('*10x2', '*36x2')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)


            if '*10x3' in phased_dup1:
                count3 = phased_dup1.count('*10x3')
                c_ind1 = phased_dup1.index('*10x3')
                c_ind2 = 1 - c_ind1

                if count3 == 1:
                    test_36 = hybrid_test_36_mod(sv_dup, cn, av_cov, cn_ex9_3pr)

                    if test_36 == 'norm_mt':
                        pass

                    elif test_36 == 'hyb_36_10': 
                        phased_dup = phased_dup.replace('*10x3', '*36+*10x2')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)


                    elif test_36 == 'hyb_36_36':
                        phased_dup = phased_dup.replace('*10x3', '*36x2+*10')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)


            if '*1x3' in phased_dup1:
                count2 = phased_dup1.count('*1x3')
                b_ind1 = phased_dup1.index('*1x3')
                b_ind2 = 1 - b_ind1


                if count2 == 1:
                    test_83 = hybrid_test_83(sv_dup, cn, av_cov, cn_ex9_3pr)


                    if test_83 == 'norm_star39':
                        pass

                    elif test_83 == 'hyb_83':
                        phased_dup = phased_dup.replace('*1x3', '*1x2+*83')
                        cn_recalib = cn_hyb_83
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)


            if '*2x2' in phased_dup1:
                count2 = phased_dup1.count('*2x2')
                b_ind1 = phased_dup1.index('*2x2')
                b_ind2 = 1 - b_ind1
                lst = ['*4', '*4x2', '*4x3']

                if count2 == 1:
                    test_13_2_v1 = hybrid_13_2_v1(cov_in4_3pr, cov_5pr_in4)
                    test_13_2_v2 = hybrid_13_2_v2(cov_2d7_ex2_in8, cov_2d7_5pr_in1)

                    if test_13_2_v1 == 'norm_var':
                        pass

                    elif test_13_2_v2 == 'norm_var':
                        pass

                    elif test_13_2_v1 == 'hyb_13_2':
                        phased_dup = phased_dup1[b_ind2] + "/" + '*13+*2'
                        cn_recalib = cn_hyb_13
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

                    elif test_13_2_v2 == 'hyb_13_2_v2':
                        phased_dup = phased_dup1[b_ind2] + "/" + '*13+*2'
                        cn_recalib = cn_hyb_13
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass
                    elif test_68 == 'hyb_68':
                        if any([x in phased_dup1 for x in lst]):
                            pass

                        elif int(cn_hyb_68) >= 4 and 'x' not in other_hap:
                            phased_dup = phased_dup.replace('*2x2', '*68+*4', 1)
                            phased_dup = phased_dup.replace(other_hap, (other_hap + 'x2'))
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                        else:
                            phased_dup = phased_dup.replace('*2x2', '*68+*4', 1)
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                elif count1 == 2:
                    pass


            if '*2x3' in phased_dup1:
                count1 = phased_dup1.count('*2x3')
                a_ind1 = phased_dup1.index('*2x3')
                a_ind2 = 1 - a_ind1
                other_hap = phased_dup1[a_ind2]
                lst = ['*4', '*4x2', '*4x3']

                if count1 == 1:

                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass
                    elif test_68 == 'hyb_68':
                        if any([x in phased_dup1 for x in lst]):
                            pass

                        elif int(cn_hyb_68) >= 4 and 'x' not in other_hap:
                            phased_dup = phased_dup.replace('*2x3', '*68+*4', 1)
                            phased_dup = phased_dup.replace(other_hap, (other_hap + 'x2'))
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                        else:
                            phased_dup = phased_dup.replace('*2x3', '*68+*4', 1)
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                elif count1 == 2:
                    pass


            gene_alleles = phased_dup
            print(phased_dup)


        elif snv_def_alleles[0] == snv_def_alleles[1]:


            rt_2 = int(cn) - 1

            phased_dup = (snv_def_alleles[0] + "/" + snv_def_alleles[1] + "x" + str(rt_2))

            phased_dup1 = phased_dup.split("/")


            if '*4x2' in phased_dup1:
                count1 = phased_dup1.count('*4x2')
                a_ind1 = phased_dup1.index('*4x2')
                a_ind2 = 1 - a_ind1


                if count1 == 1:
                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass
                    
                    elif test_68 == 'hyb_68':

                        if int(cn_hyb_68) >= 4 and 'x' not in other_hap:
                            phased_dup = phased_dup.replace('*4x2', '*68+*4')
                            phased_dup = phased_dup.replace(other_hap, (other_hap + 'x2'))
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                        else:
                            phased_dup = phased_dup.replace('*4x2', '*68+*4')
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)


                elif count1 == 2:

                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass
                    
                    elif test_68 == 'hyb_68':

                        if int(cn_hyb_68) >= 4 and 'x' not in other_hap:
                            phased_dup = phased_dup.replace('*4x2', '*68+*4')
                            phased_dup = phased_dup.replace(other_hap, (other_hap + 'x2'))
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                        else:
                            phased_dup = phased_dup.replace('*4x2', '*68+*4')
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

            if '*4x3' in phased_dup1:
                count1 = phased_dup1.count('*4x3')
                a_ind1 = phased_dup1.index('*4x3')
                a_ind2 = 1 - a_ind1
                other_hap = phased_dup1[a_ind2]

                if count1 == 1:

                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass
                    
                    elif test_68 == 'hyb_68':
                        
                        if int(cn_hyb_68) >= 4 and 'x' not in other_hap:
                            phased_dup = phased_dup.replace('*4x3', '*68+*4')
                            phased_dup = phased_dup.replace(other_hap, (other_hap + 'x2'))
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                        else:
                            phased_dup = phased_dup.replace('*4x3', '*68+*4')
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                elif count1 == 2:
                    pass

            if '*10' in phased_dup1:
                count2 = phased_dup1.count('*10')
                b_ind1 = phased_dup1.index('*10')
                b_ind2 = 1 - b_ind1


                if count2 == 1:
                    test_36 = hybrid_test_36(sv_dup, cn, av_cov, cn_ex9_3pr)


                    if test_36 == 'norm_dup':
                        pass

                    elif test_36 == 'hyb_36_10' and cn == '3' and cn_hyb_36 == '4':
                        phased_dup = phased_dup.replace('*10', '*36+*10')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

                    elif test_36 == 'hyb_36_36' and cn == '3' and cn_hyb_36 == '4':
                        phased_dup = phased_dup.replace('*10', '*36x2')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

            if '*10x2' in phased_dup1:
                count2 = phased_dup1.count('*10x2')
                b_ind1 = phased_dup1.index('*10x2')
                b_ind2 = 1 - b_ind1

                if count2 == 1:
                    test_36 = hybrid_test_36(sv_dup, cn, av_cov, cn_ex9_3pr)

                    if test_36 == 'norm_dup':
                        pass

                    elif test_36 == 'hyb_36_10':
                        phased_dup = phased_dup.replace('*10x2', '*36+*10')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

                    elif test_36 == 'hyb_36_36':
                        phased_dup = phased_dup.replace('*10x2', '*36x2')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

            if '*10x3' in phased_dup1:
                count3 = phased_dup1.count('*10x3')
                c_ind1 = phased_dup1.index('*10x3')
                c_ind2 = 1 - c_ind1

                if count3 == 1:
                    test_36 = hybrid_test_36_mod(sv_dup, cn, av_cov, cn_ex9_3pr)


                    if test_36 == 'norm_mt':
                        pass

                    elif test_36 == 'hyb_36_10':
                        phased_dup = phased_dup.replace('*10x3', '*36+*10x2')
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)


                    elif test_36 == 'hyb_36_36':
                        phased_dup = '*36+*10/*36+*10'      
                        cn_recalib = cn_hyb_36
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)


            if '*1x3' in phased_dup1:
                count2 = phased_dup1.count('*1x3')
                b_ind1 = phased_dup1.index('*1x3')
                b_ind2 = 1 - b_ind1


                if count2 == 1:
                    test_83 = hybrid_test_83(sv_dup, cn, av_cov, cn_ex9_3pr)


                    if test_83 == 'norm_star39':
                        pass

                    elif test_83 == 'hyb_83':
                        phased_dup = phased_dup.replace('*1x3', '*1x2+*83')
                        cn_recalib = cn_hyb_83
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)


            if '*2x2' in phased_dup1:
                count2 = phased_dup1.count('*2x2')
                b_ind1 = phased_dup1.index('*2x2')
                b_ind2 = 1 - b_ind1
                lst = ['*4', '*4x2', '*4x3']

                if count2 == 1:
                    test_13_2_v1 = hybrid_13_2_v1(cov_in4_3pr, cov_5pr_in4)
                    test_13_2_v2 = hybrid_13_2_v2(cov_2d7_ex2_in8, cov_2d7_5pr_in1)

                    if test_13_2_v1 == 'norm_var':
                        pass

                    elif test_13_2_v2 == 'norm_var':
                        pass

                    elif test_13_2_v1 == 'hyb_13_2':
                        phased_dup = phased_dup1[b_ind2] + "/" + '*13+*2'
                        cn_recalib = cn_hyb_13
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

                    elif test_13_2_v2 == 'hyb_13_2_v2':
                        phased_dup = phased_dup1[b_ind2] + "/" + '*13+*2'
                        cn_recalib = cn_hyb_13
                        if cn_recalib == int(0):
                            print("\nRecalibrated CN = ", cn)
                        else:
                            print("\nRecalibrated CN = ", cn_recalib)

                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass

                    elif test_68 == 'hyb_68':
                        if any([x in phased_dup1 for x in lst]):
                            pass
                        elif int(cn_hyb_68) >= 4 and 'x' not in other_hap:
                            phased_dup = phased_dup.replace('*2x2', '*68+*4', 1)
                            phased_dup = phased_dup.replace(other_hap, (other_hap + 'x2'))
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                        else:
                            phased_dup = phased_dup.replace('*2x2', '*68+*4', 1)
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                elif count1 == 2:

                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass
                    elif test_68 == 'hyb_68':
                        if any([x in phased_dup1 for x in lst]):
                            pass
                        elif int(cn_hyb_68) >= 4 and 'x' not in other_hap:
                            phased_dup = phased_dup.replace('*2x2', '*68+*4', 1)
                            phased_dup = phased_dup.replace(other_hap, (other_hap + 'x2'))
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                        else:
                            phased_dup = phased_dup.replace('*2x2', '*68+*4', 1)
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)


            if '*2x3' in phased_dup1:
                count1 = phased_dup1.count('*2x3')
                a_ind1 = phased_dup1.index('*2x3')
                a_ind2 = 1 - a_ind1
                other_hap = phased_dup1[a_ind2]
                lst = ['*4', '*4x2', '*4x3']

                if count1 == 1:

                    test_68 = hybrid_test_68(cn, av_cov, in1_3pr_float, cov_2d7_hyb, cn_hyb_68)

                    if test_68 == 'norm_dup':
                        pass
                    elif test_68 == 'hyb_68':
                        if any([x in phased_dup1 for x in lst]):
                            pass
                        else:
                            phased_dup = phased_dup.replace('*2x3', '*68+*4', 1)
                            cn_recalib = cn_hyb_68
                            if cn_recalib == int(0):
                                print("\nRecalibrated CN = ", cn)
                            else:
                                print("\nRecalibrated CN = ", cn_recalib)

                elif count1 == 2:
                    pass

            gene_alleles = phased_dup
            print(phased_dup)


elif int(cn) > 4 and snv_def_alleles != None:

    if "or" in snv_def_alleles:
        print (snv_def_alleles + "\t" + "Duplication present")

    else:
        snv_def_alleles = snv_def_alleles.split("/")
        snv_cand_alleles = "".join(snv_cand_alleles)
        snv_cand_alleles = snv_cand_alleles.split("_")

        if snv_def_alleles[0] != snv_def_alleles[1]:

            phased_dup = dup_test_cn_n(sv_dup, hap_dbs, snv_cand_alleles[0], snv_cand_alleles[1], snv_def_alleles[0], snv_def_alleles[1], cn, av_cov, in_list)

            if phased_dup == 'check':
                phased_dup = 'No_call'

            else:
                pass

            phased_dup1 = phased_dup.split("/")

            if '*10x4' in phased_dup1:
                count3 = phased_dup1.count('*10x4')
                c_ind1 = phased_dup1.index('*10x4')
                c_ind2 = 1 - c_ind1

                if count3 == 1:
                    test_36 = hybrid_test_36_multi(sv_dup, cn, av_cov, cn_ex9_3pr)


                    if test_36 == 'norm_mt':
                        pass

                    elif test_36 == 'hyb_36_10':
                        phased_dup = phased_dup.replace('*10x4', '*36+*10x3')


                    elif test_36 == 'hyb_36_36':
                        phased_dup = phased_dup.replace('*10x4', '*36x2+*10x2')

                    elif test_36 == 'hyb_36_36_36':
                        phased_dup = phased_dup.replace('*10x4','*36x3+*10')

                    else:
                        phased_dup = "No_call"


            elif '*10x3' in phased_dup1:
                count3 = phased_dup1.count('*10x3')
                c_ind1 = phased_dup1.index('*10x3')
                c_ind2 = 1 - c_ind1

                if count3 == 1:
                    test_36 = hybrid_test_36_multi(sv_dup, cn, av_cov, cn_ex9_3pr)

                    if test_36 == 'norm_mt':
                        pass

                    elif test_36 == 'hyb_36_10':
                        phased_dup = phased_dup.replace('*10x3', '*36+*10x2')

                    elif test_36 == 'hyb_36_36':
                        phased_dup = phased_dup.replace('*10x3', '*36x2+*10')

                    elif test_36 == 'hyb_36_36_36':
                        phased_dup = phased_dup.replace('*10x3','*36x3')

            elif '*10x' in phased_dup1:
                phased_dup = 'No_call'



        elif snv_def_alleles[0] == snv_def_alleles[1]:
            rt_2 = int(cn) - 1
            phased_dup = (snv_def_alleles[0] + "/" + snv_def_alleles[1] + "x" + str(rt_2))

            if phased_dup == 'check':
                phased_dup = 'No_call'

            else:
                pass

            phased_dup1 = phased_dup.split("/")

            if '*10x4' in phased_dup1:
                count3 = phased_dup1.count('*10x4')
                c_ind1 = phased_dup1.index('*10x4')
                c_ind2 = 1 - c_ind1

                if count3 == 1:
                    test_36 = hybrid_test_36_multi(sv_dup, cn, av_cov, cn_ex9_3pr)


                    if test_36 == 'norm_mt':
                        pass

                    elif test_36 == 'hyb_36_10':
                        phased_dup = phased_dup.replace('*10x4', '*36+*10x3')


                    elif test_36 == 'hyb_36_36':
                        phased_dup = '*36+*10/*36+*10x2'

                    elif test_36 == 'hyb_36_36_36':
                        phased_dup = '*36+*10/*36x2+*10'

                    else:
                        phased_dup = "No_call"

            elif phased_dup1[0].startswith('*10x') or phased_dup1[1].startswith('*10x'):

                if phased_dup1[0].startswith('*10x'):
                    dup_10_hyb = phased_dup1[0]

                elif phased_dup1[1].startswith('*10x'):
                    dup_10_hyb = phased_dup1[1]

                cn_star10 = dup_10_hyb[(dup_10_hyb.find('x') + 1 ):]

                test_36 = hybrid_test_36_multi_10(sv_dup, cn, av_cov, cn_ex9_3pr, cn_star10)

                if test_36 == 'norm_mt':
                    pass

                elif test_36 == 'check':
                    phased_dup = 'No_call'

                else:
                    c_ind1 = phased_dup1.index(dup_10_hyb)
                    c_ind2 = 1 - c_ind1
                    phased_dup = str(phased_dup1[c_ind2]) + "/" + test_36


        gene_alleles = phased_dup
        print(phased_dup)


elif int(cn) > 2 and snv_def_alleles == None:

    print("Possible rare CYP2D6/2D7 hybrid present")


print("\nActivity score:")

score_list = []

score_list1 = []
score_list2 = []
score_list3 = []

allele_dict = {}

def get_ac_score(act_score, star_alleles):
    for line in open(act_score, "r"):
        line = line.strip().split()
        score_list.append(line)

    for i in score_list:
        allele_dict[i[0]] = i[1]

    star_alleles = star_alleles.replace("/", "+")
    star_alleles = star_alleles.split("+")

    for elem in star_alleles:
        if "x" not in elem:
            m_allele = elem
            n_allele = "1"
        elif "x" in elem:
            index1 = elem.find("x")
            m_allele = elem[:index1]
            n_allele = elem[index1+1:]

        p_allele = allele_dict[m_allele] + "_" + n_allele
        p_allele = p_allele.split("_")
        score_list1.append(p_allele)

    for i in score_list1:
        score_list2.append(i[0])

    if "n" in score_list2:
        return "Indeterminate"

    else:
        for i in score_list1:
            score_list3.append(float(i[0])*float(i[1]))

        total_a_score = sum(score_list3)
        return total_a_score


if gene_alleles in ["",'No_call','check']:
    ac_score = "Indeterminate"
    print(ac_score)


elif gene_alleles != "":
    ac_score = get_ac_score(act_score, gene_alleles)
    print(ac_score)


print("\nMetaboliser status:")

if ac_score == "Indeterminate":
    print ("Indeterminate")

elif ac_score == 0:
    print("Poor metaboliser (PM)")

elif 0 < ac_score < 1.25:
    print("Intermediate metaboliser (IM)")

elif 1.25 <= ac_score <= 2.25:
     print("Normal metaboliser (NM)")

elif ac_score > 2.25:
    print("Ultrarapid metaboliser (UM)")
