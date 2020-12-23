#!/usr/bin/env python3

import os
import sys
import subprocess
from snv_def_modules import *
from bkg_modules import *


print("--------------------------------------------\n")

print("CYP2C9 Star Allele Calling with StellarPGx\n")

print("--------------------------------------------\n")



database = sys.argv[1]
infile = sys.argv[2]
infile_full = sys.argv[3]
infile_full_gt = sys.argv[4]
infile_spec = sys.argv[5]


cn = 2


supp_core_vars = get_core_variants(infile, cn)

print("\nSample core variants:")
print(supp_core_vars)


snv_def_calls = cand_snv_allele_calling(database, infile, infile_full, infile_full_gt, infile_spec, cn)

if snv_def_calls == None:

    bac_alleles = get_backgroud_alleles(database, supp_core_vars)

    if bac_alleles == None:
        print("\nResult:")
        print("Possible novel allele or suballele present: interpret with caution")


    else:
        print("\nResult:")
        print("Possible novel allele or suballele present: interpret with caution; experimental validation and expert review through PharmVar is recommended")
        print("\nLikely background alleles:")
        print("[" + bac_alleles + "]")

    sys.exit()


snv_cand_alleles = snv_def_calls[0]

print("\nCandidate alleles:")
print(snv_cand_alleles)


snv_def_alleles = snv_def_calls[-1]

dip_variants = get_all_vars_gt(infile_full_gt)


print("\nResult:")

print(snv_def_alleles)
