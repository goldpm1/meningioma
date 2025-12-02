import pandas as pd
import numpy as np
import vcf


DIR="/data/project/Meningioma/04.mutect"

#for Sample_ID in [ "221102", "221202", "230127", "230303", "230323_2", "230405_2", "230419", "230526", "230822", "230920", "231006"] :
 for Sample_ID in [ "241016", "241211", "250425", "250212", "250502", "250509" ]:
    for TISSUE in ["Tumor"]:
        print ("{}_{}".format (Sample_ID, TISSUE))
        vcf_reader = vcf.Reader(open(DIR + "/02.PASS/{}_{}.MT2.FMC.HF.RMBLACK.vcf".format (Sample_ID, TISSUE), "r"))

        MBQ_list = []
        VAF_list = []
        for record in vcf_reader:        # record.CHROM, recrod.POS ,record.ALTa
            MBQ_list.append ( int ( record.INFO["MBQ"][1] ) )
            VAF_list.append ( float ( record.samples[1].data.AF[0] ) )
#            print ("{}\t{}\t{}\t{}".format ( record.CHROM, record.POS, record.ALT, record.INFO["MBQ"] ))
        #print ("\t{}\t{}개\tBQ>25:{}개\tVAF > 0.2:{}개\tBQ>25 & VAF > 0.2:{}개".format ( np.round ( np.mean (MBQ_list), 1), len (MBQ_list), len ( [i  for i in MBQ_list if i > 25  ]  ),  len ( [ i for i in range(len(VAF_list))  if  (VAF_list[i] > 0.2)  ]),  len ( [ i for i in range(len(VAF_list))  if  (VAF_list[i] > 0.2) & ( MBQ_list[i] > 25) ])   ) )
        print ("\t{}\t{}개\tBQ>25:{}개\tVAF > 0.2:{}개".format ( np.round ( np.mean (MBQ_list), 1), len (MBQ_list), len ( [i  for i in MBQ_list if i > 25  ]  ),  len ( [ i for i in range(len(VAF_list))  if  (VAF_list[i] > 0.2)  ]),  len ( [ i for i in range(len(VAF_list))  if  (VAF_list[i] > 0.2) & ( MBQ_list[i] > 25) ])   ) )


print ( "\n\n")


for Sample_ID in [ "241016", "241211", "250425", "250212", "250502", "250509" ]:
    for TISSUE in ["Tumor"]:
        print ("{}_{}".format (Sample_ID, TISSUE))
        vcf_reader = vcf.Reader(open(DIR + "/02.PASS/{}_{}.MT2.FMC.HF.RMBLACK.vcf".format (Sample_ID, TISSUE), "r"))

        MBQ_list = []
        for record in vcf_reader:        # record.CHROM, recrod.POS ,record.ALTa
            MBQ_list.append ( int ( record.INFO["MBQ"][1] ) )
            #print ("\t{}\t{}\t{}\t{}\t{}".format ( record.CHROM, record.POS, record.ALT, record.INFO["MBQ"], record.samples[1].data.AF ))
