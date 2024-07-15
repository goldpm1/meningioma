#!/bin/bash
#$ -cwd
#$ -S /bin/bash

logPath="/data/project/Meningioma/script/51.Amplicon/log"
for sublog in 51.pysam; do
    if [ $logPath"/"$sublog ] ; then
        rm -rf $logPath"/"$sublog
    fi
    if [ ! -d $logPath"/"$sublog ] ; then
        mkdir -p $logPath"/"$sublog
    fi
done


# ', type=str, default="Dura_KLF4")
# parser.add_argument('--BAMTYPE', type=str, default="Amplicon_single")  #Amplicon, Amplicon_single

############################################ 1. SINGLE ##########################################################

BAMTYPE="Amplicon_single"

# TISSUE="221202_AKT1"
# qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE}"_"${TISSUE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}
# TISSUE="221202_TRAF7"
# qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE}"_"${TISSUE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}
# TISSUE="230405_KLF4"
# qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE}"_"${TISSUE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}
# TISSUE="230822_NF2"
# qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE}"_"${TISSUE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}

TISSUE="190426_NF2"
#qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE}"_"${TISSUE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}
TISSUE="220930_NF2"
#qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE}"_"${TISSUE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}
TISSUE="221026_NF2"
#qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE}"_"${TISSUE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}
TISSUE="230323_NF2"
#qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE}"_"${TISSUE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}
TISSUE="230920_NF2"
#qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE}"_"${TISSUE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}
TISSUE="230405_TRAF7"
#qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE}"_"${TISSUE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE ${TISSUE} --BAMTYPE ${BAMTYPE}



############################################ 2. MULTIPLEX ##########################################################

BAMTYPE="Amplicon_multiplex"
qsub -pe smp 1 -e $logPath"/51.pysam" -o $logPath"/51.pysam" -N ${BAMTYPE} "03.Amplicon_matrix_by_pysam.sh"  --TISSUE "multiplex" --BAMTYPE ${BAMTYPE}


    