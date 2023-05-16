#!/bin/bash
#$ -S /bin/bash
#$ -cwd


if ! options=$(getopt -o h --long ID:,SEQUENZA_OUTPUT_DIR:,SEQUENZA_LIFTOVER_OUTPUT_DIR:, -- "$@")
then
    echo "ERROR: invalid options"
    exit 1
fi

eval set -- $options

while true; do
    case "$1" in
        -h|--help)
            echo "Usage"
        shift ;;
        --ID)
            ID=$2
        shift 2 ;;
        --SEQUENZA_OUTPUT_DIR)
            SEQUENZA_OUTPUT_DIR=$2
        shift 2 ;;
        --SEQUENZA_LIFTOVER_OUTPUT_DIR)
            SEQUENZA_LIFTOVER_OUTPUT_DIR=$2
        shift 2 ;;
        --)
            shift
            break
    esac
done



#1.segments.txt
CrossMap.py region -r 0.7 "/home/goldpm1/resources/hg19ToHg38.over.chain.gz" ${SEQUENZA_OUTPUT_DIR}"/"${ID}"_segments.txt" ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_segments.body.txt" 
# header를 뽑아주기 (맨 뒤에 map_ratio를 하나 추가)
head -1 ${SEQUENZA_OUTPUT_DIR}"/"${ID}"_segments.txt" | sed 's/$/\tmap_ratio/'  > ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_segments.txt"
# header + body
cat ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_segments.body.txt"  >> ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_segments.txt" 

rm -rf ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_segments.body.txt.unmap" ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_segments.body.txt" 



# #2. mutation.txt
awk 'BEGIN{FS=OFS="\t"} {$3=$2+1; print}' ${SEQUENZA_OUTPUT_DIR}"/"${ID}"_mutations.txt" > ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.temp1.txt" 
CrossMap.py region -r 0.7 "/home/goldpm1/resources/hg19ToHg38.over.chain.gz" ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.temp1.txt"  ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.body.txt" 

# header를 뽑고 header와 body를 합쳐주기
head -1 ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.temp1.txt"  | sed 's/$/\tmap_ratio/' > ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.temp2.txt" 
cat ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.body.txt"  >> ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.temp2.txt" 

# # 3번째 column 을 GC.percent라고 바꿔주기. 아무 의미 없다.
awk 'BEGIN{FS=OFS="\t"} { $3="GC.percent"; print }' ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.temp2.txt"  > ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.txt" 

rm -rf ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.body.txt" ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.body.txt.unmap" ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.temp1.txt" ${SEQUENZA_LIFTOVER_OUTPUT_DIR}"/"${ID}"_mutations.temp2.txt"  