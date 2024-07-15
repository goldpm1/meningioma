import sys
import os
import glob
import datetime
from subprocess import CalledProcessError, check_output
import subprocess
from scipy import stats
import scipy as sp
from numpy import *
from collections import defaultdict
class DelInfo:
    def __init__(self):
        self.chr = ''
        self.start = ''
        self.ref = ''
        self.alt = ''
        self.totCnt = []
        self.altCnt = []
        self.freqList = []
    def setPos(self, chr, start):
        self.chr = chr
        self.start = int(start)
    def setRef(self, seq):
        self.ref = seq
        self.alt = self.ref
    def addSeq(self, seq):
        self.ref += seq
    def addCnt(self, totCnt, altCnt):
        self.totCnt.append(totCnt)
        self.altCnt.append(altCnt)
        freq = int(altCnt)/float(totCnt)
        self.freqList.append(freq)
    def getRef(self):
        return self.ref
    def getAlt(self):
        return self.alt
    def getPos(self):
        return self.chr, self.start
    def getTotCnt(self):
        statsTot = array(self.totCnt)
        return int(sp.median(statsTot))
    def getAltCnt(self):
        statsAlt = array(self.altCnt)
        return int(sp.median(statsAlt))
    def getMedianFreq(self):
        statsFreq = array(self.freqList)
        return sp.median(statsFreq)
#DIR = '/data/project/ctDNA/151112_panelData/Late/4.analysis/mpileup'
DIR = sys.argv[1]
filename = sys.argv[2]
sample = sys.argv[3]
REF=sys.argv[4]
idx=sys.argv[5]
#files = glob.glob(DIR + '/*.mpileup.call')
#print len(files)
#REF = '/data/resource/reference/human/UCSC/hg19/WholeGenomeFasta/genome.fa'
#idx = '/data/resource/reference/human/UCSC/hg19/WholeGenomeFasta/genome.fa.idx'
###############GRCH38
# REF = '/data/resource/reference/human/NCBI/GRCh38_GATK/BWAIndex/genome.fa'
# idx = '/data/resource/reference/human/NCBI/GRCh38_GATK/BWAIndex/genome.fa.idx'
##b37
# REF='/data/public/LongRanger/b37_2.1.0/fasta/genome.fa'
# idx='/home/hswbulls/download/Longranger_b37/genome.fa.fai'
file = open(REF, 'r')
REFlines = file.readlines()
file = open(idx, 'r')
lines = file.readlines()
######################b37 lines[1:0] -> lines[0:0]
# REFidx = {}
# for line in lines[0:]:
#   row = line.split('\t')
#   REFidx[row[0]] = [int(row[1]), int(row[2])]
###############GRCH38
REFidx = {}
for line in lines[1:]:
    row = line.split('\t')
    REFidx[row[0]] = [int(row[1]), int(row[2])]
# print REFidx
INPUT = DIR + '/' + filename
#for filename in files:
file = open(INPUT, 'r')
#tmp = filename.split('/')
#sampleID = tmp[len(tmp)-1].split('_')[0] + '_' + tmp[len(tmp)-1].split('_')[1]
#sampleID = tmp[len(tmp)-1].split('_')[0]
#DIR = ''
#if tmp[0] == '.':
#   DIR = tmp[0]
#else:
#   DIR = '/' + tmp[0]
#for item in tmp[1:-1]:
#   DIR += '/' + item
out_filename = DIR + '/' + filename[:-5] + '.vcf'
outfile = open(out_filename, 'w')
outfile.write('##fileformat=VCFv4.1\n')
outfile.write('##fileDate=' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + '\n')
outfile.write('##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotyping">\n')
outfile.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read count">\n')
outfile.write('##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Alt read count">\n')
outfile.write('##FORMAT=<ID=INS,Number=1,Type=String,Description="Insert sequence list">\n')
outfile.write('##FILTER=<Point mutation=altCnt>5&altFreq>2%, InDel=indelCnt>0&medianFreq>\n')
outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
preChr = ''
prePos = ''
delInfo = DelInfo()
REFSeq = ''
for idx, line in enumerate(file):
    if idx > 0:
        token = line.rstrip().split('\t')
        chr = token[1]
        if chr != preChr:
            try:
                REFSeq = ''
                for idx in range(REFidx[chr][0], REFidx[chr][1]+1):
                    REFSeq += REFlines[idx].rstrip()
                sys.stderr.write(str(len(REFSeq)) + '\n')
            except IndexError:
                print chr, idx, len(REFlines)
                sys.exit()
        pos = token[2]
        tDepth = int(token[3])
        hDepth = int(token[4])
        ref = token[5]
        delCnt = int(token[21])
        insCnt = int(token[22])
        maxAltCnt = int(token[10])
        maxAlt = token[9]
        if delCnt >= 5:
            if delInfo.getAlt() == '':
                delInfo.setPos(chr, int(pos)-1)
                delInfo.setRef(REFSeq[int(pos)-2:int(pos)-1].upper())
                delInfo.addSeq(ref)
                delInfo.addCnt(tDepth, delCnt)
                ##print chr, pos, delInfo.getRef(), delInfo.getAlt(), REFSeq[int(pos)-2:int(pos)+2]
            elif preChr == chr and int(prePos)+1 == int(pos):
                delInfo.addSeq(ref)
                delInfo.addCnt(tDepth, delCnt)
            else:
                outfile.writelines(delInfo.getPos()[0] + '\t' + str(delInfo.getPos()[1]) + '\t.\t' + delInfo.getRef() + '\t' + delInfo.getAlt())
                outfile.writelines('\t' + str(delInfo.getMedianFreq()) + '\t.')
                outfile.writelines('\tINS=.')
                outfile.writelines('\tGT:AD:DP')
                outfile.writelines('\t./.:' + str(delInfo.getTotCnt()-delInfo.getAltCnt()) + ',' + str(delInfo.getAltCnt()) + ':' + str(delInfo.getTotCnt()) + '\n')
                delInfo = DelInfo()
                delInfo.setPos(chr, int(pos)-1)
                delInfo.setRef(REFSeq[int(pos)-2:int(pos)-1].upper())
                delInfo.addSeq(ref)
                delInfo.addCnt(tDepth, delCnt)
                print chr, pos, REFSeq[int(pos)-2:int(pos)+2]
        elif delInfo.getAlt() != '':
            outfile.writelines(delInfo.getPos()[0] + '\t' + str(delInfo.getPos()[1]) + '\t.\t' + delInfo.getRef() + '\t' + delInfo.getAlt())
            outfile.writelines('\t' + str(delInfo.getMedianFreq()) + '\t.')
            outfile.writelines('\tINS=.')
            outfile.writelines('\tGT:AD:DP')
            outfile.writelines('\t./.:' + str(delInfo.getTotCnt()-delInfo.getAltCnt()) + ',' + str(delInfo.getAltCnt()) + ':' + str(delInfo.getTotCnt()) + '\n')
            delInfo = DelInfo()
        if insCnt >= 5:
            insSeq = token[23].upper().split(',')
            occurrence = defaultdict(int)
            for item in insSeq:
                occurrence[item] += 1
            maxIns = max(occurrence.iteritems(), key=lambda x:x[1])
            outfile.write(chr + '\t' + pos + '\t.\t' + ref)
            outfile.write('\t' + ref + maxIns[0] + '\t' + str(maxIns[1]/float(tDepth)) + '\t.')
            outfile.write('\tINS=' + ','.join(insSeq))
            outfile.write('\tGT:AD:DP')
            outfile.write('\t./.:' + str(tDepth-maxIns[1]) + ',' + str(maxIns[1]) + ':' + str(tDepth) + '\n')
        if maxAltCnt >= 1:
            outfile.write(chr + '\t' + pos + '\t.\t' + ref + '\t' + maxAlt + '\t' + str(maxAltCnt/float(hDepth)) + '\t.\tINS=.')
            outfile.writelines('\tGT:AD:DP')
            outfile.write('\t./.:' + str(hDepth-maxAltCnt) + ',' + str(maxAltCnt) + ':' + str(hDepth) + '\n')
        preChr = chr
        prePos = pos
outfile.close()
file.close()