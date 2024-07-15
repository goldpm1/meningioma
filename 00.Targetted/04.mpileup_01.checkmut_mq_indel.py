import os
import sys
from optparse import OptionParser
def trimseq(seq):
    newseq = ""
    insflag = False
    delflag = False
    skipsize = ""
    skipping = 0
    insCnt = 0
    insSizeStr = ""
    insSize = 0
    insSeq = ""
    for s in seq:
        if skipping > 0:
            skipping-=1
            continue
        if insSize > 0:
            insSize-=1
            insSeq+=s
            if insSize==0:
                insSeq+=","
                insCnt+=1
            continue
        if delflag:
            if s >= "0" and s <= "9":
                skipsize += s
            else:
                skipping = int(skipsize) - 1
                skipsize = ""
                delflag = False
            continue
        if insflag:
            if s >= "0" and s <= "9":
                insSizeStr += s
            else:
                insSize = int(insSizeStr) - 1
                insSizeStr = ""
                insSeq += s
                if insSize==0:
                    insSeq+=","
                    insCnt+=1
                insflag = False
            continue
        if s=="$":
            continue
        if s=="^":
            skipping = 1
            continue
        if s=="+":
            insflag = True
            continue
        if s=="-":
            delflag = True
            continue
        newseq += s
    insSeq = insSeq[:-1]
    return newseq, insCnt, insSeq
def phredToInt(phred):
    return ord(phred) - 33
def intToPhred(real):
    return str(real + 33)
def getMaxIndex(mylist):
    maxvalue = -1
    maxindex = 0
    for i in range(0, len(mylist)):
        l = mylist[i]
        if l > maxvalue:
            maxvalue = l
            maxindex = i
    return maxindex
def getcount1(ref, seq, qual, mqual, bq, mq):
    cntref = 0
    cntalt = 0
    cntATGC = [0, 0, 0, 0]
    cntlowqual = 0
    cntdel = 0
    seqlen = len(seq)
    quallen = len(qual)
    #mquallen = len(mqual)
    #if seqlen != quallen or seqlen != mquallen or quallen != mquallen:
        #print "ERROR: length does not match"
        #print seqlen, quallen, mquallen
        #return 0, 0, [], 0, 0
    if seqlen != quallen:
        print "ERROR: length does not match"
        print seqlen, quallen
        return 0, 0, [], 0, 0
    for i in range(0, seqlen):
        qi = qual[i]
        #mqi = mqual[i]
        si = seq[i]
        qualint = phredToInt(qi)
        #mqualint = phredToInt(mqi)
#       print si, qi, qualint, bq
        if si == "*":
            cntdel += 1
            continue
        #if qualint < bq or mqualint < mq:
        if qualint < bq:
            cntlowqual += 1
            continue
        if si == "." or si == ",":
            cntref += 1
        else:
            cntalt += 1
            ind = "ATGC".index(si.upper())
            cntATGC[ind] += 1
    return cntref, cntalt, cntATGC, cntlowqual, cntdel
def runpileup(pileup, freq, Freq, bq, mq, mcnt, ic, dc, outfile):
    f = open(pileup)
    w = None
    if outfile != "":
        w = open(outfile, 'w')
    i = f.readline().strip()
    if w:
        w.write("Pileup File\tChrom\tPos\tDepth\tHiQualDepth\tRef\tRefCnt\tToTAltCnt\tTotAltFreq\tMaxAllele\tMaxAlleleCnt\tMaxAlleleFreq\tCnt_A\tFreq_A\tCnt_T\tFreq_T\tCnt_G\tFreq_G\tCnt_C\tFreq_C\tLowQualFreq\tCnt_del\tCnt_ins\tIns_seq\n")
    else:
        print "Pileup File\tChrom\tPos\tDepth\tHiQualDepth\tRef\tRefCnt\tTotAltCnt\tTotAltFreq\tMaxAllele\tMaxAlleleCnt\tMaxAlleleFreq\tCnt_A\tFreq_A\tCnt_T\tFreq_T\tCnt_G\tFreq_G\tCnt_C\tFreq_C\tLowQualFreq\tCnt_del\tCnt_ins\tIns_seq"
    while i:
        values = i.split("\t")
        if len(values) < 4:
            i = f.readline().strip()
            continue
        chrom = values[0]
        pos = values[1]
        ref = values[2].upper()
        depth = values[3]
        seq, insCnt, insSeq = trimseq(values[4])
        qual = values[5]
        mapqual = ""
        if len(values) > 6:
            mapqual = values[6]
        #if len(seq) != int(depth):
        #   print "ERROR", len(seq), depth
        #   print seq, depth
        cntref, cntalt, cntATGC, cntlowqual, cntdel = getcount1(ref, seq, qual, mapqual, bq, mq)
        if cntref + cntalt == 0:
            i = f.readline().strip()
            continue
        alta = cntATGC[0]
        altt = cntATGC[1]
        altg = cntATGC[2]
        altc = cntATGC[3]
        altfreqa = round(float(alta) / float(cntref + cntalt) * 100, 2)
        altfreqt = round(float(altt) / float(cntref + cntalt) * 100, 2)
        altfreqg = round(float(altg) / float(cntref + cntalt) * 100, 2)
        altfreqc = round(float(altc) / float(cntref + cntalt) * 100, 2)
        maxalt = "ATGC"[getMaxIndex(cntATGC)]
        maxaltcnt = cntATGC[getMaxIndex(cntATGC)]
        maxaltfreq = round(float(maxaltcnt) / float(cntref + cntalt) * 100, 2)
        altfreq = round(float(cntalt) / float(cntref + cntalt) * 100, 2)
        lowfreq = round(float(cntlowqual) / float(int(depth)) * 100, 2)
        hiqual = int(depth) - cntlowqual
        values = [pileup, chrom, pos, depth, hiqual, ref, cntref, cntalt, altfreq, maxalt, maxaltcnt, maxaltfreq, alta, altfreqa, altt, altfreqt, altg, altfreqg, altc, altfreqc, lowfreq, cntdel, insCnt, insSeq]
        svalues = []
        for v in values:
            svalues.append(str(v))
        if (altfreq >= freq*100 and altfreq <= Freq*100 and cntalt >= mcnt) or (insCnt >= ic or cntdel >=dc) :
            if w:
                w.write("\t".join(svalues) + "\n")
            else:
                print "\t".join(svalues)
        i = f.readline().strip()
    if w:
        w.close()
def main(argv):
    usage = "usage: %prog [options] pileup"
    parser = OptionParser(usage)
    parser.add_option("-f", "--min_frequency", dest="freq", help="minimum allele frequeny threshold [0.01]", default = "0.01")
    parser.add_option("-F", "--max_frequency", dest="Freq", help="maximum allele frequency threshold [0.5]", default = "0.5")
    parser.add_option("-c", "--min_count", dest="mcnt", help="minimum number of mismatches [2]", default="2")
    parser.add_option("-q", "--base_quality_filter", dest="bq", help="minimum base call quality [30]", default = "30")
    parser.add_option("-Q", "--mapping_quality_filter", dest="mq", help="minimum mapping quality [30]", default= "30")
    parser.add_option("-i", "--insertion_count", dest="ic", help="minimum number of insertion [2]", default = "2")
    parser.add_option("-d", "--deletion_count", dest="dc", help="minimum number of deletion [2]", default= "2")
    parser.add_option("-o", "--output_file_name", dest="out", help="Output file name [stdout]", default="")
    options, args = parser.parse_args()
    if len(args) < 1:
        parser.error("ERROR: must provide a pileup file")
    pileup = args[0]
    if not os.path.exists(pileup):
        parser.error("ERROR: " + pileup + " does not exist.")
    freq, Freq, bq, mq, ic, dc = float(options.freq), float(options.Freq), int(options.bq), int(options.mq), int(options.ic), int(options.dc)
    outfile = str(options.out)
    mcnt = int(options.mcnt)
#   print "Working on " + pileup + "..."
    runpileup(pileup, freq, Freq, bq, mq, mcnt, ic, dc, outfile)
if __name__=="__main__":
    main(sys.argv[1:])


    