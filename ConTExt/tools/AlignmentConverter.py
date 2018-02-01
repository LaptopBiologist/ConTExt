#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Michael
#
# Created:     07/10/2015
# Copyright:   (c) Michael Peter McGurk 2015
# Licence:
#The MIT License (MIT)

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#-------------------------------------------------------------------------------
import csv
import os
import numpy
import numpy as np
import shutil
import subprocess
import sys

import time
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
from Bio import Seq

from itertools import chain

class SamLine():
    def __init__(self, row):
        self.id=row[0]
        self.flag=int(row[1])
        self.contig=row[2]
        self.pos=int(row[3])
        self.MapQ=int(row[4])
        self.cigar=row[5]
        self.rnext=row[6]
        self.pnext=int(row[7])
        self.tlen=int(row[8])
        self.seq=row[9]
        self.qual=row[10]
        self.optional=row[11:]
    def row(self):
        return [self.id, self.flag, self.contig, self.pos, int(self.MapQ), self.cigar,\
    self.rnext, self.pnext, self.tlen, self.seq, self.qual]+self.optional

class BowtieAlignment():
    def __init__(self, row):
        self.queryName=row[0]
        self.subName=row[1]
        self.pIdent=float(row[2])/100
        self.length=float(row[3])
        self.MatchingPositions=self.pIdent*self.length
        self.qStart=int( row[6])
        self.qEnd=int( row[7])
        self.sStart=int( row[8])
        self.sEnd=int( row[9])
        self.eValue=float(row[10])
        self.BTOP=row[11]


class ConTExtLine():
    strand={'0':'+', '16':'-', '4':'*'}
    def __init__(self, row, phred=33):

        self.Contig1=row[0]
        self.Strand1=ConTExtLine.strand[row[1]]
        self.Start1=int(row[2])
        self.End1=int(row[3])
        self.Seq1=row[4]
        self.Qual1=row[5]
        self.Contig2=row[6]
        self.Strand2=ConTExtLine.strand[row[7]]
        self.Start2=int(row[8])
        self.End2=int(row[9])
        self.Seq2=row[10]
        self.Qual2=row[11]
        self.Cigar1=row[12]
        self.Cigar2=row[13]
        self.MapQ1=int(row[14])
        self.MapQ2=int(row[15])
        #self.Mismatches1=row[16]
        #self.Mismatches2=row[17]


class Traceback():
    def __init__(self, alignment):
        self.BTOP=alignment.BTOP
        self.splitTrace=InterpretBTOP(alignment.BTOP)
        self.CigarParts=TraceToCigar(self.splitTrace)
        self.SubjectName=alignment.subName
        self.QueryName=alignment.queryName
        self.identity=alignment.pIdent
        self.subStart=alignment.sStart
        self.subEnd=alignment.sEnd
        self.queryStart=alignment.qStart
        self.queryEnd=alignment.qEnd
        parts=alignment.queryName.split('_')
        length=int(parts[-2])-int(parts[-3])
        self.RefStart, self.RefEnd, self.RefStrand=parts[-3],parts[-2],parts[-1]
        self.YIntercepts=DetermineYIntercepts(self.queryStart, self.subStart, self.CigarParts)
        leftInsert=[numpy.nan]*(alignment.qStart-1)
        rightInsert=[numpy.nan]*(length-alignment.qEnd)
        self.YIntercepts=numpy.hstack([leftInsert, self.YIntercepts, rightInsert])
        del self.splitTrace
        del self.CigarParts

class ConversionInfo():
    def __init__(self, row):
        header=['InsertionName', 'ConsensusName','Consensus?', 'RefStart',  'RefEnd', 'RefStrand',\
    'InsertionStart', 'InsertionEnd', 'ConsensusStart','ConsensusEnd', 'Identity', 'TracebackOperations']

        self.insName, self.consName=row[0], row[1]
        self.consensus=bool(row[2])
        self.refStart, self.refEnd,self.refStrand=int(row[3]), int(row[4]), row[5]
        self.insStart, self.insEnd= int(row[6]), int(row[7])
        self.consStart, self.consEnd= int(row[8]), int(row[9])
        self.pIdent= float(row[10])
        self.BTOP=row[11]

        self.splitTrace=InterpretBTOP(self.BTOP)
        self.CigarParts=TraceToCigar(self.splitTrace)

        length=self.refEnd-self.refStart
##        length=self.insEnd-self.insStart


        if self.consStart<self.consEnd:
            self.YIntercepts=DetermineYIntercepts(self.insStart, self.consStart, self.CigarParts)
            leftInsert=[numpy.nan]*(self.insStart-1)
            rightInsert=[numpy.nan]*(length-self.insEnd)
            self.slope=1
        elif self.consStart>=self.consEnd:
            self.YIntercepts= DetermineYIntercepts(self.insStart, self.consStart, self.CigarParts, slope=-1)
            leftInsert=[numpy.nan]*(self.insStart-1)
            rightInsert=[numpy.nan]*(length-self.insEnd)
            self.slope=-1


        self.YIntercepts=numpy.hstack([leftInsert, self.YIntercepts, rightInsert])
        del self.splitTrace
        if len( self.YIntercepts)!= length:
##            print "Warning!! {0} != {1}".format(len(self.YIntercepts), length)
            print self.consStart>=self.consEnd
        del self.CigarParts


def InterpretBTOP(btop):
    hold=''
    partList=[]
    index=0
    flag=''
    for i in range(len(btop)):
        character=btop[index]

        #If it's a number
        if 48<=ord(character)<=57:
            if hold!='' and flag!='M':
                partList.append((flag, hold))
                hold=''
            hold+=character
            flag='M'
        else:
            if hold!='' and flag=='M':
                partList.append((flag, int(hold)))
                hold=''
            elif hold!='' and flag=='N':
                partList.append((flag, hold))
                hold=''
            hold=btop[index:index+2]
            flag='N'
            index+=1
        index+=1
        if index>=len(btop): break
    if flag=='M':
        partList.append((flag, int(hold)))
    else:
        partList.append((flag, hold))

    return partList

def TraceToCigar(trace):
    cigar_rough=[]
    for part in trace:
        if part[0]=='M':
            cigar_rough.append(part)
        if part[0]=='N':
            if part[1][0]=='-':
                cigar_rough.append(('D', 1))
            elif part[1][1]=='-':
                cigar_rough.append(('I', 1))
            else: cigar_rough.append(('M',1))
    hold=''
    hold_count=0
    cigar=[]
    for oper in cigar_rough:
        oper_type=oper[0]
        oper_count=oper[1]
        if oper_type==hold:
            hold_count+=oper_count
        else:
            if hold_count!=0:
                cigar.append((hold, hold_count))
            hold=oper_type
            hold_count=oper_count
    cigar.append((hold, hold_count))
    return(cigar)

def DetermineYIntercepts(subStart, queryStart, cigar, slope=1):
    #This is the y-intercept before modification by indels: pos(Seq2)=slope* pos(Seq1) +b
    # b=pos( Seq2) -slope* pos(Seq1)
    # So at the first position of Seq2 (query ) and Seq1 (subject) the intercept is:
    # b= QueryStart- slope* SubjectStart

    currentOffset= queryStart- slope* subStart
    offsetList=[]
    for part in cigar:
        if part[0]=='M':
            offsetList+=[currentOffset]*part[1]
        if part[0]=='I':
            offsetList+=[numpy.nan]*part[1]
            currentOffset-= slope*part[1]
        if part[0]=='D':
            currentOffset+= slope*part[1]
    return numpy.array(offsetList)

##############################################################################
##############################################################################
##############################################################################

"""Functions necessary to prepare the conversion table between insertions and
consensus sequences."""

def PrepareConversionTable(consFile, refFile, gffFile, outfile, BlastDir):
    outDir='/'.join( outfile.split('/')[:-1])
    outRoot='.'.join(outfile.split('.')[:-1])
    #Create a file for each repeat family and add the sequence of each instance
    #of that family in the reference to the file
    tempInsertionDir=outDir+'/insertions_temp'

    print "Extracting insertions from reference...."
    ExtractSeqFromGFF(refFile,gffFile,tempInsertionDir)

    #Split the consensus Fasta into a separate file for each entry. This step is entirely
    #unnecessary, and all of the files it creates are subsequently deleted.
    #"Efficiency" was invented by The Man to emotionally oppress the individual:
    #ANARCHY!!!! ANARCHY!!!! ANARCHY!!!! ANARCHY!!!! ANARCHY!!!! ANARCHY!!!!
    tempConsensusDir=outDir+'/consensus_temp'
    SplitFASTA(consFile, tempConsensusDir)

    #Create the index contain both consensus sequences and all insertions:
    JoinFasta(tempConsensusDir, outRoot+'.fa')
    JoinFasta(tempInsertionDir, outRoot+'.fa')

    #Align each instance of a repeat family to its consensus
    SplitFilesInDir(tempInsertionDir) #If there are many insertions, split them into separate files to limit memory use when aligning
    tempAlignedDir=outDir+'/aligned_temp'
    BlastDirectory(tempInsertionDir, consFile, tempAlignedDir, BlastDir)

##    CheckAlignedNames(tempAlignedDir)

    #Choose the best alignment for each insertion
    #Best is define as follows:
    #The alignment with the lowest e-value
    #in case of tie:
    #The alignment with the greatest number of matching positions:
    #       =Aligned Length x Percent Identity
    print "Getting best aligments..."
    bestAlignments=GetAllAlignments(tempAlignedDir)

    WriteConversionTable(bestAlignments, consFile, outRoot+'.tsv')

    #Delete temporary files
    shutil.rmtree (tempAlignedDir)
    shutil.rmtree(tempConsensusDir)
    shutil.rmtree(tempInsertionDir)

def ExtractSeqFromGFF(refFile, gffFile, outDir):
    MakeDir(outDir)
    refSeq=GetSeq(refFile)
    gffHandle=open(gffFile, 'r')
    bedTable= csv.reader(gffHandle, delimiter='\t')
    fileDict={}
    fileList=[]
    print refSeq.keys()
    for row in bedTable:
        if row[0][0]=='#': continue
        chrom=row[0]
        start=row[3]
        end=row[4]
        strand=row[6]
        name=( row[-1].split('=')[-1].split(' '))
        if name[0][0]=='(': continue
        seqName='_'.join([name[0] , chrom,start, end, strand])

        pulledSeq=refSeq[chrom][int(start):int(end)].seq
        if strand=='-': pulledSeq=pulledSeq.reverse_complement()
        pulledSeq=str(pulledSeq)
        pulledName='_'.join(row)
        if fileDict.has_key(name[0])==False:
            fileDict[name[0]]=open(outDir+'/'+CleanName( name[0])+'.fa', 'a')
            fileList.append(name[0])
            if len(fileList)>100:
                fileDict[fileList[0]].close()
                del fileDict[fileList[0]]
                del fileList[0]

        fileDict[name[0]].write('>'+seqName+'\n')
        fileDict[name[0]].write(pulledSeq+'\n')

    for i in fileList:
        fileDict[i].close()


    gffHandle.close()

def MakeDir(path):
    if os.path.exists(path)==False:
        os.mkdir(path)

def CleanName(name):
    """Removes illegal characters from sequence names."""
    illegal=['|', '!','>', '<', '?','/','*']
    for i in illegal:
        name='_'.join(name.split(i))
    return(name)

def SplitFASTA(infile, outDir):
    refSeq=GetSeq(infile)
    MakeDir(outDir)
    for key in refSeq:
        seq=str( refSeq[key].seq)
        name= CleanName( str(refSeq[key].name))
        outHandle=open(outDir+'/'+name+'.fa', 'w')
        outHandle.write('>'+name+'\n')
        outHandle.write(seq+'\n')
        outHandle.close()

def GetSeq(ref, parser='fasta', gzipped=False):
    if gzipped==False:
        handle=open(ref, 'r')
    else:
        handle=gzip.open(ref)
    lib=SeqIO.parse(handle, parser)
    SeqLen={}
    for rec in lib:
        #if refName.count(rec.name)==0: continue
        SeqLen[CleanName(rec.name)]=rec
    handle.close()
    return SeqLen

def JoinFasta(inDir, outfile, mode='a'):
    outHandle=open(outfile, mode)
    inFiles=os.listdir(inDir)
    for f in inFiles:
        inhandle=open(inDir+'/'+f, 'r')
        for line in inhandle:
            outHandle.write(line)
        inhandle.close()
    outHandle.close()

def SplitFilesInDir(inDir):
    fileList=os.listdir(inDir)
    for f in fileList:
        splitFile(inDir+'/'+f)

def splitFile(infile, limit=1000):
    refSeq=GetSeq(infile)
    inDir='/'.join(infile.split('/')[:-1])
    root=infile.split('/')[-1]
    rootname='.'.join(root.split('.')[:-1])
    ext=root.split('.')[-1]
    fileCount=1
    outfile="{0}/{1}_temp_{2}.{3}".format(inDir, rootname, fileCount, ext)
    outHandle=open(outfile, 'w')
    count=0
    for key in refSeq:
        if count>=limit:
            outHandle.close()
            fileCount+=1
            count=0
            outfile="{0}/{1}_temp_{2}.{3}".format(inDir, rootname, fileCount, ext)
            outHandle=open(outfile, 'w')
        seq=str( refSeq[key].seq)
        name= CleanName( str(refSeq[key].name))

        outHandle.write('>'+name+'\n')
        outHandle.write(seq+'\n')
        count+=1

    outHandle.close()
    os.remove(infile)


def BlastDirectory(insDir, consPath, outDir, BlastDir='c:/ncbi-blast'):

    insFile=os.listdir(insDir)
##    consFile=os.listdir(consDir)
    MakeDir(outDir)
    for f in insFile:
        print "BLASTing {0}...".format(f)
##        consPath=consDir+'/'+f
        insPath=insDir+'/'+f
##        if os.path.exists(consPath)==False: continue
        BlastSequences(insPath, consPath, outDir,  '.'.join(f.split('.')[:-1]) , BlastDir)


def BlastSequences(Query, Subject, Out, name, BlastDir):
    """BLASTs two fastas against each other."""
    OutPath='.'.join(Out.split('.'))
    OutFile=OutPath+'/'+name+'.csv'
    print (OutPath)
    errlog=open(OutPath+'/_err.log', 'a')
    column_spec='10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue btop'
    BLAST=subprocess.Popen([BlastDir+'/bin/blastn', '-query',Query, '-subject',Subject, '-outfmt', column_spec,  '-out', OutFile], stderr=errlog)
    BLAST.communicate()
    errlog.close()
    return OutFile

def GetAllAlignments(indir):

    fileList=os.listdir(indir)
    traceDict={}
    for f in fileList:
        if f=='_err.log':continue
        print 'Getting alignments from {0}'.format(f)
        newKeys=ChooseBestAlignments(indir+'/'+f)
        for key in newKeys.keys():
            traceDict[key]=newKeys[key]
    return traceDict

def ChooseBestAlignments(infile):
    alignments={}
    inhandle=open(infile, 'r')
    inTable=csv.reader(inhandle, delimiter=',')
    for row in inTable:
        try:
            line=BowtieAlignment(row)
        except:
            print row
            print infile
            print Jabber
        if line.pIdent<.8: continue
        if alignments.has_key(line.queryName)==True:
            if alignments[line.queryName].eValue>line.eValue:
                alignments[line.queryName]=line
            if alignments[line.queryName].eValue==line.eValue and\
            alignments[line.queryName].MatchingPositions<line.MatchingPositions:
                alignments[line.queryName]=line
        else:
            alignments[line.queryName]=line
    #Convert into useable traces
    traceDict={}
    for key in alignments.keys():
        traceDict[key]=Traceback(alignments[key])
    return traceDict

def WriteConversionTable(alignments, consensusFile, outfile):

    outHandle=open(outfile, 'wb')
    header=['InsertionName', 'ConsensusName','Consensus?', 'RefStart',  'RefEnd', 'RefStrand',\
    'InsertionStart', 'InsertionEnd', 'ConsensusStart','ConsensusEnd', 'Identity', 'TracebackOperations']
    outTable=csv.writer(outHandle, delimiter='\t')
    outTable.writerow(header)
    consSeq=GetSeq(consensusFile)
    for key in sorted(consSeq.keys()):
        line=[key, key, 'True']+['']*(len(header)-3)
        outTable.writerow(line)
    for key in sorted(alignments.keys()):
        alignment=alignments[key]
        line=[alignment.QueryName, alignment.SubjectName, 'False', alignment.RefStart, alignment.RefEnd, alignment.RefStrand\
        , alignment.queryStart, alignment.queryEnd, alignment.subStart, alignment.subEnd, alignment.identity, alignment.BTOP]
        outTable.writerow(line)
    outHandle.close()

################################################################################
################################################################################
################################################################################
"""Functions necessary to convert alignements in a sam file from insertions to
consensus sequences"""


def TroubleshootInsertionsToConsensus(inSam, conversionFile):
    #Pull the lines we would normally convert


    #Read the conversion table
    #   consDict is used to identify reads that already align to the consensus
    #   insDict is used to identify which reads are aligned to insertions
    #   and stores the information necessary to convert insertion postions to consensus positions
    consDict, insDict=ReadConversionTable(conversionFile)
    #Open the input and output files as Tab-delimited files

    tempRoot='.'.join( inSam.split('.')[:-1])
    tempFile=tempRoot+'_consensus.sam'
    tempHandle=open(tempFile, 'w')
    tempTable=csv.writer(tempHandle , delimiter='\t')

    inHandle=open(inSam, 'r')
    inTable=csv.reader(inHandle, delimiter='\t')
    convertCount=0.
    total=0.
    line_list=[]
    for row in inTable:

        #Handle the header rows
        if row[0][0]=='@':
            if row[0]=='@SQ':
                name=row[1].split(':')[-1]
                if consDict.has_key(name)==True:
                    tempTable.writerow(row)
                continue
            else:
                tempTable.writerow(row)
                continue

        #Handle the non-header rows
        line=SamLine(row) #Organize the information in this row using the SamLine class
        if consDict.has_key( line.contig)==True:    #Read aligns to consensus, don't change anything about the sam line
            tempTable.writerow(row)
        elif insDict.has_key(line.contig)==True:    #Read aligns to insertion, convert its position to consensus
            total+=1
##            newPos, newCigar= CombineAlignments(line.pos, insDict[line.contig].insStart\
##            ,insDict[line.contig].insEnd, line.cigar, insDict[line.contig].YIntercepts)
            line_list.append(line)


    inHandle.close()
    tempHandle.close()

    return line_list, insDict


def DebugConvertSAMFromInsertionsToConsensus(inSam, conversionFile):
    #Read the conversion table
    #   consDict is used to identify reads that already align to the consensus
    #   insDict is used to identify which reads are aligned to insertions
    #   and stores the information necessary to convert insertion postions to consensus positions
    consDict, insDict=ReadConversionTable(conversionFile)
    #Open the input and output files as Tab-delimited files

    tempRoot='.'.join( inSam.split('.')[:-1])
    tempFile=tempRoot+'_consensus.sam'
    tempHandle=open(tempFile, 'w')
    tempTable=csv.writer(tempHandle , delimiter='\t')

    inHandle=open(inSam, 'r')
    inTable=csv.reader(inHandle, delimiter='\t')
    convertCount=0.
    total=0.
    for row in inTable:

        #Handle the header rows
        if row[0][0]=='@':
            if row[0]=='@SQ':
                name=row[1].split(':')[-1]
                if consDict.has_key(name)==True:
                    tempTable.writerow(row)
                continue
            else:
                tempTable.writerow(row)
                continue

        #Handle the non-header rows
        line=SamLine(row) #Organize the information in this row using the SamLine class
        if consDict.has_key( line.contig)==True:    #Read aligns to consensus, don't change anything about the sam line
            tempTable.writerow(row)
        elif insDict.has_key(line.contig)==True:    #Read aligns to insertion, convert its position to consensus
            total+=1
##            newPos, newCigar= CombineAlignments(line.pos, insDict[line.contig].insStart\
##            ,insDict[line.contig].insEnd, line.cigar, insDict[line.contig].YIntercepts)
            newPos,newCigar=MapToConsensus(line.pos, insDict[line.contig].insStart, line.cigar, insDict[line.contig].YIntercepts, insDict[line.contig].slope, insDict[line.contig],  len( line.seq))
            if newCigar==False: return newPos
            #Calculate how long the alignment to the insertion was
            originalCigar=ParseCIGAR(line.cigar)
            originalLength=originalCigar['M']+originalCigar['D']

            #Calculate how long the new alignment to the consensus is
            newCigarParts=ParseCIGAR(newCigar)
            newLength=newCigarParts['M']+newCigarParts['D']
            if newCigar=='': continue
        # I am ignoring internal NaNs which represent insertions relative to the individual element!!!!
            if newCigarParts['M']+newCigarParts['I']!=len( line.seq):
                print newCigar
                print newCigarParts['M']+newCigarParts['I'] , len( line.seq)
                print jabber
            continue
            if insDict[line.contig].slope==-1: #The read is now on the opposte strand, so note that and reverse the CIGAR string
                line.flag=16-line.flag
                newCigar=ReverseCIGAR(newCigar)

            if numpy.isnan( newPos)==False and .5*originalLength<=newLength<=1.5*originalLength:
                #Provided the read contains positions represented in the consensus,
                #and its alignment to the consensus is no less that 50% of the original
                #aligned length and no greater than 150% assign it to the consensus
                #convert the sam line to the consensus alignment
                convertCount+=1
                line.pos=int(newPos)
                line.cigar=newCigar
                line.contig=insDict[line.contig].consName
            else:
                #Otherwise, convert it to an unmapped read
                line.pos=0
                line.cigar='*'
                line.contig='*'
                line.flag=4
##            tempTable.writerow(line.row())  #write the samline
##            print newPos, newCigar
        else:   #Read is unmapped
            line.pos=0
            line.cigar='*'
            line.contig='*'
            line.flag=4
            tempTable.writerow(line.row())
    inHandle.close()
    tempHandle.close()
    return convertCount, total, convertCount/total


def ConvertSAMFromInsertionsToConsensus(inSam, conversionFile, consensus_file=''):
    #Read the conversion table
    #   consDict is used to identify reads that already align to the consensus
    #   insDict is used to identify which reads are aligned to insertions
    #   and stores the information necessary to convert insertion postions to consensus positions
    consSequences=GetSeq(consensus_file)
    consDict, insDict=ReadConversionTable(conversionFile)
    #Open the input and output files as Tab-delimited files

    tempRoot='.'.join( inSam.split('.')[:-1])
    tempFile=tempRoot+'_consensus.sam'
    tempHandle=open(tempFile, 'w')
    tempTable=csv.writer(tempHandle , delimiter='\t')

    inHandle=open(inSam, 'r')
    inTable=csv.reader(inHandle, delimiter='\t')
    convertCount=0.
    total=0.
    timer_conv=0.
    timer_md=0.
    for row in inTable:
        conv_start=time.clock()
        #Handle the header rows
        if row[0][0]=='@':
            if row[0]=='@SQ':
                name=row[1].split(':')[-1]
                if consDict.has_key(name)==True:
                    tempTable.writerow(row)
                continue
            else:
                tempTable.writerow(row)
                continue

        #Handle the non-header rows
        line=SamLine(row) #Organize the information in this row using the SamLine class
        if consDict.has_key( line.contig)==True:    #Read aligns to consensus, don't change anything about the sam line
            tempTable.writerow(row)
        elif insDict.has_key(line.contig)==True:    #Read aligns to insertion, convert its position to consensus
            total+=1
##            newPos, newCigar= CombineAlignments(line.pos, insDict[line.contig].insStart\
##            ,insDict[line.contig].insEnd, line.cigar, insDict[line.contig].YIntercepts)
            newPos,newCigar=MapToConsensus(line.pos, insDict[line.contig].insStart, line.cigar, insDict[line.contig].YIntercepts, insDict[line.contig].slope, insDict[line.contig], len(line.seq))

            #Calculate how long the alignment to the insertion was
            originalCigar=ParseCIGAR(line.cigar)
            originalLength=originalCigar['M']+originalCigar['D']

            #Calculate how long the new alignment to the consensus is
            newCigarParts=ParseCIGAR(newCigar)
            newLength=newCigarParts['M']+newCigarParts['D']
            if newCigar=='':   #Read is unmapped
                line.pos=0
                line.cigar='*'
                line.contig='*'
                line.flag=4
                tempTable.writerow(line.row())
                continue
            if newCigarParts['M']+newCigarParts['I']!=len( line.seq):
                print newCigar
                print jabber
            slope=insDict[line.contig].slope
            if insDict[line.contig].slope==-1: #The read is now on the opposte strand, so note that and reverse the CIGAR string
                line.seq=str(Seq.Seq(line.seq).reverse_complement())
                line.qual=line.qual[::-1]
                line.flag=16-line.flag
                newCigar=ReverseCIGAR(newCigar)

            if numpy.isnan( newPos)==False and .5*originalLength<=newLength<=1.5*originalLength:
                #Provided the read contains positions represented in the consensus,
                #and its alignment to the consensus is no less that 50% of the original
                #aligned length and no greater than 150% assign it to the consensus
                #convert the sam line to the consensus alignment
                convertCount+=1
                line.pos=int(newPos)
                line.cigar=newCigar
                line.contig=insDict[line.contig].consName
            else:
                #Otherwise, convert it to an unmapped read
                line.pos=0
                line.cigar='*'
                line.contig='*'
                line.flag=4
            timer_conv+=time.clock()-conv_start
            #This is a place where the MD string can be updated
            #We have: The new cigar, the new position, the read sequence and the read quality string
            #The only thing we are missing is the consensus sequence, but that is easily remedied.

            #The functions that go here will need to take as input SamLine objects
            #We need to obtain two strings, one for the read and one for the consensus,
            #containing the nucleotides at matched positions. We need also to account for strandedness
            #That is, the consensus sequence needs to come from the proper strand.
            md_start=time.clock()
            QD_string=''
            if line.contig!="*":
                MD_string, QD_string=UpdateMDString(line, consSequences,slope)

                #Find the MD string in the optional flags
                optional_headers=[flag.split(':')[0] for flag in line.optional]
                MD_index=optional_headers.index('MD')

                #Replace the original MD string with the updated one
                line.optional[MD_index]=MD_string
            timer_md+=time.clock()-md_start
            #write the samline
            tempTable.writerow(line.row() + [QD_string])
##            print newPos, newCigar
        else:   #Read is unmapped
            line.pos=0
            line.cigar='*'
            line.contig='*'
            line.flag=4
            tempTable.writerow(line.row())
    inHandle.close()
    tempHandle.close()
    os.remove(inSam)
    os.rename(tempFile, inSam)
    print timer_conv, timer_md
    return convertCount, total, convertCount/total


def UpdateMDString(line, consDict,slope):

#Does not produce a proper MD string. Rather it only considers nucleotides at aligned positions.
#SAM specification MD includes deletions, that information is already present in the CIGAR
#MD string in normal SAM spec indicates the reference position, not the read's nt.
##    line=SamLine(line) #delete this. It's just to tell the IDE what I'm doing
    CIGAR_array=ExpandCIGAR(line.cigar)
##    print CIGAR_array
    cons_seq=consDict[line.contig].seq.upper()

    #Determine the interval
    align_begin=line.pos-1      #Python uses a zero-coordinate system, while SAM format uses 1-coordinate
    aligned_length=((CIGAR_array=='M')+(CIGAR_array=='D')).sum()
    align_end=align_begin+aligned_length+1

    cons_align=cons_seq[align_begin:align_end]
##    if slope!=1:
##        print line.flag
##        print str(cons_align)
##        print str(cons_align)
        #Confuse homophones and take back all of the nice things you've ever said about the consensus

##    if slope==-1: cons_align= cons_align.reverse_complement()

    #Don't need to reverse complement if flag=16 because Bowtie2 already does that automatically

    #Positions that are insertions relative to the reference are present in the read
    read_cigar=CIGAR_array[(CIGAR_array=='M')+(CIGAR_array=='I')]
##    print read_cigar
    #Positions that are deletions relative to the reference are present in the reference
    cons_cigar=CIGAR_array[(CIGAR_array=='M')+(CIGAR_array=='D')]
##    print cons_cigar

    #Convert the read and the quality string into numpy arrays for easy indexing
    read_array=numpy.array( [i for i in line.seq.upper()])
    qual_array=numpy.array([i for i in line.qual])
    #Convert the corresponding region of the consensus to an array
    cons_array=numpy.array( [i for i in cons_align])

    #Identify regions of each array where the read positions are present in the consensus (and vice versa)
    read_matches=read_array[read_cigar=='M']
    cons_matches=cons_array[cons_cigar=='M']
    qual_matches=qual_array[read_cigar=='M']

    #Determine  those aligned positions where the read matches the consensus
    mismatch_positions = read_matches==cons_matches

##    print mismatch_positions

    #Pull out the nucleotides in the read that do not math the consensus
    mismatch_values=list( read_matches[~mismatch_positions] )+[''] #Need to pad the end of this with an empty string to ensure
    #Get the corresponding quality scores
    qual_values=list( qual_matches[~mismatch_positions] )+['']

##    print mismatch_values
    mismatch_indices=list( numpy.where(mismatch_positions==False)[0])
    padded_indices= [-1]+ mismatch_indices + [len(read_matches)]
    match_lengths=numpy.diff( padded_indices)-1     # Why???


    MD_string="MD:Z:{0}".format( ''.join(chain( *zip(match_lengths.astype(str), mismatch_values))))
    qual_string="QD:Z:{0}".format( ''.join(qual_values))
##    if slope!=1:
##        print read_matches
##        print cons_matches
##        print mismatch_positions
##        print mismatch_values
##        print match_lengths
##        print line.cigar, line.flag, MD_string,  match_lengths.sum()+len(mismatch_values)-1, len(read_matches)
    return MD_string, qual_string

def ReadConversionTable(conversionFile):
    consensusDict={}
    insertionDict={}
    convHandle=open(conversionFile, 'r')
    convTable=csv.reader(convHandle, delimiter='\t')
    header=convTable.next()
    for row in convTable:
        if row[2].lower()=='true':
            consensusDict[row[0]]=True
            continue
        line=ConversionInfo(row)
        insertionDict[line.insName]=line
    return consensusDict, insertionDict

##from matplotlib import pyplot

def AlignmentToIntercepts(positions, slope):
    #Each nan indicates an insertion
    #Each insertion also decreases the y-intercept by one
    #For simplicity, I'm going remove the effect on the y-intercept, so that
    #Nans are the only indicator of an insertions
    nan_pos=numpy.isnan(positions)

##    return positions
    insertion_offset=numpy.cumsum(nan_pos)
##    return insertion_offset
    if slope==1:
        expected_aligment=numpy.arange(len(positions))
    if slope==-1:
        expected_aligment=numpy.arange(len(positions)-1,-1, -1) #Remember: Python intervals are half-closed
##        expected_aligment=numpy.arange(len(positions),0)
    intercepts= positions-positions[~nan_pos][0]-expected_aligment +slope* insertion_offset
##    return intercepts
##    return intercepts
    non_nan_positions=intercepts[~nan_pos]
    #I am about to take a difference, so I need to pad the beginning
    padded_positions=numpy.hstack(([non_nan_positions[0]], non_nan_positions))
    diff_array=numpy.diff(padded_positions)
##    return diff_array
    final_array=numpy.array([numpy.nan]*len(positions))
    final_array[~nan_pos]=diff_array
    return final_array
def ConvertInterceptToCigar(intercept):
    string=''
    for i in intercept:
        if numpy.isnan(i)==True: string+='I'
        elif i<0:
            string+='D'*abs(i)
            string+='M'
        else:
            string+='M'
    return string

def ConvertInterceptToCigarEfficient(intercept):
    string=''
##    print intercept
    nan_removed=numpy.array( [0]*len(intercept))
    ins_pos=numpy.where(numpy.isnan(intercept)==True)[0]
    del_pos=numpy.where(abs( intercept) >0)[0]
    nonmatch_pos=numpy.array( sorted(list(ins_pos)+list(del_pos)))
##    return ins_pos, del_pos, nonmatch_pos
    nan_removed[~numpy.isnan(intercept)]=intercept[~numpy.isnan(intercept)]
    number_deletions=numpy.sum(abs(nan_removed))
    pos_deletions=numpy.cumsum(abs(nan_removed))
    number_insertions=numpy.sum( numpy.isnan(intercept))

    #We want to create a string that describes all of the edit operations in the alignment

    #The total number of operations is equal to the length of the intercept array plus the number of deletions
    #Make an array of that lenght and fill it with 'M' for match
    pos_array=numpy.array(['M']*(len(intercept)+number_deletions))
##    print ins_pos.shape
##    print pos_deletions.shape
    #Deletions act as spacers, so we need to shift each operation left by the number of deletions that preceded it
    #ins pos is the index of all nan values in the intercept array
    insertion_index=ins_pos +pos_deletions[ins_pos]
    pos_array[insertion_index]='I'
##    del_pos+=pos_deletions[del_pos]
    #Now add deletions into the operation array
    del_offset=0.
##    print del_pos
    for del_index in del_pos:
        del_size=abs( intercept[del_index])

        del_array= numpy.array( ['D']* del_size)
        pos_array[del_index +del_offset:del_index+del_size+del_offset]=del_array
        del_offset+=del_size

    #We padd the initial array because we are going to identify indices where the preceding operation is different
    pad_array=numpy.hstack((pos_array[0], pos_array ))
    diff_array=(pad_array[1:]==pad_array[:-1])
    false_indices=numpy.where(diff_array==False)[0]
    false_indices=numpy.hstack(([0], false_indices,[len(pos_array)] ))
    diff_indices=numpy.diff(false_indices)
    operation=pos_array[false_indices[1:]-1 ]
    cigar= ''.join([''.join(s) for s in zip(diff_indices.astype(str), operation)])
##    return diff_array, pos_array, false_indices, diff_indices, operation

    return cigar

def MapToConsensus(readStart, insertionStart, cigar_string, yintercepts,slope,conv, read_len):
    #Maybe use the updated positions to define the new cigar, because that is a function which directly maps the
    #read to the consensus.

    #Otherwise, rewrite this bit of code to make it clearer
##    try:
    cigar=SplitCigar(cigar_string)

    #build a mapping from read to insertion position: what insertion positions does the read cover?
    #A match in relative to the individual elements does not change the y-intercept
    #A deletion in the read relative to the individal element increases the y-intercept by the number of deleted positions
    #An insertion in the read relative to the individual element is represented as a NaN and decrements
    #the y-intercept at the next matching position by the number of inserted nucleotides
    f_prime= DetermineYIntercepts(insertionStart, readStart, cigar)

    f_prime_nans=numpy.isnan(f_prime)
    f_prime_insertions=numpy.cumsum(f_prime_nans)
    insertionCoordinates=f_prime+numpy.array(range(len(f_prime)))+insertionStart-1  #Python uses a zero-based coordinate system
##    print insertionCoordinates
    #NaNs are positions of the read that are not present in the insertion
    insertionPositions=insertionCoordinates[~numpy.isnan(insertionCoordinates)].astype(int)
##    insertionPositions=insertionCoordinates.astype(int)
##    pyplot.plot(insertionPositions)
##    pyplot.show()
##    pyplot.close()
##        print insertionPositions
    #convert to consensus positions
##    interceptsIC=yintercepts[insertionCoordinates]


    try:
        interceptsIC=numpy.array([numpy.nan]*len(insertionCoordinates))

##        interceptsIC=insertionCoordinates
        interceptsIC[~numpy.isnan(insertionCoordinates)]=yintercepts[insertionPositions]
    except:
        print len(yintercepts)
        print len(insertionPositions)
        print max(insertionPositions)
##        print interceptsIC
    #Ignores insertions relative to the individual sequence
##    newIntercepts=interceptsIC+f_prime[~numpy.isnan(f_prime)]
    newIntercepts=interceptsIC+  f_prime
##    Intercepts_with_nan= interceptsIC+f_prime
##    pyplot.plot(newIntercepts, c='b')
##    pyplot.plot(interceptsIC+f_prime, c='r', ls='--')
##    pyplot.show()
##    pyplot.close()
##        print newIntercepts
##    newPositions=interceptsIC+ slope* insertionPositions+1
    newPositions=interceptsIC+ slope* insertionCoordinates+1*slope
##    return newPositions
    nonNaNPos=newPositions[~numpy.isnan( newPositions)]
##    pyplot.plot(newPositions)
##    pyplot.show()
##    pyplot.close()
    if len(nonNaNPos)==0:
        newPos=numpy.nan
        return newPos, ''
    else:
        if slope==1:
##            newPos=nonNaNPos[0]
            newPos=nonNaNPos[0]
            newPos_1=min(nonNaNPos)
##            print slope, newPos, newPos_1
        else:
            newPos=nonNaNPos[-1]
            newPos_1=min(nonNaNPos)
##            print slope, newPos, newPos_1
##    print newPositions
##    print newIntercepts
##    print interceptsIC
##    print jabber
##    if slope==1:
    temp_intercepts=AlignmentToIntercepts(newPositions, slope)
    newCigar=ConvertInterceptToCigarEfficient(temp_intercepts)
##    if slope==-1:
##        newCigar=InterceptToCIGAR(reversed( -1* newIntercepts))
##    print newCigar
    newCigarParts=ParseCIGAR(newCigar)
##    if newCigar!='' and newCigarParts['M']+newCigarParts['I']!=read_len:
##        pyplot.plot(interceptsIC, c='r')
##        pyplot.plot(newIntercepts, c='g', ls='--')
##        pyplot.plot(f_prime, c='b', ls='.')
##        pyplot.show()
##        pyplot.close()
##        pyplot.plot(newPositions, c='r')
##        pyplot.show()
##        pyplot.close()
##        print newCigar, read_len
##        return newPositions, False
##        print jabbere
##    newCigar=InterceptToCIGAR(Intercepts_with_nan)
##    print newPos
    if newPos<0:
        print conv.insStart, conv.insEnd, conv.consStart, conv.consEnd
        print  conv.consStart - conv.slope* conv.insStart
        print conv.slope

        print readStart, insertionStart, cigar_string
        print f_prime
        print insertionCoordinates
        print insertionPositions
        print "Y intercepts"

        print yintercepts
        print newPositions

        print newIntercepts
        print jabber

        print newPositions

    return newPos, newCigar

def InterceptToCIGAR(intercept):

    nonNaNs=numpy.isnan(intercept)
    if nonNaNs.sum()==len(intercept): return ''
    #If the intercept begins with NaNs, it's been softclipped; add an insertion to the beginning
    begin_pad_length=numpy.where(nonNaNs==False)[0][0]
    end_pad_length=len(intercept)- numpy.where(nonNaNs==False)[0][-1] -1
    nanRemoved=intercept[~nonNaNs]
    diffArray=numpy.hstack(([0], numpy.diff(nanRemoved),[0]))
##    print diffArray
    nonzeroes=list(numpy.where(diffArray!=0)[0])
    nonzeroes.insert(0,0)
    nonzeroes.append(len(diffArray)-1)
    cigar=''
    if begin_pad_length>0: cigar+=str(begin_pad_length)+'I'
    for x in range(len(nonzeroes)-1):
        if diffArray[nonzeroes[x]]==0:
            cigar+=str(nonzeroes[x+1]-nonzeroes[x])+'M'
        if diffArray[nonzeroes[x]]>0:
            cigar+=str(int(abs( diffArray[nonzeroes[x]])))+'D'
            cigar+=str(nonzeroes[x+1]-nonzeroes[x])+'M'
        if diffArray[nonzeroes[x]]<0:
            cigar+=str(int(abs( diffArray[nonzeroes[x]])))+'I'
            cigar+=str(nonzeroes[x+1]-nonzeroes[x])+'M'
    if end_pad_length>0: cigar+=str(end_pad_length)+'I'
    return cigar



def SplitCigar(CIGAR):
    """Reads a sam CIGAR string to count deletions and insertions and matches"""

    parts={'M':0, 'I':0,'D':0}
    cigar=numpy.array(list(CIGAR))

    m=numpy.where(cigar=='M')
    i=numpy.where(cigar=='I')
    d=numpy.where(cigar=='D')
    M=['M']*len(m[0])
    I=['I']*len(i[0])
    D=['D']*len(d[0])
    indices=numpy.concatenate((m,i,d),-1)
    values=numpy.concatenate((M,I,D),-1)
    pairs=[]
    for w in range(len(indices[0])):
        pairs.append((indices[0][w], values[w]))
    pairs=sorted(pairs, key=lambda x: x[0])
    last=0
    cigarParts=[]
    for p in pairs:
        val=int(CIGAR[last:p[0]])
        cigarParts.append ((p[1], val))
        last=1+p[0]
    return(cigarParts)



def ExpandCIGAR(CIGAR):
    """Reads a sam CIGAR string to count deletions and insertions and matches"""

    parts={'M':0, 'I':0,'D':0}
    cigar=np.array(list(CIGAR))

    m=np.where(cigar=='M')[0]
    i=np.where(cigar=='I')[0]
    d=np.where(cigar=='D')[0]
##    M=['M']*len(m[0])
##    I=['I']*len(i[0])
##    D=['D']*len(d[0])
    indices=sorted(list( np.concatenate((m,i,d),-1)))
    indices=[-1]+indices
##    print indices
    CIGAR_list=[]
    for i in range(len(indices)-1):
        chunk= CIGAR[indices[ i]+1 :indices[ i+1]+1]
        count, val=int(chunk[:-1]),chunk[-1]
        CIGAR_list+=[val]*count
    return numpy.array( CIGAR_list)



def ParseCIGAR(CIGAR):
    """Reads a sam CIGAR string to count deletions and insertions and matches"""

    parts={'M':0, 'I':0,'D':0}
    cigar=np.array(list(CIGAR))

    m=np.where(cigar=='M')
    i=np.where(cigar=='I')
    d=np.where(cigar=='D')
    M=['M']*len(m[0])
    I=['I']*len(i[0])
    D=['D']*len(d[0])
    indices=np.concatenate((m,i,d),-1)
    values=np.concatenate((M,I,D),-1)
    pairs=[]
    for w in range(len(indices[0])):
        pairs.append((indices[0][w], values[w]))
    pairs=sorted(pairs, key=lambda x: x[0])
    last=0
    for p in pairs:
        val=int(CIGAR[last:p[0]])
        parts[p[1]]+=val
        last=1+p[0]
    return(parts)


def ReverseCIGAR(cigar):
    """Reverse CIGAR string"""
    cigarArray=numpy.fromstring(cigar, '|S1')
    mLoc=cigarArray=='M'
    dLoc=cigarArray=='D'
    iLoc=cigarArray=='I'
    signifierLoc= [0]+list( numpy.where( mLoc+dLoc+iLoc==True)[0]+1)
##    for i in range (len(signifierLoc)-2,-1,-1 ):
##        print cigar[signifierLoc[i]: signifierLoc[i+1]]
    reversed_CigarList=[cigar[signifierLoc[i]: signifierLoc[i+1]]   for i in range (len(signifierLoc)-2,-1,-1 )]
    return ''.join(reversed_CigarList)



def main(argv):
    param={}
    print argv
    if len(argv)==1: return
    if argv[1]=='-help':
        print """I am currently too lazy to write instructions."""
        return
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param

    if param=={}:
        return

    try:


        if param['-fxn'].lower()=='build':
            consFile=param['-cons']
            refFile=param['-ref']
            gffFile=param['-gff']
            outfile=param['-out']
            blastdir=param['-BLAST']
            PrepareConversionTable(consFile, refFile, gffFile, outfile, blastdir)
        elif param['-fxn'].lower()=='convert':
            samFile=param['-i']
            convFile=param['-conv']
            convertCount, total, percent=ConvertSAMFromInsertionsToConsensus(samFile, convFile)
            print "Converted {0} out of {1} insertion alignments to consensus alignments ({2}%)".format(convertCount, total, percent*100.)

    except:
        print """python AlignmentConverter.py -help for instructions."""


    #inDir, outDir, AssemblyIndex, TEIndex, threads, length=70, ReadIDPosition=1, PositionRange=0, phred=33, shorten=True, Trim=True



if __name__ == '__main__':
    main(sys.argv)
