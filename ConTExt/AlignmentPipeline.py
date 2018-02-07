#-------------------------------------------------------------------------------
# Name:        ConTExt
# Purpose: Align and organize paired-end reads in a way that makes thinking about repeat-derived reads tractable.
#
# Author:      Michael
#
# Created:     24/04/2015
# Copyright:   (c) Michael Peter McGurk 2015
# Licence:

#The MIT License (MIT)

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Name:        conTExt
# Purpose:  To interpret paired end data aligned to the assembled genome and
#   to TE databases and identify the genomic contexts of TEs
#
# Author:      Michael
#
# Created:     02/04/2014
# Copyright:   (c) Michael 2014
# Licence:     <your licence>
#Version 1.1
#-------------------------------------------------------------------------------

import sys
import Bio
import time
import numpy as np
import scipy as sc
import os
import csv
import gzip
import subprocess
import shutil
##import ExtractSeq
import tools.AlignmentConverter as AlignmentConverter

from numpy import array
from scipy import stats
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from numpy import random

import sys


BowtieOpt=['-v']
BowtiePar=['2']

version='conTExt-3.1'

##chrlist=['2L','2LHet','2R','2RHet','3L','3LHet','3R','3RHet','4','U','Uextra','X','XHet','YHet','dmel_mitochondrion_genome']

Strand={'0': '+', '4':'*', '16': '-'}

Phred={33:'-phred33', 64:'-phred64'}

#TrimPath='c:/trimmomatic-0.32/trimmomatic-0.32/trimmomatic-0.32.jar'

AdtSettings='4:30:10'

FastQext=set(['fq', 'fastq'])

cutoff={'as': -.64, 'te': -1.0}

class Parameters:
    AsBowtieParameters=''
    TEBowtieParameters=''
    #TrimPath='c:/trimmomatic-0.32/trimmomatic-0.32/trimmomatic-0.32.jar'
    #BowtiePath=''
    #SamtoolsPath= 'c:/samtoolswin'


    def __init__(self, conFile):
        handle=open(conFile, 'r')
        varDict={}
        for line in handle:
            variable=line.split('=')
            varDict[variable[0]]=variable[1].strip('\n')
        Parameters.BowtiePath=varDict['B']
        Parameters.TrimPath=varDict['T']
        Parameters.SamtoolsPath=varDict['S']
        handle.close()

class Read():
    ID=0
    def __init__(self, mstrand, mstart, mend, mseq, useq, sequence='', quality=''):
        self.strand=mstrand
        self.start=int(mstart)
        self.end=int(mend)
        self.mapped=mseq
        self.unmapped=useq
        self.id=Read.ID
        self.seq=sequence
        self.qual= quality
        Read.ID+=1


def configure (conFile):
    """Reads the specifiec configuration file that specifies the paths to Bowtie2,
    Trimmomatic, and SamTools."""
    handle=open(conFile, 'r')
    varDict={}
    for line in handle:
        variable=line.split('=')
        varDict[variable[0]]=variable[1].strip('\n')
    BowtiePath=varDict['B']
    TrimPath=varDict['T']
    SamtoolsPath=varDict['S']
    handle.close()
    return BowtiePath, TrimPath, SamtoolsPath

def DeterminePHRED(infile):
    if infile.split('.')[-1]!='.gz':
        handle=open(infile, 'r')
    else: handle=gzip.open(infile,'r')
    lib=FastqGeneralIterator(handle)

    read_count=0
    min_count=0
    max_count=0

    for n, s, q in lib:
        qual_list=[ord (qual) for qual in q]
        min_qual=min(qual_list)
        max_qual=max(qual_list)
##        print min_qual, max_qual
##        print jabber
        if max_qual>80: max_count+=1
        if min_qual<60: min_count+=1

        if max_count>100:
            phred=64
            break
        if min_count>100:
            phred=33
            break
        if read_count>100000:
            phred=0
            break
    print max_count, min_count
    return phred
def FindInput(inDir):
    """Examines all files in the input directory to identify paired fastq files.
    Recognizes files with the extensions *.fq and *.fastq, and will recognize
    gzipped files if they have the extension *.gz. Paired-end datasets must have a
    filename which follows this format:
        <sample_name>_{1/2}.{fastq/fq}(.gz)
    Returns a list of sample names, and a dictionary which stores the appropriate
    extension for each sample."""
    FList=os.listdir(inDir)
    FQList=[]       #Hold all fastq files
    TargetList=[]   #Hold the filename prior the extension
    ExtDict={}      #pair the extension to the filename

    #find all fastq files
    for f in FList:
        ext=f.split('.')
        if len(set([ext[-1]])&FastQext)!=0:
            FQList.append(f)
            TargetList.append('.'.join(ext[:-1])) #Store the name prior to .ext
            ExtDict['.'.join(ext[:-1])]='.'.join(ext[-1:])   #Key ext by name
        if ext[-1]=='gz':   #also handle '.gz'
            if len(set([ext[-2]])&FastQext)!=0:
                FQList.append(f)
                TargetList.append('.'.join(ext[:-2]))
                ExtDict['.'.join(ext[:-2])]='.'.join(ext[-2:])

    #Get the file roots--dataset name without the end number:
    RootSet=set()
    for t in TargetList:
        RootSet.add('_'.join(t.split('_')[:-1]))
    RootList=list(RootSet)

    #keep only the fastq files that are paired
    FileList=[]
    FileExt={}
    for t in RootList:
        e1name=t+'_1'
        e1count=TargetList.count(e1name)
        e2name=t+'_2'
        e2count=TargetList.count(e2name)
        if e1count==1 and e2count==1:
            if t=='': continue
            FileList.append(t)
            FileExt[e1name]=ExtDict[e1name]
            FileExt[e2name]=ExtDict[e2name]

    return( FileList, FileExt)



def TrimReads(inDir, root, fType, outDir,threads, logfile, TrimPath, phred=64):
    """Calls Trimmomatic to preprocess reads. Uses settings:

    ILLUMINACLIP:adapters/TruSeq2-PE.fa:4:30:10
    SLIDINGWINDOW:4:20
    MINLEN:40

    It then checks to see whether both processed PE files were created. If not, the process failed
    for some reason. If this happens, it returns a variable to indicate the pipeline
    should move on to the next dataset and record the reason for doing so."""
    #java -jar <path to trimmomatic.jar> PE [-threads <threads]
    #[-phred33 | -phred64] [-trimlog <logFile>] <input 1> <input 2>
    #<paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2>
    #<step 1> ...

#ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
    #TrimPath=Parameters.TrimPath
    os.chdir('/'.join(TrimPath.split('/')[:-1]))
    InRoot='/'.join([inDir, root])
    OutRoot='/'.join([outDir, root])

    file1=InRoot+'_1.'+fType[root+'_1']
    file2=InRoot+'_2.'+fType[root+'_2']
    goodfile1=OutRoot+'_1_p.fq'
    goodfile2=OutRoot+'_2_p.fq'
    badfile1=OutRoot+'_1_s.fq'
    badfile2=OutRoot+'_1_s.fq'
    #AdtFile='thornton_adt.fa'
    AdtFile='adapters/TruSeq2-PE.fa'

    TempLog=open(OutRoot+'_log.txt', 'w')

    #Run Trimmomatic and create a log file
    p=subprocess.Popen(['java', '-jar', TrimPath, 'PE', Phred[phred], '-threads', threads, file1, file2\
    , goodfile1, badfile1, goodfile2, badfile2, 'ILLUMINACLIP:'+AdtFile+':'+AdtSettings,'SLIDINGWINDOW:'+'4:20', 'MINLEN:'+'40'], \
    stderr=TempLog)
    #p=subprocess.Popen(['java', '-jar', TrimPath, 'PE', Phred[phred], file1, file2\
    #, goodfile1, badfile1, goodfile2, badfile2, 'ILLUMINACLIP:'+AdtFile+':'+AdtSettings], \
    #stderr=TempLog)
    p.communicate()

    #Check to see whether both PE files were created. If not, the process failed
    #some reason. If this happens the pipeline should move on to the next dataset
    #and record the reason for doing so.
    TempLog.close()

    handle=open(OutRoot+'_log.txt', 'r')
    for line in handle:
        logfile.write(line)


    CheckFiles=os.listdir(outDir)
    if CheckFiles.count(goodfile1.split('/')[-1])!=1 or \
    CheckFiles.count(goodfile2.split('/')[-1])!=1:
        Failed=True
    else: Failed=False

    try:
        if os.path.exists(badfile1)==True:
            os.remove(badfile1)
        if os.path.exists(badfile2)==True:
            os.remove(badfile2)
    except Exception:
        print ('Delete failed...')

    print ('Operation failed? '+str(Failed))

    return (Failed)

def RenameReads(length, inDir, root, ends, logfile):
    """For ConTExt to work, it requires that each read pair be assigned an integer
    identifier, and have a name of the format:

    <Sample_identifier>.<Read #>.<1/2>

    If this is not the case, the pipeline can be told to use this function to rename
    the reads.
    """
    total=0

    for e in ends:
        target=inDir+'/'+'_'.join([root, str(e), 'p'])+'.fq'
        temp=inDir+'/'+'_'.join([root, str(e), 'p','temp'])+'.fq'
        handle=open(target, 'r')
        lib=FastqGeneralIterator(handle)
        outtemp=open(temp, 'w')
        count=0
        for n, s, q in lib:
            count+=1
            new_name='{0}.{1}.{2}'.format(root, count, e)
            outtemp.write("@%s\n%s\n+\n%s\n" % ( n , s, q ) )
        outtemp.close()
        handle.close()
        print ('End {0} renamed.'.format(e))
        os.remove(target)
        os.rename(temp, target)
    return count, total


def ShortenReads(length, inDir, root, ends, logfile, rename=False, GemCode=False):
    """Because ConTExt aims to identify junctions, trimming sequence from the
    3'-end of reads increases the mate pair distance in each read pair, increasing its power.
    This function trims each read longer than the specified length."""

    total=0
    count=0
    for e in ends:
        target=inDir+'/'+'_'.join([root, str(e), 'p'])+'.fq'
        temp=inDir+'/'+'_'.join([root, str(e), 'p','temp'])+'.fq'
        handle=open(target, 'r')
        lib=FastqGeneralIterator(handle)
        outtemp=open(temp, 'w')
        read_count=0
        for n, s, q in lib:
            total+=1
            read_count+=1
            if rename==True:
                n='{0}.{1}.{2}'.format(root, read_count, e)
            if GemCode==True:
                barcode=n.split('.')[-1]
                n='{0}.{1}.{2}.{3}'.format(root, read_count, e, barcode)
            if len(s)>=length:
                count+=1
                ts, tq=s[:length], q[:length]
                outtemp.write("@%s\n%s\n+\n%s\n" % ( n , ts, tq ) )
            else: outtemp.write("@%s\n%s\n+\n%s\n" % ( n , s, q ) )
        outtemp.close()
        handle.close()
        print ('End {0} shortned.'.format(e))
        os.remove(target)
        os.rename(temp, target)
    return count, total


def CallAligner(inFile, outFile, Index, cutoff, seedLength, threads, runlog, phred=33, bowDir=''):
    """Calls BowTie2. Function parameters correspond to :
        threads: NUmber of processes to run in parallel
        phred: PHRED quality encoding of the data,
        cutoff: The alignment score cutoff: (--score-min L,0,cutoff)
        seddLength: the seed length ('-L seedLength')
        Index: The Bowtie2 index to used in the alignments ( '-x Index')
        inFile: input fastq  ('-U inFile')
        outFile: output *.sam ('-S outFile')]
        """


    #Align fastq file using bowtie
    #It might be prudent to only output aligned reads.
    print bowDir
    os.chdir(bowDir)
    logname='.'.join(outFile.split('.')[:-1])+'_bowtie_log.txt'
    logfile=open(logname, 'w')


    proc=subprocess.Popen([bowDir, '-p', threads , '-'+str(Phred[ phred]), '--score-min', 'L,0,'+ str(cutoff), '-L', str(seedLength),  '-x', Index, '-U', inFile, '-S', outFile], stderr=logfile)
    print ("Aligning {0} to {1}...".format(inFile, Index ))
    time.sleep(.5)
    crash=CheckForCrash(logname)
    if crash==True:
        proc.terminate()
    else: proc.communicate()
    #Check out file for the header lines:
    #run_failed=CheckBowtie2Log(logname)
    run_failed=False
    logfile.close()

    logfile=open(logname, 'r')
    for l in logfile:
        runlog.write(l)
    runlog.write("Run failed="+str(run_failed)+"/n")
    logfile.close()

    return (run_failed)


def CheckBowtieLog(LogFile):
    try:
        count=0
        handle=open(LogFile, 'r')
        for sentence in handle:
            word=sentence.split(' ')
            if word[0]=='#' and word[1]=='reads':
                count+=1
        if count>=3:
            failed=False
        else: failed=True
        handle.close()
    except Exception: failed=True
    return (failed)

def CheckBowtie2Log(LogFile):
    failed=False
    try:
        count=0
        handle=open(LogFile, 'r')
        for sentence in handle:
            hold=sentence
        words=sentence.split(' ')
        if ' '.join(words[-3])!='overall alignment rate' or ' '.join(words[-3])!='overall alignment rate\n':
            failed=True
        handle.close()
    except Exception: failed=True
    return (failed)


def CheckForCrash(LogFile):
    """Checks whether Bowtie2 crashed."""
    crashed=False
    try:
        handle=open(LogFile, 'r')
        for sentence in handle:
            word=sentence.split('.')
            for w in word:
                if w=='This application has requested the Runtime to terminate it in an unusual way':
                    crashed=True
        handle.close()
    except Exception: crashed=False
    return (crashed)

def SortAlignment (inDir, Root, E, I, runlog, threads=1, samDir=''):
    """Call Samtools to sort alignments by read number to simplify matching read pairs."""
    sortFailed=False
    #Move to the SamTools directory
    os.chdir(samDir)

    logName=inDir+'/'+Root+'_errlog.txt'
    logfile=open(logName, 'a')
    for e in E:
        if sortFailed==True: break
        for i in I:

            samRoot='_'.join([Root,e,'p',i])
            inSam=inDir+'/'+samRoot+'.sam'
            outBam=inDir+'/'+samRoot+'.bam'
            sortBam=inDir+'/'+samRoot+'_sorted'
            outSam=inDir+'/'+samRoot+'_sorted.sam'
            logfile.write('\nSorting {0}...\n'.format(samRoot))
            logfile.flush()
    #Convert the same file to an uncompressed Bam file
            try:
                logfile.write('\n\nConverting Sam to Bam...\n')
                logfile.flush()
                p=subprocess.Popen([samDir, 'view', '-@', threads, '-u', '-S', inSam, '-o', outBam ], stderr=logfile)
                print ('Converting Sam to Bam...')
                p.communicate()

                logfile.write('\n\nSorting Bam...\n')
                logfile.flush()
                p=subprocess.Popen([samDir, 'sort', '-n','-@', threads,  outBam, sortBam ], stderr=logfile)
                print ('Sorting Bam...')
                p.communicate()

                logfile.write('\n\nConverting Bam to Sam...\n')
                logfile.flush()
                p=subprocess.Popen([samDir, 'view', '-@', threads, '-h', sortBam + '.bam', '-o', outSam ], stderr=logfile)
                print ('Converting Bam to Sam...')
                p.communicate()

            except Exception:
                print('SAMtools failed')
                sortFailed=True
                logfile.close()
                break

            try:
                os.remove(outBam)
                os.remove(sortBam+'.bam')

            except Exception:
                print('Remove failed')
                sortFailed=True
                logfile.close()
                break


            try:
                inSize=int(os.stat(inSam)[6])
                outSize=int(os.stat(outSam)[6])
                print (inSize, outSize)
                #if outSize!=inSize:
                #    sortFailed=True
                #    logfile.close()
                #    break

            except Exception:
                sortFailed=True
                logfile.close()
                break

            try:
                os.remove(inSam)
            except:
                print('Remove fail')

    logfile.close()
    logfile=open(logName, 'r')
    for l in logfile:
        runlog.write(l)
    logfile.close()
    runlog.flush()

    return (sortFailed)

def HandleDualMappings(inDir, Root, option, end, IndexList, CountPosition, PosRange, runlog, suffix='sorted'):
    """Iterate through the Bowtie output and check for reads that mapped to both
    indexes. Do one of three things:

        d - discard both reads
        a - assign to the assembly
        t - assign to TE reference

    Otherwise, output the read as NAME_suffix-> NAME_filtered."""

    Choice={'d': 'discarded', 'a': 'assigned to the assembly', \
    't': 'assigned to the TE index'}

    Failed=False
    try:
    #Open the alignment to the Assembly
        AsIn=inDir+'/'+'_'.join([Root,end,'p',IndexList[0], suffix])+'.sam'
        AsHandle=open(AsIn, 'r')
        AsLib=csv.reader(AsHandle, delimiter='\t')
        AsOut=inDir+'/'+'_'.join([Root,end,'p',IndexList[0], 'filtered'])+'.sam'
        AsOutHandle=open(AsOut, 'w')
        AsCSV=csv.writer(AsOutHandle, delimiter='\t')
    except Exception:
        print ('Failed to open assembly alignment!')
        runlog.write('\nFailed to open assembly alignment!\n')
        Failed=True

        AsHandle.close()
        AsOutHandle.close()

        return (Failed)

    #Open the alignment to the TE index
    try:
        TEIn=inDir+'/'+'_'.join([Root, end,'p',IndexList[1], suffix])+'.sam'
        TEHandle=open(TEIn, 'r')
        TELib=csv.reader(TEHandle, delimiter='\t')
        TEOut=inDir+'/'+'_'.join([Root,end,'p',IndexList[1], 'filtered'])+'.sam'
        TEOutHandle=open(TEOut, 'w')
        TECSV=csv.writer(TEOutHandle, delimiter='\t')
    except Exception:
        print ('Failed to open TE alignment!')
        runlog.write('\nFailed to open TE alignment!\n')
        Failed=True

        AsHandle.close()
        AsOutHandle.close()
        TEHandle.close()
        TEOutHandle.close()

        return (Failed)
    UmOut=inDir+'/'+'_'.join([Root,end,'p','um', 'filtered'])+'.sam'
    UmOutHandle=open(UmOut, 'w')
    UmCSV=csv.writer(UmOutHandle, delimiter='\t')

    EndOfAs=False
    EndOfTE=False

    DualCount=0
    TECount=1
    AsCount=1

    if PosRange==0:
        rowAs=AsLib.next()
        rowTE=TELib.next()

        #Skip the  headers
        while rowAs[0][0]=='@':
            dump=AsCSV.writerow(rowAs)
            rowAs=AsLib.next()

        while rowTE[0][0]=='@':
            dump=TECSV.writerow(rowTE)
            rowTE=TELib.next()



        while EndOfAs==False or EndOfTE==False: #Go until EOF is true for both

            AsPos=CalculateReadNumber(rowAs, CountPosition, PosRange)

            TePos=CalculateReadNumber(rowTE, CountPosition, PosRange)

            if rowAs[2]!='*': AsCount+=1
            if rowTE[2]!='*': TECount+=1

            if TePos<AsPos and EndOfTE==False:
                dump=TECSV.writerow(rowTe)

            if AsPos<TePos and EndOfAs==False:
                dump=AsCSV.writerow(rowAs)

            if AsPos==TePos:
                if rowAs[2]!='*' or rowTE[2]!='*':
                    if rowAs[2]!='*' and rowTE[2]!='*':
                        DualCount+=1
                        if option.lower()=='a':
                            dump=AsCSV.writerow(rowAs)

                        if option.lower()=='t':
                            dump=TECSV.writerow(rowTE)
                    else:
                        if rowAs[2]!='*':
                            dump=AsCSV.writerow(rowAs)
                        if rowTE[2]!='*':
                            dump=TECSV.writerow(rowTE)
                else:
                    dump=UmCSV.writerow(rowAs)

            if EndOfAs==False:
                try:
                    rowAs=AsLib.next()
                except StopIteration:
                    EndOfAs=True

            if EndOfTE==False:
                try:
                    rowTE=TELib.next()
                except StopIteration:
                    EndOfTE=True

    #except Exception:
    #    print ('Failed to filter reads!')
    #    runlog.write('\nFailed to filter reads!\n')
    #    Failed=True

    #    AsHandle.close()
    #    AsOutHandle.close()

    #    TEHandle.close()
    #    TEOutHandle.close()

    #    UmOutHandle.close()

    #    return (Failed)

    AsHandle.close()
    AsOutHandle.close()

    TEHandle.close()
    TEOutHandle.close()

    UmOutHandle.close()


    PercentDualTE=100*float(DualCount)/TECount
    PercentDualAs=100*float(DualCount)/AsCount

    print ("{0} out of {1} assembly alignments also aligned to the TE index ({2}%).".format(DualCount, AsCount, PercentDualAs))
    print ("{0} out of {1} TE alignments also aligned to the assembly ({2}%).".format(DualCount, TECount, PercentDualTE))
    print ("Reads were {0}".format(Choice[option]))

    runlog.write('\n\nFiltering succeeded!')
    runlog.write("\n\n{0} out of {1} assembly alignments also aligned to the TE index ({2}%).".format(DualCount, AsCount, PercentDualAs))
    runlog.write("\n{0} out of {1} TE alignments also aligned the assembly ({2}%).".format(DualCount, TECount, PercentDualTE))
    runlog.write("\nReads were {0}".format(Choice[option]))

    return (Failed)

def RejoinSam (inDir, Root , end, CountPosition, PosRange , IndexList):
    ChrKey={}
    ChrKey['*']='um'
    Error=False
    try:
        AsIn=inDir+'/'+'_'.join([Root,end,'p',IndexList[0], 'filtered'])+'.sam'
        AsHandle=open(AsIn, 'r')
        AsLib=csv.reader(AsHandle, delimiter='\t')

        TEIn=inDir+'/'+'_'.join([Root, end,'p',IndexList[1], 'filtered'])+'.sam'
        TEHandle=open(TEIn, 'r')
        TELib=csv.reader(TEHandle, delimiter='\t')


        UmIn=inDir+'/'+'_'.join([Root, end,'p','um', 'filtered'])+'.sam'
        UmHandle=open(UmIn, 'r')
        UmLib=csv.reader(UmHandle, delimiter='\t')

        OutName=inDir+'/'+'_'.join([Root, end,'p','joined'])+'.sam'
        OutFile=open(OutName, 'w')
        OutWrite=csv.writer(OutFile, delimiter='\t')

    except Exception:
        print ('Open failed!')
        Error=True
        return (Error)

    stop=False
    Readers=[AsLib, TELib, UmLib]
    Rows=[AsLib.next(), TELib.next(), UmLib.next()]

    #Skip headers
    try:
        while Rows[0][0][0]=='@':
            Rows[0]=Readers[0].next()
            if Rows[0][0][1]=='S':
                ChrName=Rows[0][1].split(':')[1:]
                ChrKey[':'.join(ChrName)]='as'
    except StopIteration:
        del Rows[0]
        del Readers[0]
    try:
        while Rows[1][0][0]=='@':
            Rows[1]=Readers[1].next()
            if Rows[1][0][1]=='S':
                ChrName=Rows[1][1].split(':')[1:]
                ChrKey[':'.join(ChrName)]='te'

    except StopIteration:
        del Rows[1]
        del Readers[1]

    try:
        while Rows[2][0][0]=='@':
            Rows[2]=Readers[2].next()
            if Rows[2][0][1]=='S':
                ChrName=Rows[2][1].split(':')[1:]
                ChrKey[':'.join(ChrName)]='um'

    except StopIteration:
        del Rows[2]
        del Readers[2]

    #Get the read numbers:
    ReadNum=[]
    for w in Rows:

        num=CalculateReadNumber(w, CountPosition, PosRange)

        ReadNum.append(num)


    while stop==False:
        target=min(ReadNum)
        TarSpot=ReadNum.index(target)

        #write
        dump=OutWrite.writerow(Rows[TarSpot])
        try:
            Rows[TarSpot]=Readers[TarSpot].next()


            ReadNum[TarSpot]=CalculateReadNumber(Rows[TarSpot], CountPosition, PosRange)

        except StopIteration:
            del Readers[TarSpot]
            del Rows[TarSpot]
            del ReadNum[TarSpot]
            if len(Readers)==0: stop=True

    TEHandle.close()
    os.remove(TEIn)
    AsHandle.close()
    os.remove(AsIn)
    UmHandle.close()
    os.remove(UmIn)
    OutFile.close()
    return(ChrKey)

def CalculateReadNumber(row, CountPosition, Col=3, sep='.'):
    AsPos1=float(row[0].split(sep)[CountPosition])

    if Col==0:
        return(AsPos1)
    else:
        AsPos2=float('.'+row[0].split(sep)[CountPosition+1].split('#')[0])
        AsPos3=float('.000000000000000000000'\
            +row[0].split(':')[CountPosition+1].split('#')[1].split('/')[0])

        num=AsPos1+AsPos2+AsPos3

        return(num)

def GetOutputLine(FirstRow, SecondRow):
    """Concatenates the two information about the two aligned reads in a pair
    into a single output row."""
    #Assume that the reference spanned by the reads is the length of the read
    #unless the CIGAR string disagrees, in which case calculate the range
    #from the cigar
    length=len(FirstRow[9])
    MD_1=GetMDString(FirstRow)
    QD_1=GetQDString(FirstRow)
    if FirstRow[5]!=str(length)+'M':
        parsed=ParseCIGAR(FirstRow[5])
        length=parsed['M']+parsed['D']
    FirstHalf=[FirstRow[2], FirstRow[1], int(FirstRow[3]), int(FirstRow[3])+length, FirstRow[9], FirstRow[10]]
    MD_2=GetMDString(SecondRow)
    QD_2=GetQDString(SecondRow)
    length=len(SecondRow[9])
    if SecondRow[5]!=str(length)+'M':
        parsed=ParseCIGAR(SecondRow[5])
        length=parsed['M']+parsed['D']
    SecondHalf=[SecondRow[2], SecondRow[1], int(SecondRow[3]), int(SecondRow[3])+length, SecondRow[9], SecondRow[10], FirstRow[5], SecondRow[5],FirstRow[4], SecondRow[4], MD_1, MD_2, QD_1, QD_2]

    return (FirstHalf + SecondHalf)

def OrganizeOutput(inDir, Root, END, CountPosition, PosRange, IndexList, GemCode):

    """Converts the four *.sam files (as_1, as_2, te_1, te_2) into two *.sam files,
    and then creates a new output format where each read pair is specified with a row.
    Each sequence in the set of reference contigs and consensus sequences gets
    is given its own output file.

    A read pair which maps only to Repeat_A, will be assigned to Repeat_A.dat.

    A read pair which maps to both Repeat_A and Repeat_B, will be assigned to
    both Repeat_A.dat and Repeat_B.dat.

    A read pair which maps to both Repeat_A and Chr1, will be assigned to only
    Repeat_A.dat, but not Chr1.dat.

    A read pair which maps to only reference contigs will be assigned to those contigs.

    A read pair which has one unmapped reads will be assigned based on the mapped read.

    A read pair where both ends are unmapped is discarded."""
    ChrKey={}
    ChrKey['*']='um'


    #Prepare End 2 for the efficient sorting, by rejoining the component parts
    #in an ordered SAM

    ChrKey = RejoinSam(inDir, Root, END[1], CountPosition, 0, ['as', 'te'])
    ChrKey = RejoinSam(inDir, Root, END[0], CountPosition, 0,['as', 'te'])

    #Opent the rejoined End 2; hereafter this is refered to as Source

    Query={}

    SourceName=inDir+'/'+'_'.join([Root,END[1],'p','joined'])+'.sam'
    SourceHandle=open(SourceName, 'r')
    SourceLib=csv.reader(SourceHandle, delimiter='\t')

    QueryName=inDir+'/'+'_'.join([Root,END[0],'p','joined'])+'.sam'
    QueryHandle=open(QueryName, 'r')
    QueryLib=csv.reader(QueryHandle, delimiter='\t')



    Row={}
    rNumber={}
    written={}
    skip={}
    sCount=0
#Skip headers
    count={}
    SourceRow=SourceLib.next()
    try:
        while SourceRow[0][0]=='@':
            SourceRow=SourceLib.next()
    except:
            return
    Row=QueryLib.next()
    try:
        while Row[0][0]=='@':
            Row=QueryLib.next()
    except:
            return
    keycount=0
    #Make Output Files
    NameKey={}
    outFiles={}

    os.mkdir(inDir+'/'+Root)
    for k in ChrKey.keys():
        if k=='*':continue
        keycount1='_'.join(k.split('/'))
        keycount='_'.join(keycount1.split('|'))
        NameKey[k]=keycount
        outFiles[k]=open(inDir+'/'+Root+'/'+keycount+'.dat', 'w')
        header=[k, 'strand', 'start', 'end','seq', 'qual', 'Context' ,'strand', 'start', 'end','seq', 'qual']
        dump=csv.writer(outFiles[k], delimiter='\t').writerow(header)

    SourceStop=False
    while SourceStop==False:
        #Calculate Read numbers
        SourceNum=CalculateReadNumber(SourceRow, CountPosition, PosRange)
        rNumber=CalculateReadNumber(Row, CountPosition, PosRange)

        #Look for matches and write output if they are found:
        SourceType=ChrKey[SourceRow[2]]
        w=ChrKey[Row[2]]
        if GemCode==True: barcode=[Row[0].split('.')[-1]]
        else: barcode=[]
        if rNumber==SourceNum and SourceType=='te':
            if w=='te':
                line=GetOutputLine(SourceRow, Row)
                csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line+barcode)

                if SourceRow[2]!=Row[2]: #Write to both files if two TEs
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)

                if w=='um':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line+barcode)

                if w=='as':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line+barcode)


        if rNumber==SourceNum and SourceType=='as':
                if w=='te':
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line+barcode)

                if w=='um':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line+barcode)

                if w=='as':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line+barcode)

                    if SourceRow[2]!=Row[2]: #Write to both files if two Chr
                        line=GetOutputLine(Row, SourceRow)
                        csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line+barcode)


        if rNumber==SourceNum and SourceType=='um':
                if w=='te':
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line+barcode)

                if w=='as':
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line+barcode)

                   #Check which files to iterate
        #First: Is the source less than or equal to all files? Iterate it:

        if SourceNum<rNumber:

            try:
                sCount+=1
                SourceRow=SourceLib.next()
                continue
            except StopIteration:
                SourceStop=True

        if rNumber<=SourceRow:
            try:
                Row=QueryLib.next()
            except StopIteration:
                SourceStop=True

    for k in outFiles.keys():
        outFiles[k].close()
    SourceHandle.close()
    QueryHandle.close()

def OrganizeOutputII(inDir, Root, END, CountPosition, PosRange, IndexList):

    ChrKey={}
    ChrKey['*']='um'


    #Prepare End 2 for the efficient sorting, by rejoining the component parts
    #in an ordered SAM

    ChrKey = RejoinSam(inDir, Root, END[1], CountPosition, 0, ['as', 'te'])
    ChrKey = RejoinSam(inDir, Root, END[0], CountPosition, 0,['as', 'te'])

    #Opent the rejoined End 2; hereafter this is refered to as Source

    Query={}

    SourceName=inDir+'/'+'_'.join([Root,END[1],'p','joined'])+'.sam'
    SourceHandle=open(SourceName, 'r')
    SourceLib=csv.reader(SourceHandle, delimiter='\t')

    QueryName=inDir+'/'+'_'.join([Root,END[0],'p','joined'])+'.sam'
    QueryHandle=open(QueryName, 'r')
    QueryLib=csv.reader(QueryHandle, delimiter='\t')



    Row={}
    rNumber={}
    written={}
    skip={}
    sCount=0
#Skip headers
    count={}
    SourceRow=SourceLib.next()
    try:
        while SourceRow[0][0]=='@':
            SourceRow=SourceLib.next()
    except:
            return
    Row=QueryLib.next()
    try:
        while Row[0][0]=='@':
            Row=QueryLib.next()
    except:
            return
    keycount=0
    #Make Output Files
    NameKey={}
    outFiles={}

    os.mkdir(inDir+'/'+Root)
    for k in ChrKey.keys():
        if k=='*':continue
        keycount1='_'.join(k.split('/'))
        keycount='_'.join(keycount1.split('|'))
        NameKey[k]=keycount
        outFiles[k]=open(inDir+'/'+Root+'/'+keycount+'.dat', 'w')
        header=[k, 'strand', 'start', 'end','seq', 'qual', 'Context' ,'strand', 'start', 'end','seq', 'qual']
        dump=csv.writer(outFiles[k], delimiter='\t').writerow(header)

    SourceStop=False
    while SourceStop==False:
        #Calculate Read numbers
        SourceNum=CalculateReadNumber(SourceRow, CountPosition, PosRange)
        rNumber=CalculateReadNumber(Row, CountPosition, PosRange)

        #Look for matches and write output if they are found:
        SourceType=ChrKey[SourceRow[2]]
        w=ChrKey[Row[2]]
        if rNumber==SourceNum and SourceType=='te':
            if w=='te':
                line=GetOutputLine(SourceRow, Row)
                csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line)

                if SourceRow[2]!=Row[2]: #Write to both files if two TEs
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)

                if w=='um':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line)

                if w=='as':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line)


        if rNumber==SourceNum and SourceType=='as':
                if w=='te':
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)

                if w=='um':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line)

                if w=='as':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line)

                    if SourceRow[2]!=Row[2]: #Write to both files if two Chr
                        line=GetOutputLine(Row, SourceRow)
                        csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)


        if rNumber==SourceNum and SourceType=='um':
                if w=='te':
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)

                if w=='as':
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)

                   #Check which files to iterate
        #First: Is the source less than or equal to all files? Iterate it:

        if SourceNum<rNumber:

            try:
                sCount+=1
                SourceRow=SourceLib.next()
                continue
            except StopIteration:
                SourceStop=True

        if rNumber<=SourceRow:
            try:
                Row=QueryLib.next()
            except StopIteration:
                SourceStop=True

    for k in outFiles.keys():
        outFiles[k].close()
    SourceHandle.close()
    QueryHandle.close()


def BuildKey(TE_Index, As_Index):
    ChrKey={}
    ChrKey['*']='um'
    handle=open(TE_Index)
    lib=SeqIO.parse(handle, 'fasta')
    for book in lib:
        ChrKey[book.name]='te'
    handle.close()
    handle=open(As_Index)
    lib=SeqIO.parse(handle, 'fasta')
    for book in lib:
        ChrKey[book.name]='as'

    return ChrKey

def GetMDString(row):
    optional=row[11:]
    optional_headers=[flag.split(':')[0] for flag in optional]
    if optional_headers.count('MD')!=0:
        MD_index=optional_headers.index('MD')
        return optional[MD_index]
    else:
        return ''
def GetQDString(row):
    optional=row[11:]
    optional_headers=[flag.split(':')[0] for flag in optional]
    if optional_headers.count('QD')!=0:
        MD_index=optional_headers.index('QD')
        return optional[MD_index]
    else:
        return ''

def OrganizeOutputIII(inDir, Root, END, CountPosition, PosRange, IndexList):

    ChrKey=IndexList


    #Prepare End 2 for the efficient sorting, by rejoining the component parts
    #in an ordered SAM

    #Opent the rejoined End 2; hereafter this is refered to as Source

    Query={}

    SourceName=inDir+'/'+'_'.join([Root,END[1],'p','joined'])+'.sam'
    SourceHandle=open(SourceName, 'r')
    SourceLib=csv.reader(SourceHandle, delimiter='\t')

    QueryName=inDir+'/'+'_'.join([Root,END[0],'p','joined'])+'.sam'
    QueryHandle=open(QueryName, 'r')
    QueryLib=csv.reader(QueryHandle, delimiter='\t')



    Row={}
    rNumber={}
    written={}
    skip={}
    sCount=0
#Skip headers
    count={}
    SourceRow=SourceLib.next()
    try:
        while SourceRow[0][0]=='@':
            SourceRow=SourceLib.next()
    except:
            return
    Row=QueryLib.next()
    try:
        while Row[0][0]=='@':
            Row=QueryLib.next()
    except:
            return
    keycount=0
    #Make Output Files
    NameKey={}
    outFiles={}

    os.mkdir(inDir+'/'+Root)
    for k in ChrKey.keys():
        if k=='*':continue
        keycount1='_'.join(k.split('/'))
        keycount='_'.join(keycount1.split('|'))
        NameKey[k]=keycount
        outFiles[k]=open(inDir+'/'+Root+'/'+keycount+'.dat', 'w')
        header=[k, 'strand', 'start', 'end','seq', 'qual', 'Context' ,'strand', 'start', 'end','seq', 'qual']
        dump=csv.writer(outFiles[k], delimiter='\t').writerow(header)

    SourceStop=False
    while SourceStop==False:
        #Calculate Read numbers
        SourceNum=CalculateReadNumber(SourceRow, CountPosition, PosRange)
        rNumber=CalculateReadNumber(Row, CountPosition, PosRange)

        #Look for matches and write output if they are found:
        SourceType=ChrKey[SourceRow[2]]
        w=ChrKey[Row[2]]
        if rNumber==SourceNum and SourceType=='te':
            if w=='te':
                line=GetOutputLine(SourceRow, Row)
                csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line)

                if SourceRow[2]!=Row[2]: #Write to both files if two TEs
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)

                if w=='um':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line)

                if w=='as':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line)


        if rNumber==SourceNum and SourceType=='as':
                if w=='te':
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)

                if w=='um':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line)

                if w=='as':
                    line=GetOutputLine(SourceRow, Row)
                    csv.writer(outFiles[SourceRow[2]], delimiter='\t').writerow(line)

                    if SourceRow[2]!=Row[2]: #Write to both files if two Chr
                        line=GetOutputLine(Row, SourceRow)
                        csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)


        if rNumber==SourceNum and SourceType=='um':
                if w=='te':
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)

                if w=='as':
                    line=GetOutputLine(Row, SourceRow)
                    csv.writer(outFiles[Row[2]], delimiter='\t').writerow(line)

                   #Check which files to iterate
        #First: Is the source less than or equal to all files? Iterate it:

        if SourceNum<rNumber:

            try:
                sCount+=1
                SourceRow=SourceLib.next()
                continue
            except StopIteration:
                SourceStop=True

        if rNumber<=SourceRow:
            try:
                Row=QueryLib.next()
            except StopIteration:
                SourceStop=True

    for k in outFiles.keys():
        outFiles[k].close()
    SourceHandle.close()
    QueryHandle.close()


def pipeline (inDir, outDir, AssemblyIndex, TEIndex, threads, config, cons_file='', length=70, ReadIDPosition=1, PositionRange=0, phred=33, shorten=True, rename=False, Trim=True, Align=True, convTable='', GemCode=False):

    """This functions runs the alignment and post-alignment organization pipeline.

    inDir: The directory containing all of the FASTQ files to be aligned.

    outDir: The output directory (will be created if it does not exist)
    AssemblyIndex: The Bowtie2 index for the masked refernence genome
    TEIndex: The Bowtie2 index of individual insertions + consensus sequence)
    threads: Number of processes to run when parallel processing is allowed
    config: The file containing the paths to the external programs used
    cons_file: The *.fa file containing the set of repeat consensus sequences
    length: The length to which reads should be trimmed from their 3' ends
    ReadIDPosition: Position in read name (parsed by '.') that specifies read numbers
    PositionRange=
    phred: PHRED Quality encoding of the Fastqs (33/64)
    shorten: Whether reads should be trimmed from their 3' ends
    rename: Whether reads should be renamed to ensure the pipeline will work on the datasets
    Trim: Whether Trimmomatic should be used to preprocess data (recommended)
    Align: Whether data should be aligned (Always set to True)
    convTable: Path to the table used to convert alignments from individual insertions to consensus sequences
    GemCode: Whether the data have GemCode barcodes

"""

    Unzipped=False

    BowtiePath, TrimPath, SamtoolsPath=configure(config)

    Index=['as','te']
    IndexFiles={'as':AssemblyIndex, 'te':TEIndex}
    Ends=['1','2']

    if os.path.exists(outDir)==False:
        os.mkdir(outDir)

    RunLog=open('/'.join([outDir, 'runlog.txt']), 'a')

    Roots, Root_Ext = FindInput(inDir)

    RunLog.write('Version: {0}\n\n'.format(version))
    RunLog.write('Parameters: length = {0}, phred = {1}, shorten = {2}, Trim = {3}, Align = {4}\n\n'.format(length, phred, shorten, Trim, Align))
    RunLog.write('{0} < 100 == {1}; {2} >= 100 == {3}'.format(length, length<100, length , length>=100))
    print ('Paired Fastq found for the following datasets:')
    RunLog.write('Paired ends found for the following datasets:\n\n')
    for r in Roots:
        print (r)
        RunLog.write (r + "\n")

    Processed=[]
    Failed=[]

    #Identify already processed output directories in the output path:
    existing_outputs=[]
    for directory, subdir, files in os.walk(outDir):
        existing_outputs.append( directory.split('/')[-1])

    print "Identified existing outputs and removed them from the input list:"
    print "{0} files in input directory.",format( len(Roots))
    for i in existing_outputs[1:]:
        if Roots.count(i)!=0:
            index=Roots.index(i)
            print Roots[index]
            del Roots[index]

    print "{0} files to process.",format( len(Roots))



    for sample in Roots:
        #Trim reads

        print ('\n----------------------\nProcessing {0}...'.format(sample))
        RunLog.write('\n----------------------\nProcessing {0}...\n\n'.format(sample))
        print ('\n----------------------\nTrimming reads in {0}...'.format(sample))
        RunLog.write('\n----------------------\nTrimming reads in {0}...\n\n'.format(sample))


        if Trim==True:
            TrimError=TrimReads(inDir, sample, Root_Ext,  outDir,threads, RunLog,TrimPath, phred)
            if TrimError==True:
                RunLog.write('\nTrim failed! Skipping {0}...\n'.format(sample))
                Failed.append(sample+': Trimmomatic error')
                continue
        if shorten==True:
            ShortenReads(length, outDir, sample, Ends, RunLog, rename=rename, GemCode=GemCode)

        #Align both ends to both indexes
        AlignmentError=False    #Set Error flag as False to begin with
        for e in Ends:
            if AlignmentError==True: break
            for i in Index:

                print ("Aligning end {0} in {1} to {2}...".format(e, sample, i))
                RunLog.write('\n\nAligning end {0} in {1} to {2}...\n\n'.format(e, sample, i))

                if Trim==True: #If Trimmomatic was called, used its output
                    inName='{0}_{1}_p.fq'.format(sample, e)
                    target=outDir+'/'+inName
                    Unzipped=False
                else:
                    #If Trimmomatic was not called, determine the type of input
                    #files. If they are compressed, decompress before running
                    #Bowtie.

                    if Root_Ext[sample+'_'+e].split('.')[-1]=='gz':

                        try:
                            print ('File compressed. Decompressing...')

                            zipped=gzip.open(inDir+'/'+sample+"_"+e+'.'+Root_Ext[sample+'_'+e], 'r')
                            inName=sample+"_"+e+'.fq'
                            unzipped=open(outDir+'/'+inName, 'w')

                            for l in zipped:
                                unzipped.write(l)

                            zipped.close()
                            unzipped.close()
                            target=outDir+'/'+inName
                            Unzipped=True

                        except Exception:     #If there is an error decompressing, skip.
                            print ("Decompression failed!!!")
                            RunLog.write('\nDecompression failed!!! Skipping {0}...\n'.format(sample))
                            AlignmentError=True
                            break

                    else:
                        inName=sample+"_"+e+'.'+Root_Ext[sample+'_'+e]
                        target=inDir+'/'+inName
                        Unzipped=False

                outName='{0}_{1}_p_{2}.sam'.format(sample, e, i)
                SeedLength={'as': 22, 'te': 22}
                if Align==True:
                   AlignmentError=CallAligner(target,outDir+'/'+outName, IndexFiles[i], cutoff[i], SeedLength[i], str(threads), RunLog, phred, BowtiePath)

                if Unzipped==True:
                    unzipped.close()
                    os.remove(target)   #Remove the unzipped copy to save space

                if AlignmentError==True:
                    RunLog.write('\nAlignment failed! Skipping {0}...\n'.format(sample))
                    Failed.append(sample+': Alignment error')
                    break

            #Match up reads in this end:

        #remove fastq files
        if Trim==True:
          os.remove(outDir+'/'+sample+"_1_p.fq")
          os.remove(outDir+'/'+sample+"_2_p.fq")
        if AlignmentError==True: continue


        #Sort the Sam files by read number

        SortError=SortAlignment(outDir, sample, Ends, Index, RunLog, threads, SamtoolsPath)
        if SortError==True:
            RunLog.write('\nSort failed! Skipping {0}...\n'.format(sample))
            Failed.append(sample+': Sort error')
            continue

        #Convert the insertion alignments to consensus alignments
        if convTable!='':
            for e in Ends:
                RunLog.write('\nConverting insertion alignments to consensus alignments')
                alignedTEs='{0}/{1}_{2}_p_{3}_sorted.sam'.format(outDir,sample, e, 'te')
                converted, total, ratio=AlignmentConverter.ConvertSAMFromInsertionsToConsensus(alignedTEs, convTable, cons_file)
                RunLog.write('\{0} out of {1} insertions converted to consensus ({2}%)\n'.format(converted, total, ratio))

##        print jabber
        #Filter output for reads that map to both, as well as reads that map
        #to neither index
##        print jaebberwocke
        for e in Ends:
            FilterError=HandleDualMappings(outDir, sample, 't', e, Index, ReadIDPosition, 0, RunLog)
            if FilterError==True:
                RunLog.write('\nFiltering failed! Skipping {0}...\n'.format(sample))
                Failed.append(sample+': Filtering error')
                break
        if FilterError==True:
            continue
        for n in['1','2']:
            for i in ['as', 'te']:
                os.remove(outDir+'/'+sample+ '_'.join(['',n,'p',i, 'sorted' ])+'.sam')
        OrganizeOutput(outDir,sample, Ends, ReadIDPosition, 0, ['as', 'te'], GemCode)

        Processed.append(sample)
        for n in['1','2']:
            os.remove(outDir+'/'+sample+ '_'.join(['',n,'p', 'joined' ])+'.sam')
    #Note which files were processed and which failed

    print ('Samples successfully processed:')
    RunLog.write('\nSamples successfully processed:\n\n')
    for r in Processed:
        print (r)
        RunLog.write (r + "\n")

    if len(Failed)>0:
        print ('Samples failed to be processed:')
        RunLog.write('\nSamples failed to be processed:\n\n')
        for r in Failed:
            print (r)
            RunLog.write (r + "\n")

    RunLog.close()

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

def main(argv):
    param={}
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if param=={}: return()
    #inDir, outDir, AssemblyIndex, TEIndex, threads, length=70, ReadIDPosition=1, PositionRange=0, phred=33, shorten=True, Trim=True
    inDir=param['-i']
    outDir=param['-o']
    AssemblyIndex=param['-a']
    TEIndex=param['-t']
    threads=param['-p']
    cons_file=param['-cons']
    convTable=''
    if param.has_key('-conv')==True:
        convTable=param['-conv']
    if param.has_key('--phred')==True:
        phred= int(param['--phred'])
    else: phred=33
    if param.has_key('-s')==True:
        shorten=param['-s']
    if param.has_key('-config')==True:
        Parameters(param['-config'])
    if param.has_key('-B')==True:
        Parameters.BowtiePath= param['-B']
    if param.has_key('-T')==True:
        Parameters.TrimPath= param['-T']
    if param.has_key('-S')==True:
        Parameters.SamtoolsPath= param['-S']
    if param.has_key('-L')==True:
        length= int (param['-L'])
    if param.has_key('-R')==True:
        rename= True
    else: rename=False
    GemCode=False
    if param['-G']=='True': GemCode=True
    print Parameters.BowtiePath
    pipeline(inDir, outDir, AssemblyIndex, TEIndex, str(threads), param['-config'], cons_file=cons_file, phred=phred, length=length, convTable=convTable, rename=rename, GemCode=GemCode)

    pass

if __name__ == '__main__':
    main(sys.argv)

