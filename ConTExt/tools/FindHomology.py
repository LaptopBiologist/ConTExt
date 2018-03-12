#-------------------------------------------------------------------------------
# Name:        InternalRptFinder
# Purpose:
#
# Author:      Michael mcGurk
#
# Created:     20/05/2014
# Copyright:   (c) Michael McGurk 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import subprocess
import shutil
import os
import sys
import csv
import networkx
import community
import numpy

import seaborn
import matplotlib
from matplotlib import pyplot

import Bio
from Bio import SearchIO
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
import scipy
import scipy.stats


class alignment():
    def __init__(self, row):
        self.query=row[0]
        self.subject=row[1]
        self.identity=float(row[2])
        self.length=int(row[3])
        self.covs=float(row[-2])
        self.qstart=int(row[6])
        self.qend=int(row[7])
        self.pair=tuple(sorted((self.query, self.subject)))

def MakeDir(newdir):
    if os.path.exists(newdir)==False:
        os.mkdir(newdir)



def GetSeq(ref, upper=False, rename=False, clean_name=True, string=False):

    """Reads a fasta, returns of a dictionary of strings keyed by entry name."""
    if ref.split('.')[-1]=='gz':
        handle=gzip.open(ref )
    else:
        handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    SeqLen={}
    for rec in lib:
##        if refName.count(rec.name)==0: continue
        if string==True: sequence=str( rec.seq)
        else: sequence=rec
        if rename==True:
            name=rec.description.split(' ')[1].split('=')[1]
            print name
            SeqLen[CleanName(name)]=sequence
        else:
            if clean_name==True:
                SeqLen[CleanName(rec.name)]=sequence
            else:
                SeqLen[rec.name]=sequence
        if upper==True: SeqLen[CleanName(rec.name)]=SeqLen[CleanName(rec.name)].upper()
    handle.close()
    return SeqLen


def CleanName(name):
    illegal=['|', '!','>', '<', '?','/','*']
    for i in illegal:
        name='_'.join(name.split(i))
    return(name)


def ClusterIndex(index, outDir, BlastDir, cutoff=85, coverage_req=0):

    #Split the index into files containing 300 seq.
    inDir='/'.join(index.split('/')[:-1])
    tempDir=inDir+'/temp'
    splitFile(index, tempDir, 300)

    outfiles=[]

    fileList=os.listdir(tempDir)

    for j in range(len(fileList)):
        if coverage_req==0:
            I=j+1
        else:
            I=len(fileList)
        for i in range(I):
            print fileList[j], fileList[i]
            outname='alignments_{0}_{1}.tsv'.format(str(i), str(j))
            outfiles.append( BlastSeq_part(tempDir+'/'+ fileList[j],tempDir+'/'+ fileList[i], outDir, outname, BlastDir))
    AlignmentFile=outDir+'/alignment_output.tsv'
    JoinFiles(outfiles, outDir, AlignmentFile)


    aligned=ParseSequences(AlignmentFile)
    graph=BuildGraph(aligned, cutoff, coverage_req)
    clusters=GraphToCommunities(graph)
    groups, singletons=RemoveSingletons(clusters)


    WriteOutput(groups,singletons, outDir, index)
    shutil.rmtree(tempDir)
    return groups


def JoinFiles(infiles, inDir, outfile):
    outhandle=open(outfile, 'w')
    for f in infiles:
        inhandle=open(f,'r')

        for line in inhandle:
            outhandle.write(line)
        inhandle.close()
    outhandle.close()


def splitFile(infile, tempDir, limit=500):
    MakeDir(tempDir)
    refSeq=GetSeq(infile)
    inDir='/'.join(infile.split('/')[:-1])
    root=infile.split('/')[-1]
    rootname='.'.join(root.split('.')[:-1])
    ext=root.split('.')[-1]
    fileCount=1
    outfile="{0}/{1}_temp_{2}.{3}".format(tempDir, rootname, fileCount, ext)
    outHandle=open(outfile, 'w')
    count=0
    print len(refSeq)
    for key in refSeq:
        if count>=limit:
            print count
            outHandle.close()
            fileCount+=1
            count=0
            outfile="{0}/{1}_temp_{2}.{3}".format(tempDir, rootname, fileCount, ext)
            print outfile
            outHandle=open(outfile, 'w')
        seq=str( refSeq[key].seq)
        name= CleanName( str(refSeq[key].name))

        outHandle.write('>'+name+'\n')
        outHandle.write(seq+'\n')
        count+=1

    outHandle.close()

def BlastSeq(Query, Subject, Out, BlastDir):
    """BLASTs two fastas against each other."""
    print Out
    print Out.split('.')
    if len(Out.split('.'))==1:
        MakeDir(Out)
        OutPath='.'.join(Out.split('.'))
        print (OutPath)
        OutFile=OutPath+'/output.csv'
        errlog=open(OutPath+'/_err.log', 'a')
    else:
        OutFile=Out
        errfile='.'.join( Out.split('.')[:1])+'_err.log'
        errlog=open(errfile, 'a')


##    column_spec='10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue btop'
    column_spec='10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qcovs btop'
    BLAST=subprocess.Popen([BlastDir, '-query',Query, '-subject',Subject, '-outfmt', column_spec,  '-out', OutFile], stderr=errlog)
    BLAST.communicate()
    errlog.close()
    return OutFile

def BlastSeq_part(Query, Subject, OutPath, outname, BlastDir):
    """BLASTs two fastas against each other."""
    MakeDir(OutPath)
    OutFile=OutPath+'/'+outname
    print (OutPath)
    errlog=open(OutPath+'/_err.log', 'a')
    column_spec='10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue qcovs btop'
    BLAST=subprocess.Popen([BlastDir, '-query',Query, '-subject',Subject, '-outfmt', column_spec,  '-out', OutFile], stderr=errlog)
    BLAST.communicate()
    errlog.close()
    return OutFile


def ParseSequences(InFile):
    """Organizes the BLAST output and stores the information in python's working
     memory."""
    handle=open(InFile)
    table=csv.reader(handle )
    pairDict={}
    for row in table:
        aligned=alignment(row)
        if pairDict.has_key(aligned.pair)==False:
            pairDict[aligned.pair]=aligned
        if aligned.length>pairDict[aligned.pair]:pairDict[ aligned.pair]=aligned
    return(pairDict)


def BuildGraph(aligned, cutoff, coverage_req=0):
    """Takes the BLAST output and represents it as a graph using the following rules:
        Each sequence is represented as a node.
        Each alignment between sequences is represented as an unweighted,
        undirected edge."""
    #hold=set()
    #hold2=set()
    graph=networkx.Graph()
    for pair in aligned.keys():
        query,subject=pair
        if graph.has_node(subject)==False:
            graph.add_node(subject)
        if graph.has_node(query)==False:
            graph.add_node(query)
        if aligned[pair].identity>=cutoff and aligned[pair].covs>coverage_req and aligned[pair].length>50: #aligned[pair].length>100:
            graph.add_edge(subject, query)

    return graph

def GraphToCommunities(Network):
    """Identifies communities in the graph using the Louvain method as implemented in community."""
    comm=community.best_partition(Network)
    clusters={}
    for k in comm.keys():
        if clusters.has_key(comm[k])==False:
            clusters[comm[k]]=[]
        clusters[comm[k]].append(k)
    return (clusters)

def RemoveSingletons(clusters):
    """Removes clusters containing a single node. That is, entries that align
    only to themselves."""
    groups=[]
    singletons=[]
    for key in clusters.keys():
        if len(clusters[key])>1: groups.append(clusters[key])
        else: singletons+=clusters[key]
    return groups, singletons

def WriteOutput(clusters, singletons, outDir, index):
    MakeDir(outDir)
    inhandle=open(index,'r')
    lib=SeqIO.parse(inhandle, 'fasta')
    hold={}

    for book in lib:
        name=book.id.strip()
        name=name.split('\t')[0]
##        if name[-1]=='.': name=name[:-1]
        hold[CleanName( name)]=book

    for h in hold.keys():
        print h

    outfile=outDir+'/communities.txt'
    outHandle=open(outfile, 'w')
    count=0

    for community in clusters:
##        if len(community)==1: continue
        dump=outHandle.write('Community {0}, Members = {1}:\n\n'.format(count, len(community)))
        outFasta=outDir+'/community_{0}.fa'.format(count)
        FASTAHandle=open(outFasta, 'w')
        sequences=[]
        for member in community:
            sequences.append(hold[member])
            dump=outHandle.write('\t{0}, '.format(member))
        dump=outHandle.write('\n\n')
        dump=SeqIO.write(sequences, FASTAHandle, 'fasta')
        FASTAHandle.close()
        count+=1

    #Write singletons deal with later
    dump=outHandle.write('Singletons, Members = {0}:\n\n'.format( len(singletons)))
    outFasta=outDir+'/singletons.fa'.format(count)
    FASTAHandle=open(outFasta, 'w')
    sequences=[]
    for member in singletons:
        sequences.append(hold[member])
        dump=outHandle.write('\t{0}, '.format(member))
    dump=outHandle.write('\n\n')
    dump=SeqIO.write(sequences, FASTAHandle, 'fasta')
    FASTAHandle.close()

    outHandle.close()
    inhandle.close()
    return(hold)


def main(argv):
    param={}
    print argv
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if param=={}: return()
    #inDir, outDir, AssemblyIndex, TEIndex, threads, length=70, ReadIDPosition=1, PositionRange=0, phred=33, shorten=True, Trim=True
    if param.has_key('-mask')==True:

        MaskRepeats(param['-i'], param['-mask'], param['-B'], float( param['-C']))
        return()

    inDir=param['-i']
    outDir=param['-o']
    blast_path=param['-B']
    if param.has_key('-C')==True:
        cutoff=float( param['-C'])
    else:
        cutoff=80
    if param.has_key('-cvg')==True:
        cvg_req=float( param['-cvg'])
    else:
        cvg_req=0
    ClusterIndex(inDir, outDir, blast_path, cutoff, cvg_req)
    pass

if __name__ == '__main__':
    main(sys.argv)

