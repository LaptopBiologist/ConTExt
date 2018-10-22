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
##import SeqManipulations
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


def ClusterIndexByLength(index, outDir, BlastDir, cutoff=85):
    #Split file by lengths
    inDir='/'.join(index.split('/')[:-1])
    tempDir=inDir+'/temp'
    SplitFileByLength(index, tempDir)
    file_list=os.listdir(tempDir)
    for f in file_list:
        infile=tempDir+'/'+f
        outfile=tempDir+'/'+f.split('.')[0]+'_out'
        ClusterIndex(infile, outfile, BlastDir, cutoff=80)


def ClusterIndex(index, outDir, BlastDir, cutoff=85, coverage_req=0):
    MakeDir(outDir)
    #Split the index into files containing 300 seq.
    inDir='/'.join(index.split('/')[:-1])
    baseDir='/'.join(outDir.split('/')[:-1])
    tempDir=outDir+'/temp'
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


def SplitFileByLength(infile, tempdir):
    MakeDir(tempdir)
    refSeq=GetSeq(infile, clean_name=False)
    #
    length_dict={}

    for read_name in refSeq.keys():
        length=float( read_name.split('_')[1].split('=')[1])
        if length_dict.has_key(length)==False:
            length_dict[length]={}
        length_dict[length][read_name]=refSeq[read_name]

    inDir='/'.join(infile.split('/')[:-1])
    root=infile.split('/')[-1]
    rootname='.'.join(root.split('.')[:-1])
    ext=root.split('.')[-1]
    for length in length_dict.keys():
        outfile="{0}/{1}_temp_{2}.{3}".format(tempdir, rootname,int( length), ext)
        outHandle=open(outfile, 'w')
        for read_name in length_dict[length].keys():
            outHandle.write('>'+read_name+'\n')
            outHandle.write(str(length_dict[length][read_name].seq)+'\n')

        outHandle.close()

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
    BLAST=subprocess.Popen([BlastDir, '-query',Query, '-subject',Subject, \
    '-outfmt', column_spec,  '-word_size', '11','-dust', 'no', '-soft_masking','false' ,\
      '-out', OutFile], stderr=errlog)  #
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
    BLAST=subprocess.Popen([BlastDir, '-query',Query, '-subject',Subject,
                        '-outfmt', column_spec, '-max_hsps', '1',
                        '-word_size', '11' ,  '-out', OutFile], stderr=errlog)  #'-dust', 'no', '-soft_masking', 'false',
    BLAST.communicate()
    errlog.close()
    return OutFile
def BlastDir(insDir, consDir, outDir):

    insFile=os.listdir(insDir)
    consFile=os.listdir(consDir)
    MakeDir(outDir)
    for f in insFile:
        print "BLASTing {0}...".format(insFile)
        consPath=consDir+'/'+f
        insPath=insDir+'/'+f
        if os.path.exists(consPath)==False: continue
        BlastSeqII(insPath, consPath, outDir,  '.'.join(f.split('.')[:-1]) , 'c:/ncbi-blast')

def BlastDirII(insDir, consPath, outDir):

    insFile=os.listdir(insDir)
    MakeDir(outDir)
    for f in insFile:
        print "BLASTing {0}...".format(f)
        insPath=insDir+'/'+f
        BlastSeqII(insPath, consPath, outDir,  '.'.join(f.split('.')[:-1]) , 'c:/ncbi-blast')


def BlastSeqII(Query, Subject, Out, name, BlastDir):
    """BLASTs two fastas against each other."""
    OutPath='.'.join(Out.split('.'))
    OutFile=OutPath+'/'+name+'.csv'
    print (OutPath)
    errlog=open(OutPath+'/_err.log', 'a')
    column_spec='10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue btop'
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


def ClusterByLengths(indir, outfile):
    file_list=os.listdir(indir)
    outhandle=open(outfile, 'w')
    for f in file_list:
        if f=='singletons.fa':
            sequences=GetSeq(indir+'/'+f, clean_name=False)
            for key in sequences.keys():
                outhandle.write('>{0}\n'.format(key))
                outhandle.write('{0}\n'.format(str(sequences[key])))
            continue
        if f.split('_')[0]!='community':continue
        sequences=GetSeq(indir+'/'+f, clean_name=False)
        hc, clusters=HierarchicalCluster(sequences)
        for key in clusters:
            outhandle.write('>{0}\n'.format(key))
            outhandle.write('{0}\n'.format(str(sequences[key])))
    outhandle.close()

def FilterLowComplexityRepeatsFromFasta(indir, outfile, ignore=''):
    file_list=os.listdir(indir)

    if ignore!='':
        ignore_list=ReadIgnoreFile(ignore)
    else:
        ignore_list=[]
    outhandle=open(outfile, 'w')
    for f in file_list:
        file_root='.'.join( f.split('.')[:-1])
        #Add all singletons to the output
        if f=='singletons.fa':
            sequences=GetSeq(indir+'/'+f, clean_name=False)
            for key in sequences.keys():
                outhandle.write('>{0}\n'.format(key))
                outhandle.write('{0}\n'.format(sequences[key]))
            continue
        if f.split('_')[0]!='community':continue
        sequences=GetSeq(indir+'/'+f, clean_name=False)

        if ignore_list.count(  file_root.split('_')[1])!=0:
            for key in sequences.keys():
                outhandle.write('>{0}\n'.format(key))
                outhandle.write('{0}\n'.format(sequences[key]))
            continue
        clusters=RemoveLowComplexityRepeats(sequences)
        for key in clusters:
            outhandle.write('>{0}\n'.format(key))
            outhandle.write('{0}\n'.format(sequences[key]))
    outhandle.close()

def RemoveLowComplexityRepeats(sequences):
    #Set the threshold equal twice the information of the most complex repeat
    information_threshold=min( CheckSequenceInformation(sequences))*2
    complex_sequences={}
    for key in sequences.keys():
        information=(float( key.split('_')[-1].split('=')[-1]))
        if information<=information_threshold:
            complex_sequences[key]=sequences[key]
    return complex_sequences


def CheckSequenceInformation(sequences):
    #First identify the lowest information sequence:
    information_list=[]
    for key in sequences.keys():
        information_list.append(float( key.split('_')[-1].split('=')[-1]))
    return information_list


def ReadIgnoreFile(infile):
    """Parses a file indicating which communities should be skipped by
    RemoveLowComplexityRepeats.

    The numbers of communities to skip should be separated by commas
    Lines beginning with # will be treated as comments and skipped
    (eg providing justification for skipping these communities).
    """
    ignore_list=[]
    inhandle=open(infile, 'r')
    for line in inhandle:
        #Skip comments
        if line[0]=='#': continue
        ignore_list+=line.split(',')

    inhandle.close()
    return ignore_list
def FilterByLengths(infile, outfile, cutoff=100):
    sequences=GetSeq(infile, clean_name=False)
    outhandle=open(outfile, 'w')
    for key in sequences:
        length=len(sequences[key])
        if length<cutoff: continue
        outhandle.write('>{0}\n'.format(key))
        outhandle.write('{0}\n'.format(str(sequences[key].seq)))
    outhandle.close()
def HierarchicalCluster(sequences):
    """This takes a list of lengths and clusters them as follows:
        For each adjacent pair, it computes the  """
    seq_keys=numpy.array(sequences.keys())
    lengths=numpy.array( [len(sequences[k]) for k in seq_keys], float)


    #Sort the lengths and sequences by length in ascending order
    sort_ind=numpy.argsort(lengths)
    lengths=lengths[sort_ind]
    seq_keys=seq_keys[sort_ind]




    distances=abs(numpy.diff(lengths))
    midpoints=(lengths[1:]+lengths[:-1])/2
    score=distances/(2*midpoints)
    clusters=[[seq_keys[0]]]
    for i,s in enumerate(score):
        if s<.02: clusters[-1].append(seq_keys[i+1])
        else:clusters.append([seq_keys[i+1]])

    retained_keys=[]
    for c in clusters:
        counts=numpy.array( [float( SeqManipulations.ParseSeqName(k)['counts'] ) for k in c])

        best_ind=numpy.argmax(counts)
        retained_keys.append(c[best_ind])

    return clusters, retained_keys

def PlotClusters(clusters):
    count=0
    colors=matplotlib.colors.cnames.keys()
    numpy.random.shuffle(colors)
    for c in clusters:

        counts=numpy.array( [float(k.split('_')[2].split('=')[-1]) for k in c])
        lengths=numpy.array( [float(k.split('_')[1].split('=')[-1]) for k in c])
        t=pyplot.scatter(lengths, numpy.log10( counts),c=colors[count%len(colors)] )

        count+=1
    pyplot.show()

def ComputeKmerCompositionEntropy(sequence,k=5):
    """Summarize the complexity of a repeat unit by decomposing it into kmers
    and then modeling the expected counts of each observed kmer with a multinomial distribution.
    The multinomial probability is multiplied by the number of possible sets of kmers
    the same size as the observed set (the binomial coefficient), because we don't care at all about the
    identity of the kmers. This kept in log-space to avoid precision errors.
    The log-probability is divided by the number of kmers the sequence was decomposed into.

    Note: the binom coefficient becomes imprecise for k>5; 5-mers though provides
    a reasonable summary of sequece complexity."""

    #Replace ambiguous nucleotides with one of the nucleotides they represent
    sequence=SeqManipulations.ReplaceAmbiguousNucleotides(sequence)

    #Decompose the sequence into kmer counts

    kmer_counts=CountKMERS(sequence, k)

    #Number of kmers possible
    num_poss_kmers=(4.**k)

    #Assume each kmer is equally likely
    pr=1./num_poss_kmers

    #The number of kmers contained in the sequence
    num_kmers=sum(kmer_counts.values())

    obs=kmer_counts.values()+[0]* (int(num_poss_kmers-len(kmer_counts.keys())))

    #Multiply the probablity by the binomial coefficient: Don't care which kmers are enriched!
    #In log space, this means add the logs

    try:
        logprob=numpy.log( scipy.special.binom(num_poss_kmers, len(kmer_counts.keys()))) + scipy.stats.multinomial.logpmf(obs, num_kmers, [pr]*int( num_poss_kmers))
    except:
        print kmer_counts
        print len(kmer_counts.values())
        print jabber
    information=0

    if logprob!=0:
        information=-1*logprob
    else:
        information=numpy.inf

    #Output average information per kmer
    return information/ (num_kmers)#len(kmer_counts.keys())#@, len(kmer_counts.keys())/num_kmers


def CountKMERS(sequence, k=10 ):
    """Decomposes sequence into kmers."""
    kmerSet={}
    for x in range(0,len(sequence)-k+1):
        kmer=str(sequence[x:x+k])
        if kmerSet.has_key(kmer)==False:
            kmerSet[kmer]=0.
        kmerSet[kmer]+=1
    return kmerSet
##        kmerSet.add(str(complement[x:x+k]))
    return list(kmerSet)

def SplitFastaByEntropy(infile, outfile):
    low_outhandle=open(outfile+'_low.fa', 'w')
    high_outhandle=open(outfile+'_high.fa', 'w')
    seqs=GetSeq(infile, upper=True)
    for key in seqs.keys():
        entropy=ComputeKmerCompositionEntropy(seqs[key],5)
        if entropy<2:
            low_outhandle.write('>{0}_entropy={1}\n'.format(key, entropy))
            low_outhandle.write('{0}\n'.format(seqs[key]))
        else:
            high_outhandle.write('>{0}_entropy={1}\n'.format(key, entropy))
            high_outhandle.write('{0}\n'.format(seqs[key]))
    high_outhandle.close()
    low_outhandle.close()

def LabelFastaWithEntropy(infile, outfile):
    outhandle=open(outfile, 'w')

    seqs=GetSeq(infile, upper=True)
    for key in seqs.keys():
        entropy=ComputeKmerCompositionEntropy(seqs[key],5)

        outhandle.write('>{0}_information={1}\n'.format(key, entropy))
        outhandle.write('{0}\n'.format(seqs[key]))
    outhandle.close()

def LoadGraph(infile, blastdir,cutoff, cvg):
    outpath='/'.join(infile.split('/')[:-1])
    outroot='.'.join( infile.split('/')[-1].split('.')[:-1])+'_temp.csv'
    BlastSeq_part(infile, infile, outpath, outroot, blastdir)
    seq=ParseSequences(outpath+'/'+outroot)
    graph=BuildGraph(seq, cutoff, cvg)
    return graph


def MaskRepeats(infile, maskfile, blastdir, cutoff=85):
    outfile=''.join(infile.split('.')[:-1])+'_aligned.csv'

    BlastSeq(infile, maskfile, outfile, blastdir)

    aligned_handle=open(outfile, 'r')
    aligned_table=csv.reader(aligned_handle, delimiter=',')

    sequences=GetSeq(infile)
    #Make the sequences numpyarrays for easy manipulation
    for key in sequences.keys():
        sequences[key]=numpy.fromstring(str( sequences[key].seq), '|S1')
    for row in aligned_table:
        line=alignment(row)
        if line.identity<cutoff: continue
        aligned_len=abs( line.qend-line.qstart)
        if aligned_len<50: continue
        sequences[line.query][line.qstart:line.qend]='N'
    aligned_handle.close()
    outfile=''.join(infile.split('.')[:-1])+'_masked.fa'
    outhandle=open(outfile, 'w')
    for key in sorted(sequences.keys()):
        outhandle.write('>{0}\n'.format(key))
        outhandle.write('{0}\n'.format(''.join(sequences[ key])))
    outhandle.close()


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

