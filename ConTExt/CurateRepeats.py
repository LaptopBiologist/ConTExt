#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      I am
#
# Created:     16/07/2018
# Copyright:   (c) I am 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import matplotlib

import sys
if sys.platform[0]=='l':
    matplotlib.use('Agg')

from matplotlib import pyplot
pyplot.ioff()

import numpy
##import pandas
import subprocess
import scipy
from Bio import SeqIO
import scipy.sparse

import seaborn
import csv
import pickle
import shutil
import os
from collections import Iterable
try:
    from AlignmentPipeline import CallAligner
    from tools import FindHomology
##try:
    from tools import General_Tools

except:
    print "Couldn't import General_Tools"
import json
import scipy.sparse
from Bio import SeqIO
from Bio.Seq import Seq
import Bio.Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import sklearn.neighbors
from matplotlib import pyplot
import seaborn
from itertools import product
import dict_trie
from collections import Iterable
from collections import Counter
import pickle
import shutil
import os
import time


class SamLine():

    def __init__(self, row):
        self.id=row[0]
        self.flag=int(row[1])
        self.contig=row[2]
        self.pos=int(row[3])
        self.MapQ=float(row[4])
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

#------------------------------------------------------------------------------
#           Full automatic curation pipeline
#------------------------------------------------------------------------------

def AutomaticCuration(fasta_file, fastq_file, outdir, specification_file):


    spec_dict=General_Tools.ReadSpecificationFile(specification_file)

    #Want to organize the output in a directory
    General_Tools.MakeDir(outdir)
    #Increase the lengths of the sequences too short for reads to align by
    #concatenating monomers, and then trimming the ends until reads should align
    #unambiguously
    fasta_root=fasta_file.split('/')[-1]
    modified_fasta="{0}/{1}".format(outdir, fasta_root)
##    LengthenShortRepeats(fasta_file, modified_fasta)
    fasta_file=modified_fasta
    #Run FindHomology
##    FindHomology.ClusterIndex(fasta_file, outdir+'/00_Communities',spec_dict['BLAST'], 80)
    #Phase communities
##    PhaseCommunities(outdir+'/00_Communities', True)
    #Build Bowtie2 index
    phased_path="{0}/00_Communities_phased/phased_sequences.fa".format(outdir)
##    BuildIndex(spec_dict['Bowtie2'],phased_path, outdir+"/phased_index" )
    #Run Bowtie2

    aligned_dir=outdir+'/01_aligned'
    print aligned_dir
    sam_file="{0}/aligned.sam".format(aligned_dir)
    General_Tools.MakeDir(aligned_dir)
    err_file=aligned_dir+'/errlog.txt'
##    CallAligner(fastq_file,sam_file,outdir+"/phased_index", spec_dict['RepeatAlignmentParameters']+' --no-unal' ,str(int( spec_dict['threads'])),None, bowDir=spec_dict['Bowtie2'])
##    err_handle.close()
    #S
    FASTA_file="{0}/aligned.fasta".format(aligned_dir)
##    Sam2FASTA(sam_file, FASTA_file, 10,50, False, phased_path)
    #Run Blastn
    blast_file="{0}/aligned.csv".format(aligned_dir)
##    RunBLAST(FASTA_file, phased_path, blast_file)
    #Build BLAST matrix
    matrix_file="{0}/aligned_matrix".format(aligned_dir)
##    BlastToScoreMatrix(blast_file, matrix_file)
    matrix_file="{0}/aligned_matrix.npz".format(aligned_dir)
    header_file="{0}/aligned_matrix_headers.pickle".format(aligned_dir)
    comm_file="{0}/00_Communities_phased/communities.txt".format(outdir)
    cur_file="{0}/02_curated".format(outdir)
    #Curate the BLAST matrix
##    CurateRepeats(fasta_file,comm_file,  matrix_file, header_file,cur_file)
    #Return summaries17
    curated_fasta="{0}/02_curated/curated.fa".format(outdir)
    curated_index= outdir+"/curated_index"
##    BuildIndex(spec_dict['Bowtie2'],curated_fasta, curated_index )
    curated_sam_file="{0}/curated.sam".format(aligned_dir)
##    CallAligner(fastq_file,curated_sam_file,curated_index, spec_dict['RepeatAlignmentParameters']+' --no-unal' ,str(int( spec_dict['threads'])),None, bowDir=spec_dict['Bowtie2'])
    summary_file=outdir+"/summary.tsv"
    UpdateMapq(curated_sam_file)
    SummarizeSam(curated_sam_file,summary_file, comm_file,fasta_file)
def BuildIndex(bowtie_path, inpath, outpath):
    process=subprocess.Popen([bowtie_path+"-build", inpath, outpath])
    process.communicate()

def RunBLAST(query, subject,outfile):
    parameters="blastn -dust no -soft_masking false -max_hsps 1 -word_size 11 -evalue .001 -perc_identity 80 -task blastn-short -outfmt".split(' ')
    parameters+=['10 qseqid sseqid qstart qend pident evalue nident mismatch gaps gapopen bitscore']
    parameters+=["-out", '{0}'.format(outfile),"-query", '{0}'.format(query), "-subject", '{0}'.format(subject)]
    print parameters
    proc=subprocess.Popen(parameters)
##    print ("Aligning {0} to {1}...".format(inFile, Index ))

    proc.communicate()


#------------------------------------------------------------------------------
#           Functions for preparing and building the alignment index
#------------------------------------------------------------------------------


def PhaseCommunities(indir, filter_by_internal_rpt=True):
    in_root="/".join( indir.split('/')[:-1])
    phased_folder=indir+'_phased'
    General_Tools.MakeDir(phased_folder)
    file_list=os.listdir(indir)

    #Phase each community
    for f in file_list:
        ext=f.split('.')[-1]
        if f=='communities.txt':
            shutil.copyfile(indir+'/'+f,phased_folder+'/'+f)
        if ext!='fa': continue
        if f=='singletons.fa': continue
        PhaseFASTA(indir+'/'+f, phased_folder+'/'+f)
    phased_files=os.listdir(phased_folder)
    discarded_fasta=phased_folder+'/discarded_sequences.fa'
    discarded_record=phased_folder+'/discarded_sequences.txt'
    discarded_fasta_handle=open(discarded_fasta, 'w')
    discarded_record_handle=open(discarded_record, 'w')
    #Join the phased communities into a single fasta
    phased_fasta=phased_folder+'/phased_sequences.fa'
    phased_handle=open(phased_fasta, 'w')
    for f in sorted(phased_files):
        ext=f.split('.')[-1]
        if ext!='fa': continue
        if f=='singletons.fa': continue
        if f=='phased_sequences.fa': continue
        sequences=GetSeq(phased_folder+'/'+f, True)
        rpt_scores=[ComputeNeighborDistance(sequences[ key]) for key in sequences]
        med_score=numpy.nanmedian(rpt_scores)
        print f
        print "\t", med_score*1.25
        for i,seq in enumerate( sequences):

            if filter_by_internal_rpt==True and rpt_scores[i]>=1.25*med_score and med_score<=2 and len(seq.split("_"))>3:
                print "Skipped: ", rpt_scores[i],seq
                discarded_fasta_handle.write('>{0}\n'.format( seq))
                discarded_fasta_handle.write('{0}\n'.format(sequences[ seq]))
                discarded_record_handle.write('>{0}\t{1}\n'.format( seq,rpt_scores[i]))
                continue
            phased_handle.write('>{0}\n'.format( seq))
            phased_handle.write('{0}\n'.format(sequences[ seq]))

    phased_handle.close()
    discarded_fasta_handle.close()
    discarded_record_handle.close()



def CountKMERS(sequence, k=10, unbiased=False, skip_Ns=True ):
    """Decomposes sequence into kmers."""
    kmerSet={}
    seq_len=len(sequence)
    if unbiased==False:
        indices= range(0,len(sequence)-k+1)
    else:
        sequence=sequence*2
        indices= range(0, seq_len)
    for x in indices:
        kmer=str(sequence[x:x+k])
        if kmer.lower().count('n')>0:continue
        if kmerSet.has_key(kmer)==False:
            kmerSet[kmer]=0.
        kmerSet[kmer]+=1.
    return kmerSet
##        kmerSet.add(str(complement[x:x+k]))
    return list(kmerSet)

def SequenceToInt(seq):
    seq_array=numpy.fromstring(seq, '|S1')
    int_array=numpy.ndarray((len(seq_array),))
    int_array.fill(0)
    nt_array=['T', 'C', 'G']
    for i, nt in enumerate( nt_array):
        int_array[seq_array==nt]=i+1
    return int_array

def ComputeNeighborDistance(sequence,k=20, mismatches=.15, max_sample=400, unbiased=False):
    kmer_dict=CountKMERS(sequence,k, unbiased)
    kmer_list, kmer_counts=zip(* kmer_dict.items())

    radius=mismatches
    start=time.clock()
    encoded_kmers=numpy.array( [SequenceToInt(kmer) for kmer in kmer_list ])
    ball_tree=sklearn.neighbors.BallTree(encoded_kmers, metric='hamming')

    count_list=[]
    start=time.clock()
    indices=numpy.arange(len(encoded_kmers))
    if len(encoded_kmers)>max_sample:
        pr=numpy.array(kmer_counts, float)
        pr/=pr.sum()

        indices=numpy.random.choice(numpy.arange(len(pr)), size=max_sample, p=pr)
    for kmer in encoded_kmers[indices]:
        neighbors=ball_tree.query_radius(kmer[None,:], radius)[0]

        counts=numpy.array( [kmer_dict[kmer_list[n] ] for n in neighbors])

        distance_mat=(kmer[None,:] !=encoded_kmers[neighbors,:])

        distance=(1.-distance_mat.sum(1)/k )
        count_list.append(numpy.sum( distance*counts) )

    S= numpy.mean( count_list)

    return S


def ReplaceAmbiguousNucleotides(sequence):
    """Replaces ambiguous nucleotides with one of the nt that could be represented."""
    seq_array=numpy.fromstring(sequence.upper(), '|S1')
    for ambiguous_nt in IUPAC_DICT:
        amb_ind=numpy.where(seq_array==ambiguous_nt)[0]
        replacement_nt=numpy.random.choice(IUPAC_DICT[ambiguous_nt], size=len(amb_ind), replace=True)
##        print replacement_nt
        seq_array[amb_ind]=replacement_nt
    #Replace all remaining non ATGC nucleotides randomly:
    ATGC=(seq_array=='A')+(seq_array=='T')+(seq_array=='C')+(seq_array=='G')

    amb_ind=numpy.where(ATGC==0)[0]
##    print amb_ind
    replacement_nt=numpy.random.choice(IUPAC_DICT['N'], size=len(amb_ind), replace=True)
    seq_array[amb_ind]=replacement_nt
    return ''.join(seq_array)

def GetSeq(ref, string=True):
    handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    SeqLen={}
    for rec in lib:
##        if refName.count(rec.name)==0: continue
        if string==True:
            SeqLen[CleanName(rec.name)]=str( rec.seq).upper()
        else:
            SeqLen[CleanName(rec.name)]=rec.seq
    handle.close()
    return SeqLen

def GetLengths(ref):
    handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    SeqLen={}
    for rec in lib:
        SeqLen[CleanName(rec.name)] = len(rec.seq)
    handle.close()
    return SeqLen

def CleanName(name):
    illegal=['|', '!','>', '<', '?','/','*']
    for i in illegal:
        name='_'.join(name.split(i))
    return(name)

def UnambiguousFASTA(infile, outfile):
    seq=GetSeq(infile)
    outhandle=open(outfile, 'w')
    for key in seq:
        new_seq=ReplaceAmbiguousNucleotides(seq[key])
        outhandle.write('>{0}\n'.format(key))
        outhandle.write('{0}\n'.format(new_seq))
    outhandle.close()


def PhaseFASTA(infile, outfile):
    """Takes a FASTA and brings all sequences in phase with each other"""
    seqs=GetSeq(infile)
    outhandle=open(outfile, 'w')
    query_seq=seqs.values()[0]
    for i,key in enumerate( seqs.keys()):
        phased_seq, query_seq=QuickMatch(seqs[key],query_seq)
        outhandle.write('>{0}\n'.format(key))
        outhandle.write('{0}\n'.format(phased_seq))
    outhandle.close()

def QuickMatch(query, subject, seed_count=100, kmer_size=8):

    """Uses a heuristic to match the phase between two sequences."""

    #Get the reverse complement of the query
    query_rc=str( Bio.Seq.Seq(query).reverse_complement())

    #Decompose the sequences into sets of 10-mers
    kmer_set_query=set(GetKMERS(query,k=kmer_size))
    kmer_set_rc=set(GetKMERS(query_rc, k=kmer_size))
    kmer_set_subject=set(GetKMERS(subject, k=kmer_size))


    #Choose a random set of 10-mers from the subject sequence
##    seeds=numpy.random.choice(list(kmer_set_subject), size=seed_count, replace=True)
    seeds=kmer_set_subject

    #Get the intersections of the subject and query 10-mers from both the forward
    #and reverse strands
    seed_intersection=list( set(seeds)&kmer_set_query)
    seed_intersection_rc=list( set(seeds)&kmer_set_rc)
    if len(seed_intersection)<len(seed_intersection_rc):
        query=query_rc

    displacements=[]
    displacements_rc=[]

    #Check the displacement between the first occurence of each kmer in the subject
    #versus the query sequence
    for s in seed_intersection:
        distance=query.find(s)-subject.find(s)
        if distance<0:
            distance=len(subject)+distance
        displacements.append(distance)
##    print displacements
##    print displacements.count(149)
##    print displacements.count(0)
    #Check the displacements in the reverse complement
##    for s in seed_intersection_rc:
##        distance=query_rc.find(s)-subject.find(s)
##        if distance<0:
##            distance=len(query_rc)+distance
##        displacements_rc.append(distance)
##    print sorted( displacements)
    if len (displacements)>0:
        mode=scipy.stats.mode(displacements)
        print mode
        mode, count=mode.mode[0], mode.count[0]

    else: mode,count=0,0

##    if len(displacements_rc)>0:
##        mode_rc=scipy.stats.mode(displacements_rc)
##        mode_rc, count_rc=mode_rc.mode[0], mode_rc.count[0]
##    else:
##        mode_rc,count_rc=0,0

    seed_count=len(set(seeds))
    for_prop=float(count)/seed_count
##    rev_prop=float(count_rc)/seed_count
##    print len(seed_intersection)>len(seed_intersection_rc)
##    print len(seed_intersection),len(seed_intersection_rc)
##    if len(seed_intersection)>len(seed_intersection_rc):
##    print mode
    seq1_array, seq2_array=numpy.fromstring( query[mode:]+query[:mode], '|S1'),numpy.fromstring( subject, '|S1')
##    else:
##        seq1_array, seq2_array=numpy.fromstring( query_rc[mode_rc:]+query_rc[:mode_rc], '|S1'),numpy.fromstring( subject, '|S1')

    return ''.join(seq1_array),''.join( seq2_array)# , numpy.mean( seq1_array== seq2_array)

def GetKMERS(sequence, k=10, unbiased =True):
##    complement=str(Seq.Seq( sequence).reverse_complement())
    seq_len=len(sequence)
    if unbiased==False:
        indices= range(0,len(sequence)-k+1)
    else:
        sequence=sequence*2
        indices= range(0, seq_len)
    kmerSet=set()
    for x in indices:
        if sequence[x:x+k]=='N'*k: continue
        kmerSet.add(str(sequence[x:x+k]))
##        kmerSet.add(str(complement[x:x+k]))
    return list(kmerSet)

#------------------------------------------------------------------------------
#            Functions for curating repeats based on read alignments
#------------------------------------------------------------------------------

class SparseArray():

    def __init__(self):
        self.row_indices=[]
        self.col_indices=[]
        self.data=[]
        self.row_count=0
        self.col_count=0
        self.shape=(self.row_count,self.col_count)

    def add_row(self, row):
        nonzero_ind=numpy.where(row!=0)[0]
        self.col_indices+=list(nonzero_ind)
        self.row_indices+=[row_count]
        if isinstance(row, numpy.array) is True:
            self.data+=row[nonzero_ind]
        else:
            self.data+=list(numpy.array(row)[nonzero_ind])
        self.col_count=max(col_count, numpy.max(nonzero_ind))
        self.row_count+=1

    def add_col(self, col):
        nonzero_ind=numpy.where(col!=0)[0]
        self.row_indices+=list(nonzero_ind)
        self.col_indices+=[col_count]
        if isinstance(col, numpy.array) is True:
            self.data+=col[nonzero_ind]
        else:
            self.data+=list(numpy.array(col)[nonzero_ind])
        self.row_count=max(row_count, numpy.max(nonzero_ind))
        self.col_count+=1

    def add_cell(self, data,row, col):

        if isinstance(data, Iterable) is True:
            self.data+=list(data)
            self.row_indices+=list(row)
            self.col_indices+=list(col)
        if isinstance(data, Iterable) is False:
            if data!=0:
                self.data.append(data)
                self.row_indices.append(row)
                self.col_indices.append(col)

    def construct_coo(self):
        self.matrix=scipy.sparse.coo_matrix((self.data,(self.row_indices, self.col_indices)))
    def construct_csr(self):
        self.matrix=scipy.sparse.csr_matrix((self.data,(self.row_indices, self.col_indices)))
    def construct_csc(self):
        self.matrix=scipy.sparse.csc_matrix((self.data,(self.row_indices, self.col_indices)))


def GetSeq(ref, string=True):
    handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    SeqLen={}
    for rec in lib:
##        if refName.count(rec.name)==0: continue
        if string==True:
            SeqLen[CleanName( rec.name)]=str( rec.seq).upper()
        else:
            SeqLen[CleanName( rec.name)]=rec
    handle.close()
    return SeqLen

def MatrixMapq(matrix, full=False, passed=True):

    if matrix.shape[1]>2:
        AS=numpy.partition(matrix,-1,1)[:,-1]
        XS=numpy.partition(matrix,-2,1)[:,-2]
        extra_counts=( matrix==AS[:,None]).sum(1)
##        extra_counts=1
        max_score=-6.*140/5
        score=MapConf(AS,XS,70, passed)/extra_counts
    elif matrix.shape[1]==2:
        AS=numpy.partition(matrix,-1,1)[:,-1]
        XS=numpy.partition(matrix,-2,1)[:,-2]
        max_score=-6.*140/5
        score=MapConf(AS,XS,70, passed)
##        score[((XS==0)*(AS!=0)).astype(bool)]=40

    else:
        score=40*(matrix>0)
    score=1-10**(score/-10.)
    return score
##    if matrix.shape[1]>1:
##        if full==False:
##            score=numpy.array(matrix.max(1)/(matrix.sum(1)))
##        else:
##            score=numpy.array(matrix/(matrix.sum(1)[:,None]))
##        score[numpy.isnan(score)]=0.
##        one_ind= score==1
##        score=-10*numpy.log10(1-score)*.33
##        score/=52.61444578/40.
##        score[one_ind]=40
##        score[score>40]=40
##        score= 1- numpy.exp(-1.* score/10.)
##
##        score[numpy.isnan(score)]=0.
##    else:
##        score=(matrix>0)*.9999
####    print score
##    return score

def PrioritizeRepeats(mat, uniquely_mapped):
    if mat.shape[1]==1:
        return numpy.array( [0])
    best_hit=numpy.array( numpy.argmax(mat[:,:],1))
    AS=numpy.max(mat[:,:],1)
    score=numpy.array( MatrixMapq(mat[:,:]))
##    scores=numpy.array( CountReads(mat,0))
##    enrichment=scores+uniquely_mapped
##    max_score=numpy.max(score)
##    print max_score
##    print score
##    counts=Counter(best_hit[0])
    enrichment=[]
    for k in range( mat.shape[1]):
        good_ind= (score>=.9)
        bad_ind=(score<.9)*(AS>70)
        good_count=(best_hit[good_ind]==k).sum()+uniquely_mapped[k]
        bad_count=(best_hit[bad_ind]==k).sum()
        ratio=bad_count/(good_count+1.)
        enrichment.append(ratio)
#     print enrichment
    ranked_indices=numpy.argsort(enrichment)[::-1]
    return ranked_indices

def PlotMatrix(mat):
    noncluster_ind=numpy.array( (mat[:,:].sum(1)==0).T)[0]
    print noncluster_ind
    mat=mat[~noncluster_ind,:]
    mat/=mat.sum(1)
    print mat.shape
    try:
        seaborn.heatmap(mat.todense())
    except: seaborn.heatmap(mat)
    pyplot.show()


##def RemoveRepeats(mat, iterations=312):
##    print mat.shape
##    unique_reads=numpy.array( mat[0,:].todense())[0]
##    bool_array=numpy.array( [True]*mat.shape[1])
##    noncluster_ind=numpy.array( (mat[:,:].sum(1)==0).T)[0]
##    print noncluster_ind.shape
##
##    mat=mat[~noncluster_ind,:]
##    mat=numpy.array(mat.todense())
##    if mat.shape[0]>5000:
##        mat=mat[::5,:]
####    print mat
##    print mat.shape
##    height,width=mat.shape
###     if (height*width)**.5<2200:
###         mat=mat.todense()
##    scores=ProbToScore(mat[1:,:])
##    print scores.shape
##    checked_list=set([])
##    score_list=[]
##    removed_list=[]
##    ranked_indices= PrioritizeRepeats(mat[:,bool_array], unique_reads[bool_array])
##    for i in range(min( iterations, len(bool_array))):
##        max_score=0
##        max_j=0
####        ranked_indices= PrioritizeRepeats(mat[:,bool_array], unique_reads[bool_array])
####        ranked_indices=numpy.where(bool_array==True)[0][ranked_indices]
##        counter=0
##        improved=False
##        candidates=[]
##        candidate_scores=[]
##        for j in ranked_indices:
##            if bool_array[j]==False: continue
##            orig_val=numpy.copy( bool_array[j])
##            original_score=sum(CountReads(mat[:,bool_array],scores[:,bool_array]))+unique_reads[bool_array].sum() #+ mat[0,bool_array].sum()
##            if score_list==[]: score_list.append(original_score)
##            print sum(CountReads(mat[:,bool_array],scores[:,bool_array])), unique_reads[bool_array].sum()
##            bool_array[j]=False
##            if sum(bool_array)!=0:
##
##                new_score=sum(CountReads(mat[:,bool_array],scores[:,bool_array]))+unique_reads[bool_array].sum()
##            else:
##                new_score=0
##
###             print new_score
##            print "\t", sum(CountReads(mat[:,bool_array],scores[:,bool_array])), unique_reads[bool_array].sum()
##
##            if new_score>=max_score:
##
##                max_score=new_score
##                max_ind=j
##                improved=True
####            if new_score>=original_score and counter>10:
####                terminate=True
####                bool_array[j]=True
####                if improved==False:
####                    max_score=new_score
####                    max_ind=j
####                break
##            counter+=1
##            bool_array[j]=True
##
##        print i, max_score, original_score
##        removed_list.append(max_ind)
##        bool_array[max_ind]=False
##        if sum(bool_array)==0: break
####        seaborn.heatmap (mat[1:,bool_array])
####        print mat[::9,bool_array]
####        pyplot.savefig("/data/mpm289/Dmel_reannotation/{0}.png".format(i))
####        pyplot.show()
####        pyplot.close()
##        score_list.append(max_score)
##    return removed_list, score_list

def RemoveRepeats(mat,cutoff=.85, iterations=312):
    print mat.shape
    unique_reads=numpy.array( mat[0,:].todense())[0]
    bool_array=numpy.array( [True]*mat.shape[1])
    noncluster_ind=numpy.array( (mat[:,:].sum(1)==0).T)[0]
    print noncluster_ind.shape

    mat=mat[~noncluster_ind,:]
    mat=numpy.array(mat.todense())
    best_score=numpy.max(Score2PercIdent(mat),1 )
    mat=mat[best_score>=cutoff]
##    print best_score.shape
##    print mat.shape
##    print jabber
    if mat.shape[0]>5000:
        step=int (mat.shape[0]/10000.)+1
        step=5
        mat=mat[::step,:]

##    print mat
    print mat.shape
    height,width=mat.shape
#     if (height*width)**.5<2200:
#         mat=mat.todense()
    scores=ProbToScore(mat[1:,:])
    print scores.shape
    checked_list=set([])
    score_list=[]
    removed_list=[]
    ranked_indices= PrioritizeRepeats(mat[:,bool_array], unique_reads[bool_array])
    for i in range(min( iterations, len(bool_array))):
        max_score=0
        max_j=0
##        ranked_indices= PrioritizeRepeats(mat[:,bool_array], unique_reads[bool_array])
##        ranked_indices=numpy.where(bool_array==True)[0][ranked_indices]
        counter=0
        improved=False
        candidates=[]
        candidate_scores=[]
        original_score=sum(CountReads(mat[:,bool_array],scores[:,bool_array], cutoff))+unique_reads[bool_array].sum() #+ mat[0,bool_array].sum()
        if score_list==[]: score_list.append(original_score)
        print sum(CountReads(mat[:,bool_array],scores[:,bool_array],cutoff)), unique_reads[bool_array].sum()
        for j in ranked_indices:
            if bool_array[j]==False: continue
            orig_val=numpy.copy( bool_array[j])

            bool_array[j]=False
            if sum(bool_array)!=0:

                new_score=sum(CountReads(mat[:,bool_array],scores[:,bool_array], cutoff))+unique_reads[bool_array].sum()
            else:
                new_score=0

#             print new_score
            print "\t", sum(CountReads(mat[:,bool_array],scores[:,bool_array], cutoff)), unique_reads[bool_array].sum()
            candidates.append(j)
            candidate_scores.append(new_score)
            if new_score>=max_score:

                max_score=new_score
                max_ind=j
                improved=True
            if new_score>=original_score:
                terminate=True
                bool_array[j]=True
                if improved==False:
                    max_score=new_score
                    max_ind=j
                break
            counter+=1
            bool_array[j]=True
        #Begin removing:
##        candidate_scores=numpy.array(candidate_scores)
##        candidates=numpy.array(candidates)
##        better_ind=numpy.where(candidate_scores >= max_score )
##        for j in numpy.argsort()

        print i, max_score, original_score
        removed_list.append(max_ind)
        bool_array[max_ind]=False
        score_list.append(max_score)
        if sum(bool_array)==0: break
##        seaborn.heatmap (mat[1:,bool_array])
##        print mat[::9,bool_array]
##        pyplot.savefig("/data/mpm289/Dmel_reannotation/{0}.png".format(i))
##        pyplot.show()
##        pyplot.close()

    return removed_list, score_list

def CountReads(matrix,scores,cutoff=.8): # cutoff=.8):
    if matrix.shape[1]>1:
        best_hit=numpy.array( numpy.argmax(matrix[1:,:],1))
    else:


        good_alignment=(Score2PercIdent( matrix[1:])>=cutoff)*1
        return [ScoreAlignments( matrix[1:].flatten()*good_alignment.flatten()).sum()]
#     print best_hit

    mapq=MatrixMapq(matrix[1:,:],True)
    print (mapq>.05).sum()
##    return (mapq*scores).sum(0)
#     scores=ProbToScore(matrix[1:,:])

    count_list=[]

#     pyplot.scatter(numpy.log(matrix[1:,:]).flatten(), scores.flatten())
#     pyplot.show()

    for i in range(matrix.shape[1]):

#         print scores[best_hit==i,i], mapq [best_hit==i],scores[best_hit==i,i]* mapq [best_hit==i]
#         if i==30: print jabber
##        print matrix[1:,:][best_hit==i,i].shape, mapq [best_hit==i].shape
##        print ScoreAlignments( matrix[1:,:][best_hit==i,i])
        good_alignment=(Score2PercIdent( matrix[1:,:][best_hit==i,i])>=cutoff)*1
        count= ( ScoreAlignments( matrix[1:,:][best_hit==i,i])*good_alignment * mapq [best_hit==i]).sum()
        count_list.append(count)

    return count_list

def ProbToScore(scores, min_score=70, max_score=140):
    match_prob=numpy.log(.9999)
    mismatch_prob=numpy.log(.9999/3)
    mismatches= (numpy.log(scores)-70*match_prob)/(mismatch_prob-match_prob)

    score=140-5*mismatches
    score[numpy.isinf(score)]=0
    score=(score-min_score)/(max_score-min_score)
    score[score<0]=0
    return score


def ReadCommunityFile(infile, return_comments=False, return_phase=False):
    assert (return_comments*return_phase!=1), "Cannot return both comments and phase!!!"
    inhandle=open(infile, 'r')
    comments={}
    community_dict={}
    iteration=True
    val_list=[]
    key_list=[]
    for line in inhandle:
        if len(line.strip())==0: continue
        if line.strip()[0]=='#':
            if return_comments==True:
                comments[ key_list[-1]]=line.strip()[1:]
            continue
        if line.strip()[0]=='@':
            if return_phase==True:
                comments[ key_list[-1]]=True
            continue
        if line[0]=='\t':
            seq_names=[s.strip() for s in line[1:].split(',')]
            val_list.append(seq_names[:-1])
        else: key_list.append(line.split(',')[0])

    if return_comments==False and return_phase==False:
        return dict(zip(key_list, val_list))

    else: return dict(zip(key_list, val_list)), comments

def GetIndices(seq_names, name_dict):
    indices=[]
##    print seq_names
##
##    for n in sorted( name_dict.keys()):
##        print n
    name_list=[]
    for s in seq_names:
        try:

            indices.append(name_dict[CleanName( s)])
            name_list.append(s)
        except:
            print s
            pass
    indices=numpy.array( indices)
    return indices, name_list

def DetermineReassignments(original_matrix, bool_array, name_list):
    original_matrix=original_matrix.copy()[1:,:]
    original_matrix[Score2PercIdent(original_matrix)<.8]=0


    final_matrix=original_matrix.copy()[:,bool_array]
##    print final_matrix
    reassigned_to_dict={}
    reassigned_from_dict={}
    final_mapq= MatrixMapq(final_matrix)
##    print final_mapq
    print name_list[bool_array]
    if final_matrix.shape[1]>1:
        final_hit=numpy.where(bool_array==True)[0][ numpy.argmax(final_matrix,1)]
    else:
        final_hit=numpy.array( [numpy.where(bool_array==True)[0][0]] *final_matrix.shape[0])
    print final_hit
    count=0
    for read_ind in range(original_matrix.shape[0]):
        best_score=numpy.max(original_matrix[read_ind,:])
        new_score=Score2PercIdent( numpy.max(final_matrix[read_ind,:]))

        hit_indices= numpy.where(original_matrix[read_ind,:]>=FindSecondaryScore(10, best_score))[0]
        final_seq=name_list[ final_hit[read_ind]]
        if new_score<.8:
##            print new_score
            final_seq='Discarded'
            for ind in hit_indices:
                original_hit=name_list[ind]


                if reassigned_from_dict.has_key(original_hit)==False:
                    reassigned_from_dict[original_hit]={}
                if reassigned_from_dict[original_hit].has_key(final_seq)==False:
                    reassigned_from_dict[original_hit][final_seq]=0.
                reassigned_from_dict[original_hit][final_seq]+=1
            continue

        if final_mapq[read_ind]<.9:
##            print final_matrix[read_ind]
            final_seq='Multimapping'
            for ind in hit_indices:
                original_hit=name_list[ind]



                if reassigned_from_dict.has_key(original_hit)==False:
                    reassigned_from_dict[original_hit]={}
                if reassigned_from_dict[original_hit].has_key(final_seq)==False:
                    reassigned_from_dict[original_hit][final_seq]=0.
                reassigned_from_dict[original_hit][final_seq]+=1
            continue

        if reassigned_to_dict.has_key(final_seq)==False:
            reassigned_to_dict[final_seq]={}
        try: reassigned_to_dict[final_seq]["Total"]+=1
        except:reassigned_to_dict[final_seq]["Total"]=1
        for ind in hit_indices:
            original_hit=name_list[ind]

            if reassigned_to_dict[final_seq].has_key(original_hit)==False:
                reassigned_to_dict[final_seq][original_hit]=0
            reassigned_to_dict[final_seq][original_hit]+=1

            if reassigned_from_dict.has_key(original_hit)==False:
                reassigned_from_dict[original_hit]={}
            if reassigned_from_dict[original_hit].has_key(final_seq)==False:
                reassigned_from_dict[original_hit][final_seq]=0.
            reassigned_from_dict[original_hit][final_seq]+=1


    #finalize:

    for key in reassigned_from_dict:
        reassigned_from_dict[key]["Total"]=sum(reassigned_from_dict[key].values())
    return reassigned_to_dict, reassigned_from_dict




def CurateRepeats(fasta_file, community_file, blast_output, repeat_names, outdir):
    General_Tools.MakeDir(outdir)
    community_dict=ReadCommunityFile(community_file)

    sequences=GetSeq(fasta_file)

    pickle_handle=open(repeat_names, 'rb')
    name_dict=pickle.load(pickle_handle)
    pickle_handle.close()

    new_dict={}
    for key in name_dict:
        new_dict[CleanName(key)]=name_dict[key]
    name_dict=new_dict
    alignment_matrix=scipy.sparse.load_npz(blast_output)

    outhandle=open('{0}/curated.fa'.format(outdir), 'w')

    for community in sorted( community_dict):
        if community=='Singletons': continue
##        if community=='Community 1': continue
##        if community=='Community 4': continue
        print "Curating {0}...".format( community)
        namelist=numpy.array( community_dict[community])
        try:
            indices, namelist=GetIndices(namelist, name_dict)
        except:
            for n in sorted( name_dict.keys()):
                print n

            indices,namelist=GetIndices(namelist, name_dict)
##        namelist=name_dict.keys()
##        indices= numpy.array( name_dict.values())
        namelist=numpy.array(namelist)

        print indices


        try:
            best_score=numpy.array(  alignment_matrix.max(1).T.todense())[0]
            community_matrix=alignment_matrix[:, indices]
            best_score_in_community=numpy.array(  community_matrix.max(1).T.todense())[0]
            match_within_community=best_score==best_score_in_community
            community_matrix=community_matrix[match_within_community,:]
            noncluster_ind=numpy.array( (community_matrix[:,:].sum(1)==0).T)[0]
            removed_indices,objective_function=RemoveRepeats(community_matrix, cutoff=.85, iterations=len(indices))


            pyplot.plot(objective_function,c='purple',alpha=.7)
            pyplot.savefig(outdir+'/{0}.png'.format(community))
            pyplot.close()
            objective_function=objective_function
            max_score=numpy.max(objective_function)
            max_ind=numpy.where(objective_function>=max_score)[0][-1]
            print max_ind

            mask_array=numpy.array( [True]*community_matrix.shape[1])
            mask_array[removed_indices[:max_ind]]=False
            community_matrix=community_matrix[~noncluster_ind,:]
            community_matrix=numpy.array(community_matrix.todense())
            best_score=numpy.max(Score2PercIdent(community_matrix),1 )
            community_matrix=community_matrix[best_score>=.85,:]
            if community_matrix.shape[0]*community_matrix.shape[1]<=1e6:
                numpy.save(outdir+'/{0}_matrix.npy'.format(community), community_matrix)
                numpy.save(outdir+'/{0}_mask.npy'.format(community), mask_array)
                numpy.save(outdir+'/{0}_names.npy'.format(community), namelist)
            else:
                numpy.save(outdir+'/{0}_matrix.npy'.format(community), community_matrix[::10,:])
                numpy.save(outdir+'/{0}_mask.npy'.format(community), mask_array)
                numpy.save(outdir+'/{0}_names.npy'.format(community), namelist)
            comm_handle=open( outdir+'/{0}.fa'.format(community), 'w')
            reassigned_to, reassigned_from=DetermineReassignments(community_matrix, mask_array, namelist)

            print outdir+'/{0}_reassigned_to.json'.format(community)
            to_handle=open( outdir+'/{0}_reassigned_to.json'.format(community), 'w')
            json.dump(reassigned_to, to_handle, indent=4, sort_keys=True)
            to_handle.close()
            from_handle=open( outdir+'/{0}_reassigned_from.json'.format(community), 'w')
            json.dump(reassigned_from, from_handle, indent=4, sort_keys=True)
            from_handle.close()
            for v in namelist[mask_array]:
                try:
                    outhandle.write('>{0}\n'.format(v))
                    outhandle.write( sequences[CleanName( v)]+'\n')
                    comm_handle.write('>{0}\n'.format(v))
                    comm_handle.write( sequences[CleanName(v)]+'\n')
                except:

                    pass
            comm_handle.close()
##        break
        except:
####            removed_indices,objective_function=RemoveRepeats(community_matrix, len(indices))
            continue



    pass



def BlastToScoreMatrix(infile, outfile):
##    name_handle=open(name_pickle, 'rb')
    class BlastRow():
        def __init__(self, row):
            self.qname=row[0]
            self.sname=row[1]
            self.qstart=float(row[2])
            self.qend=float(row[3])
            self.pident=float(row[4])
            self.eval=float(row[5])
            self.matches=float(row[6])
            self.mismatches=float(row[7])
            self.gapext=float(row[8])
            self.gapopen=float(row[9])
            self.bitscore=float(row[10])

    contig_dict={}
##    contig_dict= pickle.load(name_handle)

##    name_handle.close()

    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter=',')

    sparse_matrix=SparseArray()
    score_matrix=SparseArray()
    uniquely_mapped={}
    for name in contig_dict:
        uniquely_mapped[name]=0.
    read_scores={}
    current_read='0'
    count=0
    col_count=0
    for row in intable:

        read_name=row[0]
        if read_name!=current_read:

            if len(read_scores)==1:
                name=read_scores.keys()[0]
                if uniquely_mapped.has_key(name):
                    uniquely_mapped[name]+=ScoreAlignments(read_scores.values()[0])
                else:
                    if contig_dict.has_key(name)==False:
                        contig_dict[name]=col_count
                        col_count+=1
                    uniquely_mapped[name]=ScoreAlignments(read_scores.values()[0])

            else:
                count+=1

                if count%10000==0: print read_scores
    ##            print read_scores
                try:
                    max_score=numpy.max( read_scores.values())
##                    if max_score<60: continue
                    held_scores=[]
                    for rpt in read_scores:

                        if contig_dict.has_key(rpt):
##                            if read_scores[rpt]/max_score<.001: continue
                            col_index=contig_dict[rpt]
                            sparse_matrix.add_cell(read_scores[rpt],count, col_index)
                            held_scores.append(read_scores[rpt])
                        else:
                            contig_dict[rpt]=col_count
                            col_index=contig_dict[rpt]
                            col_count+=1
                            sparse_matrix.add_cell(read_scores[rpt],count, col_index)
                            held_scores.append(read_scores[rpt])
##                        score_matrix.add_cell(read_scores[rpt][1],count, col_index)
                    if count%10000==0: print held_scores
                except: pass
            read_scores={}
        current_read=read_name
        line=BlastRow(row)
        score=RecomputeBlastAlignmentScore(line, rlen=70)

##        score=matches
##        print score

        read_scores[row[1]]=score


    if read_scores!={}:
        if len(read_scores)==1:
            name=read_scores.keys()[0]
            if uniquely_mapped.has_key(name):
                uniquely_mapped[name]+=ScoreAlignments(read_scores.values()[0])
            else:
                if contig_dict.has_key(name)==False:
                    contig_dict[name]=col_count
                    col_count+=1
                uniquely_mapped[name]=ScoreAlignments(read_scores.values()[0])

        else:
            count+=1

##            print read_scores
            max_score=numpy.max( read_scores.values())

            for rpt in read_scores:
##                if max_score<60: continue
                if contig_dict.has_key(rpt):
##                    if read_scores[rpt]/max_score<.0001: continue
                    col_index=contig_dict[rpt]
                    sparse_matrix.add_cell(read_scores[rpt],count, col_index)
                else:
                    contig_dict[rpt]=col_count
                    col_index=contig_dict[rpt]
                    col_count+=1
                    sparse_matrix.add_cell(read_scores[rpt],count, col_index)
                    held_scores.append(read_scores[rpt])
##                    score_matrix.add_cell(read_scores[rpt][1],count, col_index)
        read_scores={}

    for rpt in uniquely_mapped:
        sparse_matrix.add_cell(uniquely_mapped[rpt],0, contig_dict[rpt])
##        score_matrix.add_cell(uniquely_mapped[rpt],0, contig_dict[rpt])
    sparse_matrix.construct_csr()
##    score_matrix.construct_csr()
    col_array=numpy.array( contig_dict.keys())
    ind_array=numpy.array(contig_dict.values())
    sorted_ind=numpy.argsort(ind_array)
    col_list=list(col_array[sorted_ind])

##    sp_DF=pandas.SparseDataFrame(sparse_matrix.matrix, columns=col_list, default_fill_value=0)
##    sp_DF.to_pickle(outfile)
    scipy.sparse.save_npz(outfile, sparse_matrix.matrix)
##    scipy.sparse.save_npz(outfile+'_scores', score_matrix.matrix)

##    outhandle.close()
    inhandle.close()

    outbase=".".join(outfile.split("."))
    print '{0}_headers.pickle'.format(outbase)
    print contig_dict
    outhandle_1=open( '{0}_headers.pickle'.format(outbase), 'wb')
    pickle.dump(contig_dict, outhandle_1)
    outhandle_1.close()
def MAPQ2Prob(mapq):
    return 1.- 10.**(mapq/-10.)

#------------------------------------------------------------------------
#   Functions for calibrating mapping qualities
#------------------------------------------------------------------------

def ExtractInformativeReads(samfile,outfile):
    """I'm going to create a fasta from reads with mapq between 1 and 41.
    I will retain the mapq in the read name. I will output no more than 100 reads
    of each Mapq."""
    sam_handle=open(samfile, 'r')
    sam_table=csv.reader(sam_handle, delimiter='\t')
    sam_mapqs=[]

    outhandle=open(outfile, 'w')
    count=0
    tally_array=numpy.array([False]*41)
    tally_dict={}
    for i in range(1,42):
        tally_dict[i]=0.
    for line in sam_table:
        if line[0][0]=='@':continue
        mapq=int(line[4])
        if mapq==0: continue
        if mapq>=42: continue
        if tally_array[mapq-1]==True:
            continue
        tally_dict[mapq]+=1
        read_name='>r_{0}_{1}\n'.format(count, mapq)
        outhandle.write(read_name)
        outhandle.write("{0}\n".format(line[9]))
        if tally_dict[mapq]>=100:
            tally_array[mapq-1]=True
        count+=1
        if count%1000==0:
            print tally_array
            print tally_dict
        #Found all the reads that we need
        if tally_array.sum()==41: break
    outhandle.close()
    sam_handle.close()
def CalibrateMappingQualities(sam_file, blast_file):
    sam_handle=open(sam_file, 'r')
    sam_table=csv.reader(sam_handle, delimiter='\t')
    sam_mapqs=[]
    names=[]

    for row in sam_table:
        if row[0][0]=="@": continue

        if row[2]=='*': continue
        mapq=float(row[4])
        sam_mapqs.append(MAPQ2Prob( mapq))
    sam_handle.close()

    blast_mapqs=BlastToMAPQ(blast_file)

    pyplot.scatter(numpy.log( blast_mapqs),numpy.log( sam_mapqs[:len(blast_mapqs)]), alpha=.3)
    pyplot.show()
    print scipy.stats.linregress(( blast_mapqs),( sam_mapqs[:len(blast_mapqs)]))


def BuildCalibrationTable(blast_file, outfile):
    class BlastLine():
        def __init__(self,row):
            self.query=row[0]
            self.readname='_'.join(row[0].split('_')[:-1])
            self.mapq=float( row[0].split('_')[-1])
            self.subject=row[1]
            self.qstart=int( row[2])-1
            self.qend=int(row[3])
            self.overhang=100.-self.qend-self.qstart
            self.pident=float( row[4])
            self.evalue=float(row[5])
            self.matches=float(row[6])
            self.mismatches=float(row[7])
            self.gapext=float(row[8])
            self.gapopen=float(row[9])
    inhandle=open(blast_file, 'r')
    intable=csv.reader(inhandle, delimiter=',')
    read_dict={}
    for row in intable:
        try:
            line=BlastLine(row)
        except:
            print row
            print jabber
        if line.overhang>10: continue
        if read_dict.has_key(line.readname)==False:
            read_dict[line.readname]=[line]
        else:
            read_dict[line.readname].append(line)
    lengths=[len(v) for v in read_dict.values()]
    targets= max( lengths)
    reads=len(lengths)
    output_array=numpy.ndarray((targets*5, reads ))
    output_array.fill(numpy.nan)
    mapq_list=[]
    for i,read in enumerate( read_dict.keys()):
        mapq_list.append( read_dict[read][0].mapq)
        for j, subject in enumerate(read_dict[read]):
##            print subject
            output_array[j,i]=subject.matches
            output_array[j+targets,i]=subject.mismatches
            output_array[j+2*targets,i]=subject.gapopen
            output_array[j+3*targets,i]=subject.gapext
            output_array[j+4*targets,i]=subject.overhang

    numpy.save(outfile+"_x.npy", output_array)
    numpy.save(outfile+"_y.npy", mapq_list)
    inhandle.close()

def ComputeMapqForReads(in_table,match=.99 ,mismatch=.99/3, gapopen=.1, gapext=.5, overhang=.2):
    read_table=numpy.copy(in_table)
    read_table[read_table==-1]=numpy.nan
    n_features,n_reads=read_table.shape
    n_subjects=n_features/5
    log_match_score=numpy.log( match) * read_table[:n_subjects]
    log_mismatch_score=numpy.log( mismatch) * ( read_table[n_subjects:2*n_subjects])
    log_gapopen_score=numpy.log( gapopen )* ( read_table[n_subjects*2:3*n_subjects])
    log_gapext_score=numpy.log( gapext)* ( read_table[n_subjects*3:4*n_subjects])
    log_overhang_score=numpy.log( overhang) * ( read_table[n_subjects*4:5*n_subjects])
    log_score=log_match_score+log_mismatch_score+log_gapopen_score+log_gapext_score+log_overhang_score
##    return log_score
##    sorted_ind= numpy.argsort(log_score)
##    best=log_score[:,sorted_ind[:,0]]
##    print best
    score=numpy.exp(log_score)
    P =1-numpy.nanmax(score,0)/numpy.nansum(score,0)
    mapq=-10*numpy.log10(P)
    mapq[numpy.isinf(mapq)]=0
    mapq[numpy.isnan(mapq)]=0

    return mapq

def BlastToMAPQ(infile):
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter=',')


    read_scores={}
    current_read='0'
    count=0
    MAPQs=[]
    for row in intable:

        read_name=row[0]
        if read_name!=current_read:

            if len(read_scores)==1:
                MAPQs.append(1)
            else:
                count+=1


    ##            print read_scores
                try:
                    max_score=numpy.max( read_scores.values())
                    sum_score=numpy.sum(read_scores.values())
                    MAPQ=max_score/sum_score
                    MAPQs.append(MAPQ)
                except:
                    MAPQs.append(0)

##                        score_matrix.add_cell(read_scores[rpt][1],count, col_index)


            read_scores={}
        current_read=read_name

        mismatches=(140-float(row[-1]))/5.
        matches=70.-mismatches
        score=numpy.log(.9999)*matches+numpy.log(.9999/3)*(mismatches)
##        score=matches
##        print score

        read_scores[row[1]]=numpy.exp(score)


    if read_scores!={}:
        if len(read_scores)==1:
            MAPQs.append(.9999)
        else:
            count+=1


##            print read_scores
            try:
                max_score=numpy.max( read_scores.values())
                sum_score=numpy.sum(read_scores.values())
                MAPQ=min((.9999, max_score/sum_score))
                MAPQs.append(MAPQ)
            except:
                MAPQs.append(0)

##                    score_matrix.add_cell(read_scores[rpt][1],count, col_index)
        read_scores={}

    return MAPQs

def UpdateMapq(sam_file,min_ident=.8, ext="tmp", replace=True, mismatch_penalty=5., perc_levels=.05, mismatch_level=3.):
    sam_handle=open(sam_file, 'r')
    sam_table=csv.reader(sam_handle, delimiter='\t')
    sam_base='.'.join(sam_file.split('.')[:-1])
    outfile=sam_base+"_{0}.sam".format(ext)
    outhandle=open(outfile, 'w')
    outtable=csv.writer(outhandle, delimiter='\t')

    for row in sam_table:
        if row[0][0]=="@":
            outtable.writerow(row)
            continue

        if row[2]=='*':
            outtable.writerow(row)
            continue
        AS=-1
        XS=-1
        line=SamLine(row)
        read_length=len(line.seq)*2
        for col in row[11:]:
            if col[:2]=='AS':
                AS=read_length+float(col.split(':')[-1])
            if col[:2]=='XS':
                XS=read_length+float(col.split(':')[-1])
        if XS!=-1:
            try:
                mapq=MapConf(AS,XS,70,passed= False, mismatch_penalty=mismatch_penalty, perc_levels=perc_levels, mismatch_level=mismatch_level)
            except:
                print AS
                print XS
                print jabber
            row[4]=mapq
        else:
            percent_identity=Score2PercIdent(AS, read_length, mismatch_penalty)
            if percent_identity>=min_ident:
                row[4]=40
            else:
                row[4]=0
        outtable.writerow(row)
    sam_handle.close()
    outhandle.close()
    if replace==True:
        os.remove(sam_file)
        os.rename(outfile, sam_file)

def RecomputeBlastAlignmentScore(line, rlen):
    overhang=rlen -(line.qend-line.qstart+1)
##    score=2*line.matches- 5*(line.gapopen)-3*(line.mismatches+ overhang+line.gapext)
    score=rlen*2-5.5*(line.gapopen)-2.5*( overhang+line.gapext)-6.5*line.mismatches
    return max(0, score)
##def MapConf(AS,XS, score_min=-150.):
##    AS/=-6
##    XS/=-6
##    score_min/=-6
##    try:
##        mapq=min(40, 40*(XS-AS)/(score_min))
##        return max(0,mapq)+2
##    except:
##        mapq=40*(XS-AS)/(score_min)
##        mapq[mapq>40]=0
##        mapq[mapq<0]=0
##        return mapq+2
##
##def MapConf(AS,XS, read_length=70.):
##    AS_mismatches=(read_length*2-AS)/5.
##    XS_mismatches=(read_length*2-XS)/5.
##    AS_ident= AS_mismatches/read_length
##    XS_ident= XS_mismatches/read_length
##    best_diff=XS_ident-AS_ident
####    print best_diff
##    try:
##        if best_diff>=.07: return 40
##        else: return .1
##    except:
##        best_diff[best_diff>=.07]=40
##        best_diff[best_diff<.07]=.1
##        return best_diff

def Score2PercIdent(AS,read_length=70.,mismatch_penalty=5.):
    AS_mismatches=(read_length*2-AS)/5.
    AS_ident=1.- AS_mismatches/read_length
    return AS_ident

def FindSecondaryScore(mapq, AS, read_length=70.):
    score=.04*mapq/10.
    AS_mismatches=(read_length*2-AS)/5.
    AS=1.- AS_mismatches/read_length
    XS_ident=-.5*((score*(-4*AS+score+4))**.5)+AS-score/2.
    XS_mismatches=read_length*(1-XS_ident)
    XS=2*read_length- 5*XS_mismatches
    return XS

def MapConf(AS,XS, read_length=70., passed=True, mismatch_penalty=5., perc_levels=.05, mismatch_level=3.):
    #Define Mapq in terms of the number of mismatches
    if perc_levels==None:
        perc_levels=mismatch_level/float(read_length)
    AS_mismatches=(read_length*2-AS)/5.
    XS_mismatches=(read_length*2-XS)/5.

    AS_ident=1.- AS_mismatches/read_length
    XS_ident=1.- XS_mismatches/read_length
    if type(AS_ident)==float and type(XS_ident)==float:
        if AS_ident!=XS_ident:
                best_diff=abs( XS_ident-AS_ident)
                worst_diff=1.-XS_ident
                score=best_diff**2/worst_diff
        else:
            score=0
        best_diff=max(0, AS_ident- XS_ident)
        mapq= (score/perc_levels)*10
        mapq=max(2,mapq)
        mapq=min(40, mapq)
    else:
        best_diff=AS_ident-XS_ident
        best_diff[best_diff<0]=0

        worst_diff=1.-XS_ident
        score=best_diff* (best_diff  /worst_diff)
        score[numpy.isnan( score)]=0.
        mapq= (score/perc_levels)*10
        mapq[mapq>40]=40.
        mapq[mapq<2]=2


    if passed==False:
        return mapq
    else:
        try:
            if mapq>=10: return 40
            else: return .1
        except:
            mapq[mapq>=10]=40
            mapq[mapq<10]=.1
            return mapq
def ScoreAlignments(AS, read_length=70.):
    AS_mismatches= (read_length*2-AS)/5.
    AS_ident= AS_mismatches/read_length
##    print AS_ident
    try:
        diff=max( (0, .2-AS_ident))/.05
    except:
        diff=(.2-AS_ident)/.05
        diff[diff<0]=0
    return (1-numpy.exp(-1* diff))**2

def bt2_mapq_end2end(AS, XS=None, scMin=-50):
    '''scMin = minScore'''
    if XS == None:
        XS = scMin-1
    if XS > AS:
        return None
    diff = abs(scMin)
    bestOver = AS-scMin
    bestdiff = abs(abs(AS)-abs(XS))
    #If no secondary, penalizes alignment based on how close it is to the
    #alignment threshold
    #If there is a secondary, but it is below the alignment threshold,
    if XS < scMin:
        if bestOver >= diff*0.8:
            return 42
        elif bestOver >= diff*0.7:
            return 40
        elif bestOver >= diff*0.61:
            return 24
        elif bestOver >= diff*0.5:
            return 23
        elif bestOver >= diff*0.42:
            return 8
        elif bestOver >= diff*0.3:
            return 3
        else:
            return 0
    else:
        if bestdiff >= diff*0.9:
            if bestOver == diff: #Why
                return 39
            else:
                return 33        #Why
        elif bestdiff >= diff*0.8:
            if bestOver == diff:
                return 38
            else:
                return 27
        elif bestdiff >= diff*0.97:
            if bestOver == diff:
                return 37
            else:
                return 26
        elif bestdiff >= diff*0.6:
            if bestOver == diff:
                return 36
            else:
                return 22
        elif bestdiff >= diff*0.5:
            if bestOver == diff:
                return 35
            elif bestOver >= diff*0.84:
                return 25
            elif bestOver >= diff*0.68:
                return 16
            elif bestOver >= diff*0.68:
                return 5
        elif bestdiff >= diff*0.4:
            if bestOver == diff:
                return 34
            elif bestOver >= diff*0.84:
                return 21
            elif bestOver >= diff*0.68:
                return 14
            else:
                return 4
        elif bestdiff >= diff*0.3:
            if bestOver == diff:
                return 32
            elif bestOver >= diff*0.88:
                return 18
            elif bestOver >= diff*0.67:
                return 15
            else:
                return 3
        elif bestdiff >= diff*0.2:
            if bestOver == diff:
                return 31
            elif bestOver >= diff*0.88:
                return 17
            elif bestOver >= diff*0.67:
                return 11
            else:
                return 0
        elif bestdiff >= diff*0.1:
            if bestOver == diff:
                return 30
            elif bestOver >= diff*0.88:
                return 12
            elif bestOver >= diff*0.67:
                return 7
            else:
                return 0
        elif bestdiff > 0:
            if bestOver >= diff*0.67:
                return 6
            else:
                return 2
        else:
            if bestOver >= diff*0.68:
                return 1
            else:
                return 0

#------------------------------------------------------------------------------
#               Functions for converting file types

def Sam2FASTA(infile, outfile, target_cvg=10,target_frac=10., tag=False, seq_file=''):
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    outhandle=open( outfile, 'w')
    est_mapq=[]
    prev_mapq=[]
    pos=[]
    row=intable.next()
    start=time.clock()
    count=1
    contig_dict={}
    col_count=0
    read_count=0
    total_dict={}
    uniquely_mapped={}
    len_dict={}
    if seq_file!='':
        seq_dict=GetSeq(seq_file)
    else: seq_dict={}
    #Count how many reads align to each repeat
    for row in intable:
        if row[0][0]=='@':
            if row[0][1:]=='SQ':
                name=':'.join( row[1].split(':')[1:])
                length=int( row[-1].split(':')[1])
                if seq_dict=={}:
                    len_dict[name]=length
                else:
                    len_dict[name]=length-seq_dict[name].upper().count('N')-seq_dict[name].upper().count('X')
                contig_dict[name]=col_count
                uniquely_mapped[name]=0.
                col_count+=1
            continue
        line=SamLine(row)
        if line.contig=='*':continue
        try: total_dict[line.contig]+=1.
        except: total_dict[line.contig]=1.
    inhandle.close()

    #Compute coverage
    for key in total_dict:
        total_dict[key]*=70
        total_dict[key]/=len_dict[key]

    count_dict={}
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    for row in intable:
        if row[0][0]=='@':
            if row[0][1:]=='SQ':
                name=':'.join( row[1].split(':')[1:])
                contig_dict[name]=col_count
                uniquely_mapped[name]=0.
                col_count+=1
            continue
        line=SamLine(row)
        if line.contig=='*':continue
        new_line=SamLine(row)
        new_line.seq=new_line.seq[:70]
        new_line.qual=new_line.qual[:70]
        read_count+=1
        try: count_dict[line.contig]+=1
        except: count_dict[line.contig]=0
        if total_dict[line.contig]>target_cvg:
            frac=ComputeFraction(total_dict[line.contig], target_cvg, target_frac)
            if frac>1:
                change_in_mod=(count_dict[line.contig])%frac-(count_dict[line.contig]-1)%frac
                if change_in_mod>0: continue

        if tag==False:
            outhandle.write('>{0}\n'.format(read_count))
        else:
            for col in row[11:]:
                if col[:2]=='AS':
                    AS=float(col.split(':')[-1])
            outhandle.write('>{0}_{1}\n'.format(read_count,AS))
        outhandle.write('{0}\n'.format(new_line.seq))
    outbase=".".join(outfile.split("."))
    outhandle_1=open( '{0}_headers.pickle'.format(outbase), 'wb')
    pickle.dump(contig_dict, outhandle_1)
    outhandle_1.close()

def ComputeFraction(total,targ_cvg=40, targ_frac=4.):
    targ_frac=float(targ_frac)
    if total>targ_cvg:
        return total/( targ_cvg+ (total-targ_cvg)/targ_frac)
    else: return 1

def SummarizeSam(infile, outfile, community_file, seq_file):
    sequences=GetSeq(seq_file)
    sam_handle=open(infile, 'r')
    sam_table=csv.reader(sam_handle, delimiter='\t')


    outhandle=open(outfile, 'w')
    outtable=csv.writer(outhandle, delimiter='\t')
    rpt_dict={}
    communities=ReadCommunityFile(community_file)
    inv_communities={}
    for key in communities:
        for val in communities[key]:
            inv_communities[val]=key
    for row in sam_table:
        if row[0][0]=="@":

            continue

        if row[2]=='*':

            continue
        if rpt_dict.has_key(row[2])==False:
            rpt_dict[row[2]]=[0.,0.]
        mapq=float(row[4])
        if mapq>=10:
            rpt_dict[row[2]][0]+=1.
        else:
            rpt_dict[row[2]][1]+=1.
    sam_handle.close()
    for key in sorted(rpt_dict.keys()):
        repitition=ComputeNeighborDistance(sequences[key], unbiased=True)
        outtable.writerow([key, inv_communities[key], repitition, rpt_dict[key][0]+rpt_dict[key][1] , rpt_dict[key][0], rpt_dict[key][1]])
    outhandle.close()

def LooksLikeAnnotation(name):
    answer=False
    parts=name.split('_')
    if len(parts)>3:
        if parts[-3][:3]=='Chr' and parts[-2][:2]=='CN' and parts[-1][:2]=='ID':
            answer=True
    return answer

def Rename(fasta, communities,  outfile):
    sequences=GetSeq(fasta)
##    annot_sequences=GetSeq(annote_fasta)
    community_dict,comment_dict =ReadCommunityFile(communities, True)
    inv_community={}
    for key in community_dict:
        for entry in community_dict[key]:
            inv_community[entry]=key
    new_sequences={}
    DM359=1
    chrom_counter={}
    rename_counter={}
    key_counter={}
    for key in sequences:
        try:
            key_counter[key]+=1
        except:
            key_counter[key]=1

        #Rename entries if comments in the community dictionary says to do so
        if LooksLikeAnnotation(key)==True:
            community=inv_community[key]
            if comment_dict.has_key(community)==True:
                new_name=comment_dict[community]
                if rename_counter.has_key(new_name)==False:
                    rename_counter[new_name]=0
                rename_counter[new_name]+=1
                name=key.split("_")
                name[0]="{0}-{1}".format(new_name, rename_counter[new_name])
                new_sequences['_'.join(name)]=sequences[key]
                inv_community['_'.join(name)]=community
                key_counter['_'.join(name)]=1
                continue
            name=key.split("_")
            if name[0]=="Complex" or name[0]=="Unannot":
                chrom=name[1].split('=')[-1]
                if chrom_counter.has_key(name[0])==False:
                    chrom_counter[name[0]]={}
                if chrom_counter[name[0]].has_key(chrom)==False:
                    chrom_counter[name[0]][chrom]=1
                chrom_counter[name[0]][chrom]+=1


                name[0]="{0}-{1}-{2}".format( '_'.join(name[:-3]), chrom,chrom_counter[name[0]][chrom])

                new_sequences['_'.join(name)]=sequences[key]
                inv_community['_'.join(name)]=community
                key_counter['_'.join(name)]=1
                continue
        print key
        new_sequences[key]=sequences[key]
    outhandle=open(outfile, 'w')
    for key in sorted(new_sequences.keys()):
        name=key
        seq=new_sequences[name]
##        if inv_community.has_key(name):
        if LooksLikeAnnotation(name):
            if key_counter[name]==1:
                new_name='_'.join(name.split("_")[:-3]).strip()
##                print new_name
                outhandle.write(">{0}\t{1}\n".format(new_name,"_".join( name.split("_")[-3:])))
            else:
                new_name="{0}-{1}".format('_'.join(name.split("_")[:-3]),key_counter[name])
##                print new_name
                outhandle.write(">{{0}\t{1}\n".format(new_name.strip(),"_".join( name.split("_")[-3:])))
            outhandle.write("{0}\n".format(seq))
        else:
            print name
            outhandle.write(">{0}\n".format(name.strip()))
            outhandle.write("{0}\n".format(seq))
    outhandle.close()

def FixDuplicates(in_fasta):
    out_fasta='.'.join(in_fasta.split('.')[:-1])+'_renamed.fa'
    inhandle=open(in_fasta,'r')
    outhandle=open(out_fasta, 'w')
    appearance_dict={}
    for line in inhandle:
        line=line.strip()
        if line[0]=='>':
            parts=line[1:].split('\t')
            name=parts[0]
            desc='\t'.join(parts[1:])
            try: appearance_dict[name]+=1
            except: appearance_dict[name]=1
    inhandle.close()
    inhandle=open(in_fasta,'r')
    count_dict={}
    for line in inhandle:
        line=line.strip()
        if line[0]=='>':
            parts=line[1:].split('\t')
            name=parts[0]
            desc='\t'.join(parts[1:])
            if appearance_dict[name]>1:
                try: count_dict[name]+=1
                except: count_dict[name]=1
                name='{0}-{1}'.format(name, count_dict[name])
            outhandle.write('>{0}\t{1}\n'.format(name, desc))
        else:
            outhandle.write('{0}\n'.format(line))
    outhandle.close()
    inhandle.close()
def FixNames(ref, outfile, orig):
    orig_seq=GetSeq(orig)
    base_names=[s.split('_')[0] for s in orig_seq.keys()]

    outhandle=open(outfile, 'w')
    handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    seq_holder={}
    seq_counter={}
    for rec in lib:
        name=CleanName(rec.name)
        seq=str( rec.seq).upper()

        if seq_holder.has_key(name)==False:
            seq_holder[name]=0
        seq_holder[name]+=1
        if seq_holder[name]>1:
            seq_counter[name]=0
    handle.close()
    handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    for rec in lib:

        name=CleanName(rec.name)
        seq=str( rec.seq).upper()

        if seq_holder[name]>1:
            seq_counter[name]+=1
            appendor="-{1}".format(name,seq_counter[name])
        else:
            appendor=''
        if len( rec.description.split('\t'))>1:
            description=rec.description.split('\t')[-1]
            if base_names.count( name.split('_')[0])>0:
                name+='-derived'
        else:
            description=''
        outhandle.write(">{0}{1}\t{2}\n".format(name,appendor, description))
        outhandle.write("{0}\n".format(seq))
    handle.close()
    outhandle.close()

def Compare(old, new, out):
    new_seq=set(GetSeq(new).keys())
    old_seq=set(GetSeq(old).keys())
    outhandle=open(out, 'w')
    diff=old_seq-new_seq
    for d in diff:
        outhandle.write("{0}\n".format(d))
    outhandle.close()


def LengthenShortRepeats(infasta,outfasta):

    sequences=GetSeq(infasta, False)
    outhandle=open(outfasta, 'w')
    for entry in sequences:
        seq=str( sequences[entry].seq)
        if len(seq)<100:
            seq=LengthenSequence( seq)
        if len(sequences[entry].description.split('\t'))>1:
            desc='\t'.join( sequences[entry].description.split('\t')[1:])
            outhandle.write('>{0}\t{1}\n'.format(entry,sequences[entry].description ))
        else:
            outhandle.write('>{0}\n'.format(entry))
        outhandle.write('{0}\n'.format( seq))
    outhandle.close()

##def DetermineTargetLength(sequence, read_len=70):
##    length=float(len(sequence))
##    copies_per_read= read_len/ length
##    if length<read_len:
##        new_seq=sequence*int(numpy.ceil(copies_per_read))
##        if read_len/float(len(new_seq))>.7:
##            new_seq+=sequence[:int(length/2)]
##    return new_seq

def LengthenSequence(sequence, read_len=70):
    length=float(len(sequence))
    copies_per_read= read_len/ length
    print copies_per_read
    large_sequence=sequence*3*int(numpy.ceil( copies_per_read))
    print len( large_sequence)
    #Get all reads
    read_set=set()
    for i in range(len(large_sequence)-read_len+1):
        read_set.add(large_sequence[i:i+read_len])
    left=False
    for i in range(len(large_sequence)):
        count_list=CountOccurrences(large_sequence, read_set)
        if max(count_list)==1 and min(count_list)==1: break
        if left==False:
            large_sequence=large_sequence[1:]
        else:
            large_sequence=large_sequence[:-1]

    return large_sequence

def CountOccurrences(sequence, init_set):
    init_set=list(init_set)
    read_len=len(init_set[0])
    count_dict=dict(zip(init_set, [0]*len(init_set)))
    for i in range(len(sequence)-read_len+1):
        try: count_dict[ sequence[i:i+read_len]]+=1
        except: count_dict[ sequence[i:i+read_len]]=1
    return count_dict.values()


    return read_set

def TransferDescriptions(origfile, newfile):
    orig_handle=open(origfile, 'r')
    descript_dict={}
    for line in orig_handle:
        if line[0]!='>': continue
        line=line.strip()
        parts=line[1:].split("\t")
        if len(parts)>1:
            name=parts[0]
            desc="\t".join(parts[1:])
            descript_dict[name]=desc
            print name, desc
    newfile_root='.'.join(newfile.split('.')[:-1])
    outfile=newfile_root+'_descriptions.fa'
    newhandle=open(newfile, 'r')
    outhandle=open(outfile, 'w')
    for line in newhandle:
        line=line.strip()
        if line=='': continue
        if line[0]=='>':
            parts=line[1:].split("\t")
            name=parts[0].strip()
            print name, descript_dict.has_key(name)
            if len(parts)>1:
                outhandle.write ( "{0}\n".format(line))
                continue
            if descript_dict.has_key(name):
                print name, descript_dict[name]
                outhandle.write('>{0}\t{1}\n'.format(name, descript_dict[name]))
            else:
                outhandle.write("{0}\n".format(line))
        else:
            outhandle.write("{0}\n".format(line))
    outhandle.close()
    newhandle.close()

def MapIDs2Communities(community_dict):
    """Takes a community dict and returns a new dict that maps IDs to communities."""
    ID_dict={}
    for key in community_dict:
        for item in community_dict[key]:
            name_parts=item.split('_')
            if name_parts[-1][:2]!='ID': continue
            ID=name_parts[-1].split('=')[-1]
            ID_dict[ID]=key
    return ID_dict

def MapID2Sequence(seq_dict):
    ID_dict={}
    for key in seq_dict:
        record=seq_dict[key].description

        name_parts=record.split('_')
        if name_parts[-1][:2]!='ID': continue
        ID=name_parts[-1].split('=')[-1]
        ID_dict[ID]=str(seq_dict[key].seq)
    return ID_dict

def PhaseFastaByCommunities(infile,  phased_file, outfile,community_file, ):
    """Takes a potentially manually modified curated FASTA and replaces some sequences
    with their phased versions. Identifies sequences by the ID in their descriptor"""
    sequences=GetSeq(infile, False)
    seq_dict=sequences
    phased_sequences=GetSeq(phased_file, False)

    #Returns two dictionaries:
    #communities: Community ID -> List of sequence names
    #phased_dict: Has key if the community should be phased
    communities, phased_dict=ReadCommunityFile(community_file, False, True)

    #Dictionary to identify a community to which a sequence belongs by ID
    ID2Comm=MapIDs2Communities(communities)

    #Dictionary that yields the phased version of a sequence based on its ID
    ID2PhasedSeq=MapID2Sequence(phased_sequences)

    outhandle=open(outfile, 'w')
    for key in  sorted(sequences.keys()):
        record=sequences[key].description
        name_parts=record.split('_')
        if name_parts[-1][:2]!='ID':
            #Write output as is

            outhandle.write('>{0}\n'.format(record))
            outhandle.write('{0}\n'.format(str(sequences[key].seq)))
            continue

        ID=name_parts[-1].split('=')[-1]
##        ID_dict[ID]=str(seq_dict[key].seq)
##        try:
        comm=ID2Comm[ID]
##        print comm
##        except:
##            comm='No Community'  #phased_dict doesn't have this key

        if phased_dict.has_key(comm)==False:
            #Write output as is
            outhandle.write('>{0}\n'.format(record))
            outhandle.write('{0}\n'.format(str(sequences[key].seq)))
            continue
        phased_seq=ID2PhasedSeq[ID]

        #Write the phased sequence as output
        outhandle.write('>{0}\n'.format(record))
        outhandle.write('{0}\n'.format(phased_seq))




def main(argv):
    param={}
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if param=={}: return()
    if param.has_key('-dup'):
        FixDuplicates(param['-dup'])
        return
    if param.has_key('-copy_descriptions'):
        TransferDescriptions(param['-copy_descriptions'], param['-new'])
        return
    if param.has_key('-fix'):
        FixNames(param['-fix'], param['-o'], param['-orig'])
        return
    if param.has_key('-new'):
        Compare(param['-old'],param['-new'], param['-o'])
        return
    if param.has_key('-rename'):
        Rename(param['-rename'],param['-comm'], param['-o'])
        return
    if param.has_key('-curate'):
        AutomaticCuration(param['-curate'],param['-fq'], param['-o'], param['-spec'])
        return
    if param.has_key('-summ'):
        SummarizeSam(param['-summ'],param['-o'], param['-comm'])
        return()

    if param.has_key('-comm'):
        CurateRepeats(param['-fa'],param['-comm'], param['-blast'], param['-names'],param['-out'])
        return()

    if param.has_key('-fa')==True:
        if param.has_key('-frac')==True: frac=int(param['-frac'])
        else: frac=1
        BuildScoreMatrix(param['-sam'], param['-fa'],frac, True)
        return()
    if param.has_key('-phase')==True:
        BuildIndexFromFindHomology(param['-phase'], param['-o'])
        return ()

    if param.has_key('-mapq'):
        if param.has_key('-mismatches')==True:
            perc_levels=None
            mismatch_level=float( param['-mismatches'])
        else:
            mismatch_level=3.
            perc_levels=float( param['-perc_id'])
        min_ident=float(param['min_ident'])/100.
        ext=param['-tag']
        mismatch_penalty=float(param['-penalty'])
        UpdateMapq(param['-mapq'], min_ident, ext, False, mismatch_penalty, perc_levels, mismatch_level)
        return()
    if param.has_key('-extract'):
        ExtractInformativeReads(param['-extract'], param['-fa'])
        return()
    if param.has_key('-build'):
        BuildCalibrationTable(param['-build'], param['-npy'])
        return()

##    if param.has_key('-blast')==True:
##
##        CalibrateMappingQualities(param['-sam'], param['-blast'])
##        return()
    if param.has_key('-blast')==True:

        BlastToScoreMatrix(param['-blast'], param['-matrix'], param['-name'])
        return()

if __name__ == '__main__':
    main(sys.argv)
