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


import sys
import numpy
##import pandas
import scipy
from Bio import SeqIO
import scipy.sparse
from matplotlib import pyplot
import seaborn
import csv
import pickle
import shutil
import os
from collections import Iterable
try:
    import General_Tools
except:
    pass

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



#------------------------------------------------------------------------------
#           Functions for preparing and building the alignment index
#------------------------------------------------------------------------------


def BuildIndexFromFindHomology(indir, outfile, filter_by_internal_rpt=True):
    in_root="/".join( indir.split('/')[:-1])
    phased_folder=indir+'_phased_communities'
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
    out_root='.'.join(outfile.split('.')[:-1])
    phased_files=os.listdir(phased_folder)

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

            if filter_by_internal_rpt==True and rpt_scores[i]>=1.25*med_score:
                print "Skipped: ", rpt_scores[i],seq
                continue
            phased_handle.write('>{0}\n'.format( seq))
            phased_handle.write('{0}\n'.format(sequences[ seq]))

    phased_handle.close()



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
            distance=len(query)+distance
        displacements.append(distance)

    #Check the displacements in the reverse complement
##    for s in seed_intersection_rc:
##        distance=query_rc.find(s)-subject.find(s)
##        if distance<0:
##            distance=len(query_rc)+distance
##        displacements_rc.append(distance)
##    print sorted( displacements)
    if len (displacements)>0:
        mode=scipy.stats.mode(displacements)
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
            SeqLen[rec.name]=str( rec.seq).upper()
        else:
            SeqLen[rec.name]=rec.seq
    handle.close()
    return SeqLen

def MatrixMapq(matrix, full=False):

    if matrix.shape[1]>2:
        AS=numpy.partition(matrix,-1,1)[:,-1]
        XS=numpy.partition(matrix,-2,1)[:,-2]
        extra_counts=( matrix==AS[:,None]).sum(1)
        extra_counts=1
        max_score=-6.*140/5
        score=MapConf(AS,XS,70)/extra_counts
    elif matrix.shape[1]==2:
        AS=numpy.partition(matrix,-1,1)[:,-1]
        XS=numpy.partition(matrix,-2,1)[:,-2]
        max_score=-6.*140/5
        score=MapConf(AS,XS,70)
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

def PrioritizeRepeats(mat):
    best_hit=numpy.array( numpy.argmax(mat[1:,:],1))
    score=numpy.array( MatrixMapq(mat[1:,:]))

##    max_score=numpy.max(score)
##    print max_score
##    print score
##    counts=Counter(best_hit[0])
    enrichment=[]
    for k in range( mat.shape[1]):
        good_ind= (score>=.9)
        bad_ind=(score<.9)*(score!=0)
        good_count=(best_hit[good_ind]==k).sum()+mat[0,k]
        bad_count=(best_hit[bad_ind]==k).sum()
        ratio=bad_count-good_count
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


def RemoveRepeats(mat, iterations=312):
    print mat.shape
    unique_reads=numpy.array( mat[0,:].todense())[0]
    bool_array=numpy.array( [True]*mat.shape[1])
    noncluster_ind=numpy.array( (mat[:,:].sum(1)==0).T)[0]
    print noncluster_ind.shape

    mat=mat[~noncluster_ind,:]
    mat=numpy.array(mat.todense())
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
    ranked_indices= PrioritizeRepeats(mat[:,bool_array])
    for i in range(min( iterations, len(bool_array))):
        max_score=0
        max_j=0
        for j in ranked_indices:
            if bool_array[j]==False: continue
            orig_val=numpy.copy( bool_array[j])
            original_score=sum(CountReads(mat[:,bool_array],scores[:,bool_array]))+unique_reads[bool_array].sum() #+ mat[0,bool_array].sum()
            if score_list==[]: score_list.append(original_score)
            print sum(CountReads(mat[:,bool_array],scores[:,bool_array]))#, mat[0,bool_array].sum()
            bool_array[j]=False
            if sum(bool_array)!=0:

                new_score=sum(CountReads(mat[:,bool_array],scores[:,bool_array]))+unique_reads[bool_array].sum()
            else:
                new_score=0

#             print new_score
            print "\t", sum(CountReads(mat[:,bool_array],scores[:,bool_array]))#, mat[0,bool_array].sum()
            if new_score>=original_score:
                terminate=True
                max_score=new_score
                max_ind=j
                break
            if new_score>=max_score:

                max_score=new_score
                max_ind=j

            bool_array[j]=True
        print i, max_score, original_score
        removed_list.append(max_ind)
        bool_array[max_ind]=False
        score_list.append(max_score)
    return removed_list, score_list
def CountReads(matrix,scores):
    if matrix.shape[1]>1:
        best_hit=numpy.array( numpy.argmax(matrix[1:,:],1))
    else:



        return [ScoreAlignments( matrix[1:].flatten()).sum()]
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
        count= ( ScoreAlignments( matrix[1:,:][best_hit==i,i]) * mapq [best_hit==i]).sum()
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


def ReadCommunityFile(infile):
    inhandle=open(infile, 'r')
    community_dict={}
    iteration=True
    val_list=[]
    key_list=[]
    for line in inhandle:
        if len(line.strip())==0: continue
        if line[0]=='\t':
            seq_names=[s.strip() for s in line[1:].split(',')]
            val_list.append(seq_names[:-1])
        else: key_list.append(line.split(',')[0])


    return dict(zip(key_list, val_list))

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
        if community=='Community 1': continue
        if community!='Community 37': continue
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
        community_matrix=alignment_matrix[:, indices]
        try:
            removed_indices,objective_function=RemoveRepeats(community_matrix, len(indices))
        except:
            continue

        pyplot.plot(objective_function,c='purple',alpha=.7)
        pyplot.savefig(outdir+'/{0}.png'.format(community))
        pyplot.close()
        objective_function=objective_function[1:]
        max_score=numpy.max(objective_function)
        max_ind=numpy.where(objective_function>=max_score)[0][-1]+1
        print max_ind

        mask_array=numpy.array( [True]*community_matrix.shape[1])
        mask_array[removed_indices[:max_ind]]=False

        comm_handle=open( outdir+'/{0}.fa'.format(community), 'w')

        for v in namelist[mask_array]:
            try:
                outhandle.write('>{0}\n'.format(v))
                outhandle.write( sequences[CleanName( v)]+'\n')
                comm_handle.write('>{0}\n'.format(v))
                comm_handle.write( sequences[CleanName(v)]+'\n')
            except:

                pass
        comm_handle.close()




    pass



def BlastToScoreMatrix(infile, outfile, name_pickle):
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
                if contig_dict.has_key(rpt):
                    if read_scores[rpt]/max_score<.0001: continue
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

def UpdateMapq(sam_file):
    sam_handle=open(sam_file, 'r')
    sam_table=csv.reader(sam_handle, delimiter='\t')
    sam_base='.'.join(sam_file.split('.')[:-1])
    outfile=sam_base+"_tmp.sam"
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
        for col in row[11:]:
            if col[:2]=='AS':
                AS=140.+float(col.split(':')[-1])*5/6.
            if col[:2]=='XS':
                XS=140.+float(col.split(':')[-1])*5/6.
        if XS!=-1:
            try:
                mapq=MapConf(AS,XS,70,passed= False)
            except:
                print AS
                print XS
                print jabber
            row[4]=mapq
        else:
            row[4]=1.
        outtable.writerow(row)
    sam_handle.close()
    outhandle.close()
    os.remove(sam_file)
    os.rename(outfile, sam_file)

def RecomputeBlastAlignmentScore(line, rlen):
    overhang=rlen -(line.qend-line.qstart+1)
    score=2*line.matches- 5*(line.gapopen)-3*(line.mismatches+ overhang+line.gapext)
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

def MapConf(AS,XS, read_length=70., passed=True):
    AS_mismatches=(read_length*2-AS)/5.
    XS_mismatches=(read_length*2-XS)/5.

    AS_ident=1.- AS_mismatches/read_length
    XS_ident=1.- XS_mismatches/read_length
    try:
        if AS_ident!=XS_ident:
            best_diff=abs( XS_ident-AS_ident)
            worst_diff=1.-XS_ident
            score=best_diff**2/worst_diff
        else:
            score=0
    except:
        best_diff=abs( XS_ident-AS_ident)
        worst_diff=1.-XS_ident
        score=best_diff**2/worst_diff
        score[numpy.isnan( score)]=0.
##    return score
##    print best_diff
    if passed==False:
        return score
    else:
        try:
            if score>=.07: return 40
            else: return .1
        except:
            score[score>=.07]=40
            score[score<.07]=.1
            return score
def ScoreAlignments(AS, read_length=70.):
    AS_mismatches= (read_length*2-AS)/5.
    AS_ident= AS_mismatches/read_length
    try:
        diff=max( (0, .25-AS_ident))/.1
    except:
        diff=(.2-AS_ident)/.1
        diff[diff<0]=0
    return 1-numpy.exp(-1* diff)

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

def BuildScoreMatrix(infile, outfile,frac=1):
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    outhandle=open( outfile+'.fa', 'w')
    est_mapq=[]
    prev_mapq=[]
    pos=[]
    row=intable.next()
    start=time.clock()
    count=1
    contig_dict={}
    col_count=0
    read_count=0

    uniquely_mapped={}
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
        if read_count%frac!=0: continue
        outhandle.write('>{0}\n'.format(read_count))
        outhandle.write('{0}\n'.format(new_line.seq))
    outbase=".".join(outfile.split(".")[:-1])
    outhandle_1=open( '{0}_headers.pickle'.format(outbase), 'wb')
    pickle.dump(contig_dict, outhandle_1)
    outhandle_1.close()

def SummarizeSam(infile, outfile):

    sam_handle=open(infile, 'r')
    sam_table=csv.reader(sam_handle, delimiter='\t')


    outhandle=open(outfile, 'w')
    outtable=csv.writer(outhandle, delimiter='\t')
    rpt_dict={}

    for row in sam_table:
        if row[0][0]=="@":

            continue

        if row[2]=='*':

            continue
        if rpt_dict.has_key(row[2])==False:
            rpt_dict[row[2]]=[0.,0.]
        mapq=float(row[4])
        if mapq>=.07:
            rpt_dict[row[2]][0]+=1.
        else:
            rpt_dict[row[2]][1]+=1.
    sam_handle.close()
    for key in sorted(rpt_dict.keys()):
        outtable.writerow([key, rpt_dict[key][0]+rpt_dict[key][1] , rpt_dict[key][0], rpt_dict[key][1]])
    outhandle.close()
def main(argv):
    param={}
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if param=={}: return()

    if param.has_key('-comm'):
        CurateRepeats(param['-fa'],param['-comm'], param['-blast'], param['-names'],param['-out'])
        return()

    if param.has_key('-fa')==True:
        if param.has_key('-frac')==True: frac=int(param['-frac'])
        else: frac=1
        BuildScoreMatrix(param['-sam'], param['-fa'],frac)
        return()
    if param.has_key('-phase')==True:
        BuildIndexFromFindHomology(param['-phase'], param['-o'])
        return ()
    if param.has_key('-summ'):
        SummarizeSam(param['-summ'],param['-o'])
        return()
    if param.has_key('-mapq'):
        UpdateMapq(param['-mapq'])
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
