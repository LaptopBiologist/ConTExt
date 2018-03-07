#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:  Compute a summary statistic related to the probability that a sequencing
#           did not originate from a different consensus sequence
#
# Author:      I am
#
# Created:     22/02/2018
# Copyright:   (c) I am 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy
from Bio import SeqIO
from Bio.Seq import Seq
import sklearn.neighbors
from matplotlib import pyplot
import seaborn
import vptree
import time

IUPAC_DICT={'W':['A','T'], 'S':['C','G'], 'M':['A','C'],\
 'K':['G','T'], 'R':['A','G'], 'Y':['C','T'],\
 'B':['C','G','T'], 'D':['A','G','T'], 'H':['A','C', 'T'],\
 'V':['A','C','G'], 'N':['A','G','C','T'], 'U':'T' }

def ReplaceAmbiguousNucleotides(sequence):
    """Replaces ambiguous nucleotides with one of the nt that could be represented."""
    seq_array=numpy.fromstring(sequence.upper(), '|S1')
    for ambiguous_nt in IUPAC_DICT.keys():
        amb_ind=numpy.where(seq_array==ambiguous_nt)[0]
        replacement_nt=numpy.random.choice(IUPAC_DICT[ambiguous_nt], size=len(amb_ind), replace=True)
        seq_array[amb_ind]=replacement_nt
    seq_array
    return ''.join(seq_array)

def UnambiguousFASTA(infile, outfile):
    seq=GetSeq(infile)
    outhandle=open(outfile, 'w')
    for key in seq:
        new_seq=ReplaceAmbiguousNucleotides(seq[key])
        outhandle.write('>{0}\n'.format(key))
        outhandle.write('{0}\n'.format(new_seq))
    outhandle.close()

class SeedIndex():
    def __init__(self,infile, k=15, prefix_length=6):
        """Create a structure for fast KNN-searches of kmers in the set of consensus
    sequences"""

        self.k=k
        self.count=0
        self.prefix_length=prefix_length
        self.consensus_sequences=self.GetSeq(infile)

        self.kmer_location_dict={}
        self.prefix_dictionary={}

        self.index_kmers(self.consensus_sequences, k)

        kmers=[K.upper() for K in self.kmer_location_dict.keys()]

        for K in kmers:
            prefix=K[:prefix_length]
            try: self.prefix_dictionary[prefix]+=[K[prefix_length:]]
            except: self.prefix_dictionary[prefix]=[ K[prefix_length:]]
##        return
        self.BallTrees={}
##
        for prefix in self.prefix_dictionary.keys():
            prefixed_kmers_as_vectors=[SequenceToInt(prefix+kmer) for kmer in self.prefix_dictionary[prefix]]
            self.BallTrees[prefix]=sklearn.neighbors.BallTree(prefixed_kmers_as_vectors, metric='hamming')
##            self.BallTrees[prefix]=vptree.VPTree(prefixed_kmers_as_vectors, metric=scipy.spatial.distance.hamming)

##        return tree_dict

    def query_read(self,read, seed_spacing=.8, qual=''):
        #Decompose the read into seeds
        #If seed spacing is less than 1, seeds overlap and every nucleotide
        #in the read is used
        read_len=len(read)
        seeds, step_size=self.decompose_read(read, seed_spacing)

        #Wrt the Hamming metric ind the 10 nearest neighbors of each seed in the
        #consensus index

        query_hits=[self.query_seed(seed) for seed in seeds ]
        hits,sequence_hits=zip(*query_hits)

        #Determine the location of each neighbor kmer
        location_dictionary={'pos':[], 'dist':[], 'seed':[], 'sequence':[]}
        for i,hit in enumerate( hits):

                location_dictionary['pos']+=hit['pos']
                location_dictionary['dist']+=hit['dist']
                location_dictionary['seed']+=[step_size[ i]]*len(hit['dist'])
                location_dictionary['sequence']+=[i]*len(hit['dist'])



        sort_ind=numpy.argsort(numpy.array(location_dictionary['pos']))

        location_dictionary['array']=numpy.array(\
                        [location_dictionary['pos'],\
                        location_dictionary['dist'],\
                        location_dictionary['seed'],\
                        location_dictionary['sequence']], )[:,sort_ind]
##                location_dictionary[seq][strand]['dist']=numpy.array(location_dictionary[seq][strand]['dist'])[sort_ind]
##                location_dictionary[seq][strand]['seed']=numpy.array(location_dictionary[seq][strand]['seed'])[sort_ind]
        if len(location_dictionary['pos'])>1:
            split_ind=self.ClusterByDistances(location_dictionary['array'][0,:], read_len)
            splits=numpy.split(location_dictionary['array'], split_ind, 1)

        #Organize the locations of these kmers into alignments

        return location_dictionary, splits, sequence_hits


    def ClusterByDistances(self, posList, distance_cutoff=30 ):
        """Taking a list of read positions sorted by chromosome position as input,
        this carries out single-linkage agglomerative clustering with a distance
        cutoff."""
        diffArray=numpy.diff(posList)   #Distance between adjacent read pairs on the chromosome arm
        indices=numpy.where(diffArray>distance_cutoff)[0]+1 #Indices where the distnace exceeds the distance cutoff
##        splitIndices=indices #Indices that will define cluster boundaries
        return indices

    def FindAlignments(self):
        pass
    def query_seed(self,query, neighbors=10):
        query_prefix=query[:self.prefix_length]
        query_as_vector=SequenceToInt(query)[None,:]
        try:
            edit_distances, indexes=self.BallTrees[query_prefix].query(query_as_vector, neighbors)
        except:
            dim, tree_size=self.BallTrees[query_prefix].data.shape
            edit_distances, indexes=self.BallTrees[query_prefix].query(query_as_vector, tree_size)

        edit_distances=edit_distances[0]


        #score matches based on quality scores

        #return locations
        hits=[query_prefix+self.prefix_dictionary[query_prefix][i] for i in indexes[0]]

        location_dictionary={'pos':[], 'dist':[]}
        for i,hit in enumerate( hits):
            location_dictionary['pos']+=self.kmer_location_dict[hit]
            location_dictionary['dist']+=[edit_distances[i]]*len(self.kmer_location_dict[hit])
##                    print edit_distances[i]


        return location_dictionary, hits

    def decompose_read(self, read, seed_spacing=.8):
        seed_length=self.k
        step_size=max(1, int(seed_spacing*seed_length))
        seed_starts=list( set( range(0, len(read)-seed_length, step_size)+[len(read)-seed_length]))
        seeds=[read[x:x+seed_length] for x in seed_starts]
        return seeds, seed_starts

    def GetSeq(self,ref):
        """Different from Tools.General_Tools.GetSeq function in that this returns
        Bio.Seq objects."""
        handle=open(ref, 'r')
        lib=SeqIO.parse(handle, 'fasta')
        SeqLen=''
        SeqLen_rc=''
        for rec in lib:
            SeqLen+=str(rec.seq).upper()

        handle.close()
        return  Seq (SeqLen)

    def index_kmers(self,sequence, k=15):

        """Decomposes sequence into kmers."""

        seq_len=len(sequence)
        indices= range(0,len(sequence)-k+1)

        for x in indices:
            self.count+=1
##            if self.count>100000: break
            kmer=str(sequence[x:x+k])
            try:
                self.kmer_location_dict[kmer].append(x)
            except:
                self.kmer_location_dict[kmer]=[x]


            #Add the reverse complement
            try:
                kmer_rc=str(sequence[x:x+k].reverse_complement())
                try:
                    self.kmer_location_dict[kmer_rc].append(-1*x)
                except:
                    self.kmer_location_dict[kmer_rc]=[-1*x]




            except:
                p=0

def FindAlignments( seed_alignments, read,quality, index, seed_len=20, phred=33):

    sort_ind=numpy.argsort(seed_alignments[1])
    if len(sort_ind)==1: return 100
    qual=numpy.array([ord(s)-phred for s in quality])

    read_array=numpy.fromstring(read, '|S1')
    mismatches=0
    for i in sort_ind:
        pos=seed_alignments[0, i]
        hit=index.consensus_sequences[abs(int(pos)):abs(int(pos))+20]
        if pos<0: hit=hit.reverse_complement()
        hit=numpy.fromstring(str(hit),'|S1')
        read_start=int(seed_alignments[2,i])
        read_slice=read_array[read_start:read_start+seed_len]
        qual_slice=qual[read_start:read_start+seed_len]
##        print qual_slice
##        print read_slice
##        print read_slice!=hit
        mismatches+=(qual_slice*(read_slice!=hit)) [read_slice!='N'].sum()
        read_array[read_start:read_start+seed_len]='N'
##    indel= (abs((numpy.diff(seed_alignments[0,:])-numpy.diff(seed_alignments[2,:]))))
##    indel[indel!=0]=1
##    indel=indel.sum()


    return mismatches+40*(read_array!='N').sum()/20.


    pass

def PlotLocDict(dictionary):
    pyplot.scatter(dictionary['pos'],dictionary['seed'], cmap='jet', c='r', zorder=2)

    for i in range(len(dictionary['pos'])):
        pyplot.plot([dictionary['pos'][i],dictionary['pos'][i]+20], [dictionary['seed'][i]]*2, c='blue',lw=5.-10*dictionary['dist'][i], zorder=1)
    pyplot.show()

def BuildIndex(infile, k=15):
    """Create a structure for fast KNN-searches of kmers in the set of consensus
    sequences"""
    seq=GetSeq(infile)
    k_dict={}
    for key in seq:
        k_dict.update(CountKMERS(seq[key], k))
    kmers=[K.upper() for K in k_dict.keys()]
    kmer_dict={}
    for K in kmers:
        try: kmer_dict[K[:6]]+=[SequenceToInt( K)]
        except: kmer_dict[K[:6]]=[SequenceToInt( K)]
    tree_dict={}
    for prefix in kmer_dict.keys():
        tree_dict[prefix]=sklearn.neighbors.BallTree(kmer_dict[prefix], metric='hamming')

    return tree_dict


def SequenceToInt(seq):
    seq_array=numpy.fromstring(seq, '|S1')
    int_array=numpy.ndarray((len(seq_array),))
    int_array.fill(0)
    nt_array=[ 'T', 'C', 'G']
    for i, nt in enumerate( nt_array):
        int_array[seq_array==nt]=i+1
    return int_array




def CountKMERS(sequence, k=10, unbiased=False ):
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
        if kmerSet.has_key(kmer)==False:
            kmerSet[kmer]=0
        kmerSet[kmer]+=1.

    return kmerSet
##        kmerSet.add(str(complement[x:x+k]))
    return list(kmerSet)

def RandomSequence(length, gc=.5):
    AT_prob=(1-gc)/2.
    GC_prob=gc/2.
    nt=['A','T','G', 'C']
    probs=[AT_prob, AT_prob, GC_prob,GC_prob]
    return ''.join( numpy.random.choice(nt, size=length, p=probs, replace=True))

def CompareBalltoVP(length=20,reps=500):
    data=[SequenceToInt( RandomSequence(20)) for i in range(reps) ]
    return data
    start=time.clock()
    BallTree=sklearn.neighbors.BallTree(data, metric='hamming')
    build_BallTree=time.clock()-start
    start=time.clock()
    VPTree=vptree.VPTree(numpy.array( data), scipy.spatial.distance.hamming)
    build_VPTree=time.clock()-start
    start=time.clock()
    for d in data:
        test=BallTree.query(d[None,:], 10)
    BallTree_QueryTime=(time.clock()-start)/reps
    start=time.clock()
    for d in data:
        test=VPTree.get_n_nearest_neighbors(d, 10)
    VPTree_QueryTime=(time.clock()-start)/reps
    print 'BallTree: Build time = {0}; Query time = {1}'.format(build_BallTree, BallTree_QueryTime)
    print 'VPTree: Build time = {0}; Query time = {1}'.format(build_VPTree, VPTree_QueryTime)



def main():
    pass

if __name__ == '__main__':
    main()
