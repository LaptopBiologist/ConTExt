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
import scipy
from Bio import SeqIO
from Bio.Seq import Seq
import sklearn.neighbors
from matplotlib import pyplot
import seaborn
import vptree
import time
try:
    import redis
    import persistentdict
    import redis_collections
except:
    print "Failed to import Redis and persistentdict; the SeedIndex must be stored internally."
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
    def __init__(self,infile, k=15, prefix_length=6, store_externally=True):
        """Create a structure for fast KNN-searches of kmers in the set of consensus
    sequences"""
        self.redis_persistence=store_externally
        self.k=k
        self.count=0
        self.prefix_length=prefix_length
        self.consensus_sequences=self.GetSeq(infile)
        if self.redis_persistence==True:
            self.kmer_location_dict=redis_collections.dicts.Dict(writeback=True)
            self.prefix_dictionary=redis_collections.dicts.Dict(writeback=True)

        else:
            self.kmer_location_dict={}
            self.prefix_dictionary={}

        for seqname in self.consensus_sequences:
##            if self.count>100000: print jabber
            self.index_kmers(seqname, k)

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

    def query_read(self,read, seed_spacing=.8):
        #Decompose the read into seeds
        #If seed spacing is less than 1, seeds overlap and every nucleotide
        #in the read is used
        seeds, step_size=self.decompose_read(read, seed_spacing)

        #Wrt the Hamming metric ind the 10 nearest neighbors of each seed in the
        #consensus index
        hits=[self.query_seed(seed) for seed in seeds ]

        #Determine the location of each neighbor kmer
        location_dictionary={}
        for i,hit in enumerate( hits):
            for seq_name in hit.keys():
                if location_dictionary.has_key(seq_name)==False:
                    location_dictionary[seq_name]={}
                for strand in hit[seq_name].keys():
                    if location_dictionary[seq_name].has_key(strand)==False:
                        location_dictionary[seq_name][strand]={'pos':[], 'dist':[], 'seed':[], 'sequence':[]}
                    location_dictionary[seq_name][strand]['pos']+=hit[seq_name][strand]['pos']
                    location_dictionary[seq_name][strand]['dist']+=hit[seq_name][strand]['dist']
                    location_dictionary[seq_name][strand]['seed']+=[step_size[ i]]*len(hit[seq_name][strand]['dist'])
                    location_dictionary[seq_name][strand]['sequence']+=[i]*len(hit[seq_name][strand]['dist'])
        for seq in location_dictionary.keys():
            for strand in location_dictionary[seq].keys():
                sort_ind=numpy.argsort(numpy.array(location_dictionary[seq][strand]['pos']))

                location_dictionary[seq][strand]['array']=numpy.array(\
                                [location_dictionary[seq][strand]['pos'],\
                                location_dictionary[seq][strand]['dist'],\
                                location_dictionary[seq][strand]['seed'],\
                                location_dictionary[seq][strand]['sequence']], )[:,sort_ind]
##                location_dictionary[seq][strand]['dist']=numpy.array(location_dictionary[seq][strand]['dist'])[sort_ind]
##                location_dictionary[seq][strand]['seed']=numpy.array(location_dictionary[seq][strand]['seed'])[sort_ind]
                if len(location_dictionary[seq][strand]['pos'])>1:
                    split_ind=self.ClusterByDistances(location_dictionary[seq][strand]['array'][0,:])
                    splits=numpy.split(location_dictionary[seq][strand]['array'], split_ind, 1)
        #Organize the locations of these kmers into alignments

        return location_dictionary

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
        hits=[self.prefix_dictionary[query_prefix][i] for i in indexes[0]]
        location_dictionary={}
        for i,hit in enumerate( hits):

            for seq_name in self.kmer_location_dict[hit].keys():
                if location_dictionary.has_key(seq_name)==False:
                    location_dictionary[seq_name]={}
                for strand in self.kmer_location_dict[hit][seq_name].keys():
                    if location_dictionary[seq_name].has_key(strand)==False:
                        location_dictionary[seq_name][strand]={'pos':[], 'dist':[]}
                    location_dictionary[seq_name][strand]['pos']+=self.kmer_location_dict[hit][seq_name][strand]
                    location_dictionary[seq_name][strand]['dist']+=[edit_distances[i]]*len(self.kmer_location_dict[hit][seq_name][strand])
##                    print edit_distances[i]


        return location_dictionary

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
        SeqLen={}
        for rec in lib:
            SeqLen[rec.name]=rec.seq
        handle.close()
        return SeqLen

    def index_kmers(self,seqname, k=15):

        """Decomposes sequence into kmers."""
        sequence=self.consensus_sequences[seqname]
        seq_len=len(sequence)
        indices= range(0,len(sequence)-k+1)

        for x in indices:
            self.count+=1
##            if self.count>100000: break
            kmer=str(sequence[x:x+k])
            try:
                test_query=self.kmer_location_dict[kmer]
            except:
                self.kmer_location_dict[kmer]={}
            try: self.kmer_location_dict[kmer][seqname]['+'].append(x)
            except:
                    self.kmer_location_dict[kmer][seqname]={}
                    self.kmer_location_dict[kmer][seqname]['+']=[x]


            #Add the reverse complement
##            try:
##                kmer_rc=str(sequence[x:x+k].reverse_complement())
##                try:
##                    test_query=self.kmer_location_dict[kmer_rc]
##                except:
##                    self.kmer_location_dict[kmer_rc]={}
##
##                try: self.kmer_location_dict[kmer_rc][seqname]['-'].append(x)
##                except:
##                        self.kmer_location_dict[kmer_rc][seqname]={}
##                        self.kmer_location_dict[kmer_rc][seqname]['-']=[x]
##
##
##            except:
##                p=0

def FindAlignments( seed_alignments,read, qual):


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
