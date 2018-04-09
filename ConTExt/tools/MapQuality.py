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

import sys
import ParseCigar
import numpy
from Bio import SeqIO
from Bio.Seq import Seq
import sklearn.neighbors
from matplotlib import pyplot
import seaborn
try:
    import edlib
except:
    print ("Failed to load edlib. Alignments will not be possible.")
import time
import csv

import bz2file

IUPAC_DICT={'W':['A','T'], 'S':['C','G'], 'M':['A','C'],\
 'K':['G','T'], 'R':['A','G'], 'Y':['C','T'],\
 'B':['C','G','T'], 'D':['A','G','T'], 'H':['A','C', 'T'],\
 'V':['A','C','G'], 'N':['A','G','C','T'], 'U':['T'] }

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

def GetSeq(ref):
    handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    SeqLen={}
    for rec in lib:
##        if refName.count(rec.name)==0: continue
        SeqLen[CleanName(rec.name)]=str( rec.seq).upper()
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

class SeedIndex():
    def __init__(self,infile, k=15, prefix_length=6, phred=33, phred_base=10.,method='BallTree'):
        """Create a structure for fast KNN-searches of kmers in the set of consensus
    sequences"""
        self.method=method
        self.phred=phred
        self.k=k
        self.distance_cutoff=60
        self.neighbors=10
        self.count=0
        self.prefix_length=prefix_length
        self.consensus_sequences=self.GetSeq(infile)

        self.phred_const=phred_base/numpy.log(phred_base)
        self.kmer_location_dict={}
        self.prefix_dictionary={}

        self.index_kmers(self.consensus_sequences, k)

        kmers=[K.upper() for K in self.kmer_location_dict]

        self.seed_query_time=0.
        self.search_query_time=0.
        self.organize_time=0.
        self.read_count=0.
        self.BallTree_time=0.
        self.organize_seeds_time=0.
        self.recall_hits_time=0.
        self.zip_time=0

        for K in kmers:
            prefix=K[:prefix_length]
            try: self.prefix_dictionary[prefix]+=[K[prefix_length:]]
            except: self.prefix_dictionary[prefix]=[ K[prefix_length:]]
##        return
        self.BallTrees={}

        for prefix in self.prefix_dictionary:
            prefixed_kmers_as_vectors=[SequenceToInt(kmer) for kmer in self.prefix_dictionary[prefix]]
            if self.method=='BallTree':
                self.BallTrees[prefix]=sklearn.neighbors.BallTree(prefixed_kmers_as_vectors, metric='hamming')
            else:
                self.BallTrees[prefix]=numpy.array( prefixed_kmers_as_vectors)
##            self.BallTrees[prefix]=vptree.VPTree(prefixed_kmers_as_vectors, metric=scipy.spatial.distance.hamming)

##        return tree_dict

    def query_read(self,read, seed_spacing=.8, qual=''):
        #Decompose the read into seeds
        #If seed spacing is less than 1, seeds overlap and every nucleotide
        #in the read is used
        read_len=len(read)
        seeds, step_size=self.decompose_read(read, seed_spacing)

        if qual!='':
            qual_as_array=numpy.array([ord(s)-self.phred for s in qual])
            qual_as_array[qual_as_array<0]=0
        else: qual_as_array=''

        qual_slices=[qual_as_array[ step_size[i]:step_size[i]+self.k] for i in range(len(step_size))]
        #Wrt the Hamming metric ind the 10 nearest neighbors of each seed in the
        #consensus index
        query_hits=[]
        start=time.clock()
        query_hits=[self.query_seed(seeds[i], qual=qual_slices[i]) for i in range(len(seeds))]
        self.seed_query_time+=time.clock()-start
        start=time.clock()
        hits,sequence_hits=zip(*query_hits)
        self.zip_time=time.clock()-start
        #Determine the location of each neighbor kmer
        location_dictionary={'pos':[], 'dist':[], 'seed':[], 'sequence':[]}
        start=time.clock()
        for i,hit in enumerate( hits):

                location_dictionary['pos']+=hit['pos']
                location_dictionary['dist']+=hit['dist']
                location_dictionary['seed']+=[step_size[ i]]*len(hit['dist'])
                location_dictionary['sequence']+=[i]*len(hit['dist'])



        sort_ind=numpy.argsort(numpy.array(location_dictionary['pos']))

        location_dictionary['array']=numpy.array(\
                        [location_dictionary['pos' ],\
                        location_dictionary['dist'],\
                        location_dictionary['seed'],\
                        location_dictionary['sequence']], )[:,sort_ind]
##                location_dictionary[seq][strand]['dist']=numpy.array(location_dictionary[seq][strand]['dist'])[sort_ind]
##                location_dictionary[seq][strand]['seed']=numpy.array(location_dictionary[seq][strand]['seed'])[sort_ind]
        self.organize_time+=time.clock()-start
        start=time.clock()
        if len(location_dictionary['pos'])>1:
            split_ind=self.ClusterByDistances(location_dictionary['array'][0,:], read_len)
            splits=numpy.split(location_dictionary['array'], split_ind, 1)
        self.search_query_time+=time.clock()-start
        #Organize the locations of these kmers into alignments
        self.read_count+=1.
        return location_dictionary, splits, sequence_hits


    def ClusterByDistances(self, posList, distance_cutoff=60 ):
        """Taking a list of read positions sorted by chromosome position as input,
        this carries out single-linkage agglomerative clustering with a distance
        cutoff."""

        diffArray=numpy.diff(posList)   #Distance between adjacent read pairs on the chromosome arm
        indices=numpy.where(diffArray>self.distance_cutoff)[0]+1 #Indices where the distnace exceeds the distance cutoff
##        splitIndices=indices #Indices that will define cluster boundaries
        return indices

    def FindAlignments(self):
        pass
    def query_seed(self,query,qual=''):
        neighbors=self.neighbors
        start=time.clock()
        query_prefix=query[:self.prefix_length]
        query_as_vector=SequenceToInt(query[self.prefix_length:])[None,:]
        if self.method=="BallTree":
            try:
                edit_distances, indexes=self.BallTrees[query_prefix].query(query_as_vector, neighbors)
            except:
                dim, tree_size=self.BallTrees[query_prefix].data.shape
                edit_distances, indexes=self.BallTrees[query_prefix].query(query_as_vector, tree_size)
            hits=[query_prefix+self.prefix_dictionary[query_prefix][i] for i in indexes[0]]
            edit_distances=edit_distances[0]
        else:
            edit_distances=1.-((query_as_vector==self.BallTrees[query_prefix]).sum(1)/float( self.k))

            sorted_ind=numpy.argsort(edit_distances)
            if len(sorted_ind)>neighbors:
                indexes=sorted_ind[:10]
            else:
                indexes=sorted_ind
            hits=[query_prefix+self.prefix_dictionary[query_prefix][i] for i in indexes]
        self.BallTree_time+=time.clock()- start
        #return locations

        start=time.clock()

        self.recall_hits_time=time.clock()-start

        start=time.clock()
        if qual!='':
        #score matches based on quality scores
            query_as_array=numpy.fromstring(query, '|S1')

            weight=numpy.exp(-1* qual/self.phred_const)

            error_prob=numpy.array([0.]*len(query), float)
            for index, hit in enumerate(hits):
                hit_as_array=numpy.fromstring(hit, '|S1')
                #Computes a weighted distance
##                edit_distances[index]= ((query_as_array!=hit_as_array)*(1.-weight)).sum()/self.k
                #Could instead be expressed as a probability
                matching_positions=query_as_array==hit_as_array
                edit_distances[index]=numpy.log( 1.-weight[matching_positions]).sum()\
                                      +numpy.log( weight[~matching_positions]/3.).sum()


        location_dictionary={'pos':[], 'dist':[]}
        for i,hit in enumerate( hits):
            location_dictionary['pos']+=self.kmer_location_dict[hit]
            location_dictionary['dist']+=[edit_distances[i]]*len(self.kmer_location_dict[hit])
##                    print edit_distances[i]

        self.organize_seeds_time=time.clock()-start
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
    def reset_time(self):
        self.search_query_time=0
        self.seed_query_time=0
        self.organize_time=0
    def report_performance(self):
        print ("Spent {0} s on seed searches ({1})".format(self.seed_query_time, self.seed_query_time/self.read_count))
        print ("\tSpent {0} s on BallTrees ({1})".format(self.BallTree_time, self.BallTree_time/self.read_count))
        print ("\tSpent {0} s determine hits ({1})".format(self.recall_hits_time, self.recall_hits_time/self.read_count))
        print ("\tSpent {0} s organizing seeds ({1})".format(self.organize_seeds_time, self.organize_seeds_time/self.read_count))
        print ("Spent {0} s unzipping seeds ({1})".format(self.zip_time, self.zip_time/self.read_count))
        print ("Spent {0} s on organizing searches ({1})".format(self.organize_time, self.organize_time/self.read_count))
        print ("Spent {0} s on alignment building ({1})".format(self.search_query_time, self.search_query_time/self.read_count))


class HashedIndex():
    def __init__(self, k=10,prefix_length=5, phred=33, phred_base=10.,method='BallTree'):
        """Create a structure for fast KNN-searches of kmers in the set of consensus
    sequences"""
        self.method=method
        self.phred=phred
        self.prefix_length=prefix_length
        self.k=k
        self.distance_cutoff=60
        self.neighbors=10
        self.count=0
        self.phred_base=phred_base
        self.kmer_location_dict={}
        self.phred_base=phred_base

    def BuildIndex(self, infile):

        self.consensus_sequences=self.GetSeq(infile)
        phred_base=self.phred_base
        self.phred_const=phred_base/numpy.log(phred_base)

        self.prefix_dictionary={}
        self.suffix_dictionary={}

        self.index_kmers(self.consensus_sequences, self.k)

        kmers=[K.upper() for K in self.kmer_location_dict]

        self.seed_query_time=0.
        self.search_query_time=0.
        self.organize_time=0.
        self.read_count=0.
        self.BallTree_time=0.
        self.organize_seeds_time=0.
        self.recall_hits_time=0.
        self.zip_time=0
        prefix_length=self.prefix_length
        for K in kmers:
            prefix=K[:prefix_length]
            try: self.prefix_dictionary[prefix].append( K[prefix_length:])
            except: self.prefix_dictionary[prefix]=[ K[prefix_length:]]
            suffix=K[prefix_length:]
            try: self.suffix_dictionary[suffix].append( K[:prefix_length])
            except: self.suffix_dictionary[suffix]=[ K[:prefix_length]]

        self.BallTrees_prefix={}

        for prefix in self.prefix_dictionary:
            prefixed_kmers_as_vectors=[SequenceToInt(kmer) for kmer in self.prefix_dictionary[prefix]]
            self.BallTrees_prefix[prefix]=numpy.array( prefixed_kmers_as_vectors)

        self.fuzzy_location_dict={}
        start=time.clock()
        #Add all kmers within 1-edit distance disallowing mismatches the in the prefix
        for kmer in self.kmer_location_dict.keys():
            fuzzy_matches=self.GetHitsWithinHammingDistance(kmer)
            if self.fuzzy_location_dict.has_key(kmer)==False:
                self.fuzzy_location_dict[kmer]=[]
            for hit in fuzzy_matches:
                self.fuzzy_location_dict[kmer]+=self.kmer_location_dict[hit]
        del self.BallTrees_prefix
        self.BallTrees_suffix={}

        for suffix in self.suffix_dictionary:
            suffixed_kmers_as_vectors=[SequenceToInt(kmer) for kmer in self.suffix_dictionary[suffix]]
            self.BallTrees_suffix[suffix]=numpy.array( suffixed_kmers_as_vectors)

        #Add all kmers within 1-edit distance disallowing mismatches the in the suffix
        for kmer in self.kmer_location_dict.keys():
            fuzzy_matches=self.GetHitsWithinHammingDistanceSuffix(kmer)
            if self.fuzzy_location_dict.has_key(kmer)==False:
                self.fuzzy_location_dict[kmer]=[]
            for hit in fuzzy_matches:
                if hit==kmer: continue
                self.fuzzy_location_dict[kmer]+=self.kmer_location_dict[hit]

        self.kmer_location_dict=self.fuzzy_location_dict
        print time.clock()-start


    def GetHitsWithinHammingDistance(self,query):
        query_prefix=query[:self.prefix_length]
        query_as_vector=SequenceToInt(query[self.prefix_length:])[None,:]
        edit_distances=(query_as_vector!=self.BallTrees_prefix[query_prefix]).sum(1)
        matches=numpy.where( edit_distances<=1)[0]

        hits=[query_prefix+self.prefix_dictionary[query_prefix][i] for i in matches]
        return hits

    def GetHitsWithinHammingDistanceSuffix(self,query):
        query_suffix=query[self.prefix_length:]
        query_as_vector=SequenceToInt(query[:self.prefix_length])[None,:]
        edit_distances=(query_as_vector!=self.BallTrees_suffix[query_suffix]).sum(1)
        matches=numpy.where( edit_distances<=1)[0]

        hits=[self.suffix_dictionary[query_suffix][i]+query_suffix for i in matches]
        return hits


    def query_read(self,read, seed_spacing=.8, qual=''):
        #Decompose the read into seeds
        #If seed spacing is less than 1, seeds overlap and every nucleotide
        #in the read is used
        read_len=len(read)
        seeds, step_size=self.decompose_read(read, seed_spacing)

        if qual!='':
            qual_as_array=numpy.array([ord(s)-self.phred for s in qual])
            qual_as_array[qual_as_array<0]=0
        else: qual_as_array=''

        qual_slices=[qual_as_array[ step_size[i]:step_size[i]+self.k] for i in range(len(step_size))]
        #Wrt the Hamming metric ind the 10 nearest neighbors of each seed in the
        #consensus index
        query_hits=[]
        start=time.clock()
        query_hits=[self.query_seed(seeds[i], qual=qual_slices[i]) for i in range(len(seeds))]
        self.seed_query_time+=time.clock()-start
        start=time.clock()
        hits,sequence_hits=zip(*query_hits)
        self.zip_time=time.clock()-start
        #Determine the location of each neighbor kmer
        location_dictionary={'pos':[], 'dist':[], 'seed':[], 'sequence':[]}
        start=time.clock()
        for i,hit in enumerate( hits):

                location_dictionary['pos']+=hit['pos']
                location_dictionary['dist']+=hit['dist']
                location_dictionary['seed']+=[step_size[ i]]*len(hit['dist'])
                location_dictionary['sequence']+=[i]*len(hit['dist'])



        sort_ind=numpy.argsort(numpy.array(location_dictionary['pos']))

        location_dictionary['array']=numpy.array(\
                        [location_dictionary['pos' ],\
                        location_dictionary['dist'],\
                        location_dictionary['seed'],\
                        location_dictionary['sequence']], )[:,sort_ind]
##                location_dictionary[seq][strand]['dist']=numpy.array(location_dictionary[seq][strand]['dist'])[sort_ind]
##                location_dictionary[seq][strand]['seed']=numpy.array(location_dictionary[seq][strand]['seed'])[sort_ind]
        self.organize_time+=time.clock()-start
        start=time.clock()
        if len(location_dictionary['pos'])>1:
            split_ind=self.ClusterByDistances(location_dictionary['array'][0,:], read_len)
            splits=numpy.split(location_dictionary['array'], split_ind, 1)
        self.search_query_time+=time.clock()-start
        #Organize the locations of these kmers into alignments
        self.read_count+=1.
        return location_dictionary, splits, sequence_hits


    def ClusterByDistances(self, posList, distance_cutoff=60 ):
        """Taking a list of read positions sorted by chromosome position as input,
        this carries out single-linkage agglomerative clustering with a distance
        cutoff."""

        diffArray=numpy.diff(posList)   #Distance between adjacent read pairs on the chromosome arm
        indices=numpy.where(diffArray>self.distance_cutoff)[0]+1 #Indices where the distnace exceeds the distance cutoff
##        splitIndices=indices #Indices that will define cluster boundaries
        return indices

    def FindAlignments(self):
        pass
    def query_seed(self,query,qual=''):
        neighbors=self.neighbors
        start=time.clock()
        query_prefix=query[:self.prefix_length]
        query_as_vector=SequenceToInt(query[self.prefix_length:])[None,:]
        if self.method=="BallTree":
            try:
                edit_distances, indexes=self.BallTrees[query_prefix].query(query_as_vector, neighbors)
            except:
                dim, tree_size=self.BallTrees[query_prefix].data.shape
                edit_distances, indexes=self.BallTrees[query_prefix].query(query_as_vector, tree_size)
            hits=[query_prefix+self.prefix_dictionary[query_prefix][i] for i in indexes[0]]
            edit_distances=edit_distances[0]
        else:
            edit_distances=1.-((query_as_vector==self.BallTrees[query_prefix]).sum(1)/float( self.k))

            sorted_ind=numpy.argsort(edit_distances)
            if len(sorted_ind)>neighbors:
                indexes=sorted_ind[:10]
            else:
                indexes=sorted_ind
            hits=[query_prefix+self.prefix_dictionary[query_prefix][i] for i in indexes]
        self.BallTree_time+=time.clock()- start
        #return locations

        start=time.clock()

        self.recall_hits_time=time.clock()-start

        start=time.clock()
        if qual!='':
        #score matches based on quality scores
            query_as_array=numpy.fromstring(query, '|S1')

            weight=numpy.exp(-1* qual/self.phred_const)

            error_prob=numpy.array([0.]*len(query), float)
            for index, hit in enumerate(hits):
                hit_as_array=numpy.fromstring(hit, '|S1')
                #Computes a weighted distance
##                edit_distances[index]= ((query_as_array!=hit_as_array)*(1.-weight)).sum()/self.k
                #Could instead be expressed as a probability
                matching_positions=query_as_array==hit_as_array
                edit_distances[index]=numpy.log( 1.-weight[matching_positions]).sum()\
                                      +numpy.log( weight[~matching_positions]/3.).sum()


        location_dictionary={'pos':[], 'dist':[]}
        for i,hit in enumerate( hits):
            location_dictionary['pos']+=self.kmer_location_dict[hit]
            location_dictionary['dist']+=[edit_distances[i]]*len(self.kmer_location_dict[hit])
##                    print edit_distances[i]

        self.organize_seeds_time=time.clock()-start
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
        name_list=[]
        left_list=[]
        right_list=[]
        offset=0
        position=0
        offset=200
        for rec in lib:
            length=len(rec)
            left_list.append(position)
            right_list.append(position+length)
            name_list.append(rec.name)
            new_seq=str(rec.seq).upper()+'N'*offset
            SeqLen+=new_seq
            position+=len(new_seq)

        handle.close()
        self.left=numpy.array(left_list)
        self.right=numpy.array(right_list)
        self.name_list=numpy.array(name_list)
        return  Seq (SeqLen)

    def QueryReference(self,pos):
        strand=numpy.sign(pos)
        pos=abs(pos)
        hit_ind=(pos>=self.left)*(pos<self.right)
        name=self.name_list[hit_ind]
        hit_pos=pos-self.left[hit_ind]
        if len(name)!=0:
            name=name[0]
            hit=hit_pos[0]
        else:
            name=''
            hit_pos=-1
        return name, hit_pos


    def index_kmers(self,sequence, k=15):

        """Decomposes sequence into kmers."""

        seq_len=len(sequence)
        indices= range(0,len(sequence)-k+1)

        for x in indices:
            self.count+=1

            kmer=str(sequence[x:x+k])
            if kmer.count('N')!=0: continue
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
    def reset_time(self):
        self.search_query_time=0
        self.seed_query_time=0
        self.organize_time=0
    def report_performance(self):
        print ("Spent {0} s on seed searches ({1})".format(self.seed_query_time, self.seed_query_time/self.read_count))
        print ("\tSpent {0} s on BallTrees ({1})".format(self.BallTree_time, self.BallTree_time/self.read_count))
        print ("\tSpent {0} s determine hits ({1})".format(self.recall_hits_time, self.recall_hits_time/self.read_count))
        print ("\tSpent {0} s organizing seeds ({1})".format(self.organize_seeds_time, self.organize_seeds_time/self.read_count))
        print ("Spent {0} s unzipping seeds ({1})".format(self.zip_time, self.zip_time/self.read_count))
        print ("Spent {0} s on organizing searches ({1})".format(self.organize_time, self.organize_time/self.read_count))
        print ("Spent {0} s on alignment building ({1})".format(self.search_query_time, self.search_query_time/self.read_count))


def CleanKmerLocations(hi, cutoff=50):
    for key in hi.kmer_location_dict.keys():
        locations=hi.kmer_location_dict[key]
        sorted_loc=sorted(locations)
        diff=numpy.diff([min(locations)-1]+  hi.kmer_location_dict[key])
        good_seeds=diff>cutoff
##    return sorted_loc
        hi.kmer_location_dict[key]=list(numpy.array(sorted_loc)[good_seeds])



def AlignRead(read, subject, k):
    alignments=[]
    all_found=False
    while all_found==False:
        alignment=edlib.align(read, subject, 'HW', 'path', k=k)
        locations=alignment['locations']
        alignments.append(alignment)
        if len(locations)>1:
            l,r=locations[0]
            subject=subject[:l]+'N'*(r-l)+subject[r:]
        else:
            all_found=True
        if float(subject.count('N'))/len(subject)>.7:
            all_found=True
    return alignments

def ScoreRead(read, qual, has_index,exp_seq, exp_pos, phred=33):
    read=read.upper()
    start=time.clock()
    subject_sequences=GetPotentialTargetSequences(read, has_index)
    if subject_sequences==[]: return 0

    seed_time=time.clock()-start
    quality_array=numpy.array([ord( q)-phred for q in qual])
    #Align the read to each target sequence and  save the alignment score
    alignment_scores=[]
##    print subject_sequences
    align_time=0
    parse_time=0
    score_time=0
    start=time.clock()
    aligned_subjs=[]
    alignments=[]

    for subject in subject_sequences:

        new_alignments=AlignRead(read, subject[2],k=20)
        alignments+=new_alignments
        aligned_subjs+=[subject]*len(new_alignments)
    distances=numpy.array( [a['editDistance'] for a in alignments])
##    print distances
##    print alignments
##    if len(distances[distances>0])>3:
##        max_dist=sorted(distances[distances>0])[2]
    max_dist=max(distances)
    align_time+=time.clock()-start
    exp_score=-1
    best_score=-200
    best_cigar=''
    best_seq, best_loc='',0
    for i, alignment in enumerate( alignments):
        l,r,subject=aligned_subjs[i]
        if alignment['editDistance']>max_dist:continue
        if alignment['editDistance']==-1: continue

        cigar_string=alignment['cigar']
##        print cigar_string
        loc=alignment['locations'][0]

        start=time.clock()
        parse_time+=time.clock()-start

        start=time.clock()
##        try:
##        print read
##        print subject
        alignment_score=ScoreAlignment(quality_array, cigar_string)


##        print alignment_score
##        print alignment_score
        alignment_scores.append(alignment_score)
        subj_seq= has_index.QueryReference(l)[0]
        if subj_seq=='': subj_seq= has_index.QueryReference(r)[0]
        subj_l= has_index.QueryReference(abs(l))[1]
        if subj_l==-1:subj_l=0
        subj_r= has_index.QueryReference(abs(r))[1]
        if subj_r==-1: subj_r=subj_l+200
        if subj_seq==exp_seq and exp_pos<subj_r and subj_l<=exp_pos:
            exp_score=alignment_score
##            best_cigar=cigar_string
##        except:
##            continue
        if alignment_score>best_score:
            best_score=alignment_score
            best_cigar=cigar_string
            best_loc=subj_l
            best_seq=subj_seq
        score_time+=time.clock()-start


    if exp_score!=-1:
        if len(alignment_scores)>1:
            linear_probs= numpy.exp( alignment_scores)
            very_small_probs=numpy.isnan(linear_probs)+numpy.isinf(linear_probs)
            mapq=-4.34*numpy.log(1.-( numpy.exp( exp_score))/numpy.nansum( linear_probs[~very_small_probs]))*.33
            mapq/=52.61444578/40.
            if numpy.isnan(mapq) or numpy.isinf(mapq):
                mapq=40
        else:
            mapq=40
    else:
        mapq=0

##        print 1.-( numpy.exp( exp_score)/(numpy.exp( alignment_scores).sum()))
##        print alignment_scores
##    print alignment_scores
##    print mapq
    return mapq, best_seq, best_loc, best_cigar


def GetPotentialTargetSequences(read, hash_index):
    begin=time.clock()
    start=time.clock()
    read_length=len(read)
    seeds=hash_index.decompose_read(read, .4)
##    print seeds
    pos=[]
    start=time.clock()
    seed_count=float( len(seeds)+1.)

    for i in range(len(seeds[0])):
        try:  #If the seed is present in the reference
            hit_pos= hash_index.kmer_location_dict[seeds[0][i]]
            pos+=hit_pos
        except:  #If it is absent
            continue
    print time.clock()-start, 'Build'

    start=time.clock()

    pos=numpy.sort(pos)
##    pos=numpy.array(pos)[sort_ind]
    singletons=FindSingletons(pos,len(read))
    pos=pos[singletons]
    print time.clock()-start, 'Sort'
    start=time.clock()
##    split_indices=numpy.concatenate(([0], ClusterByDistances(pos, len(read)), [len(pos)]))
    split_indices=ParseCigar.ClusterByDistances(pos, 60)
    print time.clock()-start, 'Cluster'
    start=time.clock()

    grouped_positions=[]
    diff=numpy.diff(split_indices)

    size_set=sorted(diff)[::-1]

    if len(size_set)>2:
        min_size=max(2,min( (size_set[1]-1,size_set[1]*.5)))
    else: min_size=2

    accept_ind=numpy.where(diff>=min_size)[0]
    accept_ind=sorted(accept_ind)
    if len(accept_ind)>5:
        accept_ind=accept_ind[-5:]

    for ind in accept_ind:

        grouped_positions.append(pos[split_indices[ind]:split_indices[ind]+diff[ind]])

    mean_pos=[int(numpy.mean(g)) for g in grouped_positions]
    strand=numpy.sign(mean_pos)
    subject_sequences=[]
    print time.clock()-start, 'Split'
    start=time.clock()
    for i in range(len(mean_pos)):
        left,right=abs(mean_pos[i])-read_length,abs(mean_pos[i])+read_length
        if strand[i]==1:
            subject=str(hash_index.consensus_sequences[left-50:right+50])
        else: #Hit is on the reverse strand
            subject=str(hash_index.consensus_sequences[left-50:right+50].reverse_complement())
        subject_sequences.append((left,right,subject))
    print time.clock()-start, 'Extract'
    print time.clock()-begin, 'Total'
    return subject_sequences


def GetPotentialTargetSequences_svaed(read, hash_index):
    start=time.clock()
    read_length=len(read)
    seeds=hash_index.decompose_read(read, .4)
##    print seeds
    pos,seed_id=[],[]
    start=time.clock()

    for i in range(len(seeds[0])):
        try:
            hit_pos= hash_index.kmer_location_dict[seeds[0][i]]
            pos+=hit_pos
            seed_id+=[i]*len(hit_pos)
        except:
##            print seeds[0][i]
            continue
##    print time.clock()-start

    start=time.clock()
    sort_ind=numpy.argsort(pos)
    pos=numpy.array(sorted( pos))
    print len(pos)
    sorted_seeds=numpy.array(seed_id)[sort_ind]
##    print time.clock()-start
##    return pos

    singletons=FindSingletons(pos,len(read))
    pos=pos[singletons]
    print len(pos)
    sorted_seeds=sorted_seeds[singletons]
##    print pos


##    print test
##    print len(pos)


    start=time.clock()
    split_indices=numpy.concatenate(([0], ClusterByDistances(pos, len(read)), [len(pos)]))
    print split_indices
##    print time.clock()-start
##    split_count=
##    return split_indices, pos, sorted_seeds
    start=time.clock()
    grouped_positions=numpy.split(pos, split_indices)
    print grouped_positions
    grouped_seeds=numpy.split(sorted_seeds, split_indices)
##    print time.clock()-start
    group_sizes=numpy.array( [len(set( g)) for g in grouped_seeds])


    #exclude clusters with few counts
    start=time.clock()
    size_set=sorted(group_sizes)[::-1]

    if len(size_set)>2:
        min_size=max(2,min( (size_set[1]-1,size_set[1]*.5)))
    else: min_size=2
##    min_size=2
    start=time.clock()
    keep_indices=numpy.where(group_sizes>=min_size)[0]
##    if len(keep_indices)>5:
##        keep_indices=numpy.argsort(group_sizes)[-5:]

    grouped_positions=[grouped_positions[i] for i in keep_indices ]
    grouped_seeds=[grouped_seeds[i] for i in keep_indices ]

    mean_pos=[int(numpy.mean(g)) for g in grouped_positions]
    strand=numpy.sign(mean_pos)
    subject_sequences=[]

##    start=time.clock()
    for i in range(len(mean_pos)):
        left,right=abs(mean_pos[i])-read_length,abs(mean_pos[i])+read_length
        if strand[i]==1:
            subject=str(hash_index.consensus_sequences[left-50:right+50])
        else: #Hit is on the reverse strand
            subject=str(hash_index.consensus_sequences[left-50:right+50].reverse_complement())
        subject_sequences.append((left,right,subject))
##    print time.clock()-start
    return subject_sequences


def CythonScore(cigar, qual, phred=33):
    t=ParseCigar.ParseEdlibCigar(cigar)
##    print t
    s=ParseCigar.ComputeAlignmentScore(qual, numpy.array( t[0]),numpy.array( t[1], int), t[2], t[3], t[4])
    return s


def FindSingletons(pos, cutoff=60):
    test=abs(pos[1:]-pos[:-1])

    test=numpy.concatenate(([1000], test, [1000]))

    return numpy.minimum(test[1:], test[:-1])<cutoff

def ClusterByDistances( posList, distance_cutoff=60 ):
    """Taking a list of read positions sorted by chromosome position as input,
    this carries out single-linkage agglomerative clustering with a distance
    cutoff."""

    diffArray=numpy.diff(posList)   #Distance between adjacent read pairs on the chromosome arm
    indices=numpy.where(diffArray>distance_cutoff)[0]+1 #Indices where the distnace exceeds the distance cutoff

##    posList=posList[indices[:-1][ numpy.diff(indices)>1]]
##    diffArray=numpy.diff(posList)   #Distance between adjacent read pairs on the chromosome arm
##    indices=numpy.where(diffArray>distance_cutoff)[0]+1 #Indices where the distnace exceeds the distance cutoff
##        splitIndices=indices #Indices that will define cluster boundaries
    return indices
def StoreIndex(index, outfile):
    outhandle=open(outfile, 'w')
    outtable=csv.writer(outhandle, delimiter='\t')
    #Write the general parameters
    method=index.method
    phred=index.phred
    k=index.k
    prefix_length=index.prefix_length
    row=['@a:method', method]
    outtable.writerow(row)
    row=['@a:phred', phred]
    outtable.writerow(row)
    row=['@a:k', k]
    outtable.writerow(row)
    row=['@a:prefix', prefix_length]
    outtable.writerow(row)
    #Write the location dictionary
    for key in index.kmer_location_dict.keys():
        sorted_loc=sorted([str(k) for k in index.kmer_location_dict[key]])
        row=['@d:{0}'.format(key), ','.join(sorted_loc)]
        outtable.writerow(row)


    #Write consensus sequences
    cons_seq=index.consensus_sequences

    for i in range(0, len(cons_seq)-10000, 10000):
        seq_slice=str(cons_seq[i:i+10000])
        row=['@C:{0}'.format(i),seq_slice ]
        outtable.writerow(row)

    seq_slice=str(cons_seq[i+10000:])
    row=['@C:{0}'.format(i),seq_slice ]
    outtable.writerow(row)

    left_list=[str(l) for l in index.left]
    row=['@Q:L',','.join(left_list)]
    outtable.writerow(row)

    right_list=[str(r) for r in index.right]
    row=['@Q:R',','.join(right_list)]
    outtable.writerow(row)

    name_list=[str(n) for n in index.name_list]
    row=['@Q:N',','.join(name_list)]
    outtable.writerow(row)

    outhandle.close()

def LoadIndex(infile, index):
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')

    #Write the general parameters
    sequence=''
    for row in intable:
##        print row
        rtype,attribute= row[0].split(':')
        if rtype=='@a': #Attribute
            value=row[1]
            if attribute=='method': index.method=value
            if attribute=='phred': index.phred=int(value)
            if attribute=='k': index.k=int(value)
            if attribute=='prefix': index.prefix_length=int(value)
        if rtype=='@d': #Kmer dictionary
            key=attribute

            values=[int(v) for v in row[1].split(',')]

            index.kmer_location_dict[key]=values
        if rtype=='@C':#Consensus sequence
            sequence+=row[1]
        if rtype=='@Q':#Query informartion
            value=row[1]
            if attribute=='L':
                index.left=numpy.array([int(v) for v in value.split(',')])
            if attribute=='R':
                index.right=numpy.array([int(v) for v in value.split(',')])
            if attribute=='N':
                index.name_list=numpy.array(value.split(','), str)
    sequence=Seq(sequence)
    index.consensus_sequences=sequence

    inhandle.close()

def ScoreAlignment(quality_string, cigar,phred_base=33, phred_const=4.34):

    quality=quality_string
    start=time.clock()
    matching_positions, mismatching_positions, indel= ParseEdlibCigar(cigar)
##    print time.clock()-start
    start=time.clock()
    matching_positions= numpy.array( matching_positions[0])
    mismatching_positions= numpy.array(mismatching_positions[0])
    weight=numpy.exp(-1* quality/phred_const)
##    print weight
    match_error_prob=numpy.log( 1.-weight[matching_positions]).sum()
##    print '\t', weight
##    print '\t', numpy.log( 1.-weight[matching_positions])
    if len(mismatching_positions)>0:
        mismatch_error_prob=numpy.log( weight[mismatching_positions]/3.).sum()
    else:
        mismatch_error_prob=0

##    print '\t', numpy.log( weight[mismatching_positions]/3.)
##    print '\t', match_error_prob, mismatch_error_prob
    score=match_error_prob+mismatch_error_prob - indel
##    print time.clock()-start
##    print cigar, len(mismatching_positions), score
    return score

def ParseEdlibCigar(cigar):
    cigar_array=numpy.fromstring(cigar, '|S1')
    matches=cigar_array=='='
    mismatches=cigar_array=='X'
    insertions=cigar_array=='I'
    deletions=cigar_array=='D'
    operations=matches+mismatches+insertions+deletions
    operation_indices= numpy.where(operations==True)[0]
    last_index=0
    read_mismatch_locations=[]
    read_match_locations=[]
    subj_mismatch_locations=[]
    subj_match_locations=[]
    read_pointer=0
    subj_pointer=0
    indel=0
    operation_array=[]
    for current_index in operation_indices:
        number_of_edits=int(cigar[last_index:current_index])
        edit_type=cigar[current_index]
        operation_array+=[edit_type]*number_of_edits
##        if edit_type=='=':
##            read_match_locations+=range(read_pointer, read_pointer+number_of_edits)
##            subj_match_locations+=range(subj_pointer, subj_pointer+number_of_edits)
##            read_pointer+=number_of_edits
##            subj_pointer+=number_of_edits
##        if edit_type=='I':
##            indel+=2
##            read_pointer+=number_of_edits
##        if edit_type=='D':
##            indel+=2
##            subj_pointer+=number_of_edits
##        if edit_type=='X':
##            read_mismatch_locations+= range(read_pointer, read_pointer+number_of_edits)
##            subj_mismatch_locations+=range(subj_pointer, subj_pointer+number_of_edits)
##            read_pointer+=number_of_edits
##            subj_pointer+=number_of_edits


        last_index=current_index+1
##
##    return(numpy.array( read_match_locations) ,numpy.array(  subj_match_locations) ,numpy.array( read_mismatch_locations), numpy.array(subj_mismatch_locations), indel)
##
    operation_array=numpy.array(operation_array)
##
    read_mismatch_locations=numpy.where( operation_array[operation_array!='D']=='X')
    read_match_locations= numpy.where( operation_array[operation_array!='D']=='=')
    indel=(operation_array=='D').sum()+(operation_array=='I').sum()
    return read_match_locations, read_mismatch_locations, indel

def ScoreSingleAlignment( alignment, read_len=50, seed_len=20, expected_seed=3):
    seed_array=alignment[2]
    seed_set=set(list( seed_array))
    score=0
    penalty_array=numpy.array([1]*read_len)
    for seed in seed_set:
##        print seed_array
##        print seed
        matching_seeds=seed_array==seed
##        print matching_seeds
        score+=numpy.max(alignment[1,matching_seeds])
        penalty_array[int( seed):int(seed)+seed_len]=0
    penalty=penalty_array[15:].sum()*numpy.log(10)/30.
    return score-penalty

def ScoreAllAlignments(alignments, expected_seeds=9):
    lengths=numpy.array( [len(a[1,:]) for a  in alignments])
    max_length=3
##    print max_length
    scores=[]
    for align in alignments:
        if len(align[1,:])==1: continue
        scores.append(ScoreSingleAlignment(align))
    mapq=min(40, -4.34*numpy.log(1-( numpy.exp( max(scores))/sum(numpy.exp( scores)))))
    return mapq

def TestMapQ(infile, outfile, index):
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    outhandle=open(outfile, 'w')
    outtable=csv.writer(outhandle, delimiter='\t')
    est_mapq=[]
    prev_mapq=[]
    pos=[]
    row=intable.next()
    start=time.clock()
    count=0
    score_tracker=[]
    for row in intable:

        if row[0]!=row[6]: continue
        if row[1]==row[7]: continue
        count+=1
        r1, q1=row[4], row[5]
        exp_seq=row[0]
        exp_pos=(float(row[2])+float(row[3]))/2

        r2, q2=row[10], row[11]

        exp_seq_2=row[0]
        exp_pos_2=(float(row[8])+float(row[9]))/2

        score1, best_seq_1, best_loc_1, best_cigar_1=ScoreRead(r1,q1,index, exp_seq, exp_pos)
        score2, best_seq_2, best_loc_2, best_cigar_2=ScoreRead(r2,q2,index, exp_seq_2, exp_pos_2)
        row=[int(row[2]),int(row[8]), score1, score2, float(row[14]), float(row[15]),row[1], row[4],r2, best_seq_1, best_loc_1, best_cigar_1, best_seq_2, best_loc_2, best_cigar_2]
        outtable.writerow(row)
##        except:
##            continue

    inhandle.close()
    outhandle.close()
    print (time.clock()-start)
    print (count)
    return numpy.array( est_mapq), numpy.array( prev_mapq), numpy.array(pos), score_tracker


def PlotLocations(alignments, color):
    for x in alignments:
        pyplot.scatter(x[0,:], x[2,:], c=color, alpha=.3)

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
    kmers=[K.upper() for K in k_dict]
    kmer_dict={}
    for K in kmers:
        try: kmer_dict[K[:6]]+=[SequenceToInt( K)]
        except: kmer_dict[K[:6]]=[SequenceToInt( K)]
    tree_dict={}
    for prefix in kmer_dict:
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




def main(argv):
    param={}
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if param=={}: return()
    index=HashedIndex()

    LoadIndex(param['-ind'], index)
    print len (index.consensus_sequences)
    print len(index.kmer_location_dict.keys())
    if param.has_key('-o'):
        TestMapQ(param['-i'], param['-o'], index)
        return()
    r=param['-r']
    q=param['-q']
    count=0
    S,A,P,SC=0,0,0,0
    for t in range(1):
        start=time.clock()
        s,a,p,sc=ScoreRead(r,q,index)
        S+=s
        A+=a
        P+=p
        SC+=sc
        count+=time.clock()-start
    print count
    print S,A,P,SC
if __name__ == '__main__':
    main(sys.argv)
