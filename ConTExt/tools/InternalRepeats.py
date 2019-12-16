import os
import numpy
import scipy
import sklearn.neighbors
from Bio import SeqIO
from Bio import Seq
import csv

from annoy import AnnoyIndex




def SeqToExpanded(seq,ignore_N=False):
    """Converts each nucleotide in a sequences to a series of three zeroes and a one:
    A - 1000
    T - 0100
    G - 0010
    C - 0001

    The dot product of two sequences so represented is equal to the number of
    matching nucleotides in an ungapped end-to-end alignment between them.

    For example:

    [A=A]=1

    and

    1*1 + 0*0 + 0*0 + 0*0 = 1

    But

    [A=T]=0

    and

    1*0 + 0*1 + 0*0 + 0*0 = 0
    """
    char_array=numpy.fromstring(seq.upper(), '|S1')
    if ignore_N==True:
        char_array=char_array[char_array!='N']

    num_array=numpy.ndarray(((len(char_array))*4))
    num_array.fill(0.)

    nt=['A', 'T', 'G', 'C']
    for beg in range(4):

        num_array[beg::4][char_array==nt[beg]]=1.

    return num_array

class KmerIndex():

  def __init__(self,training_set, f=70):
      #Store original kmers
      self.training_set=[]
      f = len(training_set[0])
      for value in training_set:
        assert len(value)==f, "All kmers in the training set must be the same length"
        self.training_set.append(SequenceToInt(value))
      self.training_set=numpy.array(self.training_set)
      #Build the Annoy index

      f = len(value)

      self.index = AnnoyIndex(f*4, metric='angular')  # Length of item vector that will be indexed
      for i,v in enumerate(training_set):
          self.index.add_item(i, SeqToExpanded( v))

      self.index.build(20) # 20 trees

#       self.index=t


  def query_radius(self, query,rad=.05, k=50, use_annoy_distances=False):
      """Identify the k nearest neighbors and return those within the specified radius (rad)
      of the query."""
      int_query=SequenceToInt(query)
      if use_annoy_distances==True: #Convert the angular distance into the Hamming distance
        result_indices,distances=self.index.get_nns_by_vector(SeqToExpanded(query), k, include_distances=True)
        result_indices=numpy.array(result_indices)
        distances=numpy.array(distances)
        cossim=numpy.cos(numpy.pi*(distances)/2)
        hamming_distance=(1.-cossim)/2

      else: #Exactly compute the Hamming distance
        result_indices=self.index.get_nns_by_vector(SeqToExpanded(query), k, include_distances=False)
        result_indices=numpy.array(result_indices)
        results=self.training_set[result_indices,:]
        hamming_distance=1-numpy.mean ((results-int_query)==0,1)


      within_radius=hamming_distance<=rad


      return result_indices[within_radius]

  def compare_distances(self, query,rad=.05, k=50, distance='Angular'):
      """Identify the k nearest neighbors and return those within the specified radius (rad)
      of the query."""
      int_query=SequenceToInt(query)

      result_indices,distances=self.index.get_nns_by_vector(SeqToExpanded(query), 30, include_distances=True)
      distances=numpy.array(distances)

      result_indices=numpy.array(result_indices)


      #pull results
      results=self.training_set[result_indices,:]

      hamming_distance=1-numpy.mean ((results-int_query)==0,1)
      within_radius=hamming_distance<=rad
      if distance=="Angular":
        pyplot.scatter(hamming_distance, distances, c='black', s=10)
#       pyplot.xlabel('Hamming Distance')
#       pyplot.ylabel('Angular Distance')
#       pyplot.show()
      if distance=="Hamming":
        cossim=numpy.cos(numpy.pi*(distances)/2)
        approx_hamming=(1.-cossim)/2
        pyplot.scatter(hamming_distance, approx_hamming, c='black', s=10 )
#       pyplot.xlabel('Hamming Distance')
#       pyplot.ylabel('Approximate Hamming')
#       pyplot.show()

      return result_indices[within_radius]


def CleanName(name):
    illegal=['|', '!','>', '<', '?','/','*']
    for i in illegal:
        name='_'.join(name.split(i))
    return(name)

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

def CountKMERS(sequence, k=10, unbiased=False, skip_Ns=True ):
    """Decomposes sequence into kmers. If unbiased, repeat the first k nucleotides
    at the end of the sequence (assumes the true sequence is tandemly repeated)."""
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
    """Converts a DNA sequence to a sequence of integers.
  The hamming distance used by the Ball Tree requires integers."""
    seq_array=numpy.fromstring(seq, '|S1')
    int_array=numpy.ndarray((len(seq_array),))
    int_array.fill(0)
    nt_array=['T', 'C', 'G']
    for i, nt in enumerate( nt_array):
        int_array[seq_array==nt]=i+1
    return int_array



def SummarizeRepetiveness(sequence,k=20, mismatches=.15, max_sample=400,N=10, unbiased=False, method='BallTree'):

    assert mismatches<k, "Number of mismatches must be less than kmer size"
    assert mismatches>=0, "Number of mismatches cannot be negative. That's not how counting works."

    #Assume values between 0 and 1 reflect percent divergence
    #Assume values great than one reflect the absolute number of mismatches
    if mismatches>1:
      mismatches/=float( k)

    kmer_dict=CountKMERS(sequence,k, unbiased)
    kmer_list, kmer_counts=zip(* kmer_dict.items())

    radius=mismatches

    encoded_kmers=numpy.array( [SequenceToInt(kmer) for kmer in kmer_list ])
    if method=='BallTree':
      ball_tree=sklearn.neighbors.BallTree(encoded_kmers, metric='hamming')
      count_list=[]

      indices=numpy.arange(len(encoded_kmers))
      if len(encoded_kmers)>max_sample:
          pr=numpy.array(kmer_counts, float)
          pr/=pr.sum()

          indices=numpy.random.choice(numpy.arange(len(pr)), size=max_sample, p=pr)
      for index,kmer in enumerate( encoded_kmers):
          if numpy.isin(index, indices)==False: continue
          neighbors=ball_tree.query_radius(kmer[None,:], radius)[0]

          counts=numpy.array( [kmer_dict[kmer_list[n] ] for n in neighbors])

          distance_mat=(kmer[None,:] !=encoded_kmers[neighbors,:])

          distance=(1.-distance_mat.sum(1)/k )
          count_list.append(numpy.sum( distance*counts) )

      S= numpy.mean( count_list)
      return S
    else:

      ball_tree=KmerIndex(kmer_list)

      count_list=[]

      indices=numpy.arange(len(encoded_kmers))
      if len(encoded_kmers)>max_sample:
          pr=numpy.array(kmer_counts, float)
          pr/=pr.sum()

          indices=numpy.random.choice(numpy.arange(len(pr)), size=max_sample, p=pr)
      for i,kmer in enumerate(kmer_list):
          neighbors=ball_tree.query_radius(kmer, radius,N)

          counts=numpy.array( [kmer_dict[kmer_list[n] ] for n in neighbors])

          distance_mat=(encoded_kmers[i][None,:] !=encoded_kmers[neighbors,:])

          distance=(1.-distance_mat.sum(1)/k )
          count_list.append(numpy.sum( distance*counts) )

      S= numpy.mean( count_list)

      return S

def GetKMERLocations(sequence, k=10, unbiased=False, skip_Ns=True ):
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
            kmerSet[kmer]=[]
        kmerSet[kmer].append(x)
    return kmerSet
##        kmerSet.add(str(complement[x:x+k]))
    return list(kmerSet)
def IntToSequence(seq):
    seq_array=numpy.array(['']*len(seq), '|S1')

    nt_array=['T', 'C', 'G']
    seq_array[seq==0]="A"
    seq_array[seq==1]="T"
    seq_array[seq==2]="C"
    seq_array[seq==3]="G"

    return ''.join(seq_array)


def InitializeDotPlot(sequence,k=20, mismatches=.15, unbiased=True,N=30):

    assert mismatches<k, "Number of mismatches must be less than kmer size"
    assert mismatches>=0, "Number of mismatches cannot be negative. That's not how counting works."

    #Assume values between 0 and 1 reflect percent divergence
    #Assume values great than one reflect the absolute number of mismatches
    if mismatches>1:
      mismatches/=float( k)

    kmer_dict=CountKMERS(sequence,k, unbiased)
    kmer_locations=GetKMERLocations(sequence,k, unbiased)
    kmer_list, kmer_counts=zip(* kmer_dict.items())

    radius=mismatches

    encoded_kmers=numpy.array( [SequenceToInt(kmer) for kmer in kmer_list ])
    ball_tree=KmerIndex(kmer_list)



    indices=numpy.arange(len(encoded_kmers))
#     if len(encoded_kmers)>max_sample:
#         pr=numpy.array(kmer_counts, float)
#         pr/=pr.sum()

#         indices=numpy.random.choice(numpy.arange(len(pr)), size=max_sample, p=pr)
    dot_dict={}
    visited_locations=set()
    for index,kmer in enumerate( encoded_kmers[indices]):
        skip=False

        neighbors=ball_tree.query_radius(kmer_list[index], radius,N)
#         print neighbors
        kmer=IntToSequence( kmer)
        if dot_dict.has_key(IntToSequence( kmer))==False:
          dot_dict[kmer]=set()

        for location in kmer_locations[kmer]:
            visited_locations.add(location)
            dot_dict[kmer].add(location)
        for n in neighbors:
          target_kmer=IntToSequence(encoded_kmers[ n])
          for location in kmer_locations[target_kmer]:
            visited_locations.add(location)
            dot_dict[kmer].add(location)



    for key in dot_dict:
      dot_dict[key]=numpy.array(sorted(list( dot_dict[key])))

    return dot_dict


def PlotDots(dot_dict, method='half', max_sample=None):
  indices=numpy.arange(len(dot_dict.keys()))
  if max_sample is not None:

    if len(indices)>max_sample:
        indices=numpy.random.choice(indices, size=max_sample)
  for i in indices:
      key= dot_dict.keys()[i]

      if method=='half':
        for loc in dot_dict[key]:
          pyplot.scatter(dot_dict[key][0],loc, c='black', alpha=.5, s=8)

      elif method=='diag':
        pyplot.scatter(dot_dict[key][0],dot_dict[key][0], c='black', alpha=.5, s=8)
      else:
        for loc in dot_dict[key]:
          for loc1 in dot_dict[key]:
            pyplot.scatter(loc1,loc, c='black', alpha=.5, s=8)
#   pyplot.show()


def GetFirstInstance(dot_dict):
  pos_list=[]
  for key in dot_dict:

    pos_list.append( dot_dict[key][0])
  return numpy.array(sorted( pos_list))


def MaskSequenceByKmers(sequence,k=70, mismatches=5, unbiased=True, N=30, method="unmask"):

    assert mismatches<k, "Number of mismatches must be less than kmer size"
    assert mismatches>=0, "Number of mismatches cannot be negative. That's not how counting works."

    #Assume values between 0 and 1 reflect percent divergence
    #Assume values great than one reflect the absolute number of mismatches
    if mismatches>1:
      mismatches/=float( k)

    kmer_dict=CountKMERS(sequence,k, unbiased)
    kmer_locations=GetKMERLocations(sequence,k, unbiased)
    kmer_list, kmer_counts=zip(* kmer_dict.items())

    radius=mismatches

##    encoded_kmers=numpy.array( [SequenceToInt(kmer) for kmer in kmer_list ])
    encoded_kmers=[kmer for kmer in kmer_list ]
    ball_tree=KmerIndex(kmer_list)



    indices=numpy.arange(len(encoded_kmers))

    dot_dict={}
    visited_locations=set()

    seq_len=len(sequence)

    if unbiased==False:
        indices= range(0,len(sequence)-k+1)
    else:
        sequence=sequence*2
        indices= range(0, seq_len)
    seq_array=numpy.fromstring(sequence, '|S1')
    masked_sequence=numpy.fromstring(sequence, '|S1')
    if method=='unmask':
      masked_sequence.fill('N')
    for kmer in kmer_list:
#         kmer=''.join(sequence[x:x+k])
#         kmer=sequence[x:x+k]

        x=min(kmer_locations[kmer])
        if kmer.count("N")>0: continue
        encoded_kmer=kmer
        skip=False

        neighbors=ball_tree.query_radius(kmer, radius,N)
#         print neighbors

        if method=='elimination':
          for i,location in enumerate( sorted(kmer_locations[kmer])):
              if location<=x: continue

              masked_sequence[location:location+k]='N'
          for n in neighbors:
            target_kmer=encoded_kmers[ n]
            for i,location in enumerate( sorted(kmer_locations[target_kmer])):
              if location<=x: continue

              masked_sequence[location:location+k]='N'
        elif method=='unmask':
          locations=[]
          locations+=kmer_locations[kmer]
          for n in neighbors:
            target_kmer=encoded_kmers[ n]
            locations+=kmer_locations[target_kmer]
          location=min(locations)
##          masked_sequence[location:location+k]=seq_array[location:location+k]
          masked_sequence[location]=seq_array[location]
    if method=='unmask':
        #Corrects for the uniqueness of the junction: Want to exclude the junction between
        #two repeats and only retain the first monomer
        new_masked_sequence=numpy.copy(masked_sequence)
        window_masked_count=0.
        for i,nt in enumerate( masked_sequence[indices]):
          if nt=='N': window_masked_count+=1
          else: window_masked_count-=1
          if window_masked_count>20: window_masked_count=20
          if window_masked_count<0: window_masked_count=0
          if window_masked_count/20.>=.75:
            new_masked_sequence[i:i+k-10]='N'
        masked_sequence=new_masked_sequence


    return ''.join(masked_sequence[indices])

def MaskInternalRepeats(infile, outfile, k=70, mismatches=5, unbiased=True):
    """Takes a FASTA and brings all sequences in phase with each other"""
    seqs=GetSeq(infile)
    outhandle=open(outfile, 'w')
    report_file='.'.join(outfile.split('.')[:-1])+'_report.tsv'

    report_handle=open(report_file, 'w')
    report_table=csv.writer(report_handle, delimiter='\t')
    report_line=['Name',"Removed",  "Length","Masked (init)", "Masked (final)"]
    report_table.writerow(report_line)
    query_seq=seqs.values()[0]
    for i,key in enumerate( sorted( seqs.keys())):
        masked_count=seqs[key].upper().count('N')
        seq_len=float(len(seqs[key]))
        unmasked_count=seq_len-masked_count
        if unmasked_count<100:
            print key
            report_line=[key,False,  seq_len,masked_count, -1]
            report_table.writerow(report_line)
            continue
        try:
            masked_seq=MaskSequenceByKmers(seqs[key], k=k, mismatches=mismatches, unbiased=True, method='elimination')
        except:
            print key
            report_line=[key,False,  seq_len,masked_count, -1]
            report_table.writerow(report_line)
            continue
        masked_count_new=masked_seq.upper().count('N')
        outhandle.write('>{0}\n'.format(key))
        outhandle.write('{0}\n'.format(masked_seq))
        report_line=[key,True,  seq_len,masked_count, masked_count_new]
        report_table.writerow(report_line)
        report_handle.flush()
    outhandle.close()