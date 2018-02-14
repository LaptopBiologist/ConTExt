
#-------------------------------------------------------------------------------
# Name:        CoverageModeller
# Purpose:     Determine, model the biases of, and interpret the coverage of datasets
#               processed using context.
#
# Author:      Michael Peter McGurk
#
# Created:     10/12/2015
# Copyright:   (c) Michael 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy
import sys
import scipy
import scipy.stats
import os
import csv
import Bio
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
from scipy import optimize
from matplotlib import pyplot
from General_Tools import *

class ConTExtLine():
    """Class to organize ConTExt lines in a readable way."""
    strand={'0':'+', '16':'-', '4':'*'}
    def __init__(self, row):

        self.Seq1=row[0]
        self.Strand1=ConTExtLine.strand[row[1]]
        self.Start1=int(row[2])
        self.End1=int(row[3])
        self.sequence1=row[4]
        self.sequence2=row[10]
        self.Seq2=row[6]
        self.Strand2=ConTExtLine.strand[row[7]]
        self.Start2=int(row[8])
        self.End2=int(row[9])
        self.Cigar1=row[12]
        self.Cigar2=row[13]
        self.MapQ1=int(row[14])
        self.MapQ2=int(row[15])
        #self.Mismatches1=row[16]
        #self.Mismatches2=row[17]
        if self.Strand1=='+':
            self.threePrime1=self.End1
            self.fivePrime1=self.Start1
        else:
            self.threePrime1=self.Start1
            self.fivePrime1=self.End1
        if self.Strand2=='+':
            self.threePrime2=self.End2
            self.fivePrime2=self.Start2
        else:
            self.threePrime2=self.Start2
            self.fivePrime2=self.End2


def ReadDistributions(infile, cutoff=.005):
    """Read the summary of size distributions output by the structure calling pipeline."""
    handle=open(infile,'r')
    table=csv.reader(handle, delimiter='\t')
    mpd=table.next()
    kde=table.next()
    pr=table.next()
    try:
        domain=table.next()
    except:
        return numpy.array( [float(m) for m in mpd]),numpy.array( [float(k) for k in kde]),numpy.array( [float(p) for p in pr])
    m,k,p,d=numpy.array( [float(m) for m in mpd]),numpy.array( [float(k) for k in kde]),numpy.array( [float(p) for p in pr]), numpy.array( [float(d) for d in domain])
    handle.close()
    p=p[d>=0]
    cumul_k=numpy.cumsum(k)
    max_cutoff=d[ numpy.where(cumul_k>=(1.-cutoff))[0][0]]
    min_cutoff=d[ numpy.where(cumul_k>=(cutoff))[0][0]]

    return m, (min_cutoff, max_cutoff), p

def ReadGCContent(refFile, autosomes):
##    autosomes=['4']
    """Convert the major arms of a Drosophila reference genome into an array
    valued 1 where the nucleotide is a G or a C and valued 0 otherwise."""
    if refFile[-2:]=='gz':
        inhandle=gzip.open(refFile, 'r')
    else:
        inhandle=open(refFile, 'r')
    inParser=SeqIO.parse(inhandle, 'fasta')     #Read file as a *.fasta
    GC_dict={}
    for record in inParser:
        if autosomes.count( record.name)==0: continue
        sequence=numpy.fromstring(str(record.seq).lower(), '|S1', sep='')
        GC=(sequence=='g')+(sequence=='c')
        GC_dict[record.name]=GC
    return GC_dict

def ReadMaskedContent(refFile, autosomes):
##    autosomes=['4']
    """Convert the major arms of a Drosophila reference genome into an array
    valued True where the nucleotide is masked and valued False otherwise."""
    if refFile[-2:]=='gz':
        inhandle=gzip.open(refFile, 'r')
    else:
        inhandle=open(refFile, 'r')
    inParser=SeqIO.parse(inhandle, 'fasta')
    masked_dict={}
    for record in inParser:
        if autosomes.count( record.name)==0: continue
        sequence=numpy.fromstring(str(record.seq).lower(), '|S1', sep='')
        masked=(sequence=='n')
        masked_dict[record.name]=masked
    return masked_dict

def ConstructReads(len_dist, gap_size  ,  max_len=100):

    read=[sum(len_dist[:max_len]) ]
    for l in range(2,max_len):
        read.append(read[-1]-len_dist[l-1])
    rev_array= numpy.array(read)
    for_array=numpy.flipud(rev_array)
    max_prob=numpy.max(rev_array)
    gap=numpy.array([max_prob]*gap_size)
    read_array=numpy.hstack((for_array,gap, rev_array ))

    return read_array/read_array.sum()

def ConstructKernel(I, R):
    """Construct the kernel representing the probability that a nt is contained in
an insert spanning some position X nucleotides away. This kernel depends on the
insert size distribution (I) and the trimmed read length distribution (R). Note
that for simplicity this construction treats the two distributions as independent,
when in reality more heavily trimmed reads (r~R is small) should be negatively correlated
with i~I."""

    if type(I)==str:
        m,k,I=ReadDistributions(I)
    if type(R)==str:
        R,k,p=ReadDistributions(R)
    count=0.
    print I
    print R
    for i in range (len(I)):
        gap_size=i+1
        count+=1
        curr=I[i] * numpy.convolve(ConstructReads(R, gap_size), [1./gap_size]*gap_size)
        if count==1:
            last=curr

        else:
            curr[1:-1]+=last
            last=curr
    CDF=Cumulative(last)    #Stupid. Should have used numpy.cumsum.... alas!
    good_indices=numpy.where((( CDF>.0001)*(CDF<.9999))==True)[0]
    min_index=good_indices[0]
    max_index=good_indices[-1]
    trimmed_kernel=last[min_index:max_index]
    return trimmed_kernel/trimmed_kernel.sum()

def Cumulative(dist):
    CDF=[dist[0]]
    for i in range(1,len(dist)):
        CDF.append( CDF[i-1]+dist[i])
    return numpy.array( CDF)

def BuildPECoverageSignalAsArray(infile,chrom_length, ins_dist):
##    refLen=GetLengths(refFile)
    #Initialize a coverage array
    cvg=numpy.array([0.]*chrom_length)
    if len(ins_dist)==2:
        min_gap, max_gap=ins_dist
    else:
        cumul_dist=numpy.cumsum(ins_dist)
        max_gap=numpy.where(cumul_dist>=.999)[0][-1]
    print "Concordant definition:{0}<gap<{1}".format( min_gap, max_gap)
    handle=open(infile, 'r')
    table=csv.reader(handle, delimiter='\t')
    table.next()
    total_reads=0
    count=0
    length_list=[]
    q_list=[]
    for line in table:
        readpair=ConTExtLine(line)
        total_reads+=1    #Record that we've examined a read

        #Concordant read pairs must
        #1.) ...both align to the same sequence
        if readpair.Seq1!=readpair.Seq2: continue
##        if refNames.count( data.Seq1)==0: continue
##            length=data.Start2-data.Start1
            #Maybe base it on the coverage dist instead?
        #Only count a read pair if:

        #2.) ...map to opposite strands. If they map to the same strand, skip line.
        if readpair.Strand1==readpair.Strand2: continue


        #We n
        if readpair.Strand1=='+':
            length=readpair.threePrime2-readpair.threePrime1
            r,l=readpair.threePrime2,readpair.threePrime1
        elif readpair.Strand1=='-':
            length=readpair.threePrime1-readpair.threePrime2
            r,l=readpair.threePrime1,readpair.threePrime2

            #2.) Distance between reads is what's expected given the insert size distribution. If not, skip line.
        length_list.append(length)
        q_list.append(readpair.MapQ1+readpair.MapQ2)
        if length<min_gap or length>max_gap: continue

            #3.) The mapping quality isn't complete crap.
        if readpair.MapQ1+readpair.MapQ2<=30: continue
        count+=1
        #Increment by one the sequence coverer by the read pair
        cvg[l:r+1]+=1.
    print "{0} concordant reads, {1} total read:".format( count, total_reads)
    print "Q1,Q2,Q3 for gap distribution:{0}, {1}, {2}". format(  scipy.stats.scoreatpercentile(length_list, 25), scipy.stats.scoreatpercentile(length_list, 50), scipy.stats.scoreatpercentile(length_list, 75))
    print "Q1,Q2,Q3 for MapQ1+MapQ2 distribution:{0}, {1}, {2}".format(  scipy.stats.scoreatpercentile(q_list, 25), scipy.stats.scoreatpercentile(q_list, 50), scipy.stats.scoreatpercentile(q_list, 75))
    return cvg


def ConstructWeights(obs_gc, exp_gc):
    """Construct two histograms, and """
    obs_hist=numpy.histogram(obs_gc, bins=numpy.arange(0,1.01,.01))[0].astype(float)
    exp_hist=numpy.histogram(exp_gc, bins=numpy.arange(0,1.01,.01))[0].astype(float)
    return(exp_hist/obs_hist)

def DetermineGCofREads(infile,chrom_length, chrom_GC, ins_dist):
##    refLen=GetLengths(refFile)
    cvg=numpy.array([0.]*chrom_length)
    if len(ins_dist)==2:
        min_gap, max_gap=ins_dist
    else:
        cumul_dist=numpy.cumsum(ins_dist)
        max_gap=numpy.where(cumul_dist>=.999)[0][-1]
    print "Concordant definition:{0}<gap<{1}".format( min_gap, max_gap)
    handle=open(infile, 'r')
    table=csv.reader(handle, delimiter='\t')
    table.next()
    total=0
    count=0
    length_list=[]
    q_list=[]
    gc_list_l=[]
    gc_list_r=[]
    for line in table:
        data=ConTExtLine(line)
        if data.Seq1!=data.Seq2: continue
            #Maybe base it on the coverage dist instead?
        #Only count a read pair if:
            #1.) Reads map to opposite strands. If they map to the same strand, skip line.
        if data.Strand1==data.Strand2: continue
        total+=1
        if data.Strand1=='+':
            length=data.threePrime2-data.threePrime1
            r3,l3=data.threePrime2,data.threePrime1
            r5,l5=data.fivePrime2,data.fivePrime1
        elif data.Strand1=='-':
            length=data.threePrime1-data.threePrime2
            r3,l3=data.threePrime1,data.threePrime2
            r5,l5=data.fivePrime1,data.fivePrime2
##        if length>=len(cumul_dist):continue
            #2.) Distance between reads is what's expected given the insert size distribution. If not, skip line.
        length_list.append(length)
        q_list.append(data.MapQ1+data.MapQ2)
        if length<min_gap or length>max_gap: continue
            #3.) The mapping quality isn't complete crap.
        if data.MapQ1+data.MapQ2<=30: continue
        count+=1
        gc_list_r.append(read_gc_r)
        gc_list_l.append(read_gc_l)
    print "{0} concordant reads, {1} total read:".format( count, total)
    print "Q1,Q2,Q3 for gap distribution:{0}, {1}, {2}". format(  scipy.stats.scoreatpercentile(length_list, 25), scipy.stats.scoreatpercentile(length_list, 50), scipy.stats.scoreatpercentile(length_list, 75))
    print "Q1,Q2,Q3 for MapQ1+MapQ2 distribution:{0}, {1}, {2}".format(  scipy.stats.scoreatpercentile(q_list, 25), scipy.stats.scoreatpercentile(q_list, 50), scipy.stats.scoreatpercentile(q_list, 75))
    return numpy.array(gc_list_r), numpy.array(gc_list_l)

def StackArrays(array_list, pad_dim):
    stacked_array=numpy.array( array_list[0])
    for array in array_list[1:]:
        print array.shape
        pulled_array=numpy.array(array)
        stacked_dim=list( stacked_array.shape)
        array_dim=list( array.shape)
        pad_size= stacked_dim[pad_dim]-array_dim[pad_dim]
        print stacked_array.shape, pulled_array.shape
        if pad_size<0:
            pad_width=[[0,0]]*len(stacked_dim)
            pad_width[pad_dim]=[0,abs(pad_size)]
            print pad_width
            stacked_array=numpy.pad( stacked_array, pad_width, 'constant')

        elif pad_size>0:
            pad_width=[[0,0]]*len(array_dim)
            pad_width[pad_dim]=[0, abs(pad_size)]
            print pad_width
            pulled_array=numpy.pad( pulled_array, pad_width, 'constant')
        print stacked_array.shape, pulled_array.shape
        stacked_array=numpy.dstack((stacked_array, pulled_array))
    return stacked_array

def WriteCvgHistogramByArray(indir, kdeFile,lenFile,  refFile, outfile, window_size=470, refNames=['2R', '2L', '3L', '3R', 'X']):
    """indir is the path to the conTExt output for this sample,
    kdeFile is the path to a file which contains the insert size distribution for this sample
    lenFile is the path to a file which contains the aligned read length distribution for this sample,
    refFile is the path to the repeat masked reference genome to which the data was aligned."""
    m,k,p=ReadDistributions(kdeFile)
    meanins=int( WeightedAverage(m, numpy.arange(len(m))))
    print k
    m_len,k_len,p_len=ReadDistributions(lenFile)
    meanread=int( WeightedAverage(m_len, numpy.arange(len(m_len))))+1
    gc_raw=ReadGCContent(refFile, refNames)
    masked_raw=ReadMaskedContent(refFile, refNames)
    kernel=ConstructKernel(p, k_len)
    window_size=meanins+2*meanread
    flat_window=numpy.array( [1./window_size]*window_size, float)
    count=0
    gc_Dict={}
    cvg_Dict={}
##    refNames=['2R']
    for chrom in refNames:
        print chrom
##        print wigHandle.ChromLength[chrom]
##        iv=HTSeq.GenomicInterval(chrom, 0, wigHandle.ChromLength[chrom])
        gc_chr=numpy.convolve(gc_raw[chrom], flat_window, 'same')
        masked_chr=numpy.convolve(masked_raw[chrom], flat_window, 'same')
        chrom_file='{0}/{1}.dat'.format(indir, chrom)
        cvg_chr=BuildPECoverageSignalAsArray(chrom_file, len(gc_raw[chrom]), k)
        poorly_mapped_cvg=numpy.convolve(cvg_chr,flat_window , 'same')
        print cvg_chr.shape
        good_indices = (masked_chr<=.01)*(poorly_mapped_cvg>=.001)
        gc_Dict[chrom] = gc_chr[good_indices]
        cvg_Dict[chrom] = cvg_chr[good_indices]

##        print cvg_array.shape
    array_list=[]
    gc_array_list=[]
    for key in refNames:
        cvg_array=cvg_Dict[key]
        gc_array=gc_Dict[key]
        gc_hist=numpy.histogram(gc_array, bins=numpy.arange(0,1.01,.01))
        gc_array_list.append(gc_hist)
##        pyplot.scatter(gc_array, cvg_array, alpha=.01)
##        pyplot.show()
##        pyplot.close()
        max_cvg=cvg_array.max()

        hist=numpy.histogram2d(gc_array,cvg_array, bins=(numpy.arange(0,1.01,.01),numpy.arange(0,max_cvg)))

        p_hist=hist[0]/hist[0].sum(1)[:,None]
        array_list.append(p_hist)
##        numpy.savetxt(outfile+'_hist_{0}.tsv'.format(key), p_hist,'%f', delimiter='\t')
        exp_array=p_hist*numpy.arange(0,max_cvg-1) [None,:]
        numpy.savetxt(outfile+'_exp_{0}.tsv'.format(key), exp_array.sum(1),'%f', delimiter='\t')
    stacked_array=StackArrays(array_list,1)
    gc_stacked=numpy.array( gc_array_list)
    numpy.save(outfile+'_hist.npy', stacked_array)
    numpy.save(outfile+'_GCcounts.npy', gc_stacked)




def WriteGCHistograms(indir, kdeFile,lenFile,  refFile, outfile, window_size=300, refNames=['X', '2L','3L', '2R', '3R']):
    """indir is the path to the conTExt output for this sample,
    kdeFile is the path to a file which contains the insert size distribution for this sample
    lenFile is the path to a file which contains the aligned read length distribution for this sample,
    refFile is the path to the repeat masked reference genome to which the data was aligned."""
    m,k,p=ReadDistributions(kdeFile)
    print k

    m_len,k_len,p_len=ReadDistributions(lenFile)
    meanins=int( WeightedAverage(m_len, numpy.arange(len(m_len))))+1
    print meanins
    gc_raw=ReadGCContent(refFile)
    masked_raw=ReadMaskedContent(refFile)
    kernel=ConstructKernel(p, k_len)
    flat_window=numpy.array( [1./window_size]*window_size, float)
    count=0
    gc_Dict={}
    mean_Dict={}
    cvg_Dict={}
    for chrom in refNames:
        print chrom
##        print wigHandle.ChromLength[chrom]
##        iv=HTSeq.GenomicInterval(chrom, 0, wigHandle.ChromLength[chrom])
        gc_chr_approx=numpy.convolve(gc_raw[chrom], kernel, 'same')
        gc_chr_mean=numpy.convolve(gc_raw[chrom], [1./meanins]*meanins, 'same')
        gc_chr_approx=numpy.hstack((gc_chr_approx, gc_chr_approx))
        gc_chr_mean=numpy.hstack((gc_chr_mean, gc_chr_mean))
        masked_chr=numpy.convolve(masked_raw[chrom], kernel, 'same')
        chrom_file='{0}/{1}.dat'.format(indir, chrom)
        cvg_chr_r, cvg_chr_l=DetermineGCofREads (chrom_file, len(gc_raw[chrom]), gc_raw[chrom],k)
        cvg_chr=numpy.hstack((cvg_chr_l, cvg_chr_r))
        poorly_mapped_cvg=numpy.convolve(cvg_chr,flat_window , 'same')
        print cvg_chr.shape
        good_indices = (masked_chr<=.01)
        bin_step=1./meanins
        bins=numpy.arange(0, 1+bin_step, bin_step)
        gc_Dict[chrom] =numpy.histogram( gc_chr_approx, bins=bins)[0]
        mean_Dict[chrom]=numpy.histogram( gc_chr_mean, bins=bins)[0]
        cvg_Dict[chrom] = numpy.histogram( cvg_chr, bins=bins)[0]

##        print cvg_array.shape
    array_list=[]
    gc_array_list=[]
    mean_array_list=[]
    for key in ['2R', '2L', '3L', '3R', 'X']:
        array_list.append(cvg_Dict[key])
        gc_array_list.append(gc_Dict[key])
        mean_array_list.append(mean_Dict[key])
    stacked_array=numpy.vstack(array_list)
    gc_stacked=numpy.vstack(gc_array_list)
    mean_stacked=numpy.vstack( mean_array_list)
    numpy.save(outfile+'_obs_reads.npy', stacked_array)
    numpy.save(outfile+'_exp_reads_approx.npy', gc_stacked)
    numpy.save(outfile+'_exp_reads_mean.npy', mean_stacked)


#Read depth distribution processing

def FastHistoKDE_BW(histo, bw, ignore_zero=True, kernel='norm', var_array=[]):
    """A quick kernel density estimator for discrete data. If W(x) is a bin in a histogram
    describing the frequency with which a given value x is a observed, the kde F(x) can be written:
    F(x)=sum_x' (W(x')*K_bw(x,x')
    where sum_x' sums over all bins in the histogram, and K_bw(x,x') is the kernel function with bandwidth bw"""

    bw=float(bw)
    if ignore_zero==True:   #If true, ignore the bin counting positions with zero read depth; may not represent single copy sequence, so more error in this estimate
        histo[0]=histo[1] #Fill in the zero with the observed frequency of positions with read depth equals 1
        histo/=sum(histo) #Renormalize
    n=sum(histo)
    x1=numpy.arange(0, len(histo)) #The domain of all possible x'
    x2=numpy.arange(0, len(histo))
    dist_array=numpy.subtract.outer(x1, x2) #Compute all pairwise distances with domains. This produces a 2-dimensional matrix D_ij = i - j
    #The above matrix will be passed into the kernel; it's easier to write if this is divided by the bandwidth
    u=dist_array/bw

    #Some options for different kernels
    if kernel=='norm':
        norm_array=numpy.exp(-(u**2)/2)/(2*numpy.pi)**.5
##    missing=.5*(1+scipy.special.erf( (0-numpy.arange(len(histo)))/(bw*(2**.5))))

    elif kernel=='epanechnikov':
        norm_array=.75*(1-u**2)*(abs(u)<=1)
    else:
        norm_array=(35./32)*(1-u**2)**3*(abs(u)<=1)
    #We have a matrix describing the probability of x' given x (or vice versa) for each pair of values in the data's domain
    norm_array/=norm_array.sum(1) #For the hell of it, let's normallize this matrix so it sums to one along one axis
                                  #actually not for the hell of it... kernels are unit distributions....
    weighted=histo[:,None]*  norm_array     #Great! Now weight these pairs by W(x)...

    u_kde=weighted.sum(0)/(n*bw)    #A KDE is the sum of unit distributions, so the integral will be equal N... or possibly N*bw.... I don't remember.
##    missing=.5*(1+scipy.special.erf( (0-numpy.arange(len(histo)))/(bw*(2**.5))))
####    print missing
##    u_kde/=(1.-missing)

    return u_kde/sum(u_kde)         #Doesn't matter if I divided by bw one too many times... it's still proportional to the correct KDE, so now Imma force it to sum to one.


def Cross_validate(tuple_data, min_bw, max_bw, steps, depth=3, fold=100, cutoff=.99):
    """Idea: We have read count distributions for each autosome arm and the X chromosome.
    So, the we can get a sense of the noise in our estimates by comparing the distributions
    of different chromosome arms. Use this to choose a good bandwidth for kernel density estimation."""

    data_count=float( tuple_data.shape[1])  #Number of chromosomes to choose from
    full_data=tuple_data.sum(1)
    max_index=numpy.where(Cumulative(numpy.asarray( full_data,float)/sum(full_data))>=cutoff)[0][0] #Only going to consider the first 99% of the read count distribution. Don't want outlines to favor undersmoothing

    kernels=['norm', 'epanechnikov', 'triweight']
    samples=[]


    bw_list=[]
    entropy_list=[]

    #Create a set of permuted read depth distribution pairs where Fr(read count=X)
    #is the freq of read count=x of a single autosome and Fr'(read count=X) is the freq for a different autosome.
    #Samples then ends up being a list of tuples
    for i in range (fold):
        samples.append(HalveData(tuple_data[:max_index,:]))

    for rep in range(depth):
##        kde_main=FastHistoKDE_BW(data, bw)
        entropy=0.
        step_size=float(max_bw-min_bw)/steps
        params=numpy.arange(min_bw, max_bw, step_size)
        bw_list, entropy_list=OptimizeInRange(samples,params)
        try:
            best_bw=bw_list[numpy.nanargmin(entropy_list)]
        except:
            best_bw=numpy.nan
            break
        max_bw=best_bw+step_size
        min_bw=best_bw-step_size

##    print bw_list[numpy.nanargmin(entropy_list)]
##    return bw_list, entropy_list
    return best_bw
def OptimizeInRange(samples, param):
    bw_list=[]
    entropy_list=[]
    fold=len(samples)
    for bw in param:
        entropy=0.
        for sample in samples:

            f1,f2=sample

            kde=FastHistoKDE_BW(f1,float( bw), False )
            entropy+=( (kde-f2)**2).sum()/len(f2)

            kde=FastHistoKDE_BW(f2,float( bw), False )
            entropy+=( (kde-f1)**2).sum()/len(f2)

        entropy_list.append(entropy/(2*fold))

        bw_list.append(bw )
    return bw_list, entropy_list

def HalveData(data):
    rows, columns=data.shape
    fold_size=float( columns)/2
    if fold_size!=int(fold_size):
        fold_size=int(fold_size)+1
##    print fold_size
    fold_size=int(fold_size)
    index_array= numpy.ndarray((rows, fold_size), int)
    for i in range(rows):
        index_array[i,:]=numpy.random.choice(range(columns), size=fold_size, replace=False)

    fold_1=numpy.array([0.]*rows)

    for i in range(fold_size):
        fold_1+=data[range(rows),index_array[:,i]]
    fold_1=data[range(rows), index_array[:,0]]
##    fold_2=data.sum(1)-fold_1
    fold_2=data[range(rows), index_array[:,1]]
##    fold_1/=fold_size
##    fold_2/=(columns- fold_size)
    return fold_1/sum(fold_1), fold_2/sum(fold_2)

def SmoothCvgByGC(table):
    autosome_list=[]
    x_list=[]
    print "%GC\tOptimal bandwidth"
    for i in range(table.shape[0]):
        table_at_gc=table[i,:,:4]
##        print table_at_gc.sum()
        if numpy.isnan( table_at_gc.sum())==True:

            autosome_list.append(table_at_gc.sum(1)/4.)
            x_list.append(table[i,:,4])
            continue
        bw=Cross_validate(table_at_gc, 0, 100, 10, 3)
        print "{0}\t{1}".format( i, bw)
        kde_autosome=FastHistoKDE_BW(table_at_gc.sum(1)/4., bw, True)
        kde_x=FastHistoKDE_BW(table[i,:,4], bw, True)
        autosome_list.append(kde_autosome)
        x_list.append(kde_x)
##        print i
    return autosome_list, x_list


import seaborn
from matplotlib import pyplot
def PlotCoverageScatters(table, read_depth, auto_kde=[], x_kde=[]):
    colors=['r','b','g','purple', 'orange']
    chrom=['2R', '2L', '3L', '3R', 'X']
    sizes=[40]*4+[80]
    R2=WeightedAverage( table[read_depth,:,0],numpy.arange(len(table[read_depth,:,0])))
    for i in range(5):
##        pyplot.scatter( range(1,len(table[read_depth,:,i])), table[read_depth,1:,i]/(sum(table[read_depth,1:,i]+table[read_depth,1,i])), c=colors[i], s=sizes[i], label=chrom[i], alpha=.5  )
        pyplot.scatter( range(len(table[read_depth,:,i])), table[read_depth,:,i], c=colors[i], s=sizes[i], label=chrom[i], alpha=.5  )
        print chrom[i],WeightedAverage( table[read_depth,:,i],numpy.arange(len(table[read_depth,:,i]))),
        print WeightedStdDev( table[read_depth,:,i],numpy.arange(len(table[read_depth,:,i])))

##    pyplot.scatter( range(len(table[read_depth,:,i])), table[read_depth,:,:4].sum(1)/4., c='teal', s=80, label='Autosomes', alpha=.5  )

    pyplot.plot(scipy.stats.poisson.pmf(range(len(table[read_depth,:,i])),28))
    if auto_kde!=[]: pyplot.plot(auto_kde[read_depth], lw=2, c='indianred')
    if x_kde!=[]: pyplot.plot(x_kde[read_depth], lw=2, c='powderblue')
    pyplot.legend()
    pyplot.ylabel('Frequency')
    pyplot.xlabel('Read depth')


def main(argv):
    param={}
    print argv
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if len( param.keys())==0: return
    #inDir, outDir, AssemblyIndex, TEIndex, threads, length=70, ReadIDPosition=1, PositionRange=0, phred=33, shorten=True, Trim=True
    if param.has_key('-join'):
        StackKDEs(param['-i'], param['-o'])
        return
    wigfile=param['-wig']
    reffile=param['-ref']
    kdefile=param['-kde']
    outroot=param['-o']
    lenFile=param['-len']

    if param.has_key('-exp'):
        ExpectedGCContent(kdefile, lenFile)
    elif param.has_key('-hist'):
        WriteGCHistograms(wigfile, kdefile,lenFile, reffile, outroot)
    else:
##        WriteGCHistograms(wigfile, kdefile,lenFile, reffile, outroot)
        WriteCvgHistogramByArray(wigfile, kdefile,lenFile, reffile, outroot)
if __name__ == '__main__':
    main(sys.argv)