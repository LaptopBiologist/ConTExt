#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      I am
#
# Created:     06/03/2018
# Copyright:   (c) I am 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import csv
import numpy
import os
import scipy
import scipy.signal
import sys
try:
    import tools.CoverageModeller as CoverageModeller
    from tools.General_Tools import *
except:
    print('Failed to import tool modules')

csv.field_size_limit(sys.maxsize)


def ReadDistributions(infile, cutoff=.005, return_full=True):
    """Reads the *.tsv file which stores the insert size or read length distribution.

    If return_full==True, it returns the observed MPD distribution, its KDE, and the
    KDE conditioned on read pairs spanning a junction.

    If return_full==False, it returns the MPD distribution, the cutoff% and (1-cutoff)%
    quantiles of the KDE, and the KDE conditioned on read pairs spanning a junction. """

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
##    print cumul_k
##    print numpy.where(cumul_k>=(1.-cutoff))
    max_cutoff=d[ numpy.where(cumul_k>=(1.-cutoff))[0][0]]
    min_cutoff=d[ numpy.where(cumul_k>=(cutoff))[0][0]]
##    median= numpy.where(cumul_k>=.5)[0][0]
##    offset=min(abs(max_cutoff-median), abs(min_cutoff-median))

##    print (median-offset,median+offset)

    if return_full==True:
        return m[d>=0]/m[d>=0].sum(),k[d>=0]/k[d>=0].sum(), p[d>=0]/p[d>=0].sum()

    return m, (min_cutoff, max_cutoff), p


def mergeDicts(dict_1, dict_2):
    new_dict={}
    key_intersect=set(dict_1)&set(dict_2)
    if len(key_intersect):
        raise()
    for key in dict_1.keys():
        new_dict[key]=dict_1[key]
    for key in dict_2.keys():
        new_dict[key]=dict_2[key]
    return new_dict

def EstimateJunctionCopyNumber(infile,outfile, kdeFile, lenFile,seqDict, kdefile ):

    cvg=numpy.load(kdefile)
##    kde=numpy.load(kdefile)

    kde, kde_x=CoverageModeller.SmoothCvgByGC(cvg)
    inDir='_'.join(infile.split('_')[:-1])
    kde_root=inDir+"_cvg_kde"
    numpy.save(kde_root+'_auto.npy', kde)
    numpy.save(kde_root+'_x.npy', kde_x)

##    table=numpy.load(cvg_fil/e)
    kde_auto=kde
    conv_auto=ConvolveAll( kde_auto)
    conv_x=ConvolveAll( kde_x)

    kernel=CoverageModeller.ConstructKernel(kdeFile, lenFile)
    kernel_len=len(kernel)
    inHandle=open(infile, 'rb')
    inTable=csv.reader(inHandle, delimiter='\t')
    outHandle=open(outfile, 'wb')
    outTable=csv.writer(outHandle, delimiter='\t')
    row=inTable.next()
    header=['Seq1', 'Seq2', 'Strand1', 'Strand2','Feature', 'X_est', 'X_range', 'Y_est', 'Y_range', 'Read Count', "X 3'-positions", "Y 3'-positions","X 5'-positions", "Y 5'-positions", "X MapQ", "Y MapQ",'ID', 'X CIGAR', 'Y CIGAR', 'X MD String', 'Y MD String', 'Copy Number (MAP)', 'Lower', 'Upper']
    outTable.writerow(header)
    out_root='.'.join(outfile.split('.')[:-1])
##    numpy.save(out_root+'_cvg_auto.npy', kde_auto)
##    numpy.save(out_root+'_cvg_x.npy', kde_x)
    for row in inTable:
        line=cluster(row)
        line.ClassifyJunction()

    #Check this and fix it if necessary
        if numpy.isnan( line.j0_y)==True or numpy.isnan( line.j0_x)==True:
            outTable.writerow(line.row()+['nan','nan','nan'])
            continue
        outHandle.flush()
        if line.Quad1=='+' and line.Quad2=='-':
            seq2_bound_right=min(seqDict[line.Seq2], int( line.j0_y))
            seq2_bound_left=max(0, seq2_bound_right-kernel_len)
            seq2=seqDict[line.Seq2][seq2_bound_left:seq2_bound_right]

            seq1_bound_left=max(0, int( line.j0_x))
            seq1_bound_right=min(len(seqDict[line.Seq1]), seq1_bound_left+kernel_len)
            seq1=seqDict[line.Seq1][seq1_bound_left:seq1_bound_right]

        if line.Quad1=='-' and line.Quad2=='+':
            seq1_bound_right=min(seqDict[line.Seq1], int( line.j0_y))
            seq1_bound_left=max(0, seq1_bound_right-kernel_len)
            seq1=seqDict[line.Seq1][seq1_bound_left:seq1_bound_right]

            seq2_bound_left=max(0, int( line.j0_x))
            seq2_bound_right=min(len(seqDict[line.Seq2]), seq2_bound_left+kernel_len)
            seq2=seqDict[line.Seq2][seq2_bound_left:seq2_bound_right]

        if line.Quad1=='+' and line.Quad2=='+':
            seq1_bound_right=min(seqDict[line.Seq1], int( line.j0_x))
            seq1_bound_left=max(0, seq1_bound_right-kernel_len)
            seq1=seqDict[line.Seq1][seq1_bound_left:seq1_bound_right]

            seq2_bound_left=min(seqDict[line.Seq2], int( line.j0_y))
            seq2_bound_right=min(len(seqDict[line.Seq2]), seq2_bound_left+kernel_len)
            seq2=seqDict[line.Seq2][seq2_bound_left:seq2_bound_right][::-1]

        if line.Quad1=='-' and line.Quad2=='-':
            seq1_bound_left=max(0, int( line.j0_x))
            seq1_bound_right=min(len(seqDict[line.Seq1]), seq1_bound_left+kernel_len)
            seq1=seqDict[line.Seq1][seq1_bound_left:seq1_bound_right][::-1]

            seq2_bound_left=max(0, int( line.j0_y))
            seq2_bound_right=min(len(seqDict[line.Seq2]), seq2_bound_left+kernel_len)
            seq2=seqDict[line.Seq2][seq2_bound_left:seq2_bound_right]


        gc=ExpectedGCContentForRow(seq1, seq2, kernel)
        if gc<0 or gc >1:
            print gc

        max_index=(len(conv_auto)-1)
        gc_bin=int(numpy.round( (gc*max_index)))
        cn_a,lcn_a,ucn_a=EstimateCopyNumber(line.count, conv_auto[gc_bin])



        new_row= line.row()+[gc]+[cn_a,lcn_a,ucn_a]
        outTable.writerow(new_row)
        outHandle.flush()
    inHandle.close()
    outHandle.close()

def ExpectedGCContentForRow(leftSeq, rightSeq, kernel):
    """Reconstructs the expected sequence of the junction and calculates the
    expected GC-content of a read pairs spanning the junction as the number of Gs
    and Cs in the flanking sequence weighted by the probability that they would
    be found in a read spanning the junction (described by the kernel)."""
    jxn_pos=len(leftSeq)
    fullseq=leftSeq.lower()+rightSeq.lower()
    seq_array=numpy.fromstring(fullseq, '|S1')
    gc_array=(seq_array=='g')+(seq_array=='c')
    kernel_mode=len(kernel)/2 + 1
    kernel_left_edge=max(0, kernel_mode-jxn_pos)
    kernel_right_edge=min(len(kernel)/2, len(rightSeq))
##    print kernel_left_edge, kernel_right_edge+kernel_mode
    kernel=kernel[kernel_left_edge:kernel_right_edge+kernel_mode]/kernel[kernel_left_edge:kernel_right_edge+kernel_mode].sum()

    gc=(gc_array[jxn_pos-kernel_mode+kernel_left_edge:jxn_pos+kernel_right_edge]*kernel).sum()

    return gc


def ConvolveSeries(distribution, length):
    series=[numpy.array( distribution)]
    for i in range(length-1):
        convolved=scipy.signal.fftconvolve(series[i], distribution)
        series.append(convolved)
    return(numpy.array(series))


def ConvolveAll(smoothed_dists, size=50):
    convolved_list=[]
    for i in range(len(smoothed_dists)):
        convolved_list.append(ConvolveSeries(smoothed_dists[i], size))
    return convolved_list

def EstimateCopyNumber(readCount, convolvedDist):

    #Create LkPlot
    p_dist=ComputeApproximatePosterior(readCount, convolvedDist, 50)
    estimate=numpy.argmax(p_dist)
    cdf=numpy.cumsum(p_dist)
    try:
        lower=numpy.where(cdf>=.01)[0][0]
        upper=numpy.where(cdf>=.99)[0][0]
    except:
        lower=numpy.nan
        upper=numpy.nan
    return estimate, lower, upper


def ComputeApproximatePosterior(x, conv,threshold=100, bias=1):
    pdf=[0.]
    avg=WeightedAverage(conv[0],numpy.arange(len(conv[0])))*bias
    sd=WeightedStdDev(conv[0],numpy.arange(len(conv[0])))
    exp=x/avg+1000
    if numpy.isnan(exp)==True:
        return numpy.array([0.])
    if numpy.isinf(exp)==True:
        print exp, avg
        return numpy.array([0.])
    for i in range(1, int( exp)):

        if i<threshold:
            if len(conv[i-1])>x:        #Convolution index 0 corresponds to CN=1
                pdf.append(conv[i-1][x])
            else: pdf.append(0.)
        else:
            avg_i=avg*(i)
            sd_i=sd*(i)**.5
            pdf.append(NormalPDF(x, avg_i, sd_i))
    pdf=numpy.array(pdf)
    return pdf/numpy.nansum(pdf)

def NormalPDF(x,mu, sigma):
    return numpy.exp(-.5* (x-mu)**2/(sigma)**2) /(sigma * (2*numpy.math.pi )**.5)

def main(argv):
    param={}
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if param=={}: return()


    spec_dict=ReadSpecificationFile(param['-spec'])
    reffile=spec_dict['Ref']
    consfile=spec_dict['Cons']
    refSeq=GetSeq(reffile)
    consSeq=GetSeq(consfile)
    seqDict=mergeDicts(refSeq, consSeq)
    file_list=os.listdir(param['-i'])
    count=0
    for file_path in file_list:

        count+=1
        if count==1: continue
##            f=directory.split('/')[-1]
        file_root=file_path.split('.')[0]
        if file_root.split('_')[-1]!='out': continue
##            if file_root.split('_')[1]=='2': continue
        if param.has_key('-head'):
            if file_root[0].lower()!=param['-head'].lower(): continue

##            print directory
        directory=param['-i']+"/"+ '_'.join(file_root.split('_')[:-1])
        infile=directory+'_out.tsv'
        outfile=directory+'__outCN.tsv'
##            directory='_'.join(file_root.split('_')[:2])
        dist_file=directory+'_kde.tsv'
        len_file= directory+'_len.tsv'
##            kde_file=directory+'_cvg_kde_auto.npy'
##            kde_file_X=directory+'_cvg_kde_x.npy'
        kde_file=directory+'_cvg_hist.npy'
        EstimateJunctionCopyNumber(infile, outfile, dist_file, len_file, seqDict, kde_file)


if __name__ == '__main__':
    main(sys.argv)
