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

# General Feature:
# Want to allow the user to input a samples.ignore file
#
#
#
#
#
#
#
#
#
#

import csv
import numpy
import numpy as np
import os
import scipy
import scipy.signal
import sys
import pickle
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

def ExpectedReadDepth(hist_table):
    exp_cvg=[]
    for i in range(hist_table.shape[0]):
        exp_cvg.append(WeightedAverage(hist_table[i,:],numpy.arange(len(hist_table[i,:]))))
    return numpy.array(exp_cvg)

def EstimateRepeatAbundance(infile, outfile, consFile, kdefile, lenFile, insertFile):
##    m,k,p=ReadDistributions(lenFile)
##    read_len=2*WeightedAverage(m, range(len(m)))
##    m,k,p=ReadDistributions(insertFile)
    kernel=CoverageModeller.ConstructKernel(insertFile, lenFile)
##    ins_len=1*WeightedAverage(k,range(len(k)))
##    read_length=read_len+ins_len
    cvg=numpy.load(kdefile)
    kde_root='_'.join(infile.split('_')[:-1]) +"_cvg_kde"
    if os.path.exists(kde_root+'_auto.npy')==False:
        kde, kde_x=CoverageModeller.SmoothCvgByGC(cvg)
##        kde_root=inDir+"_cvg_kde"
        numpy.save(kde_root+'_auto.npy', kde)
    else:
        kde=numpy.load(kde_root+'_auto.npy')
##    numpy.save(kde_root+'_x.npy', kde_x)
    kde=numpy.array(kde)
##    exp_cvg=ExpCvg( kde)

    consSeq=GetSeq(consFile)
##    cumul_dist=numpy.cumsum(k)


    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    outhandle=open(outfile, 'w')
    countDict={}
    cons_dict={}
##    refNames=['X', '2L','3L', '2R', '3R']
    row=intable.next()
    exp_cvg=ExpectedReadDepth(kde)
    abundance_dictionary={}
    for sequence in consSeq.keys():
        abundance_dictionary[sequence]=0.
    for row in intable:

        line=cluster(row)
        if line.feature!='Consensus': continue
        if consSeq.has_key(line.Seq1)==False: continue
        cons_dict[line.Seq1]=numpy.array([0.]*len(consSeq[line.Seq1]))
        cons_dict[line.Seq1].fill(0.)
        for i in range( len(line.x_list)):
            l,r=line.y_list[i], line.x_list[i]
            cons_dict[line.Seq1][l:r-1]+=1.
##        if key=='':continue
        seq_array=numpy.fromstring ( consSeq[line.Seq1].lower(), '|S1')
        gc_array=(seq_array=='g')+(seq_array=='c')
    ##        mean_gc=numpy.mean(gc_array)
        if len(cons_dict[line.Seq1])>max(600, len(kernel)):

            exp_GC=numpy.convolve(gc_array, kernel, 'same')
            rounded_gc=numpy.asarray(100*exp_GC, int)
            rounded_gc[rounded_gc>100]=100

            expected_coverage=exp_cvg[rounded_gc]
            estimated_CN=cons_dict[line.Seq1][300:-300]/expected_coverage[300:-300]

            estimated_abundance=estimated_CN.sum()/(len(estimated_CN)/float(len(cons_dict[line.Seq1])))
            cons_dict[line.Seq1]=(cons_dict[line.Seq1], estimated_CN, expected_coverage)
            abundance_dictionary[line.Seq1]=estimated_abundance
    pickle.dump(cons_dict, outhandle)
    outhandle.close()
    return abundance_dictionary
##
##        if key=='':continue
##        seq_array=numpy.fromstring ( consSeq[f].lower(), '|S1')
##        gc_array=(seq_array=='g')+(seq_array=='c')
##    ##        mean_gc=numpy.mean(gc_array)
##        if len(cvg_array)>len(kernel):
##
##            exp_GC=numpy.convolve(gc_array, kernel, 'same')
##            rounded_gc=numpy.asarray(100*exp_GC, int)
##            exp_cvg=numpy.array( ExpCvg(kde))
##            good_indices=(exp_GC<100)
##            print len( cvg_array),len( good_indices),len( rounded_gc)
##            CN_array=cvg_array[good_indices]/exp_cvg[rounded_gc[good_indices]]
##
##
##        outhandle.write('{0}\t{1}\t{2}\n'.format(key, numpy.mean (gc_array) , numpy.nansum( CN_array)))
##    else: outhandle.write('{0}\t{1}\t{2}\n'.format(key, numpy.mean(gc_array) , numpy.nan))
##    outhandle.flush()

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

def SplitConTEXtOutput(infile, outDir):
    #Create the output directory if it doesn't already exist
    MakeDir(outDir)
    ID_Header=['Subject', 'Target', 'Quad X', 'Quad Y', 'Feature Type']
    Jxn_Header=['Jxn X', 'Jxn X (range)','Jxn Y', 'Jxn X (range)']
    Position_Header=["Read Count","Reads 3' X", "Reads 3' Y"]
    #Open the input file and read as a tab-delimited table
    inhandle=open(infile, 'r')
    table=csv.reader(inhandle, delimiter='\t')

    #Create a dictionary that will hold open output files
    #all files will be opened for appending.
    #Becauswe python puts a limit on how many files will be opened at once,
    #we will have three objects:

    #PathDictionary -- Stores the paths to all output files
    #HandleDictionary -- Stores the handles for open output files
    #StackList -- Stores the keys to open output files in the order that they
                    #were opened, with the most recent at the end of the list.
    PathDictionary={}
    HandleDictionary={}
    WriterDictionary={}
    StackList=[]
    fileHeader=True

    for row in table:
        if fileHeader==True:
            fileHeader=False
            continue
        line=cluster(row)
        Entry1=line.Seq1
        Entry2=line.Seq2

        #Only true if this is the first time an entry is encountered
        #If so, write a header row.
        Header1=False
        Header2=False

        #Create a new output path if the repeat family doesn't have one

        if PathDictionary.has_key(Entry1)==False:
            newPath=outDir + '/' + Entry1 + '.tsv'
            PathDictionary[Entry1]=newPath
            Header1=True

        if HandleDictionary.has_key(Entry1)==False:
            HandleDictionary[Entry1]=open(PathDictionary[Entry1], 'a')
            WriterDictionary[Entry1]=csv.writer(HandleDictionary[Entry1], delimiter='\t')
            StackList.append(Entry1)

            if len(StackList)>100:
                #If more than 100 output files are open, close the oldest and remove it from the stacks of open files
                HandleDictionary[StackList[0]].close()
                del HandleDictionary[StackList[0]]
                del WriterDictionary[StackList[0]]
                del StackList[0]
        if Header1==True:
            WriterDictionary[Entry1].writerow(ID_Header+Jxn_Header+Position_Header)
        WriterDictionary[Entry1].writerow(line.row())

        if Entry1!=Entry2:

            if PathDictionary.has_key(Entry2)==False:
                newPath=outDir + '/' + Entry2 + '.tsv'
                PathDictionary[Entry2]=newPath
                Header=True

            if HandleDictionary.has_key(Entry2)==False:
                HandleDictionary[Entry2]=open(PathDictionary[Entry2], 'ab')
                WriterDictionary[Entry2]=csv.writer(HandleDictionary[Entry2], delimiter='\t')
                StackList.append(Entry2)

                if len(StackList)>100:
                    HandleDictionary[StackList[0]].close()
                    del HandleDictionary[StackList[0]]
                    del WriterDictionary[StackList[0]]
                    del StackList[0]
            if Header2==True:
                WriterDictionary[Entry2].writerow(ID_Header+Jxn_Header+Position_Header)
            WriterDictionary[Entry2].writerow(line.inverse_row())


    #Close all output files
    for entry in StackList:
        HandleDictionary[entry].close()

    inhandle.close()
    return StackList

#------------------------------------------------------------------------------#
#
#Functions for building SNP pileups
#
#------------------------------------------------------------------------------#


def PileupConsensusSequences(indir,  outdir,in_tail, seq_file,sample_pos=1, rep_pos=-1 ):
    MakeDir(outdir)
    files=os.listdir(indir)
    image_directory={}
    rpt_dict=GetSeq(seq_file)
    snp_dict={}
    #Fill the dictionary firs
    for rpt in sorted( rpt_dict.keys()):
        length=len(rpt_dict)+1
        snp_dict[rpt]={}
        for f in sorted( files ):

            filename=f.split('/')[-1]
            root=filename.split('.')[0]
            tail=root.split('_')[-1]
            name='_'.join(root.split('_')[:-1])
            if tail!=in_tail: continue
            root_parts= root.split('_')
            sample=root_parts[sample_pos]
            if rep_pos!=-1:
                rep=root_parts[rep_pos]
            else:
                rep='1'

            snp_dict[rpt][sample]=numpy.ndarray((5,length))
            snp_dict[rpt][sample].fill(0.)
    sample_file=outdir+'/sample.tsv'
    sample_handle=open(sample_file, 'w')
    sample_table=csv.writer(sample_handle, delimiter='\t')
    for f in sorted( files ):

        #Only process files with the appropriate tag
        filename=f.split('/')[-1]
        root=filename.split('.')[0]
        tail=root.split('_')[-1]
        name='_'.join(root.split('_')[:-1])
        if tail!=in_tail: continue
        print name

        root_parts= root.split('_')
        sample=root_parts[sample_pos]
        if rep_pos!=-1:
            rep=root_parts[rep_pos]
        else:
            rep='1'
        if rep=='2': continue
##        continue

        inhandle=open(indir+'/'+f, 'r')
        intable=csv.reader(inhandle, delimiter= '\t')
        intable.next()
        sample_table.writerow([sample])
        for row in intable:
            try:
                line=cluster(row)
            except:
                print row
                print jabber
            #Delete this tomorrow
            if line.Seq1=='L1HS' or line.Seq1=='AluYb10': continue
            if line.feature=='Consensus':
##                print line.ID
                try:
                    length=len(rpt_dict[line.Seq1])+1
                except:
                    #Not in the repeat index:
                    continue

                snp_dict[line.Seq1][sample]=(PileUp(line, length=length))
##                break
        inhandle.close()
##        print len(snp_list)
    #Write the output files
    for repeat in sorted( snp_dict.keys()):
        repeat_array=[]
        for sample in sorted(snp_dict[repeat].keys()):
            repeat_array.append(snp_dict[repeat][sample])
        outfile='{0}/{1}_snps.npy'.format(outdir, repeat)
        numpy.save(outfile, numpy.array(repeat_array))

def PileUp(clust, Seq='', length=8119):
    SNP_array=numpy.ndarray((5,length))
    SNP_array.fill(0.)
    quad_1=clust.Quad1
    quad_2=clust.Quad2
    cvg=numpy.array([0.]*length)
    count=0.
    cig_list,md_list, pos_list=[],[],[]
    for i in range(clust.count):
        if clust.MapQY[i]<20 and clust.MapQX[i]<20: continue
        x_3=clust.x_list[i]
        y_3=clust.y_list[i]
        try:
            cigX=clust.CigarX[i]
            cigY=clust.CigarY[i]
        except:
            print clust.count
            print len(clust.CigarX)
            print jabber
        qdX=clust.QD_X[i]
        qdY=clust.QD_Y[i]
        #Maybe need to reverse?

        parsed_cigX=ParseCIGAR(cigX)
        len_x=parsed_cigX['D']+parsed_cigX['M']
        cig_array_X=ExpandCIGAR(cigX)
        matched_indices_X=numpy.where(cig_array_X[cig_array_X!='I'] =='M')[0]
        md_array_X=ExpandMD(clust.MD_X[i])

        parsed_cigY=ParseCIGAR(cigY)
        len_y=parsed_cigY['D']+parsed_cigY['M']
        cig_array_Y=ExpandCIGAR(cigY)
        matched_indices_Y=numpy.where(cig_array_Y[cig_array_Y!='I'] =='M')[0]
        md_array_Y=ExpandMD(clust.MD_Y[i])


        if md_array_X=='Failed' or md_array_Y=='Failed': continue
        if quad_1=='+':
            x_5=x_3-len_x
        else:
            x_5=x_3+len_x
        if quad_2=='+':
            y_5=y_3-len_y
        else:
            y_5=y_3+len_y
        x_l=min(x_3, x_5)
        x_r=max(x_3, x_5)
        y_l=min(y_3, y_5)
        y_r=max(y_3, y_5)
        cvg[x_l:x_r]+=1.
        cvg[y_l: y_r]+=1.
        try:
            if qdY=='':

                SNP_array[:,y_l:y_r][md_array_Y[5:-5], matched_indices_Y[5:-5]]+=1

            else:
    ##            try:
                mask_ind=set(numpy.asarray(qdY.split(','),int))
                con_ind=set(range(0,len(matched_indices_Y)))
                if max(mask_ind)>max(con_ind): print max(mask_ind),max(con_ind)
                ind= numpy.array(sorted( list( con_ind-mask_ind)))
                if len(ind)>10:
##                    SNP_array[:,x_l:x_r][md_array_X[ind[5:-5]], matched_indices_X[ind[5:-5]]]+=1
                    SNP_array[:,y_l:y_r][md_array_Y[ind[5:-5]], matched_indices_Y[ind[5:-5]]]+=1

        except:

            z=1


        if qdX=='':

            try:
                SNP_array[:,x_l:x_r][md_array_X[5:-5], matched_indices_X[5:-5]]+=1
            except:
                print md_array_X.shape

                print matched_indices_X.shape
                print SNP_array[:,x_l:x_r+1].shape
                print clust.MD_X[i]
                print cigX
                print '=--------------------------------------------'

        else:
##                try:
            mask_ind=set(numpy.asarray(qdX.split(','),int))

            con_ind=set(range(0,len(matched_indices_X)))
            ind= numpy.array(sorted( list( con_ind-mask_ind)))
            if max(mask_ind)>max(con_ind): print max(mask_ind),max(con_ind)
            try:
                if len(ind)>10:
                    SNP_array[:,x_l:x_r][md_array_X[ind[5:-5]], matched_indices_X[ind[5:-5]]]+=1
            except:
                print x_l, x_r
                print clust.MD_X[i]
                print mask_ind
                print md_array_X.shape
                print matched_indices_X.shape
                print ind.shape

    print count
    print SNP_array.sum()
    return SNP_array

def ParseCIGAR(CIGAR):
    """Reads a sam CIGAR string to count deletions and insertions and matches"""

    parts={'M':0, 'I':0,'D':0}
    cigar=np.array(list(CIGAR))

    m=np.where(cigar=='M')
    i=np.where(cigar=='I')
    d=np.where(cigar=='D')
    M=['M']*len(m[0])
    I=['I']*len(i[0])
    D=['D']*len(d[0])
    indices=np.concatenate((m,i,d),-1)
    values=np.concatenate((M,I,D),-1)
    pairs=[]
    for w in range(len(indices[0])):
        pairs.append((indices[0][w], values[w]))
    pairs=sorted(pairs, key=lambda x: x[0])
    last=0
    for p in pairs:
        val=int(CIGAR[last:p[0]])
        parts[p[1]]+=val
        last=1+p[0]
    return(parts)

def ExpandMD(MD_string):
    """Consensus: 0
    A: 1
    T: 2
    C: 3
    G: 4 """
    if MD_string.count('N')!=0:  return 'Failed'
    nt={'A':1, 'T':2, 'C':3,'G':4}
    hold='0'
    null=False
    expanded_md=[]
    for char in MD_string:
##        print null
        if null==True:
            if nt.has_key(char)==False:
                null=False
            else: continue
##            continue
        if char=='^':
            null=True
            expanded_md+=[0]*int(hold)
            hold='0'
            continue
        if nt.has_key(char)==False:
            hold+=char
        else:
##            print hold

            expanded_md+=[0]*int(hold)
            expanded_md+=[nt[char]]
            hold='0'
    expanded_md+=[0]*int(hold)
    return numpy.array( expanded_md )
def ExpandCIGAR(CIGAR):
    """Reads a sam CIGAR string to count deletions and insertions and matches"""

    parts={'M':0, 'I':0,'D':0}
    cigar=np.array(list(CIGAR))

    m=np.where(cigar=='M')[0]
    i=np.where(cigar=='I')[0]
    d=np.where(cigar=='D')[0]
##    M=['M']*len(m[0])
##    I=['I']*len(i[0])
##    D=['D']*len(d[0])
    indices=sorted(list( np.concatenate((m,i,d),-1)))
    indices=[-1]+indices
##    print indices
    CIGAR_list=[]
    for i in range(len(indices)-1):
        chunk= CIGAR[indices[ i]+1 :indices[ i+1]+1]
        count, val=int(chunk[:-1]),chunk[-1]
        CIGAR_list+=[val]*count
    return numpy.array( CIGAR_list)



def main(argv):
    param={}
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param
    if param=={}: return()
    if param.has_key('-split'):
        outdir='.'.join(param['-split'].split('.')[:-1])+'_split'
        SplitConTEXtOutput(param['-split'], outdir)
        return()

    spec_dict=ReadSpecificationFile(param['-spec'])
    reffile=spec_dict['Ref']
    consfile=spec_dict['Cons']
    refSeq=GetSeq(reffile)
    consSeq=GetSeq(consfile)
    seqDict=mergeDicts(refSeq, consSeq)
    file_list=os.listdir(param['-i'])
    count=0

    if param.has_key('-snp')==True:
        PileupConsensusSequences(param['-i'], param['-o'], param['-tail'], consfile)
        return()
    sample_dict={}
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
        sample='_'.join( infile.split('/')[-1].split('_')[:-1])
        outfile=directory+'__outCN.tsv'
##            directory='_'.join(file_root.split('_')[:2])
        dist_file=directory+'_kde.tsv'
        len_file= directory+'_len.tsv'
##            kde_file=directory+'_cvg_kde_auto.npy'
##            kde_file_X=directory+'_cvg_kde_x.npy'
        kde_file=directory+'_cvg_hist.npy'
        if param.has_key('-count')==True:

            outfile=directory+'__pickle.dict'
            sample_dict[sample]= EstimateRepeatAbundance(infile, outfile, consfile, kde_file, len_file, dist_file)
        else:
            outfile=directory+'__outCN.tsv'
            EstimateJunctionCopyNumber(infile, outfile, dist_file, len_file, seqDict, kde_file)
    if sample_dict!={}:
        outhandle=open(param['-i']+'/count_table.tsv', 'w')
        outtable=csv.writer(outhandle, delimiter='\t')
        ref_sample=sample_dict.keys()[0]
        header=['Repeat Family']+sorted( sample_dict.keys())
        outtable.writerow(header)
        for repeat in sorted( sample_dict[ref_sample].keys()):
            repeat_row=[repeat]
            for sample in sorted( sample_dict.keys()):
                repeat_row.append(sample_dict[sample][repeat])
            outtable.writerow(repeat_row)

        outhandle.close()

if __name__ == '__main__':
    main(sys.argv)
