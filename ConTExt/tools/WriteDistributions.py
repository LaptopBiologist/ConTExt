#-------------------------------------------------------------------------------
# Name:        WriteDistributions
# Purpose: Estimate the read length and gap distributions from ConTExt aligned
# output.
#
# Author:      Michael Peter McGurk
#
# Created:     06/11/2015
# Copyright:   (c) Michael Peter McGurk 2015

#-------------------------------------------------------------------------------
import csv
import numpy
import sys
import CoverageModeller

def HistoSD(histo):
    n=sum(histo)
    mean=0
    for i in range(len(histo)):
        mean+=float(i*histo[i])
    mean/=n

    SD=0
    for i in range(len(histo)):
        SD+=histo[i]*(mean-i)**2
    SD/=n
    sd=numpy.sqrt(SD)
    return(sd)

def HistogramToKDE_Silverman(histo):
    kde=[0.0]*len(histo)
    SD=HistoSD(histo)
    n=sum(histo)
    bandwidth=1.144*SD*n**(-.2)
    #bandwidth=100
    #print bandwidth
    for x in range(len(histo)):
        F=0.0
        for i in range (len(histo)):
            F+=histo[i]*Kernel(float(x-i)/bandwidth)
        F/=n*bandwidth
        kde[x]=F
    return(kde)


def EstimateGapDistAll(inDir, RefNames):
    dist={}
    dist['total']=0.0
    for ref in RefNames:
        inFile=inDir+'/'+ref+'.dat'
        handle=open(inFile, 'r')
        lib=csv.reader(handle, delimiter='\t')
        l,r=0,0

        header=True
        for row in lib:
            if header==True:
                header=False
                continue
            #pick ends:
            if row[0]!=row[6]:continue
            if row[1]==row[7]:continue
            if row[1]=='0':
                l3=int(row[3])
                r3=int(row[8])
                l5=int(row[2])
                r5=int(row[9])
            if row[1]=='16':
                r5=int(row[3])
                l5=int(row[8])
                r3=int(row[2])
                l3=int(row[9])
            gap=r3-l3
            insert_len=r5-l5
            if insert_len<1000 and insert_len>0:
                if dist.has_key(gap)==False:
                    dist[gap]=0.0
                dist[gap]+=1
                dist['total']+=1
        handle.close()
    G=[0.0]*1002
    for w in dist.keys():
        if w!='total':G[w]=dist[w]#/dist['total']
    return(G)


def EstimateLenDistAll(inDir, RefNames):
    dist={}
    dist['total']=0.0
    for ref in RefNames:
        inFile=inDir+'/'+ref+'.dat'
        handle=open(inFile, 'r')
        lib=csv.reader(handle, delimiter='\t')
        l,r=0,0
        print ref

        header=True
        for row in lib:
            if header==True:
                header=False
                continue
            #pick ends:
            if row[0]!=row[6]:continue
            if row[1]==row[7]:continue
            gap=abs(int(row[3])-int(row[2]))

            if dist.has_key(gap)==False:
                dist[gap]=0.0
            dist[gap]+=1
            dist['total']+=1

            gap=abs(int(row[8])-int(row[9]))
            if dist.has_key(gap)==False:
                dist[gap]=0.0
            dist[gap]+=1
            dist['total']+=1
        handle.close()
    G=[0.0]*1002
    for w in dist.keys():
        if w!='total':G[w]=dist[w]#/dist['total']
    return(G)


def EstimateDistAll(inDir, RefNames):
    gap_list=[]
    domain_list=[]
    len_dist={}
    dist={}
    dist['total']=0.0
    for ref in RefNames:
        print ref
        gap_dist={}
        gap_dist['total']=0.0
        inFile=inDir+'/'+ref+'.dat'
        handle=open(inFile, 'r')
        lib=csv.reader(handle, delimiter='\t')
        l,r=0,0
        print ref
        header=True
        for row in lib:
            if header==True:
                header=False
                continue
            #pick ends:
            if row[0]!=row[6]:continue
            if row[1]==row[7]:continue



            gap=abs(int(row[3])-int(row[2]))

            if dist.has_key(gap)==False:
                dist[gap]=0.0
            dist[gap]+=1
            dist['total']+=1

            gap=abs(int(row[8])-int(row[9]))

            if dist.has_key(gap)==False:
                dist[gap]=0.0
            dist[gap]+=1
            dist['total']+=1
            if row[1]=='0':
                l3=int(row[3])
                r3=int(row[8])
                l5=int(row[2])
                r5=int(row[9])
            if row[1]=='16':
                r3=int(row[2])
                l3=int(row[9])
                r5=int(row[3])
                l5=int(row[8])
            gap=r3-l3
            insert_len=r5-l5
            if insert_len<1000 and insert_len>0:
                if gap_dist.has_key(gap)==False:
                    gap_dist[gap]=0.0
                gap_dist[gap]+=1
                gap_dist['total']+=1

        handle.close()
        gap_set=set(gap_dist.keys())-set(['total'])
        min_gap=min(0, min( gap_set))
        print min_gap
        G=[0.0]*(1002-min_gap)

        domain=numpy.arange(min_gap,1002)
        domain_list.append(domain)
        for w in gap_dist.keys():
            if w!='total' :
                index=numpy.where(domain==w)[0][0]
                G[index]=gap_dist[w]#/dist['total']
        gap_array=numpy.array(G, float)
        gap_array/=gap_array.sum()
        gap_list.append(gap_array)
    len_hist=[0.0]*1002
    for w in dist.keys():
        if w!='total':len_hist[w]=dist[w]#/dist['total']
    len_array=numpy.array(len_hist, float)
    len_array/=len_array.sum()
    #Align domains
    min_dom_list=[min(i) for i in domain_list ]
    full_domain=numpy.arange(min( min_dom_list), 1002)
    new_gap_list=[]
    for i in range(len( domain_list)):
        new_gap=numpy.array( [0.]*len(full_domain))
        min_index=numpy.where(full_domain==min_dom_list[i])[0][-1]
        new_gap[min_index:]=numpy.array( gap_list[i])
        new_gap_list.append(list(new_gap))
    return(new_gap_list, len_array, full_domain )

def Kernel(u):
    K=numpy.exp(-(u**2)/2)/(2*numpy.pi)**.5
    return(K)
def BayesTransform(gap_distribution):
    transformed=[]
    scale=sum([l*gap_distribution[l] for l in range(len(gap_distribution))])
    for size in range(len(gap_distribution)):
        probability=size*gap_distribution[size]/scale
        transformed.append(probability)

    return transformed


def WriteDistributions(inDir,outfile='', names=['2L', '2R', '3L', '3R']):
##    outfile='{0}.tsv'.format(inDir)
    if outfile=='': outfile='{0}.tsv'.format(inDir)
    print "Getting insert size distribution..."
    gap_list, len_dist, domain=EstimateDistAll(inDir, names)
    gap_list=numpy.array(gap_list)
    zero_indices=(gap_list==0)
    nonzero_min=numpy.min(gap_list[~zero_indices])
    gap_list[zero_indices]=nonzero_min
    print len(domain)
    gap_array=numpy.vstack(gap_list).transpose()
    bw= CoverageModeller_9.CV_4( gap_array,0,100,10, 3,fold=100, cutoff=.99)
    mpd=gap_array.sum(1)
    kde=numpy.array( CoverageModeller_9.FastHistoKDE_BW(mpd, bw))
    Pr=[0.]*(domain<0).sum()+list( BayesTransform(kde[domain>=0]))

    totalReads=sum(mpd)

    imageName='.'.join(outfile.split('.')[:-1])+'_kde.tsv'
    imageHandle=open(imageName,'w')
    imageWrite=csv.writer(imageHandle, delimiter='\t')
    imageWrite.writerow(mpd )
    imageWrite.writerow(kde)
    imageWrite.writerow(Pr)
    imageWrite.writerow(domain)
    imageHandle.close()

    imageName='.'.join(outfile.split('.')[:-1])+'_len.tsv'
    imageHandle=open(imageName,'w')
    imageWrite=csv.writer(imageHandle, delimiter='\t')
    imageWrite.writerow(len_dist )
    imageWrite.writerow(len_dist)
    imageWrite.writerow(len_dist)
    imageHandle.close()



def WriteLenDistributions(inDir, outfile='', names=['2L', '2R', '3L', '3R']):
    if outfile=='': outfile='{0}.tsv'.format(inDir)
    print "Getting insert size distribution..."

    gap_list, len_dist=EstimateLenDistAll(inDir, names)
    gap_array=numpy.vstack(gap_list).transpose()

    bw= CoverageModeller_9.CV_4(gap_array, numpy.arange(0,40,.5),cutoff=.99)
    kde=CoverageModeller_9.FastHistoKDE_BW(gap_list)

    Pr=BayesTransform(kde)
    totalReads=sum(mpd)

    imageName='.'.join(outfile.split('.')[:-1])+'_len.tsv'
    imageHandle=open(imageName,'w')
    imageWrite=csv.writer(imageHandle, delimiter='\t')
    imageWrite.writerow([m/totalReads for m in mpd])
    imageWrite.writerow(kde)
    imageWrite.writerow(Pr)
    imageHandle.close()

def main(argv):
    param={}
    print argv
##
##    if argv[1]=='-help':
##        print """I am currently too lazy to write instructions."""
##        return
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param

    if param=={}:
        return
    if param.has_key('-L')==False:
        WriteDistributions(param['-i'])
    else: WriteLenDistributions(param['-i'])

    #inDir, outDir, AssemblyIndex, TEIndex, threads, length=70, ReadIDPosition=1, PositionRange=0, phred=33, shorten=True, Trim=True



if __name__ == '__main__':
    main(sys.argv)
