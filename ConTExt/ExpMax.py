#-------------------------------------------------------------------------------
# Name:        K-Jxn clustering
# Purpose:
#
# Author:      Michael Peter McGurk
#
# Created:     25/08/2015
# Copyright:   (c) Michael Peter McGurk 2015
# License:
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#-------------------------------------------------------------------------------

"""The distribution of read pairs in an image arises from the set of junctions
present in the genome sequenced; each jxn will generate read pairs that are distributed in
a structured manner: the distance between a read pair and the junction it spans
follows the insert size distribution. So, the distribution of read pairs in an
image is generated from a mixture of many component distributions. The aim then
is to assign each read pair to a component; thus, for each read pair there exists
a latent variable that represents the component from which it arose. The aim then
is two-fold:

    1.) Estimate the location and read count of each component distribution
    2.) Assign each read pair to a component distribution

Expectation-Maximization is a successive approximation strategy often successful
when applied to this sort of problem. However, the algorithm can be extremely
operation intensive. But, for mixtures of Gaussian distributions, the equations
can be efficiently computed. While the distribution of read pairs spanning a junction
is clearly not normal (there a exists vector along which the distribution may be approximately
normal, but orthogonal to that vector, the distribution is a stove-top distribution), some
of the distribution's structure can be approximated using a bivariate Gaussian
with the same covariance matrix as the expected read pair distribution. Then
the EM algorithm only needs to approximate the mean and weight of each component
distribution. This is in essence a K-means algorithm that assumes non-round clusters.

A further advantage is that constraining the shape of the component distributions with
a fixed covariance matrix appears to avoid overfitting--choosing a model with more
components than the generative distribution doesn't generally yield spurious clusters
because the means of extraneous components still find local maxima; when it happens that
multiple components have the same mean, then only the component with the highest
weight will be assigned read pairs."""

import matplotlib
import sys
if sys.platform[0]=='l':
    matplotlib.use('Agg')
from matplotlib import pyplot
import seaborn

#Set plotting defaults

pyplot.ioff()

import time
import math
import numpy
import numpy as np
import scipy
import csv
import collections
import random
import copy
import os


from scipy import stats
from itertools import chain

#Imports for Scikit-Learn
import sklearn
import sklearn.gaussian_process
import sklearn.metrics
from sklearn import gaussian_process


from Bio import SeqIO
#imports for plotting


import tools.CoverageModeller as CoverageModeller
import tools.WriteDistributions as WriteDistributions
import tools.GeneralizedLogisticRegression as GeneralizedLogisticRegression
import tools.PlotsForContext as PlotsForContext

from tools.General_Tools import *

#Set figure defaults
matplotlib.rcParams['font.size']=14


class ConTExtLine():
    strand={'0':'+', '16':'-', '4':'*'}
    reverse_strand={'+':'0', '-':'16', '*':'4'}
    def __init__(self, row):

        self.Seq1=row[0]
        self.Strand1=ConTExtLine.strand[row[1]]
        self.Start1=int(row[2])
        self.End1=int(row[3])
        self.Seq2=row[6]

        self.Read1=row[4]
        self.Read2=row[10]

        self.Qual1=row[5]
        self.Qual2=row[11]
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
        self.MD1=row[16]
        self.MD2=row[17]
        self.QD1=row[18]
        self.QD2=row[19]
        if len(row)>19:self.optional=row[19]
        else: self.optional=''

    def row(self):
        row=[self.Seq1, ConTExtLine.reverse_strand[ self.Strand1], self.Start1, self.End1, self.Read1, self.Qual1, self.Seq2, ConTExtLine.reverse_strand[self.Strand2], self.Start2, self.End2,self.Read2,self.Qual2, self.Cigar1, self.Cigar2, self.MapQ1, self.MapQ2, self.MD1, self.MD2, self.QD1, self.QD2]+[self.optional]
        return row

    def inverse_row(self):
        row=[self.Seq2, ConTExtLine.reverse_strand[ self.Strand2], self.Start2, self.End2, self.Read2,self.Qual2, self.Seq1,ConTExtLine.reverse_strand[ self.Strand1], self.Start1, self.End1,self.Read1,self.Qual1, self.Cigar2, self.Cigar1, self.MapQ2, self.MapQ1, self.MD2, self.MD1, self.QD2, self.QD1]+[self.optional]
        return row

class ModelParam:
    def __init__(self, mean, wt, cov):
        self.weight=wt
        self.mean=numpy.asarray( mean, 'float64')
        self.cov=numpy.asarray( cov, 'float64')


    def __repr__(self):
        return "{0}. {1}, {2}".format( self.weight, self.mean, self.cov)
    def __sub__( minuend, subtrahend):
        diff_weight=minuend.weight-subtrahend.weight
        diff_mean=minuend.mean-subtrahend.mean
        diff_cov=minuend.cov-subtrahend.cov
        return ModelParam(diff_mean, diff_weight, diff_cov)

    def __add__( minuend, subtrahend):
        diff_weight=minuend.weight+subtrahend.weight
        diff_mean=minuend.mean+subtrahend.mean
        diff_cov=minuend.cov+subtrahend.cov
        return ModelParam(diff_mean, diff_weight, diff_cov)
    def __mul__(multiplier, multiplicand):
        prod_weight=multiplier.weight*multiplicand[0]
        prod_mean=multiplier.mean*multiplicand[1]
        prod_cov=multiplier.cov*multiplicand[2]
        return ModelParam(prod_mean, prod_weight, prod_cov)

class Cluster:
    def __init__(self, data_array, line_array):
        self.read_positions=data_array
        self.line_list=line_array
##        self.consensus=consensus

    def truncateData(self, percentile, quadrant):
        indices=[]
        for i in range(2):
            if quadrant[i]=='+':
                X_pos=self.read_positions[i,:]
                ##            X_cutoff=numpy.percentile(X_pos, 100-percentile, interpolation='higher')
                X_cutoff=numpy.percentile(X_pos, 100-percentile)
                if len(self.read_positions[i,:])<=10:
                    X_cutoff=numpy.max(X_pos)
                indices.append(numpy.where(self.read_positions[i,:]<=X_cutoff)[0])
            if quadrant[i]=='-':
                X_pos=self.read_positions[i,:]
                ##            X_cutoff=numpy.percentile(X_pos, percentile, interpolation='lower')
                X_cutoff=numpy.percentile(X_pos, percentile)
                if len(self.read_positions[i,:])<=10:
                    X_cutoff=numpy.min(X_pos)
                indices.append(numpy.where(self.read_positions[i,:]>=X_cutoff)[0])
        new_indices=numpy.intersect1d(indices[0], indices[1])
        self.truncatedReadPositions=self.read_positions[:,new_indices]

    def FindClosestJunction(self, quadrant):
        self.min_junction=[]
        for i in range(2):
            if quadrant[i]=='+':
                self.min_junction.append(numpy.max(self.truncatedReadPositions[i,:]))
            if quadrant[i]=='-':
                self.min_junction.append(numpy.min(self.truncatedReadPositions[i,:]))
        pass



#----------------------------------------------------------------------------#
#General Functions

def EstimateJunctionLocation(data, dist, quadrant):

    """Estimatest the location of a junction as the average position
    of the reads in the cluster shifted by one half the expected gap size
    along each axis. The direction of the offset is determined by the strands
    involved."""
    if quadrant[0]=='+':
        x_mod=1
    else:
        x_mod=-1
    if quadrant[1]=='+':
        y_mod=1
    else:
        y_mod=-1
    mean_dist=WeightedAverage(dist, numpy.arange(len(dist)))
    x,y=numpy.array(data[0]), numpy.array(data[1])
    avg_avg=numpy.mean(x+x_mod*mean_dist/2), numpy.mean(y+y_mod*mean_dist/2)
    return avg_avg

def MakeDir(newdir):
    if os.path.exists(newdir)==False:
        os.mkdir(newdir)



def CleanName(name):
    illegal=['|', '!','>', '<', '?','/','*']
    for i in illegal:
        name='_'.join(name.split(i))
    return(name)


def Cumulative(dist):
    CDF=[dist[0]]
    for i in range(1,len(dist)):
        CDF.append( CDF[i-1]+dist[i])
    return numpy.array( CDF)

def Simpson(data, is_count=True):
    count_dict=collections.Counter(data)
    data_array=numpy.array(count_dict.values(), float)
    if is_count==True:
        fr_array=data_array/data_array.sum()
    else:
        fr_array=data_array
##    diversity=(fr_array**2).sum()
    diversity=max(fr_array)
    return diversity

def HarmonicMean(data):
    data_array=numpy.array(data)
    return len(data)/(1./data_array).sum()

def DecomposeMatrix(cov):
    """Decomposes 2D covariance matrix into its eigenvectors and eigenvalues"""
    cov_matrix= numpy.matrix(cov)
    eig_val, eig_vec=numpy.linalg.eig(cov_matrix)
    eig_val_diag=[[eig_val[0],0],[0,eig_val[1]]]
    return eig_vec, eig_val_diag

def RecomposeMatrix(eig_vec, eig_val_matrix):
    """Returns a covariance matrix from eigenvectors and eigenvalues"""
    return eig_vec*eig_val_matrix*numpy.linalg.inv(eig_vec)


def PowerFunction(x,c):
    return c*(x**-1)

def ExpandCIGAR(CIGAR):
    """Reads a sam CIGAR string to count deletions and insertions and matches"""

    parts={'M':0, 'I':0,'D':0}
    cigar=np.array(list(CIGAR))

    m=np.where(cigar=='M')[0]
    i=np.where(cigar=='I')[0]
    d=np.where(cigar=='D')[0]

    indices=sorted(list( np.concatenate((m,i,d),-1)))
    indices=[-1]+indices

    CIGAR_list=[]
    for i in range(len(indices)-1):
        chunk= CIGAR[indices[ i]+1 :indices[ i+1]+1]
        count, val=int(chunk[:-1]),chunk[-1]
        CIGAR_list+=[val]*count
    return numpy.array( CIGAR_list)

def ReadDistributions(infile, cutoff=.005, return_full=False):
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



#----------------------------------------------------------------------------#
#Monte Carlo Simulations

def GenSample(insert_size_probability, read_count=5):

    """Simulates a read pair distribution of a given read count from the provided
    insert size distribution."""

    value_range=numpy.array( range( len(insert_size_probability)))
    x_list, y_list=[],[]
    Pr=numpy.array( insert_size_probability)/sum(insert_size_probability)
    insert_vals= numpy.random.choice(value_range, read_count, p=Pr)
    for val in insert_vals:
        x_val=numpy.random.random_integers(0, val)
##        if 50< x_val <80: continue
        y_val=val-x_val
        x_list.append(x_val)
        y_list.append(y_val)
    return x_list, y_list

def DrawFromReadDistanceDistribution(gap_dist, read_count, sample_num):
    samples=[]
    for i in range(sample_num):
        cluster=GenSample(gap_dist, read_count)
        marginal_x=sorted( cluster[0])
        marginal_y=sorted(cluster[1])
        samples.append(max(numpy.diff(marginal_x)))
        samples.append(max(numpy.diff(marginal_y)))
    return samples

def SimulateDataFromPoints(distribution, x_,y_, cov_min, cov_max, max_clusters=20):
    distribution=numpy.copy(distribution)
    x_list, y_list,label_list =[],[],[]
    cdf=numpy.cumsum(distribution)
    low,upp=numpy.where(cdf>.005)[0][0], numpy.where(cdf>=.995)[0][0]
    distribution[:low]=0.
    distribution[upp:]=0.
    distribution/=distribution.sum()
    for i in range(max_clusters):
        coverage=(max(cov_min, min( scipy.stats.binom.rvs(100,.20),cov_max)))
        x_init,y_init=GenSample(distribution, coverage)
        jxn_x= x_[i]
        jxn_y= y_[i]
        x_list+=[x+jxn_x for x in  x_init]
        y_list+=[x+jxn_y for x in  y_init]
        label_list+=[i]*coverage
    return numpy.vstack((x_list, y_list)).transpose(), label_list

def SimulateErrorByReadCount(dist, max_read_count=1000, rep_size=1000, regfile=''):
    """Returns the coefficients of power functions fitted to error in junction
    location estimates along the first and second eigenvectors as a function of
    read count."""

    mean_dist=WeightedAverage(dist, numpy.arange(len(dist)))
    X_list, Y_list=[],[]
    comp_1, comp_2, sample_size=[],[],[]
    for t in range (10):
        for rep in range(rep_size):
            x,y=GenSample(dist,max_read_count)
            X_list.append(x)
            Y_list.append(y)
        X_array=numpy.array(X_list)-mean_dist/2
        Y_array=numpy.array(Y_list)-mean_dist/2
        for count in range(1, max_read_count):
            x_mean=numpy.mean(X_array[:, :count],1)
            y_mean=numpy.mean(Y_array[:, :count],1)
            data_count=[x_mean, y_mean]
            cov=numpy.cov(data_count)
    ##        print cov
    ##        print  cov[0][0], cov[1][1], cov[0][1], cov[1][0]
            cov_a=numpy.mean([ cov[0][0], cov[1][1]])
            cov_b=numpy.mean([cov[0][1], cov[1][0]])
            new_cov=numpy.array( [[cov_a, cov_b], [cov_b, cov_a]])
            mat, val=DecomposeMatrix(new_cov)

            comp_1.append(val[0][0])
            comp_2.append(val[1][1])
            sample_size.append(count)
    #Fit power function
    p1= scipy.optimize.curve_fit(PowerFunction, comp_1, sample_size)[0]
    p2= scipy.optimize.curve_fit(PowerFunction, comp_2, sample_size)[0]
    if regfile!='':
        reg_handle=open(regfile, 'wb')
        reg_table=csv.writer(reg_handle, delimiter='\t')
        reg_table.writerow(['alpha_1', p1[0]])
        reg_table.writerow(['alpha_2', p2[0]])
        for i in range(len(sample_size)):
            line=[sample_size[i], comp_1[i], comp_2[i]]
            reg_table.writerow(line)
##    print p1
##    print p2
##    return comp_1, comp_2, sample_size, p1, p2
        reg_handle.close()
    return numpy.array([p1,p2])



#----------------------------------------------------------------------------#
#Expectation Maximization

def SquareEM(data, cov, components=60, iterations=3, Agglomerate=False, seed_size=100,  param=[], verbose=False,space_seeds=True):
    """This is an EM algorithm to solve a GMM using the SquareEM acceleration procedure."""


    start_time=time.clock()
    Agglomerate=True
    init_start=time.clock()
    dataLength=data.shape[0]
    if dataLength*components>=4000000:
        large_array=False
    else: large_array=False
##    print large_array
##    print dataLength
    transposedData=data.transpose()
    eig_vec, eig_val=DecomposeMatrix(cov)
    scale_1, scale_2=eig_val[0][0], eig_val [1][1]
    min_scale=min(scale_1, scale_2)
    seed_size=(scale_1**.5+scale_2**.5)/4.
##    print seed_size
    #Initialize theta:
    alpha=1./components #Choose uniform weights to begin
    lk=numpy.inf
    count =0
    seed_mod=1.
    if verbose==True: print "\t\t\tInitializing EM parameters..."
    while (numpy.isinf(lk)==True or numpy.isnan(lk)==True):
        #Initialize the component locations of the model.

        count+=1

        if space_seeds==True:
            init_theta=InitializeByFogleButSpaceSeeds(components, data, alpha, cov, seed_size*seed_mod)
        else:
            init_theta=InitializeByFogle(int (components/seed_mod), data, alpha, cov)

        lk_list=[]
        theta_old=init_theta
        cov=init_theta[0].cov
        lk= CalculateLk(data, theta_old)
        if verbose==True: print "\t\t\t\tLikelihood: {0}; Seed modifier: {1}".format(lk, seed_mod)
        seed_mod*=.9


    lk_list.append(lk)
    count=0
    init_time=time.clock()-init_start
    update_time=0.
    recalc_time=0.
    alpha_search=0.
    lk_time=0.
    error_list=[]
    weight_array_time=0
    if iterations==0: return theta_old, [], lk_list
    for i in range(iterations):

        update_start=time.clock()

        #Update parameters from the current theta
        theta_1, wt_time= UpdateEMParameters(data, theta_old, cov, large_array, transposedData)
        weight_array_time+=wt_time

        #Update parameters from the new theta
        theta_2, wt_time= UpdateEMParameters(data, theta_1, cov, large_array, transposedData)
        weight_array_time+=wt_time
        update_time+=time.clock()-update_start
        recalc_start=time.clock()

        #Compute r and v which contain information about how much the parameters
        #changed in the past two updates.

        r=[theta_1[x]-theta_old[x] for x in range(len(theta_old))]
        r_norm=ModelNorm(r)

        v=[theta_2[x]-theta_1[x]-r[x] for x in range(len(r))]
        v_norm=ModelNorm(v)

        #If r_norm = 0, alpha will equal 0. If v_norm = 0, alpha will be undefined
        #If alpha = 0, no EM update is performed. If alpha = -1, a single standard EM update
        #is performed. We always want at least one EM update, so in the worst case
        #we set alpha = -1
        if r_norm!=0 and v_norm!=0: alpha=-1* r_norm/v_norm
        else: alpha=-1
        if numpy.isnan(alpha)==True:
            alpha=-1.

        #Is the step size good???
        if alpha>-1. : alpha=-1.
        else:
            #Check to make sure the model likelihood improves with the proposed alpha
            #If not, iteratively decrease the size of alpha such that it approaches -1.
            #If after a set number of iterations, this fails to increase the model likelihood
            #set alpha = -1
##            print alpha
            alpha_search_start=time.clock()
            alpha=FindGoodAlpha(alpha, v, r, data, theta_old, lk_list[i])
            alpha_search+=time.clock()-alpha_search_start
##            print alpha


        #Update theta using alpha
        resized_r=[k*[alpha, alpha, alpha]*[2,2,2] for k in r]
        resized_v=[k*[alpha**2, alpha**2, alpha**2] for k in v]
##        print alpha
        theta_prime=[theta_old[j] -resized_r[j]+resized_v[j] for j in range(len(theta_old))]
        recalc_time+=time.clock()-recalc_start
##        theta_prime=theta_old-(r*[alpha, alpha, alpha])*[2,2,2]+v*[alpha**2, alpha**2, alpha**2]

        #Update the model parameters with the alpha-informed theta
        update_start=time.clock()
        theta_old, wt_time, weight_array= UpdateEMParameters(data, theta_prime, cov, large_array, transposedData, True)
        weight_array_time+=wt_time

        update_time+=time.clock()-update_start

        #If told to do so, merge any components that are sufficiently close to each other.
        if Agglomerate==True: theta_old=AgglomerateModels(theta_old)
        lk_start=time.clock()

        #Compute the model likelihood.
        lk_list.append(CalculateLk(data, theta_old))

        if dataLength>100000:
            print i, abs( lk_list[-1]-lk_list[-2])

        #Terminate if the change in log-likelihood is sufficiently small.
        if  abs( lk_list[-1]-lk_list[-2])<=1e-3:
                count+=1
        else: count=0
        lk_time+=time.clock()-lk_start

        if count>=1: break
    #Assign each data point to the most likely component distribution

    if large_array==False:
        labels=numpy.argmax( weight_array,0)
    else:
        labels=ChooseLabels( data, theta_old)

    return theta_old,labels, lk_list

def InitializeByFogle(components,data, alpha, cov ):
    """Construct an initial model by centering components on randomly chosen data points."""

    init_theta=[]
    dataLength=data.shape[0]

    #Choose data points as the initial means
    if components<dataLength:
        index=numpy.random.choice(range(0, dataLength), components, replace=False)
    else:
        index=range(0, dataLength)
    for i in range( len (index)):


        init_theta.append(ModelParam(data[index[i]], alpha, cov))

    return init_theta

def InitializeByFogleButSpaceSeeds(components,data, alpha, cov, buffer_size=10 ):
    """Construct an initial model by centering components on randomly chosen
    data points, but with the additional requirement that no two components be within
    a set distance of each other. Place components until no eligible data points remain."""
    init_theta=[]
    dataLength=data.shape[0]
    numpy.random.seed(int(time.clock()*(10000)))

    #Choose data points as the initial means
##    for i in range(components):
    while range(0, dataLength)!=[]:
        index=numpy.random.choice(range(0, dataLength))
        init_theta.append(ModelParam(data[index], alpha, cov))

        pulledRead= data[index]
        distance_from_seed= numpy.sum( abs( data-pulledRead),1)
        RemainingReads=numpy.where(distance_from_seed>buffer_size)

        data=data[RemainingReads]

        dataLength=data.shape[0]
##        print dataLength
##        if range(0, dataLength)==[]: break
    return init_theta

def CalculateLk(data, theta_current, large_array=False):
    """"Compute log-likelihood of the GMM."""
    r_array=CalcBivariateNormal(data, theta_current)
    p_array=numpy.log(numpy.sum(r_array,0))
##    else:
##
##        for index in range( len(theta_current)):
##            r_k=stats.multivariate_normal.pdf(data, theta_current[index].mean,theta_current[index].cov, allow_singular=False)*theta_current[index].weight
##            r_list.append(r_k)
##        r_array=numpy.array(r_list)
##        p_array=numpy.log(numpy.sum(r_array,0))

##    print p_array
    try:
        lk=sum(p_array)
    except:
        lk=p_array
    return lk

def CalcBivariateNormal(data, theta_old):

    """Compute the probability of each of N data points given each of the K components.
    Returns a NxK matrix of probabilities."""

    #Unpack the model parameters
    means=numpy.array( [x.mean for x in theta_old])
    weight=numpy.array( [x.weight for x in theta_old])
    cov=theta_old[0].cov

    #Structure of the covariance matrix allows for some simplifications
    rho=float( cov[0][1])/float( cov[0][0])
    sigma_x=abs( float( cov[0][0]))

    alpha=(2*numpy.math.pi*sigma_x  *(1-(rho)**2)**.5)**-1  #Constant depending only on covariance
    beta=-1*(2*(1-(rho)**2)*sigma_x)**-1                    #Constant depending only on covariance

    #Compute distances from the mean
    x_dist = numpy.subtract.outer( means[:,0], data[:,0]) #(x-E[X]) for all x in X and all k in K
    y_dist= numpy.subtract.outer(means[:,1], data[:,1])  #(y-E[Y]) for all y in Y and all k in K

    # F(x,y)= alpha* exp( beta/sigma^2 *[(x-mean)^2 + (y-mean)^2 - 2 * correl * (x-mean)(y-mean)])

    #Compute the probability oft the bivariate normal
    PDF=alpha*numpy.exp(beta *( (x_dist**2) + (y_dist**2) - 2*rho* x_dist*y_dist))

    #Multiple the probabilities by the mixing proportions
    PDF*=weight[:,None]
##    PDF=PDF.transpose()

    return PDF



def UpdateEMParameters(data, theta_old, cov, large_array, transposedData, return_weight=False):
    #Compute the cluster assignments given the current parameters
    start=time.clock()
    weight_array=ComputeClusterAssignments(data, theta_old, large_array=large_array)
    wt_time=time.clock()-start

    #Update the mixing proportions given the cluster assignments
    new_weights=UpdateMixingProportions(weight_array,data, theta_old, large_array=large_array)

    #Update the means given the cluster assignments
    new_means=UpdateMeans(weight_array, transposedData, theta_old, large_array=large_array)


    #Store the updated parameters as instances of the ModelParam class:
    theta_new=[0]*len(theta_old)
    for j in range(len(new_weights)):
        theta_new[j]=ModelParam(new_means[j], new_weights[j], cov)
##        theta_new[j]=ModelParam(new_means[j], new_weights[j], new_cov[j])
    if return_weight==False:
        return theta_new, wt_time
    else:
        return theta_new, wt_time, weight_array

def ComputeClusterAssignments(data, theta_old, large_array=False):

    """Theta_old is a list of model params

    w_ik = P(X_i|Z_k,Theta_k)*alpha_k / (sum(m=1 to k) P(X_i|Z_m, Theta_m)*alpha_m)
    """

    PDF=CalcBivariateNormal(data, theta_old)
    PDF/=sum(PDF)
##        PDF.transpose()
    return PDF



def CalculateProbabilities(data, parameters, model_index):
    r_k=stats.multivariate_normal.pdf(data, parameters[model_index].mean,parameters[model_index].cov, allow_singular=True)*parameters[model_index].weight
    return r_k


def UpdateMixingProportions(weight_array, data, theta, large_array=False):
    """Update the mixing proportions of the model given cluster assignments"""

    if large_array==False:  #If the data set is small enough, calculate time efficiently
        #Number of data points
        N=weight_array.shape[1]

        #Expected cluster sizes divided by the total number of data points
        weights=weight_array.sum(1)/N

    else: #Otherwise, calculate memory efficiently.
        weights=[]
        for i in range( len(theta)):
            weights.append( sum(CalculateProbabilities(data,theta, i)/ weight_array))

        weights=numpy.array(weights)/sum(weights)
    return weights

def UpdateMeans(weight_array, data, theta, large_array=False):
    """Update the means of the model given cluster assignments"""
    data_x, data_y=numpy.vsplit(data,2)

    if large_array==False:
        N=weight_array.shape[1]
        N_k=weight_array.sum(1)
        mean_x=(weight_array*data_x).sum(1)/N_k
        mean_y=(weight_array*data_y).sum(1)/N_k
        return(numpy.vstack(( mean_x, mean_y)).transpose())

    if large_array==True:
        mean_x=[]
        mean_y=[]
        for i in range(len(theta)):
            r_k=CalculateProbabilities(data.transpose(),theta, i)
            r_k/=weight_array
            mean_x.append(sum( sum( (r_k.transpose()*data_x)/sum(r_k))))
            mean_y.append(sum( sum( ( r_k.transpose()*data_y)/sum( r_k))))
        return(numpy.vstack(( mean_x, mean_y)).transpose())


def ModelNorm(model):
    weight_norm=sum( [d.weight**2 for d in model ])
    mean_norm=sum(sum([d.mean**2 for d in model ]))
    return (weight_norm+ mean_norm)**.5

def FindGoodAlpha(alpha, v, r, data, theta, Lk_zero):
    good_alpha=False
    count=0
    while good_alpha==False and count<10 and alpha<=-1:
        resized_r=[k*[alpha, alpha, alpha]*[2,2,2] for k in r]
        resized_v=[k*[alpha**2, alpha**2, alpha**2] for k in v]
        theta_alpha=[theta[i] -resized_r[i]+resized_v[i] for i in range(len(theta))]
        Lk_alpha=CalculateLk(data, theta_alpha)

        count+=1
        if Lk_alpha>=Lk_zero: good_alpha=True
        else:
            alpha=(alpha-1.)/2.
    if good_alpha==True:
        return alpha
    else:
        return -1.

def ChooseLabels(data, theta_old, large_array=False):

    """Theta_old is a list of model params

    w_ik = P(X_i|Z_k,Theta_k)*alpha_k / (sum(m=1 to k) P(X_i|Z_m, Theta_m)*alpha_m)

    """
    data_length= max( data.shape)
    last_row=numpy.array([0.]*data_length)
    labels=numpy.array([0]*data_length)
    for index in range(0, len(theta_old)):
        next_row= CalculateProbabilities(data, theta_old, index)
        better_indices=numpy.where(next_row>=last_row)[0]
        labels[better_indices]=index
        last_row=next_row
    return labels



def AgglomerateModels(theta):

    """Search through a list of ModelParam objects and merge any components with similar means.
    Merging is accomplished  by averaging the means of the two components, and summing
    the mixing proportions."""
    meanDict={}
    totalWeigh=0.
    x_list=[]
    y_list=[]

    #Unpack the parameter means into x and y components
    for t in theta:
        x_list.append(t.mean[0])
        y_list.append(t.mean[1])

    #Compute the pairwise Manhattan distance between components
    y_diff=abs( numpy.subtract.outer(y_list, y_list))
    x_diff=abs( numpy.subtract.outer(x_list,x_list))
    pw_dist=y_diff+x_diff

    #Identify which components are within one unit apart
    close_parameters=numpy.where(pw_dist<=1)
    not_diag=close_parameters[0]!=close_parameters[1]
    pairs_to_join=zip(close_parameters[0][not_diag], close_parameters[1][not_diag])

    #Remove doubles
    check_dict={}
    join_list=[]
    for pair in pairs_to_join:
        sorted_pair=sorted(pair)
        if check_dict.has_key( ','.join(numpy.array(sorted_pair,str)))==False:
            join_list.append(sorted_pair)
            check_dict[ ','.join(numpy.array(sorted_pair,str))]=True

    #Identify which components to keep as they are
    keep_unmodified= set(close_parameters[0][~not_diag])-set(close_parameters[0][not_diag])


    #Keep these components
    newTheta=[]
    for i in keep_unmodified:
        newTheta.append(theta[i])

    #Merge the remaining pairs;
    possibly_merge=set(range(len(theta)))-keep_unmodified
    merge_tracker={}
    for p in possibly_merge:
        merge_tracker[p]=False


    already_joined={}
    for j in join_list:
        theta1,theta2=j
        if merge_tracker[theta1]==False and merge_tracker[theta2]==False:
            new_mean_x=numpy.mean([ theta[j[0]].mean[0], theta[j[1]].mean[0]])
            new_mean_y=numpy.mean([ theta[j[0]].mean[1], theta[j[1]].mean[1]])
            new_weight=theta[j[0]].weight+theta[j[1]].weight
            newTheta.append(ModelParam([new_mean_x, new_mean_y], new_weight, theta[j[0]].cov))
            merge_tracker[theta1]=True
            merge_tracker[theta2]=True

    for p in merge_tracker.keys():
        if merge_tracker[p]==False:
            newTheta.append(theta[p])

    if len(newTheta)>len(theta):
        print "Theta grew from {0} to {1}".format(len( theta), len(newTheta))
        print pairs_to_join
        print keep_unmodified
        print theta
        print newTheta
        print jabber
    return newTheta


#----------------------------------------------------------------------------#
#Feeding ConTExt output throught the EM algorithm


def ClusterDirectory(indir, outfile, dist_file, cov, ref_dict, GemCode=False, ref_names=['2L', '3L', '2R', '3R', 'X', '4', 'U']):

    """Iterate through each file in the directory indir, fit GMMs to the read pair
    distributions and write the output to outfile."""

    cluster_count=0
    print ref_dict
    Pair_List={}
    start_time=time.clock()
    read_timer, cluster_timer, write_timer=0.,0.,0.

    m,k,p=ReadDistributions(dist_file)
    p/=sum(p)

    #Simulate the estimation error as a function of read count
    reg_file='{0}_varreg.tsv'.format(indir) #Needs to be redirected to the output directory
    regression_coefficients=SimulateErrorByReadCount(p, 200, 1000, reg_file)

    files=sorted( os.listdir(indir))
    outhandle=open(outfile, 'wb')
    outTable=csv.writer(outhandle, delimiter='\t')
    header=['Seq1', 'Seq2', 'Strand1', 'Strand2','Feature', 'X_est', 'Y_est', 'Err1', 'Err2', 'Read Count', "X 3'-positions", "Y 3'-positions","X 5'-positions", "Y 5'-positions", "X MapQ", "Y MapQ", "Cluster ID", "X CIGAR", "Y CIGAR", "X MD", "Y MD",]
    outTable.writerow(header)

    for f in files:
        file_root=f.split('.')[0]
##        if f[0]=='(':continue
##        if file_root!= 'Jockey-3_DSim': continue
        if ref_names.count(file_root)!=0: continue
        infile='{0}/{1}'.format(indir, f)
        Pair_List, read_time, cluster_time, write_time, cluster_count=ClusterFile(infile, outTable, Pair_List, cov, k, p, ref_dict, cluster_count, regression_coefficients, GemCode)

        print "Time spent reading: {0}\nTime spent clustering: {1}\nTime spent writing: {2}".format (read_time, cluster_time, write_time)
        read_timer+=read_time
        cluster_timer+=cluster_time
        write_timer+=write_time
    outhandle.close()
    print "It took {0} seconds to process.\n Time spent reading: {1}\nTime spent clustering: {2}\nTime spent writing: {3}".format (time.clock()-start_time, read_timer, cluster_timer, write_timer)

def ClusterFile(infile,outTable, Pair_List, cov, k, p, ref_dict, cluster_count, regression_coefficients, GemCode=False):


    print ref_dict
    read_start=time.clock()
    image, seq_name=LoadImages(infile,ref_dict)
    scatter, image=PrepareImages(image, seq_name)
    read_time=time.clock()-read_start
    cluster_start=time.clock()
    clusters=ClusterAllImages(scatter, image, seq_name,cov, k )
    cluster_time=time.clock()-cluster_start
    write_start=time.clock()
    for key in clusters.keys():
        pair=tuple( sorted((seq_name, key)))
        if Pair_List.has_key(pair)==True:
            continue
        else:
            Pair_List[pair]=True
        for quadrant in clusters[key]:
            print '.',
            cluster_count= WriteClusters(outTable, clusters[key][quadrant], seq_name, key, quadrant, p, cluster_count,  regression_coefficients, GemCode)
    write_time=time.clock()-write_start
    return Pair_List, read_time, cluster_time, write_time, cluster_count

def LoadImages(infile, refLen):
    """Reads a ConTExt file stores each row as ConTExtLine object, and organizes these into separate scatterplots.
    Reads where one end is anchored to a reference chromosome (as given by keys in teh refLen)
    are discarded if the MapQ to the reference is less than 20. As MapQ is meaningless for repeats
    in our alignment pipeline, this is only considered for reference alignments.

    Returns a dictionary keyed by sequence name and a tuple of strands, as well as the name of repeat
    represented by the file."""

    print refLen
    inHandle=open(infile,  'r')
    inTable=csv.reader(inHandle, delimiter='\t')
    row=inTable.next()
    imageDict={}
    #load the lines
    seq_name=row[0]
    for row in inTable:
        line=ConTExtLine(row)

        context=line.Seq2
        quadrant=(line.Strand1, line.Strand2)
        if refLen.has_key(context) and line.MapQ2<20: continue
        if line.Seq1==line.Seq2 and line.Strand1!=line.Strand2:
            quadrant=('-','+')
        if imageDict.has_key(context)==False:
            imageDict[context]={}
        if imageDict[context].has_key(quadrant)==False:
            imageDict[context][quadrant]=[]

        imageDict[context][quadrant].append(line)

    return imageDict, seq_name


def PrepareImages(images, anchor_sequence):

    """Takes a list of ConTExt lines and extracts read positions from it.
    Some scatterplots place the reads in a pair axes other than those indicated by the
    input file. In these case, the ConTExtLines are reorganized to reflect this."""
    scatterDict={}
    new_image={}

    for context in images.keys():
        scatterDict[context]={}
        new_image[context]={}
        self_self_image= anchor_sequence==context
        opposite_strands_processed=False

        for quadrant in images[context].keys():
            if self_self_image==False:
                scatterDict[context][quadrant], new_image[context][quadrant] =GetPointsFromLines(images[context][quadrant], quadrant,self_self_image )

            else:
                if (quadrant==('-','+') or quadrant==('+','-')) and opposite_strands_processed==False:
                    joinedList=[]

                    #Merges +- and -+ images.
                    if images[context].has_key(('-','+'))==True:
                        joinedList+=images[context][('-','+')]

                    if images[context].has_key(('+','-'))==True:
                        joinedList+=images[context][('+','-')]

                    scatterDict[context][quadrant], new_image[context][quadrant]=GetPointsFromLines(joinedList, quadrant,self_self_image )
                    opposite_strands_processed=True

                if quadrant[0]==quadrant[1]:
                    scatterDict[context][quadrant], new_image[context][quadrant]=GetPointsFromLines(images[context][quadrant], quadrant,self_self_image )

    return scatterDict, new_image

def GetPointsFromLines(lineList, quadrant, self_self_image=False):
    """Extracts read positions from the list of ConTExt lines. If reads are assigned
    to an axis different from that assumed by the ConTExt output, it also modifies
    the organization of the lines to reflect this reassignment.

    Returns
    numpy.array containing the x and y positions
    list containing updated lines.

    """

    xList, yList=[],[]
    new_linelist=[]

    for line in lineList:
        if self_self_image==False:
            xList.append(line.threePrime1)
            yList.append(line.threePrime2)
            new_linelist.append(line)
        else:
            if line.Strand1==line.Strand2:

#In self-self images the coordinates of reads mapping to the same strand are
#ambiguous. To resolve this, I assign them consistently to one of the two positions
#they could be assigned.

                if line.Strand1=='+':
                    if line.threePrime1< line.threePrime2:
                        new_linelist.append(ConTExtLine( line.inverse_row())) #The second read is on the X-axis, so we need to rearrange the ConTExt line
                        xList.append(line.threePrime2)
                        yList.append(line.threePrime1)
                    else:
                        new_linelist.append(line)
                        xList.append(line.threePrime1)
                        yList.append(line.threePrime2)

                elif line.Strand2=='-':
##                    xList.append(min( line.threePrime1, line.threePrime2))
##                    yList.append(max( line.threePrime1, line.threePrime2))
                    if line.threePrime1> line.threePrime2:
                        new_linelist.append(ConTExtLine( line.inverse_row())) #The second read is on the X-axis, so we need to rearrange the ConTExt line
                        xList.append(line.threePrime2)
                        yList.append(line.threePrime1)
                    else:
                        new_linelist.append(line)
                        xList.append(line.threePrime1)
                        yList.append(line.threePrime2)
            else:

#In self-self images, for reads mapping to the different strands, I always want
#reverse read on the x-axis.

                if line.Strand1=='+':
                    xList.append( line.threePrime2)
                    yList.append( line.threePrime1)
                    new_linelist.append(ConTExtLine( line.inverse_row())) #The second read is on the X-axis, so we need to rearrange the ConTExt line
                elif line.Strand1=='-':
                    xList.append( line.threePrime1)
                    yList.append(line.threePrime2)
                    new_linelist.append(line)

    return numpy.vstack((xList, yList)), new_linelist


def ClusterAllImages(scatter, lines, anchor, cov, dist):
    """Takes a dictionary holding the scatterplots for a file and a list of ConTExt lines
    and clusters each using SquareEM"""

    start_time=time.clock()
    skip_splitting=True
    has_consensus=False
    cluster_dict={}

    print "Cluster reads in {0}".format(anchor)
    for contig in scatter.keys():   #For each contig
        if contig=='*': continue
        cluster_dict[contig]={}
        for quadrant in scatter[contig].keys():     #for each pair of strands

            if anchor==contig and quadrant==('-','+'): #Are the two contigs in the image the same?
                self_self=True
            else:
                self_self=False

            if quadrant[0]==quadrant[1]:
                orientation=[[1,1],[1,1]]
            else:
                orientation=[[1,-1],[-1,1]]

            print "\t\tProcessing sequences {1}, {0}, strands {2}".format(anchor,  contig, quadrant)

            if self_self==True and quadrant==('-','+'): #This image contains concordant read pairs... we need to remove them
                concordant, concordant_lines, concordantIndices,  nonconcordantReads, nonconcordantIndices=SplitConcordantAndNonconcordantByThreshold(scatter[contig][quadrant], lines[contig][quadrant], dist, .999)
                has_consensus=True
            else:
                nonconcordantReads=scatter[contig][quadrant]
                nonconcordantIndices=numpy.array( range(0,len(nonconcordantReads[0,:])))

            if len( nonconcordantIndices)==0: return cluster_dict

            #This condition is impossible, routine never accessed. Let's comment it out and see what happens

##            if nonconcordantReads[1,:].max()-nonconcordantReads[1,:].min()<-20000: #If the image is very large,
##            #seeds will not be able to cross the space between clusters during EM
##            #because of limits on the precision of the bivariate PDF.
##                partitioning_indices=PartitionReads(nonconcordantReads, 500)
##                partitioned_read=[]
##                partitioned_indices=[]
##                cluster_list=[]
##
##                for i in partitioning_indices:
##                    partitioned_read.append(nonconcordantReads[:,i])
##                    partitioned_indices.append(nonconcordantIndices[i])
##                    reads_x, reads_y, labels, nonconcordantIndices_2, clusters_part=ClusterImage(scatter[contig][quadrant], quadrant, numpy.array( cov)*orientation, lines[contig][quadrant], self_self, dist,partitioned_read[-1], partitioned_indices[-1])
##                    cluster_list.append(clusters_part)
##
##                full_clusters={}
##                count=0
##
##                for dictionary in cluster_list:
##                    for key in dictionary.keys():
##                        full_clusters[count]=dictionary[key]
##                        count+=1
##
##                cluster_dict[contig][quadrant]=full_clusters

##            else:
            reads_x, reads_y, labels, nonconcordantIndices_2, clusters=ClusterImage(scatter[contig][quadrant], quadrant, numpy.array( cov)*orientation, lines[contig][quadrant], self_self, dist,nonconcordantReads, nonconcordantIndices)
            cluster_dict[contig][quadrant]=clusters

    if has_consensus==True and len(concordantIndices)!=0:
        consensus_cluster=Cluster(concordant, concordant_lines)
        cluster_dict['Consensus']={}
        cluster_dict['Consensus'][('-','+')]={0:consensus_cluster}

    print "\tTook {0} seconds to process {1}.".format( time.clock()-start_time, anchor)
    return cluster_dict


def SplitConcordantAndNonconcordantByThreshold(image, lines, dist, cutoff=.005):
    """If a scatterplot is self-self and represents two different strands,
    it contains concordant read pairs. We only want to cluster discordant read pairs.
    This functions pulls out concordant read pairs as using percentiles of the MPD
    distributions and gives them their own clusters.

    In future, it might be better to include a component in the model in that corresponds to
    the diagonal formed by concordant read pairs. It wouldnt' require any mean updates,
    just the mixing proportion updates, so it might not influence the update statistics.
    """

    insertSize= image[0,:]-image[1,:]
    min_diff, max_diff=dist

    concordancyTest=(min_diff<insertSize)*(insertSize<=max_diff)
    concordantIndices=numpy.where(concordancyTest==True)[0]
    nonconcordantIndices=numpy.where(concordancyTest==False)[0]

    concordant=image[:, concordantIndices]
    nonconcordant=image[:,nonconcordantIndices]
    concordant_lines=[lines[i] for i in concordantIndices]

    print concordant.shape, nonconcordant.shape

    return concordant, concordant_lines, concordantIndices, nonconcordant, nonconcordantIndices

def PartitionReads(reads, threshold):
    """Not used. Possibly delete."""
    indices=range(len(reads[1,:]))
    dummy_array=numpy.vstack((reads[1,:],indices)).transpose()
##    print dummy_array
    sorted_array=numpy.array( sorted(dummy_array, key=lambda x: x[0])).transpose()
    positions,  indices=numpy.vsplit(sorted_array,2)
    diff_array=numpy.diff(positions)[0]
##    print diff_array
    gap_array=list( numpy.where(diff_array>threshold)[0]+1)
##    print gap_array
    gap_array=[0]+gap_array+[len(positions[0])+1]
##    print gap_array
    partitioned_indices=[]
    for i in range(len(gap_array)-1):
        partitioned_indices.append(indices[0][gap_array[ i]:gap_array[ i+1]])
##    print partitioned_indices

    partitioned_read=[]
    for i in partitioned_indices:
        partitioned_read.append(reads[:,i])
##    print partitioned_read

    return partitioned_indices


def ClusterImage(image, quadrant, cov, lines, self_self_image, dist,nonconcordantReads, nonconcordantIndices, block_size=15):
    eig_vec, eig_val=DecomposeMatrix(cov)
    block_size=(min([eig_val[0][0], eig_val[1][1]])**.5)*2
    startTime=time.clock()
    agg_Bool=True
    xlen=numpy.max( nonconcordantReads[0,:])-numpy.min( nonconcordantReads[0,:])+1
    ylen=numpy.max( nonconcordantReads[1,:])-numpy.min( nonconcordantReads[1,:])+1
    components=max(10, int( (xlen*ylen)/(block_size**2))+1)

    if self_self_image==True and quadrant==('-','+'):
        agg_Bool=True

    clusterDict={}
    if len(nonconcordantReads[0,:])==1:
        clusterDict[0]= Cluster( image[:,nonconcordantIndices], lines)
        reads_x, reads_y, labels=1,2,3
        return reads_x, reads_y, labels, nonconcordantIndices, clusterDict

    good_clustering=numpy.nan

    #Cluster the scatterplot with SquareEM, check whether the the algorithm has converged
    #and if it has not, repeat until it converges.
    while numpy.isnan( good_clustering)==True or numpy.isinf( good_clustering)==True:
        final_theta, final_labels, final_Lk= SquareEM(nonconcordantReads.transpose(), cov, components, 200, Agglomerate=agg_Bool, seed_size=block_size, verbose=True)
        print "Likelihood is {0}  after {1} iterations".format(final_Lk[-1],len(final_Lk))
        if numpy.isnan( final_Lk[-1])==True or numpy.isinf(final_Lk[-1])==True:
            print '\t\t\tFailed to converge. Retrying.'
        good_clustering=final_Lk[-1]


    print "\t\t\tConverged."
    reads_x, reads_y, labels=nonconcordantReads[0,:], nonconcordantReads[1,:], final_labels
    #Create output clusters


    #Identify the set of unique labels
    labelSet=sorted( list(set(labels)))
##        print labelSet
    label_array=numpy.array(labels)

    #Construct clusters of read positions and ConTExt lines based on the EM labels.
    for unique_label in labelSet:
        label_indices=numpy.where(label_array==unique_label)[0]
##        original_indices=old_indices[label_indices]
        original_indices=label_indices
        read_indices=nonconcordantIndices[original_indices]
        line_list=[]
        for index in list( read_indices):
            line_list.append(lines[index])
##        try:
        clusterDict[unique_label]=Cluster( image[:,read_indices], line_list)

    print "after {0} seconds".format( time.clock()- startTime)

    return reads_x, reads_y, labels, nonconcordantIndices, clusterDict

#Functions used to write output

def WriteClusters(outTable, clusters,anchor_name, target_name, quadrant, Pr_dist, cluster_count,  regression_coefficients, GemCode, phred=33):
    """Write the clustering output to the output file."""

    #outTable should be a csv.writer object

    consensus=False
    if target_name=='Consensus':
        consensus=True
        target_name=anchor_name
    for key in clusters.keys():

        read_count=len(clusters[key].read_positions[0,:])
        if consensus==False:
##            closest_jxn, offset, CI=EstJxnDist(clusters[key], quadrant, Pr_dist)
##            best_jxn, min_jxn, max_jxn=EstJxnDist(clusters[key], quadrant, Pr_dist)
            if read_count<=1000:
                x5_positions=[str(i.fivePrime1) for i in clusters[key].line_list[:] ]
                y5_positions=[str( i.fivePrime2 ) for i in clusters[key].line_list[:] ]
                ends5prime=zip(x5_positions,y5_positions)
                unique_read_count=len(set(ends5prime))
            else: unique_read_count=read_count
            x_positions=numpy.asarray( clusters[key].read_positions[0,:],'float')
            y_positions=numpy.asarray( clusters[key].read_positions[1,:],'float')
            data=[x_positions, y_positions]
            best_jxn=EstimateJunctionLocation(data, Pr_dist, quadrant)
            err1, err2=PowerFunction(unique_read_count, regression_coefficients)
        else:
            best_jxn, err1, err2=('nan', 'nan'), 'nan', 'nan'

        image_information=[anchor_name,target_name, quadrant[0], quadrant[1] ]
        if consensus==False:

            junction_information=['Junction', best_jxn[0], best_jxn[1], err1[0]**.5, err2[0]**.5]
        else:
            junction_information=['Consensus', best_jxn[0], best_jxn[1], err1, err2]

        if consensus==False and read_count>=500:
            right_index=read_count
        else:right_index=read_count

        x_positions=numpy.asarray( clusters[key].read_positions[0,:],'int')
        y_positions=numpy.asarray( clusters[key].read_positions[1,:],'int')
        x_string=';'.join(numpy.char.mod('%i',x_positions))
        y_string=';'.join(numpy.char.mod('%i',y_positions))
        clusters[key].line_list=RepairAllMDStrings(clusters[key].line_list)
        x5_positions=[str(i.fivePrime1) for i in clusters[key].line_list[:right_index] ]
        y5_positions=[str( i.fivePrime2 ) for i in clusters[key].line_list[:right_index] ]
        x5_string=';'.join(x5_positions)
        y5_string=';'.join(y5_positions)

        x_mapq=[str(i.MapQ1) for i in clusters[key].line_list[:] ]
        y_mapq=[str(i.MapQ2 ) for i in clusters[key].line_list[:] ]
        x_mapq_string=';'.join(x_mapq)
        y_mapq_string=';'.join(y_mapq)

        x_CIGAR=[str(i.Cigar1) for i in clusters[key].line_list[:] ]
        y_CIGAR=[str(i.Cigar2) for i in clusters[key].line_list[:] ]
        x_CIGAR_string=';'.join(x_CIGAR)
        y_CIGAR_string=';'.join(y_CIGAR)

        x_MD=[str(i.MD1) for i in clusters[key].line_list[:] ]
        y_MD=[str(i.MD2 ) for i in clusters[key].line_list[:] ]
        x_MD_string=';'.join(x_MD)
        y_MD_string=';'.join(y_MD)

        x_QD, y_QD=GetQualityMasks(clusters[key].line_list[:], phred=phred)
        x_QD_string=';'.join(x_QD)
        y_QD_string=';'.join(y_QD)


        cluster_information=[ read_count, x_string, y_string, x5_string, y5_string, x_mapq_string,  y_mapq_string, str(cluster_count), x_CIGAR_string,y_CIGAR_string, x_MD_string, y_MD_string, x_QD_string, y_QD_string]
##        cluster_information=[]

        if GemCode==True:
            barcodes=[str( i.optional ) for i in clusters[key].line_list ]
            barcode_string=[';'.join(barcodes)]
            row=image_information+junction_information+cluster_information+barcode_string
        else:
            row=image_information+junction_information+cluster_information
        cluster_count+=1
        outTable.writerow(row)
    return cluster_count



def RepairAllMDStrings(context_lines):

    for i in range(len( context_lines)):

        line=context_lines[i]
        #repair first
        line.MD1=RepairMDString(line.MD1, line.Cigar1, line.Read1)
        #repair second:
        line.MD2=RepairMDString(line.MD2, line.Cigar2, line.Read2)
        context_lines[i]=line

    return context_lines


def RepairMDString(md_string, cigar, seq):
    """Use the md_string to identify positions in the read that do not match the
    consensus. The MD string display the consensus nucleotide. We are not retaining
    the actual sequencing read in the cluster output, so use this function to replace the
    consensus nucleotide with the read nucleotide."""

    header=md_string[:5]
    md=md_string[5:]
    hold='0'
    null=False
    expanded_md=[]
    CIGAR=False
    count=0
    new_md=[]
    nt={'A':1, 'T':2, 'C':3,'G':4, 'N':5}

    for char in md:
        new_md.append( char)
##        print null
        if null==True:
            if nt.has_key(char)==False:
                null=False
            else: continue
##            continue
        if char=='^':
            null=True
            count+=int(hold)
            hold='0'

            continue

        if nt.has_key(char)==False:
            hold+=char
        if nt.has_key(char)==True:
            count+=int(hold)
            if CIGAR==False:
                CIGAR_array=ExpandCIGAR(cigar)

                read_cigar=CIGAR_array[(CIGAR_array=='M')+(CIGAR_array=='I')]
##            print hold

                read_array=numpy.fromstring(seq.upper(),'|S1')
                read_matches=read_array[read_cigar=='M']
                new_md[-1]=read_matches[count]


            hold='0'
    return header+ ''.join( new_md)


def GetQualityMasks(context_lines, phred=33, cutoff=30):
    """Constructs a string that can be used to masked read positions with low Phred Quality"""
    mask_x, mask_y=[],[]
    for line in context_lines:

        cigar_x=line.Cigar1
        if cigar_x!='*':
            #Convert the quality string to an array
            qual_string_x=numpy.fromstring(line.Qual1, '|S1')
            #Convert the quality string to Phred scores
            qual_x=ConvertFromQualityString( line.Qual1, phred)

            #Expand the CIGAR and identify positions in the read which map to the reference
            cig_array_x=ExpandCIGAR(cigar_x)
            #Place the indices within the set of positions belonging in the read (not deleted relative to the reference)
            #which are considered matching.
            matched_indices_x=numpy.where(cig_array_x[cig_array_x!='D'] =='M')[0]
            qual_x_match=qual_x[matched_indices_x]

            #Identify aligned positions with poor quality scores
            bad_calls_x=numpy.where( qual_x_match<cutoff)[0]

            #Interleave the lengths of good intervals and the pr
            MD_string_x="QD:Z:{0}".format( ','.join(bad_calls_x.astype(str)))
            mask_x.append(MD_string_x)
        else: mask_x.append('')

        cigar_y=line.Cigar2
        if cigar_y!='*':
            qual_string_y=numpy.fromstring(line.Qual2, '|S1')
            qual_y=ConvertFromQualityString( line.Qual2, phred)
            cig_array_y=ExpandCIGAR(cigar_y)
            matched_indices_y=numpy.where(cig_array_y[cig_array_y!='D'] =='M')[0]
            qual_y_match=qual_y[matched_indices_y]

            bad_calls_y=numpy.where( qual_y_match<cutoff)[0]

            QD_string_y="QD:Z:{0}".format( ','.join(bad_calls_y.astype(str)))
            mask_y.append(QD_string_y)
        else: mask_y.append('')
    return mask_x, mask_y

def ConvertFromQualityString(string, phred=33):
    out=numpy.array([ord(s)-phred for s in string])
    return out

#----------------------------------------------------------------------------#
#Optimization Functions


def OptimizeOnFileWithGP(indir,dist_file, outfile, max_x, max_y, axis, steps, reps,cutoff=.99, distance_cutoff=400 ,filter_x=400, filter_y=400,assembly_list=['2L', '2R', '3L', '3R', 'X'], files=['BAGGINS1.dat', 'ROO_LTR.dat', 'HOBO.dat', 'DMCR1A.dat']):


    print axis

    kernel=gaussian_process.kernels.RationalQuadratic(length_scale_bounds=( 1e-20, 1e20))* (gaussian_process.kernels.ConstantKernel()+ gaussian_process.kernels.WhiteKernel())#+gaussian_process.kernels.ExpSineSquared()*gaussian_process.kernels.WhiteKernel()

    #Read the MPD distribution
    m,k,p=ReadDistributions(dist_file)
##    upper_size=numpy.where(numpy.array( Cumulative(p))>=.99)[0][0]
    samples=[]
    samples_opp=[]
    #Obtain a training set of insertions into  unique sequence
    for f in files:
        clusters=RoughClusters(indir+'/'+f +'.dat', distance_cutoff,3, assembly_list)
        samples+=clusters[('-','-')]
        samples+=clusters[('+','+')]
        orientation=numpy.array([1,-1])

        samples+=[i*orientation[:,None] for i in clusters[('-','+')]]
        samples+=[i*orientation[:,None] for i in clusters[('+','-')]]
##    for i in range(150):
##        samples.append(GenSample(p, 200))
    #Exclude poorly behaved clusters
    print len (samples)
    #Identify items in the training set that do not return a single cluster with the
    #largest covariance matrix tested.
    eig_vec=numpy.matrix([[ 0.70710678,  0.70710678], [-0.70710678,  0.70710678]])
    cov_max=numpy.array( RecomposeMatrix(eig_vec, [[filter_x**2,0],[0,filter_y**2]]))
    orientation=[[1,1],[1,1]]
    score, max_scores=ScoreClustering(p, cov_max*orientation, len(samples), 10, 100, samples, bootstrap=False)
    max_scores=numpy.array(max_scores)
    bad_indices= numpy.where( (max_scores==1.)==False)[0]

    #Remove these from the training set
    for i in reversed(bad_indices):
        del samples[i]

    print "First filter leaves {0}".format(len( samples))

    #Identify items in the training set that return a single cluster with the
    #smallest covariance matrix tested.
    if axis=='x':
        cov_min=RecomposeMatrix(eig_vec, [[10,0],[0,10]])
    else:
        cov_min=RecomposeMatrix(eig_vec, [[1,0],[0,1]])
    score, min_scores=ScoreClustering(p, cov_min,  len(samples), 10, 100, samples, bootstrap=False)
    max_scores=numpy.array(min_scores)
    bad_indices=numpy.where(max_scores==1.)
    bad_indices= numpy.where( (max_scores==1.)==True)[0]

    #Remove these from the training set
    for i in reversed(bad_indices):
        del samples[i]

    print "Second filter leaves {0}".format( len (samples))

    #Assess clustering performance as a function of one eigenvalue
    score, scale=FindBestCovV(p, steps,max_x,1, max_y,1,reps, axis=axis, samples=samples)

    #Estimate the true relationship between performance and the eigenvalue using
    #GP regression
    GP=gaussian_process.GaussianProcessRegressor(kernel,n_restarts_optimizer=20)
    GP.fit(scale[:,None],score) #Learn to predict score as a function of scale

    if axis=='x':
        test_scales=numpy.arange(0, max_x, .1)
    else: test_scales=numpy.arange(0, max_y, .1)
    expected_score=GP.predict(test_scales[:,None])


    best_index=numpy.min(  numpy.where(expected_score>=cutoff)[0])
    print best_index
    best_scale=test_scales[best_index]
    print best_scale

    outhandle=open(outfile, 'ab')
    outTable=csv.writer(outhandle, delimiter='\t')
    outTable.writerow(['Average Error'] +list(score) )
    outTable.writerow(['Gaussian Scale']+ list(scale))
    outTable.writerow(['GP Prediction']+list(GP.predict(scale[:,None])))
    print best_scale
    outTable.writerow(['Best Scale for dimension {0}'.format(axis)]+[best_scale])

    plotfile='{0}_{1}.jpeg'.format('.'.join(outfile.split('.')[:-1]), axis)


    pyplot.scatter(scale, score, alpha=.7, c='firebrick', zorder=1)
    pyplot.plot(test_scales,expected_score, lw=2, alpha=1., c='dodgerblue',  zorder=2)
    pyplot.plot([best_scale,best_scale], [min(score), expected_score[best_index]], c='orange', lw=1.5, alpha=1, ls='--',  zorder=3)
    pyplot.plot([0, best_scale], [expected_score[best_index], expected_score[best_index]], c='forestgreen', lw=1.5, alpha=1, ls='--',  zorder=4)
##    pyplot.scatter(best_scale,expected_score[best_index], c='black', lw=2,alpha=1 ,  zorder=5)

    pyplot.xlabel('Scale', size=14)
    pyplot.ylabel('Clustering Performance', size=14)
    pyplot.ylim(min(score), 1.05)
    pyplot.savefig(plotfile)
    pyplot.close()

    return best_scale

def RoughClusters(infile, distance_cutoff=400, read_depth=3, assembly_list=['2L', '2R', '3L', '3R', 'X']):

    """Accomplishes single-linkage agglomerative clustering with the specified
    distance cutoff on reads from the provided file to identify TE insertions
    into the specified reference chromosomes."""

    inHandle=open(infile,  'r')
    inTable=csv.reader(inHandle, delimiter='\t')
    inTable.next()
    readDict={}

    #load the read pairs
    for row in inTable:
        line=ConTExtLine(row)
        context=line.Seq2
        quadrant=(line.Strand1, line.Strand2)
        if line.MapQ2<20: continue  #Don't cluster on poor alignments
        if readDict.has_key(context)==False:
            readDict[context]={}
        if readDict[context].has_key(quadrant)==False:
            readDict[context][quadrant]=[]
        readDict[context][quadrant].append(line)

    #For all lines in assembly sort them into scatterplots
    TrainByQuadrant={}
    for key in readDict:
        if assembly_list.count(key)==0: continue
        for quad in readDict[key].keys():
##            print key, quad
            sortedPositions=    sortedLines=sorted(readDict[key][quad], key=lambda w:w.threePrime2)
            if TrainByQuadrant.has_key(quad)==False:
                TrainByQuadrant[quad]=[]

            #Cluster the scatterplot; this is a list of lists
            TrainByQuadrant[quad]+=ContextLinesByPosition(sortedPositions, distance_cutoff)
    cov={}
    data_points={}
    for key in TrainByQuadrant.keys():
##        print key
        data_points[key]=[]
        count=0.
        for item in TrainByQuadrant[key]:
            if len (item)<read_depth: continue
            point_array, new_linelist=GetPointsFromLines(item, key)
            data_points[key].append(point_array )
        if count==0: continue

    return data_points

def ContextLinesByPosition(lines, distance_cutoff):
    """Rename to ClusterLinesByPositions"""
##    sortedLines=sorted(lines, key=lambda w:w.threePrime2)
    posList=[line.threePrime2 for line in lines ]
    splitIndices=DetermineIndicesToSplitBy(posList, distance_cutoff)
    clusteredLines=SplitListByIndices( lines, splitIndices)
    return clusteredLines

def DetermineIndicesToSplitBy(posList, distance_cutoff):
    """Taking a list of read positions sorted by chromosome position as input,
    this carries out single-linkage agglomerative clustering with a distance
    cutoff."""
    diffArray=numpy.diff(posList)   #Distance between adjacent read pairs on the chromosome arm
    indices=list( numpy.where(diffArray>distance_cutoff)[0]+1) #Indices where the distnace exceeds the distance cutoff
    splitIndices=[0]+indices+[len(posList)] #Indices that will define cluster boundaries
    return splitIndices

def SplitListByIndices(inList, indices):
    """"""
    newLists=[]
    for i in range( len(indices)-1):
        newLists.append(inList[indices[ i]:indices[ i+1]])
    return newLists

def ScoreClustering(dist,cov_test, reps, components, iterations, samples=[] , bootstrap=True):
    """Assesses clustering performance on a set of real or simulated read pair
    distributions. The GMM is initialized using the standard Forgy method with
    10 components. As each training distribution is expected to represent one
    cluster, performance is measures as the proportion of read pairs in the
    largest cluster."""
    start_time=time.clock()
    simpson_list=[]
    for i in range(reps):
        if samples==[]:
            data=GenSample(dist, 5)
        else:
            if bootstrap==True:
                j=numpy.random.randint(0, len(samples))
            else:
                j=i
            data=samples[j]
        t,l,lk=SquareEM(numpy.vstack(data).transpose(), cov_test, 10, iterations, Agglomerate=True, space_seeds=False)
        SI=Simpson(l)
        simpson_list.append(SI)
    print time.clock()-start_time
    freq=HarmonicMean(simpson_list)
    return freq, simpson_list

def FindBestCovV( Pr, step_count, scale_1_max, scale_1_min, scale_2_max, scale_2_min, reps,axis='x', samples=[]):
    """Assesses clustering performance over a range of covariance eigenvalues."""

    upper_size=numpy.where(numpy.array( Cumulative(Pr))>=.999)[0][0]
    print axis
    print upper_size

    X_vals, X_scores=[],[]

    Y_vals, Y_scores=[],[]

    eig_vec=numpy.matrix([[ 0.70710678,  0.70710678], [-0.70710678,  0.70710678]])

    new_scale_1=scale_1_max
    new_scale_2=scale_2_max
    last_2=new_scale_2
    last_1=new_scale_1
    print scale_1_max, scale_2_max

    step_1=int( (scale_1_max-scale_1_min)/step_count)
    step_2=int( (scale_2_max-scale_2_min)/step_count)
    perf_list=[]
    scale_list=[]

    if axis!='x':

        for j in numpy.linspace( (scale_2_max),(scale_2_min) ,step_count):
            test_val=[[last_1**2,0],[0,j**2]]
            cov_test=numpy.array( RecomposeMatrix(eig_vec, test_val))

            if cov_test[0][0]==0 or cov_test[0][0]==cov_test[0][1]:
                print cov_test
                continue

            print cov_test
            score, score_list=ScoreClustering(Pr[:upper_size], cov_test, reps, 6, 100, samples)
            perf_list. append(numpy.mean( score_list))
            print  score
            scale_list.append(j)

    else:

        for j in numpy.linspace( (scale_1_max),(scale_1_min) ,step_count):
            print j
            test_val= [[j**2,0],[0,last_2**2]]
            cov_test=numpy.array( RecomposeMatrix(eig_vec, test_val))
            if cov_test[0][0]==0 or cov_test[0][0]==cov_test[0][1]:
                print cov_test
                continue
            print cov_test
            score, score_list=ScoreClustering(Pr[:upper_size], cov_test, reps, 6, 100, samples)
            perf_list. append(numpy.mean( score_list))

            print  score
            scale_list.append(j)

    return numpy.array( perf_list), numpy.array( scale_list)

#----------------------------------------------------------------------------#
#Functions for assessing performance on a dataset

def ValidateDirectory(indir, outfile, dist_file, cov, distance_cutoff=400, ref_names=['2L', '3L', '2R', '3R', 'X', '4', 'U']):
    """Compare the performance of the EM algorithm to the agglomerative cluster algorithm."""

    start_time=time.clock()
    m,k,p=ReadDistributions(dist_file)
    files=sorted( os.listdir(indir))
    outhandle=open(outfile, 'w')
    outTable=csv.writer(outhandle, delimiter='\t')
    for f in files:
        file_root=f.split('.')[0]
        if ref_names.count(file_root)!=0: continue
        infile='{0}/{1}'.format(indir, f)
        ValidateFile(infile, outTable, cov, k, p, distance_cutoff, ref_names=ref_names)
        outhandle.flush()
    outhandle.close()
    print "It took {0} seconds to process.".format (time.clock()-start_time)


def ValidateFile(infile,outTable, cov, k, p, distance_cutoff=400, ref_names=['2L', '3L', '2R', '3R', 'X', '4', 'U']):
    """Cluster the input file using single linkage agglomerative clustering, and
    then recluster each cluster using EM. Compare performance."""
    seq_name= infile.split('/')[-1].split('.')[0]

    distance_based_clusters=RoughClusters(infile, distance_cutoff, assembly_list=ref_names)
    scores, keys, labels, correct, incorrect=ValidateClustering( distance_based_clusters, cov)
    if len(scores)!=0:
        arith_mean=numpy.mean(scores)
        harm_mean=HarmonicMean(scores)
        fraction_consistent=float( (numpy.array( scores)==1.).sum())/len(scores)
        num_corr=sum(correct)
        num_incorr=sum(incorrect)
        total=num_corr+num_incorr
        row=[seq_name,len(scores) , arith_mean, harm_mean, fraction_consistent, num_corr, total, float(num_corr)/total]
    else:
        row=[seq_name,0 , 0, 0, 0]
    print row
    outTable.writerow(row)

def ValidateClustering(testData, cov):
    """"Recluster each input cluster using the provided covariance matrix.
    Assess performance."""
    score_list=[]
    key_list=[]
    label_list=[]
    correct_list=[]
    incorrect_list=[]
    for key in testData. keys():
        if key[0]==key[1]:
            cov_quad= numpy.array( cov)
        else: cov_quad=numpy.array(cov)*[[1,-1],[-1,1]]

        for index in range(len( testData[key])):
            t,l,lk=SquareEM(testData[key][index].transpose(), cov_quad, 100,100, Agglomerate=True)
            score=Simpson(l)
            correct_list.append(score*len(l))
            incorrect_list.append((1-score)*(len(l)))
            score_list.append(score)
            if score<1:
                key_list.append((key, index))
                label_list.append(l)
    return score_list,key_list, label_list, correct_list, incorrect_list

def TestPrecisionAlongAngle(ins, cov, max_displacement, degrees,steps=50, reps=10):
    """Examines behavior of the EM algorithm with a given covariance (cov) on a
    mixture of two distributions of data simulated from the insert size
    distribution (ins) two determine how far apart junctions must be to distinguish
    them. The relevant score is pairwise precision."""

    distribution=numpy.copy(ins)
    data_list,label_list, l_list =[],[],[]
    cdf=numpy.cumsum(distribution)
    low,upp=numpy.where(cdf>.005)[0][0], numpy.where(cdf>=.995)[0][0]
    distribution[:low]=0.
    distribution[upp:]=0.
    distribution/=distribution.sum()

    radians=math.radians(degrees)
    precision, recall=[],[]
    FM=[]
    disp=[]
    cluster_count=[]
    print numpy.sign(degrees)
##    orientation=numpy.asarray( [[1,numpy.sign(degrees)],[numpy.sign(degrees),1]],float)
##    cov=numpy.array( cov) *orientation
    print cov
    displacement=numpy.linspace(0, max_displacement, steps)
    for i in range( steps):
        for r in range(reps):
            #construct the mixture
            y=numpy.cos(radians)*displacement[i]
            x=-numpy.sin(radians)*displacement[i]
            dist1=numpy.vstack( GenSample(distribution,20)).transpose()
            dist2=numpy.vstack( GenSample(distribution,20)+numpy.array([x,y])[:,None]).transpose()
            labels=numpy.array( [0]*20+[1]*20)
            mixture=numpy.vstack((dist1, dist2))
            data_list.append(mixture)
            label_list.append(labels)
            #cluster
            t,l,lk=SquareEM(mixture,cov,100,200)
            l_list.append(l)
            mask=numpy.mask_indices(len(labels), numpy.tril,k=-1)
            real_pw=numpy.equal.outer(labels,labels)[mask].flatten().astype(int)
    ##            print real_pw
            identified_pw=numpy.equal.outer(l,l)[mask].flatten().astype(int)
            precision.append(sklearn.metrics.precision_score(real_pw, identified_pw))
            recall.append(sklearn.metrics.recall_score(real_pw, identified_pw))
            FM.append(sklearn.metrics.fowlkes_mallows_score(labels, l))
            cluster_count.append(len(set(l)))
            disp.append(displacement[i])

    return disp, precision, recall, FM


def TestOverRangeOfSizes(ins, cov):
    """Examines behavior of the EM algorithm with a given covariance (cov) on a
    mixture of data simulated from the insert size distribution (ins) for a range
    of true mixture sizes. The distributions are placed on a grid 500 units apart
    to ensure that none overlap. The performance is summarized with the number of
    identified clusters, adjusted Rand index, and Fowlkes-Mallows
    (an average of precision and recall)."""
    size,real, clusters, rand, FM=[],[],[],[],[]
    recall, precision=[],[]
    error_rate=[]
    good_sizes, bad_sizes=[],[]
    x,y=numpy.meshgrid(numpy.linspace(0, 10000,20.), numpy.linspace(0, 10000,20.))
    x=numpy.hstack(x)
    y=numpy.hstack(y)
    for j in range(2,20,1):
        for r in range(20):
            data, labels=SimulateDataFromPoints(ins,x,y,j,j,10)
            t,l,lk=SquareEM(data,cov,100,100, True)
            real.append(20)
            size.append(j)
            clusters.append(len(set(l)))
            rand.append(sklearn.metrics.adjusted_rand_score(labels, l))
            FM.append(sklearn.metrics.fowlkes_mallows_score(labels,l))
            real_pw=numpy.equal.outer(labels,labels).flatten().astype(int)
##            print real_pw
            identified_pw=numpy.equal.outer(l,l).flatten().astype(int)
            recall.append(sklearn.metrics.recall_score(real_pw, identified_pw))
            precision.append(sklearn.metrics.precision_score(real_pw, identified_pw))

            cluster_counts=collections.Counter(l)
            sorted_counts=sorted( cluster_counts.items(), reverse=True, key=lambda c:c[1])
            good_sizes+=[c[1] for c in sorted_counts[:j]]
            bad_sizes+=[c[1] for c in sorted_counts[j:]]
            error_rate.append((len(set(l))-j)/float(j))
    means=[numpy.mean(numpy.array(recall)[numpy.array(size)==i]) for i in sorted(list(set(size)))]
    return size, real, clusters, rand, FM, good_sizes, bad_sizes,error_rate, precision, recall, means


def ReadSpecificationFile(infile):
    handle=open(infile)
    spec_dict={}
    for line in handle:
        line=line.strip()
##        print line
        if len(line)==0: continue   #Empty line
        if line[0]=='#': continue #These are comments
        if line[0]=='[': continue #Section title
        param,value=line.split('=')[0], '='.join(line.split('=')[1:])
        param_type, param_name=param[0], param[1:]
        if param_type=='$':
            spec_dict[param_name]=value
        if param_type=='@':
            spec_dict[param_name]=value.split(',')
        if param_type=='%':
            spec_dict[param_name]=float( value)
        if param_type=='!':
            try:
                spec_dict[param_name]=float( value)
            except:
                spec_dict[param_name]=str(value).lower()
    return spec_dict

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

    infile= param['-i']
    outfile=param['-o']


    #Read the specification file
    spec_dict=ReadSpecificationFile(param['-spec'])

    #Check whether the specification file is properly specified.
##    executable_list=['Bowtie2', 'Samtools','Trimmomatic']
##    for executable in executable_list:
##        assert os.path.exists(spec_dict[executable]), 'The path to the {0}'

    assert spec_dict.has_key('Cvg'), 'Specification file must have @Cvg field specifying which chromosomes to use in read depth estimation.'
    assert spec_dict.has_key('CV'), 'Specification file must have @CV field specifying which chromosomes to use in the validation and parameter choice steps.'


    refFile=spec_dict['Ref']
    consFile=spec_dict['Cons']
    ref_seq=GetLengths(refFile)
    #Pipeline to cluster a directory of alignment outputs
    count=0



    GemCode=False
    if param.has_key('-head')==True:
        head_tag=param['-head']
    else:
        head_tag=''
    if param.has_key('-o')==True:
        MakeDir(param['-o'])
        summary_file=param['-o']+'/sample_summary_{0}.tsv'.format(head_tag)
        summary_handle=open( summary_file,'w')
        summary_table=csv.writer(summary_handle, delimiter='\t')
    else:
        summary_file=param['-i']+'/sample_summary_{0}.tsv'.format(head_tag)
        summary_handle=open( summary_file,'w')
        summary_table=csv.writer(summary_handle, delimiter='\t')
    header_fields=['Sample','Gap Size Distributions']+['']*3+['Agglomerative Clustering', 'Clustering Parameters','', 'Precision along eigenvector 1']+['']*6+ ['Precision along eigenvector 2']+['']*6+['Recall at read count']+['']*len(range(3,20))
    summary_header=['Sample Name', 'Mean Gap Size','Std Dev Gap Size', 'Gap Q0.5', 'Gap Q95.5','Distance Cutoff', 'Scale 1', 'Scale 2','A1','B1','Q1', 'M1', 'v1', 'Exp-phase', 'Prec>=.99','A2','B2','Q2', 'M2', 'v2', 'Exp-phase2', 'Prec>=.99 on Scale 2' ]+range(2,20)
    summary_table.writerow(header_fields)
    summary_table.writerow(summary_header)
    summary_handle.flush()
    if param.has_key('-G')==True: GemCode=True
    for directory, subdir, files in os.walk(infile):
        count+=1


        if param.has_key('-o')==False:
            sample_name=directory.split('/')[-1]
            dist_file=directory+'_kde.tsv'
            len_file= directory+'_len.tsv'
##                cvg_file=directory+'_cvg_hist.npy'
            cvg_file=directory+'_cvg'
            count_file=directory+"_counts.tsv"
            opt_file=directory+'_opt.tsv'
            val_file=directory+'_val.tsv'
            output_file=directory+'_out.tsv'
            hist_file=directory+'_cvg_hist.npy'
            eig_file1=directory+'_eigenvector1_precision.jpeg'
            eig_file2=directory+'_eigenvector2_precision.jpeg'
        else:
            sample_name=directory.split('/')[-1]
            sample_path='{0}/{1}'.format(param['-o'], sample_name)
            kde_file=sample_path+'.tsv'
            dist_file=sample_path+'_kde.tsv'
            len_file= sample_path+'_len.tsv'
##                cvg_file=directory+'_cvg_hist.npy'
            cvg_file=sample_path+'_cvg'
            hist_file=sample_path+'_cvg_hist.npy'
            count_file=sample_path+"_counts.tsv"
            opt_file=sample_path+'_opt.tsv'
            val_file=sample_path+'_val.tsv'
            output_file=sample_path+'_out.tsv'
            eig_file1=sample_path+'_eigenvector1_precision.jpeg'
            eig_file2=sample_path+'_eigenvector2_precision.jpeg'
            recall_file=sample_path+'_recall.jpeg'

        if count==1: continue

##        if param.has_key('-ref'):
        #Check whether file should be skipped if the -head
        if param.has_key('-head'):  #Only process files that begin with a certain character
            dir_root=directory.split('/')[-1]
            if dir_root[0].lower()!=param['-head'].lower(): continue

        if param.has_key('-find'): #Only process files that contain a certain string
            dir_root=directory.split('/')[-1]
            if dir_root.lower().find(param['-find'].lower())==-1 : continue

        print directory

        sample_summary=[sample_name]
        #Estimate the gap size distribution
        if param.has_key('-o')==True:
            WriteDistributions.WriteDistributions(directory, kde_file, names=spec_dict['Cvg'])
        else:
            WriteDistributions.WriteDistributions(directory,names=spec_dict['Cvg'])

        m,k,p=ReadDistributions(dist_file)
        mean_insert_size= WeightedAverage(p, numpy.arange(len(p)))
        sd_insert_size= WeightedStdDev(p, numpy.arange(len(p)))

        if spec_dict['Cov1']=='auto':
            #Along the first eigenvector read pairs should be approximately uniformly
            #distributed. The variance of read pairs with the mean gap size m is then
            #       var=((2m)^2)/12
            #So the standard deviation is:
            #       std= m/(3^.5)
            #So the mean insert size should generally be larger than the
            #optimal 1st eigenvector

            cov1=mean_insert_size
        else:
            assert type( spec_dict['Cov1'])==float, '!Cov1 must either be a float or AUTO.'
            cov1=float( spec_dict['Cov1'])

        if spec_dict['Cov1']=='auto':
            #Along the second eigenvector read pairs should be distributed according to the gap size
            #distribution. So the twice the gap size distribution's standard deviation
            #should generally be larger than the optimal second eigenvector

            cov2=sd_insert_size*2.
        else:
            assert type( spec_dict['Cov2'])==float, '!Cov1 must either be a float or AUTO.'
            cov2=float( spec_dict['Cov2'])


        #Maybe place a warning here if the cov matrix looks too small compared to the
        #gap size distribution
        cov=numpy.array([[cov1, cov2], [cov2, cov1]])
        sample_summary+=[mean_insert_size, sd_insert_size, k[0], k[1]]
        print k

        #Generate 10000 simulated distributions of 5 reads
        #Calculate the maximum distance between reads on each axis
        #Determine the 99th percentile of distances. Use this as the cutoff for the
        #simple clustering strategy
        sample_distances= DrawFromReadDistanceDistribution(p, 3, 10000)
        distance_cutoff=scipy.stats.scoreatpercentile(sample_distances, 99)
        print "Distance cutoff = {0}".format( distance_cutoff)
        sample_summary+=[distance_cutoff]

        #Construct the coverage distributions conditioned on %GC
        CoverageModeller.WriteCvgHistogramByArray(directory, dist_file,len_file, refFile, cvg_file, refNames=spec_dict['Cvg'])

        #Count the number of reads mapping to each consensus sequence
        #ReadConTExt_9.CountPEReads(directory,count_file , refFile, consFile, hist_file, len_file, dist_file)

        #Choose the covariance parameter for the library
        scale_x=OptimizeOnFileWithGP(directory, dist_file, opt_file, cov1, cov2, 'x', 200, 150, spec_dict['Cutoff1'], filter_x=cov1, filter_y=cov2, distance_cutoff=distance_cutoff,assembly_list= spec_dict['CV'],files= spec_dict['Train'])
        scale_y=OptimizeOnFileWithGP(directory, dist_file, opt_file, scale_x, cov2, 'y', 200, 150, spec_dict['Cutoff2'], filter_x=cov1, filter_y=cov2, distance_cutoff=distance_cutoff,assembly_list= spec_dict['CV'],files= spec_dict['Train'] )
        sample_summary+=[scale_x,scale_y]
        eig_vec=numpy.matrix([[ 0.70710678,  0.70710678], [-0.70710678,  0.70710678]])
        cov_opt=RecomposeMatrix(eig_vec, [[scale_x**2,0],[0,scale_y**2]])

        #Assess clustering performance on this dataset
        #Test precision along first eigenvector:
        displacement_1, precision_1, recall_1, FM_1=TestPrecisionAlongAngle(p, cov_opt, 4*scale_x, 45, 200,10)
        regression, shoulder1_point, inflection_point, shoulder2_point, prec_99=PlotsForContext.PlotPrecRecFM(displacement_1, precision_1, recall_1, FM_1,45)
        sample_summary+=[regression.A]+ list( regression.best_param)
        sample_summary+=[shoulder1_point, prec_99]
        pyplot.savefig(eig_file1)
        pyplot.close()

        #Test precision along second eigenvector:
        displacement_2, precision_2, recall_2, FM_2=TestPrecisionAlongAngle(p, cov_opt, 4*scale_y, -45, 200,10)
        regression, shoulder1_point, inflection_point, shoulder2_point, prec_99=PlotsForContext.PlotPrecRecFM(displacement_2, precision_2, recall_2, FM_2, -45)
        sample_summary+=[regression.A]+ list( regression.best_param)
        sample_summary+=[shoulder1_point, prec_99]
        pyplot.savefig(eig_file2)
        pyplot.close()

        opt=TestOverRangeOfSizes(p, cov_opt)
        mean_recall= PlotsForContext.PlotRecallAgainstSize(opt[0], opt[-2])
        sample_summary+=mean_recall
        pyplot.savefig(recall_file)
        pyplot.close()

    #inDir, outDir, AssemblyIndex, TEIndex, threads, length=70, ReadIDPosition=1, PositionRange=0, phred=33, shorten=True, Trim=True
        #Comparing
        summary_table.writerow(sample_summary)
        summary_handle.flush()
        ValidateDirectory(directory,val_file, dist_file, cov_opt, distance_cutoff)
        ClusterDirectory(directory,output_file, dist_file, cov_opt,ref_seq, GemCode)

if __name__ == '__main__':
    main(sys.argv)

