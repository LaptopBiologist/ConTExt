#-------------------------------------------------------------------------------
# Name:        PlotsForContext
# Purpose:
#
# Author:      Michael Peter McGurk
#
# Created:     07/01/2017
# Copyright:   (c) Michael Peter McGurk 2017
#-------------------------------------------------------------------------------

import numpy
import scipy
import seaborn
import csv
import matplotlib
from matplotlib import pyplot

from matplotlib.patches import Ellipse
import sklearn
import sklearn.linear_model
from sklearn import gaussian_process
from sklearn.base import BaseEstimator

import GeneralizedLogisticRegression

matplotlib.rcParams['font.size']=14

def GenLogGrad(A,K,B,t,M,Q, v):
    delta=K-A
    mu=t-M
    alpha=B*(t-M)
    dB=delta*mu*Q*numpy.exp(-1*alpha)/(v*(1+Q*numpy.exp(-1*alpha))**(1/v+1) )

def GenLog(x,K, B,C, Q, M, v):
    """Generalized logistic function."""
    A=.5
    K=1.
    alpha=B*(x-M)
    return A+ (K-A)/(C+Q*numpy.exp(-1*alpha))**(1/v)

def PlotPrecision(x,y,alpha=.05):
##    x=numpy.array(x)
##    y=numpy.array(y)
    means=[numpy.mean(numpy.array(y)[numpy.array(x)==i]) for i in sorted(list(set(x)))]
    cond_dists=[numpy.array(y)[numpy.array(x)==i] for i in sorted(list(set(x)))]
##    seaborn.violinplot(x=sorted(list(set(x))), y=cond_dists, inner=None)
    pyplot.scatter(x,y, alpha=alpha, s=40, c='firebrick')
##    pyplot.plot(sorted(list(set(x))), means, c='dodgerblue', lw=3)
##    kernel=gaussian_process.kernels.PairwiseKernel( metric='sigmoid') *(gaussian_process.kernels.ConstantKernel()+gaussian_process.kernels.WhiteKernel())
    GLR=GeneralizedLogisticRegression.GeneralizedLogisticRegressor()
##    GP=sklearn.linear_model.LogisticRegression()
    GLR.fit(x, y)
    x_domain=sorted(list( set( x)))
    predicted=GLR.predict(x_domain)
    diff_pred=numpy.diff(predicted)

    inflection_ind= numpy.argmax(numpy.diff(predicted))+1
    inflection_point=x_domain[inflection_ind]
    print inflection_point
    prec_99= numpy.where(predicted>=.99)[0][0]
    diff2_pred=numpy.diff(diff_pred)
    shoulder1_point= x_domain[numpy.argmax(diff2_pred)]
    shoulder2_point= x_domain[numpy.argmin(diff2_pred)]

    pyplot.plot(x_domain, predicted, c='goldenrod', lw=3, alpha=1)
    pyplot.scatter(shoulder1_point, GLR.predict(shoulder1_point), c='grey', marker='D', zorder=5)
    pyplot.scatter(shoulder2_point, GLR.predict(shoulder2_point), c='grey',  marker='D', zorder=5)
    pyplot.scatter(inflection_point, GLR.predict (inflection_point), c='grey',  marker='D', zorder=5)
    pyplot.scatter(x_domain[ prec_99], GLR.predict(x_domain[ prec_99]), c='black',  marker='D', zorder=5)
##    pyplot.plot(x_domain[:-2], diff2_pred, c='b')
##    pyplot.xlabel('Read Count')
    pyplot.ylabel('Pairwise Precision')

##    pyplot.xticks(numpy.arange(0, 25), numpy.arange(0, 25))
##    pyplot.show()
    return GLR,  shoulder1_point,inflection_point, shoulder2_point, x_domain[ prec_99]
def PlotClusteringChoice(data, true_labels,param,l):
    x,y=numpy.meshgrid(numpy.linspace(-200, 800,200.), numpy.linspace(0, 1000,200.))
##    x=numpy.hstack(x)
##    y=numpy.hstack(y)
    grid=numpy.dstack((x,y))
    Gau1=scipy.stats.multivariate_normal(param[0].mean, t[0].cov)
    Gau2=scipy.stats.multivariate_normal(param[1].mean, t[1].cov)
    Lk1=t[0].weight* Gau1.pdf(grid)
    Lk2=t[1].weight* Gau2.pdf(grid)
    stackedLk=numpy.dstack((Lk1,Lk2))
    dec= numpy.nanargmax(stackedLk,2)
    print stackedLk.shape
    pyplot.contourf(x,y,dec, cmap='summer')
    pyplot.scatter(data[:,0], data[:,1],s=100, c=l)
    pyplot.scatter(data[:,0], data[:,1], c=true_labels,s=60, cmap='PiYG')
    pyplot.show()
    return x,y,Lk1, Lk2, dec
    pass


def PlotRecall(x,y,alpha=.05):
##    x=numpy.array(x)
##    y=numpy.array(y)
    means=[numpy.mean(numpy.array(y)[numpy.array(x)==i]) for i in sorted(list(set(x)))]
    cond_dists=[numpy.array(y)[numpy.array(x)==i] for i in sorted(list(set(x)))]

    pyplot.scatter(x,y, alpha=alpha, s=40, c='firebrick')

    pyplot.ylabel('Pairwise Recall')

##    print x_domain[numpy.argmax(numpy.diff(predicted))]

def PlotRecallAgainstSize(x,y):
##    x=numpy.array(x)
##    y=numpy.array(y)
    means=[numpy.mean(numpy.array(y)[numpy.array(x)==i]) for i in sorted(list(set(x)))]
    cond_dists=[numpy.array(y)[numpy.array(x)==i] for i in sorted(list(set(x)))]
    x_domain=sorted(list(set(x)))
    jitter=numpy.random.random(size=len(x))/4.-.125
    pyplot.scatter(x+jitter,y, alpha=.8, s=40, c='dodgerblue',zorder=4)
##    pyplot.plot(x_domain, means, c='dodgerblue', lw=3)
    bx_plt=pyplot.boxplot(cond_dists,positions= x_domain, widths=.8, showmeans=True)
    set_med_colors=[med.set_c('firebrick') for med in    bx_plt['medians']]
    set_med_colors=[med.set_linewidth(3) for med in    bx_plt['medians']]
    set_med_colors=[med.set_zorder(7) for med in    bx_plt['medians']]
    set_med_colors=[med.set_zorder(8) for med in    bx_plt['means']]
    set_med_colors=[med.set_c ('black') for med in    bx_plt['whiskers']]
    set_med_colors=[med.set_c ('black') for med in    bx_plt['caps']]
    set_med_colors=[med.set_linewidth(3) for med in    bx_plt['boxes']]
    set_med_colors=[med.set_c('black') for med in    bx_plt['boxes']]
    pyplot.ylabel('Pairwise Recall')
    pyplot.xlabel('Read Count')

##    print x_domain[numpy.argmax(numpy.diff(predicted))]
    return means

def PlotFM(x,y,alpha=.05):
##    x=numpy.array(x)
##    y=numpy.array(y)

##    seaborn.violinplot(x=sorted(list(set(x))), y=cond_dists, inner=None)
    pyplot.scatter(x,y, alpha=alpha, s=40, c='firebrick')
    pyplot.ylabel('Fowlkes-Mallow Score')

def PlotPrecRecFM(x, prec, rec, fm,degrees=-45,alpha=.1):
    ax=pyplot.subplot(311)
    GLR, shoulder1_point,inflection_point, shoulder2_point, prec_99= PlotPrecision(x,prec,alpha)
    pyplot.subplot(312, sharex=ax)
    PlotRecall(x,rec,alpha)
    pyplot.subplot(313, sharex=ax)
    PlotFM(x,fm,alpha)
    pyplot.xlim(min(x), max(x))
    pyplot.xlabel('Distance Between Components at {0} deg'.format(degrees))
    pyplot.show

    return GLR, shoulder1_point, inflection_point, shoulder2_point, prec_99


def GenSample(insert_size_probability, size=5):
    """Given an insert size distribution, generate a distribution of read pairs
    spanning a junction:
            insert_size_probability: The insert size distribution
            size:   The number of read pairs to simulate"""

    value_range=numpy.array( range( len(insert_size_probability)))
    x_list, y_list=[],[]
    Pr=numpy.array( insert_size_probability)/sum(insert_size_probability)
    insert_vals= numpy.random.choice(value_range, size, p=Pr)
    for val in insert_vals:
        x_val=numpy.random.random_integers(0, val)
##        if 50< x_val <80: continue
        y_val=val-x_val
        x_list.append(x_val)
        y_list.append(y_val)
    return x_list, y_list

def WeightedMean(values, weights):
    avg=numpy.nansum(weights*values)
    return avg

def ReadDistributions(infile, cutoff=.005, return_full=False):
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
        return m[d>=0]/m[d>=0].sum(),k[d>=0]/k[d>=0].sum(),p/p.sum()

    return m, (min_cutoff, max_cutoff), p

def PlotInsertSizeDistribution(m,k,p):
    pyplot.scatter(numpy.arange(len(m)), m, s=20, alpha=1.)
    pyplot.plot(k, c='r', alpha=.7, lw=1.5)
    cumul_k=numpy.cumsum(k)
##    print cumul_k
##    print numpy.where(cumul_k>=(1.-cutoff))
    max_cutoff=numpy.where(cumul_k>=(1.-.005))[0][0]
    min_cutoff=numpy.where(cumul_k>=(.005))[0][0]
    pyplot.xlim(min_cutoff, max_cutoff)
    pyplot.ylim(0, 1.1*max(k))
    pyplot.show()


def RecomposeMatrix(eig_vec, eig_val_matrix):
    return eig_vec*eig_val_matrix*numpy.linalg.inv(eig_vec)

def PlotScatterWithCovAndVectors(ins, scale_1, scale_2):
    mean=WeightedMean(numpy.arange(len(ins)), ins)
    x=numpy.linspace(-700,700,5.)
    y=numpy.linspace(-700,700,5.)

    X_pos,Y_pos=numpy.meshgrid(x,y)

    pos = numpy.empty(X_pos.shape + (2,))
    pos[:, :, 0] = X_pos; pos[:, :, 1] = Y_pos

##    X=numpy.hstack(X_pos)
##    Y=numpy.hstack(Y_pos)


##    grid= numpy.vstack(( X,Y)).transpose()

    eig_vec=numpy.matrix([[ 0.70710678,  0.70710678], [-0.70710678,  0.70710678]])
    cov=numpy.array( RecomposeMatrix(eig_vec, [[scale_1**2,0],[0,scale_2**2]]))

    pdf=scipy.stats.multivariate_normal( (mean,mean), cov)
    Z=pdf.logpdf(pos)
##    seaborn.heatmap(Z)
    pyplot.contour(X_pos, Y_pos, Z)
    pyplot.show()

    ax=pyplot.gca()
    confidence_region=Ellipse((mean/2., mean/2.), 4*scale_1, 4*scale_2,-45, alpha=.5)
    confidence_region.set_edgecolor('dodgerblue')
    confidence_region.set_linewidth(1.5)
    confidence_region.set_fill(True)
    sample=GenSample(ins, 50)
    pyplot.scatter(sample[0], sample[1])
    ax.add_artist(confidence_region)
    pyplot.show()

def GenLog(x,param):
        """Generalized logistic function."""
##        A=.5
        C=1
        K=1.
        A, B, Q, M, v=param
        alpha=B*(x-M)
        return A+ (K-A)/(C+Q*numpy.exp(-1*alpha))**(1/v)


def PlotSummaryTable(path):
    ax=pyplot.gca()
    ax.yaxis. grid(True)
    #Read the table
    handle=open(path, 'r')
    intable=csv.reader(handle, delimiter='\t')
    row_holder=[]
    intable.next()
    intable.next()
    for row in intable:
        if row[0].split('_')[-1]!='2':
            row_holder.append(row[1:])

    table=numpy.asarray(row_holder, float)

    sample_num, summary_num=table.shape
    PlotGenLog(numpy.arange(800), table[:, 7:12])
    pyplot.ylabel('Precision', size=20)
    pyplot.xlabel('Distance', size=20)
    pyplot.xticks(size=14)
    pyplot.yticks(size=14)
    pyplot.show()

    ax=pyplot.gca()
    ax.yaxis. grid(True)
    PlotGenLog(numpy.arange(100), table[:, 14:19])
    pyplot.ylabel('Precision', size=20)
    pyplot.xlabel('Distance', size=20)
    pyplot.xticks(size=14)
    pyplot.yticks(size=14)
    pyplot.show()

    seaborn.distplot(table[:,5])
    pyplot.show()
    seaborn.distplot(table[:,6])
    pyplot.show()
    RecallBoxPlot(table[:,-18:])
    pyplot.show()

    return table


def PlotGenLog(x_domain ,table):
    regressions=[GenLog(x_domain, table[i,:]) for i in range(table.shape[0])]
    min_lin,max_lim= min([min(r) for r in regressions]), max([max(r) for r in regressions])
    f=[pyplot.plot(x_domain, regressions[i] , c='r', alpha=.1) for i in range(table.shape[0])]
    pyplot.ylim(min_lin-.01, max_lim+.01)
def RecallBoxPlot(table):
    print range(2,20)
    print numpy.mean(table,0)
    print numpy.median(table,0)
    seaborn.set_style('dark')
    ax=pyplot.gca()

    bx_plt=ax.boxplot(table,positions=range(2,20), widths=.8, showmeans=True, patch_artist=True)

    set_med_colors=[med.set_color('dodgerblue') for med in    bx_plt['medians']]
    set_med_colors=[med.set_alpha(.8) for med in    bx_plt['medians']]
    set_med_colors=[med.set_linewidth(3) for med in    bx_plt['medians']]
    set_med_colors=[med.set_zorder(7) for med in    bx_plt['medians']]
    set_med_colors=[med.set_zorder(8) for med in    bx_plt['means']]
    set_med_colors=[med.set_color ('dodgerblue')  for med in    bx_plt['whiskers']]
    set_med_colors=[med.set_color ('dodgerblue') for med in    bx_plt['caps']]
    set_med_colors=[med.set_linewidth(3) for med in    bx_plt['boxes']]
    set_med_colors=[med.set_linewidth(3) for med in    bx_plt['boxes']]
    set_med_colors=[med.set_edgecolor('dodgerblue') for med in    bx_plt['boxes']]
    set_med_colors=[med.set_facecolor('aliceblue') for med in    bx_plt['boxes']]
    set_med_colors=[med.set_alpha(.8) for med in    bx_plt['boxes']]

    pyplot.ylabel('Pairwise Recall', size=14)
    pyplot.xlabel('Read Count', size=14)
    pyplot.yticks(size=14)
    pyplot.xticks(range(2,20), range(2,20), size=14)
    pyplot.ylim(.795,1.005)

    ax.yaxis. grid(True)

    return bx_plt

def main():
    pass

if __name__ == '__main__':
    main()
