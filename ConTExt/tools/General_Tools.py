#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      I am
#
# Created:     31/01/2018
# Copyright:   (c) I am 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from Bio import SeqIO
import numpy



def GetSeq(ref, upper=False):
    """Reads a fasta, returns of a dictionary of strings keyed by entry name."""
    handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    SeqLen={}
    for rec in lib:
        SeqLen[CleanName(rec.name)]=str( rec.seq)
        if upper==True: SeqLen[CleanName(rec.name)]=SeqLen[CleanName(rec.name)].upper()
    handle.close()
    return SeqLen

def GetLengths(ref):
    """Reads a fasta, returns of a dictionary of strings keyed by entry name."""
    handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    SeqLen={}
    for rec in lib:
        SeqLen[CleanName(rec.name)]=len( str( rec.seq))


    handle.close()
    return SeqLen

def CleanName(name):
    illegal=['|', '!','>', '<', '?','/','*']
    for i in illegal:
        name='_'.join(name.split(i))
    return(name)

def MakeDir(newdir):
    if os.path.exists(newdir)==False:
        os.mkdir(newdir)

def RecomposeMatrix(eig_vec, eig_val_matrix):
    """Returns a covariance matrix from eigenvectors and eigenvalues"""
    return eig_vec*eig_val_matrix*numpy.linalg.inv(eig_vec)

def DecomposeMatrix(cov):
    """Decomposes 2D covariance matrix into its eigenvectors and eigenvalues"""
    cov_matrix= numpy.matrix(cov)
    eig_val, eig_vec=numpy.linalg.eig(cov_matrix)
    eig_val_diag=[[eig_val[0],0],[0,eig_val[1]]]
    return eig_vec, eig_val_diag

def ReadSpecificationFile(infile):
    handle=open(infile)
    spec_dict={}
    for line in handle:
        line=line.strip()
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

def WeightedAverage(weights, items):
    avg=numpy.nansum(weights*items)
    return avg

def WeightedStdDev(weights, items):
    avg=WeightedAverage(weights, items)
    sd=numpy.nansum(weights* (items-avg)**2)**.5
    return sd

def main():
    pass

if __name__ == '__main__':
    main()
