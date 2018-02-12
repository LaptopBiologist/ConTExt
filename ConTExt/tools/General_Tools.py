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

import csv
class FileManager():
    """A class for opening many handles at once"""
    def __init__(self, maxfiles=100, delimiter='', mode='w'):
        self.maxfiles=maxfiles
        self.__handle_dictionary={}
        self.__wrapped_handles={}
        self.delimiter=delimiter
        self.mode=mode

    def __len__(self):
        return len(self.__handle_dictionary)

    def __getitem__(self, _file):
        """If the file is already open"""

        #If the file isn't currently open, open it
        if self.__handle_dictionary.has_key(_file)==False:
            if len(self)== self.maxfiles:   #Need to close a file first
                file_to_close= self.__handle_dictionary.keys()[-1]
                self.__handle_dictionary[file_to_close].close()
                del self.__handle_dictionary[file_to_close]
                del self.__wrapped_handles[file_to_close]

            self.__handle_dictionary[_file]=open(_file, 'a')

            #If told the file is delimited, open it with the csv parser
            if self.delimiter!='':
                if self.mode=='r':
                    self.__wrapped_handles[_file]=csv.reader(self.__handle_dictionary[_file], delimiter=self.delimiter)
                if self.mode=='w':
                    self.__wrapped_handles[_file]=csv.writer(self.__handle_dictionary[_file],delimiter=self.delimiter)

        #If it's not a delimited file, return the handle
        if self.delimiter=='':
            return self.__handle_dictionary[_file]
        #If it is a delimited file, return the parser
        else:
            return self.__wrapped_handles[_file]

    def close(self):
        for _file in self.__handle_dictionary.keys():
            self.__handle_dictionary[_file].close()

        self.__handle_dictionary={}
        self.__wrapped_handles={}



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
    bool_dict={'true':True, 'false':False}
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
        if param_type=='?':
            spec_dict[param_name]=bool_dict[str(value).lower()]

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
