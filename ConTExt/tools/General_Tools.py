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
from Bio import Seq
import numpy
import gzip
import csv
import os

class cluster():
    def __init__(self, row , Seq1=0, Seq2=1, Quad1=2, Quad2=3,\
    ID=0, feature=4, jxn_x=5, jxn_y=6, range_x=7 , range_y=8, count=9, x_list=10,\
    y_list=11, x5_list=12, y5_list=13, mapQ_x=14, mapQ_y=15, counter=16, cigar_x=17, cigar_y=18, MD_x=19, MD_y=20, QD_x=21, QD_y=22):
        self.Seq1=row[Seq1]
        self.Seq2=row[Seq2]
        self.Quad1=row[Quad1]
        self.Quad2=row[Quad2]
        self.ID=int(0)
        self.feature=row[feature]
        if row[jxn_x]!='err' and row[jxn_x]!='nan' and row[jxn_x]!='nan' and row[jxn_y]!='nan':
            try:
                self.j0_x=float(row[jxn_x])
            except:
                print row
                print row[jxn_x]
                print jabber
            try:
                self.j0_y=float(row[jxn_y])
            except:
                print row
                print jabber
            try:
##                self.x_min, self.x_max=[toInt( s) for s in row[range_x].split(';')]
                self.err_1=float( row[range_x])
            except:
                self.err_1=numpy.nan()
            try:
                self.err_2=float(row[range_y])
            except:
                self.err_2=numpy.nan()
##                print row[range_x].split('-')
##                print jabber
        else: self.j0_x,self.j0_y, self.err_1, self.err_2 =numpy.nan, numpy.nan, numpy.nan,numpy.nan

        self.CI=0
        self.ID=row[ counter]
        self.exact=False
        self.count=int(row[count])
        try:

            self.x_list=[int(x) for x in row[x_list].split(';')]
            self.y_list=[int(y) for y in row[y_list].split(';')]
        except:
            self.x_list=[]
            self.y_list=[]
        try:

            self.x5_list=[int(x) for x in row[x5_list].split(';')]
            self.y5_list=[int(y) for y in row[y5_list].split(';')]
        except:
            self.x5_list=[]
            self.y5_list=[]
        try:
            self.CigarX=[x for x in row[cigar_x].split(';')]
            self.CigarY=[y for y in row[cigar_y].split(';')]
        except:
            self.CigarX=[]
            self.CigarY=[]

        self.MD_X=[x.split(':')[-1] for x in row[MD_x].split(';')]
        self.MD_Y=[y.split(':')[-1] for y in row[MD_y].split(';')]
        self.QD_X=[x.split(':')[-1] for x in row[QD_x].split(';')]
        self.QD_Y=[y.split(':')[-1] for y in row[QD_y].split(';')]


        try:
            self.MapQX=[int(x) for x in row[mapQ_x].split(';')]
            self.MapQY=[int(y) for y in row[mapQ_y].split(';')]
        except:
            self.MapQX=[]
            self.MapQY=[]
        if len(row)>21:
            self.optional=row[21:]
            self.barcodes=row[21].split(';')
        else: self.optional, self.barcodes =[],[]

        self.labels=["ID: {0}, Read Count: {1}".format(self.ID,self.count)]*count

    def __repr__(self):
        print "Seq1: {0}\tSeq2: {1}\nQuadrant: ({2}, {3}); ID: {4}".format(self.Seq1, self.Seq2, self.Quad1, self.Quad2, self.ID)
        print "J0: ({0}, {1}); Exact: {2} Offset: {3} +/- {4}".format(self.j0_x, self.j0_y, self.exact, self.offset, self.CI)
        print "Read Count: {0}".format(self.count)
        if len (self.x_list)<=20:
            print "X positions: {0}".format(','.join([str(x) for x in self.x_list]))
            print "Y positions: {0}".format(','.join([str(y) for y in self.y_list]))
        else:
            print "X positions: {0},...,{1}".format(','.join([str(x) for x in self.x_list[:10]]),','.join([str(x) for x in self.x_list[10:]]))
            print "Y positions: {0},...,{1}".format(','.join([str(y) for y in self.y_list[:10]]),','.join([str(y) for y in self.y_list[10:]]))
        return ''

    def row(self):
        Identifier=[self.Seq1, self.Seq2, self.Quad1, self.Quad2,  self.feature]
        jxn=[self.j0_x,  self.j0_y, self.err_1, self.err_2]
        positions=[self.count,';'.join([str(x) for x in self.x_list]), ';'.join([str(y) for y in self.y_list]),';'.join([str(x) for x in self.x5_list]), ';'.join([str(y) for y in self.y5_list]), ';'.join([str(x) for x in self.MapQX]), ';'.join([str(y) for y in self.MapQY]), self.ID, ';'.join([str(x) for x in self.CigarX]), ';'.join([str(y) for y in self.CigarY]), ';'.join([str(x) for x in self.MD_X]), ';'.join([str(y) for y in self.MD_Y]), ';'.join([str(x) for x in self.QD_X]), ';'.join([str(y) for y in self.QD_Y])]
        return (Identifier+jxn+positions+ self.optional)

    def inverse_row(self):
        Identifier=[self.Seq2, self.Seq1, self.Quad2, self.Quad1, self.feature]
        jxn=[self.j0_y, self.j0_x,self.err_1, self.err_2]
        positions=[self.count,';'.join([str(y) for y in self.y_list]), ';'.join([str(x) for x in self.x_list]),';'.join([str(y) for y in self.y5_list]), ';'.join([str(x) for x in self.x5_list]), ';'.join([str(y) for y in self.MapQY]), ';'.join([str(x) for x in self.MapQX]), self.ID, ';'.join([str(y) for y in self.CigarY]), ';'.join([str(x) for x in self.CigarX]), ';'.join([str(y) for y in self.MD_Y]), ';'.join([str(x) for x in self.MD_X]), ';'.join([str(y) for y in self.QD_Y]), ';'.join([str(x) for x in self.QD_X])]

        return (Identifier+jxn+positions+ self.optional)

    def ClassifyJunction(self):
        if self.feature=='Consensus': return ()
        unique_reads=len(set(self.x5_list ))*len(set(self.y5_list))
        if unique_reads==1:
            self.feature='Probable Artifact'
            return()
        if self.Seq1!=self.Seq2: return ()
        if self.Quad1=='-' and self.Quad2=='+':
            if self.j0_x<self.j0_y: self.feature='Tandem'
            else: self.feature='Deletion'
        if self.Quad1==self.Quad2:
            self.feature='Inversion'



class clustered_line():
    def __init__(self, row , Seq1=0, Seq2=1, Quad1=2, Quad2=3,\
    ID=0, feature=4, jxn_x=5, jxn_y=6, range_x=7 , range_y=8, count=9, x_list=10,\
    y_list=11, x5_list=12, y5_list=13, mapQ_x=14, mapQ_y=15, counter=16, cigar_x=17, cigar_y=18, MD_x=19, MD_y=20):
        self.Seq1=row[Seq1]
        self.Seq2=row[Seq2]
        self.Quad1=row[Quad1]
        self.Quad2=row[Quad2]
        self.ID=int(0)
        self.feature=row[feature]
        if row[jxn_x]!='err' and row[jxn_x]!='nan' and row[jxn_x]!='nan' and row[jxn_y]!='nan':
            try:
                self.j0_x=float(row[jxn_x])
            except:
                print row
                print row[jxn_x]
                print jabber
            try:
                self.j0_y=float(row[jxn_y])
            except:
                print row
                print jabber
            try:
##                self.x_min, self.x_max=[toInt( s) for s in row[range_x].split(';')]
                self.err_1=float( row[range_x])
            except:
                self.err_1=numpy.nan()
            try:
                self.err_2=float(row[range_y])
            except:
                self.err_2=numpy.nan()
##                print row[range_x].split('-')
##                print jabber
        else: self.j0_x,self.j0_y, self.err_1, self.err_2 =numpy.nan, numpy.nan, numpy.nan,numpy.nan

        self.CI=0
        self.ID=row[ counter]
        self.exact=False
        self.count=int(row[count])
        try:

            self.x_list=[int(x) for x in row[x_list].split(';')]
            self.y_list=[int(y) for y in row[y_list].split(';')]
        except:
            self.x_list=[]
            self.y_list=[]
        try:

            self.x5_list=[int(x) for x in row[x5_list].split(';')]
            self.y5_list=[int(y) for y in row[y5_list].split(';')]
        except:
            self.x5_list=[]
            self.y5_list=[]
        try:
            self.CigarX=[x for x in row[cigar_x].split(';')]
            self.CigarY=[y for y in row[cigar_y].split(';')]
        except:
            self.CigarX=[]
            self.CigarY=[]

        self.MD_X=[x.split(':')[-1] for x in row[MD_x].split(';')]
        self.MD_Y=[y.split(':')[-1] for y in row[MD_y].split(';')]

        try:
            self.MapQX=[int(x) for x in row[mapQ_x].split(';')]
            self.MapQY=[int(y) for y in row[mapQ_y].split(';')]
        except:
            self.MapQX=[]
            self.MapQY=[]
        if len(row)>21:
            self.optional=row[21:]
            self.barcodes=row[21].split(';')
        else: self.optional, self.barcodes =[],[]

        self.labels=["ID: {0}, Read Count: {1}".format(self.ID,self.count)]*count

    def __repr__(self):
        print "Seq1: {0}\tSeq2: {1}\nQuadrant: ({2}, {3}); ID: {4}".format(self.Seq1, self.Seq2, self.Quad1, self.Quad2, self.ID)
        print "J0: ({0}, {1}); Exact: {2} Offset: {3} +/- {4}".format(self.j0_x, self.j0_y, self.exact, self.offset, self.CI)
        print "Read Count: {0}".format(self.count)
        if len (self.x_list)<=20:
            print "X positions: {0}".format(','.join([str(x) for x in self.x_list]))
            print "Y positions: {0}".format(','.join([str(y) for y in self.y_list]))
        else:
            print "X positions: {0},...,{1}".format(','.join([str(x) for x in self.x_list[:10]]),','.join([str(x) for x in self.x_list[10:]]))
            print "Y positions: {0},...,{1}".format(','.join([str(y) for y in self.y_list[:10]]),','.join([str(y) for y in self.y_list[10:]]))
        return ''

    def row(self):
        Identifier=[self.Seq1, self.Seq2, self.Quad1, self.Quad2,  self.feature]
        jxn=[self.j0_x,  self.j0_y, self.err_1, self.err_2]
        positions=[self.count,';'.join([str(x) for x in self.x_list]), ';'.join([str(y) for y in self.y_list]),';'.join([str(x) for x in self.x5_list]), ';'.join([str(y) for y in self.y5_list]), ';'.join([str(x) for x in self.MapQX]), ';'.join([str(y) for y in self.MapQY]), self.ID, ';'.join([str(x) for x in self.CigarX]), ';'.join([str(y) for y in self.CigarY]), ';'.join([str(x) for x in self.MD_X]), ';'.join([str(y) for y in self.MD_Y])]
        return (Identifier+jxn+positions+ self.optional)

    def inverse_row(self):
        Identifier=[self.Seq2, self.Seq1, self.Quad2, self.Quad1, self.feature]
        jxn=[self.j0_y, self.j0_x,self.err_1, self.err_2]
        positions=[self.count,';'.join([str(y) for y in self.y_list]), ';'.join([str(x) for x in self.x_list]),';'.join([str(y) for y in self.y5_list]), ';'.join([str(x) for x in self.x5_list]), ';'.join([str(y) for y in self.MapQY]), ';'.join([str(x) for x in self.MapQX]), self.ID, ';'.join([str(y) for y in self.CigarY]), ';'.join([str(x) for x in self.CigarX]), ';'.join([str(y) for y in self.MD_Y]), ';'.join([str(x) for x in self.MD_X])]

        return (Identifier+jxn+positions+ self.optional)

    def ClassifyJunction(self):
        if self.feature=='Consensus': return ()
        unique_reads=len(set(self.x5_list ))*len(set(self.y5_list))
        if unique_reads==1:
            self.feature='Probable Artifact'
            return()
        if self.Seq1!=self.Seq2: return ()
        if self.Quad1=='-' and self.Quad2=='+':
            if self.j0_x<self.j0_y: self.feature='Tandem'
            else: self.feature='Deletion'
        if self.Quad1==self.Quad2:
            self.feature='Inversion'




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

def GetLengths(ref):
    """Reads a Fasta, returns a dictionary storing the lenghts of the sequences."""
    if ref[-2:]=='gz':
        handle=gzip.open(ref, 'r')
    else:
        handle=open(ref, 'r')
    lib=SeqIO.parse(handle, 'fasta')
    SeqLen={}
    for rec in lib:
        SeqLen[CleanName(rec.name)] = len(rec.seq)
    handle.close()
    return SeqLen

def GetSeq(ref, upper=False):
    """Reads a fasta, returns of a dictionary of strings keyed by entry name."""
    if ref[-2:]=='gz':
        handle=gzip.open(ref, 'r')
    else:
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

def ReadReferenceSummaries(infile):
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    ref_dict={}
    cons_dict={}
    entry_dict={}
    for row in intable:
        entry_type, entry_name=row[0], row[1]

        entry_dict[entry_name]=True
        if entry_type=='@ref':
            ref_dict[entry_name]=True
        if entry_type=='@cons':
            cons_dict[entry_name]=True
    inhandle.close()
    return entry_dict, cons_dict, ref_dict

def main():
    pass

if __name__ == '__main__':
    main()
