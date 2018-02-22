#-------------------------------------------------------------------------------
# Name:        curation
# Purpose:
#
# Author:      Michael P McGurk
#
# Created:     14/02/2018
# Copyright:   (c) I am 2018
# Licence:     MIT license; see accompanying license.txt in Git repository
#-------------------------------------------------------------------------------
import numpy
import scipy
import csv
import os
from General_Tools import *
from matplotlib import pyplot



def BuildLinkageMatrix(indir):
    file_list=os.listdir(indir)
    linkage_matrix=numpy.ndarray((len(file_list),len(file_list)))
    linkage_matrix.fill(0.)
    sorted_keys=sorted(file_list)
    rpt_names=['.'.join(s.split('.')[:-1]) for s in sorted_keys]
    for f in sorted_keys:
        inhandle=open(indir+'/'+f,'r')
        intable=csv.reader(inhandle, delimiter='\t')
        for row in intable:
            try:
                line=clustered_line(row)
            except:
                continue
            try:
                i=rpt_names.index(line.Seq1)
                j=rpt_names.index(line.Seq2)
            except:
                continue
            linkage_matrix[i,j]+=line.count
    return linkage_matrix


def main():
    pass

if __name__ == '__main__':
    main()
