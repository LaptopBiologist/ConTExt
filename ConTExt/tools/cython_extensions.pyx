# distutils: language = c++
#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      I am
#
# Created:     26/03/2018
# Copyright:   (c) I am 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#distutils: language=c++
#C++ imports
import math
from libc.stdio cimport sscanf
from libcpp.vector cimport vector
from libc.math cimport exp
from libc.math cimport log
#Python imports
import time
import numpy
import cython

cimport numpy
ctypedef numpy.int_t DTYPE_t
cdef list operations=['=','D','I','X']




#@cython.boundscheck(False)

@cython.boundscheck(False)	
cpdef ComputeAlignmentScore(const unsigned char [:] quality_string, long[:] read_match_locations,long[:] read_mismatch_locations, long indels, long num_matches, long num_mismatches, double phred_base=33, double phred_scale= .23):

	cdef:
		int string_length=quality_string.shape[0]
		
		
		
		
		long match_pos
		long mismatch_pos
		double prob
		double score=0
		double quality
	
	#Interpret the CIGAR
	#read_match_locations, read_mismatch_locations,indels, num_matches, num_mismatches=ParseEdlibCigarPython(cigar)
	#Handle the matches
	for match_pos in range(num_matches):
		quality=quality_string[match_pos]
		prob=exp( -1*(quality-phred_base)*phred_scale)
		score+=log( 1-prob)
	#Handle the mismatches
	for mismatch_pos in range(num_mismatches):
		quality_string[mismatch_pos]
		prob=exp( -1*(quality-phred_base)*phred_scale)
		score+=log(prob/3.)
	score-=indels
	#print score
	return score
	
	
@cython.boundscheck(False)			
cpdef ParseEdlibCigar(const unsigned char [:] cigar):
	cdef:
		vector[int] read_mismatch_locations
		vector[int] read_match_locations
		char c
		int read_pointer=0
		int pos
		int indel=0
		int string_length=cigar.shape[0]
		int number_of_edits
		int num_matches=0
		int num_mismatches=0
		
	number_of_edits=0
	for current_index in xrange(string_length):
		c=cigar[current_index]
		if c<61:  #This is a number, not an edit operation
			#Update the number of edits for the upcoming operation
			number_of_edits*=10  #Each new number increments the digit place of number_of_edits
			number_of_edits+=c-48	#Characters are encoded as ASCII; zero is position 48
			continue
		if c==61:  #"='
			for pos in range(read_pointer,read_pointer+number_of_edits):
				read_match_locations.push_back(pos)
				num_matches+=1
			read_pointer+=number_of_edits
		if c==88:
			for pos in range(read_pointer,read_pointer+number_of_edits):
				read_mismatch_locations.push_back(pos)
				num_mismatches+=1
			read_pointer+=number_of_edits
		if c==73:
			read_pointer+=number_of_edits
			indel+=number_of_edits
		if c==68:
			indel+=number_of_edits
		number_of_edits=0 	#Reset the number of_edits

	return read_match_locations[:],read_mismatch_locations[:], indel, num_matches, num_mismatches
	
@cython.boundscheck(False)	
cpdef ClusterByDistances(long [:] pos_list, long cutoff):
	cdef:
		int	list_size=pos_list.shape[0]
		int i
		long distance
		vector[int] split_indices
		 
	split_indices.push_back(0)
	for i in range(1, list_size):
		distance=pos_list[i]-pos_list[i-1]
		if distance>cutoff:
			split_indices.push_back(i)
	split_indices.push_back(list_size)
	return split_indices[:]

cpdef FindSingletons(long [:] pos_list, long cutoff):
	cdef:
		int	list_size=pos_list.shape[0]
		int i
		long distance
		vector[int] split_indices
	
	split_indices.push_back(0)
	for i in range(1, list_size):
		distance=pos_list[i]-pos_list[i-1]
		if distance>cutoff:
			split_indices.push_back(i)
	split_indices.push_back(list_size)
	return split_indices[:]	
	
	

def main():
	pass
