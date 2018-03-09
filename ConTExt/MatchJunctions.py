#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      I am
#
# Created:     30/01/2018
# Copyright:   (c) I am 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy

import scipy
import scipy.stats
import scipy.spatial

import os
import sys
import csv
import time
import shutil

from tools.General_Tools import *

csv.field_size_limit(sys.maxsize)

#------------------------------------------------------------------------------#
#
#Object classes
#
#------------------------------------------------------------------------------#

class MetaRow():
    """Class for parsing the ConTExt metatable output."""
    def __init__(self, row, include=False):
        self.sample=row[0].split('_')[0]
        self.replicate=row[0].split('_')[1]
        self.ID=row[1]
        self.seq1=row[2]
        self.seq2=row[3]
        self.quad1=row[4]
        self.quad2=row[5]
        self.feature=row[6]
        self.readcount=float( row[7])
        self.x=float( row[8])
        self.y=float(row[9])
        if self.feature!='Consensus':
            self.err1=float( row[10])
            self.err2=float(row[11])
            self.expGC=float(row[13])
            self.CN_map=float(row[14])
            self.CN_low=float(row[15])
            self.CN_upp=float(row[16])
        else:
            self.err1=numpy.nan
            self.err2=numpy.nan
            self.expGC=numpy.nan
            self.CN_map=numpy.nan
            self.CN_low=numpy.nan
            self.CN_upp=numpy.nan

class MetaJxn():
    """Class for organizing information about junctions identified by paired-end reads."""
    def __init__(self):
        self.IDs=[]
        self.Xs=[]
        self.Ys=[]
        self.XYs=[]
        self.readcount=[]
    def add_jxn(self, line):
        self.IDs.append(line.ID)
        self.readcount.append(line.readcount)
        self.Xs.append(line.x)
        self.Ys.append(line.y)
        self.XYs.append( [line.x, line.y])

class Insertion():
    def __init__(self):
        self.r_ids=[]
        self.r_ins=[] #Where on the reference chrom
        self.r_edge=[]     #Where in the element
        self.r_count=[]
        self.l_ids=[]
        self.l_ins=[]
        self.l_edge=[]
        self.l_count=[]
        self.motif=[]
        self.type='single'
        self.gene=''
        self.dist_to_TSS=''
    def AddRightJxn(self, ID, position,edge,  readcount ):
        self.r_ids.append(ID)
        self.r_ins.append(position)
        self.r_edge.append(edge)
        self.r_count.append(readcount)
    def AddLeftJxn(self, ID, position,edge,  readcount ):
        self.l_ids.append(ID)
        self.l_ins.append(position)
        self.l_edge.append(edge)
        self.l_count.append(readcount)
    def AddMotif(self, row):
        self.motif.append(row)
        self.type='tandem'
    def AddAnnotation(self, gene, distance):
        self.gene=gene
        self.TSS_los=distance
    def FinalizeInsertion(self):
        if self.r_ids!=[]:
            self.r_pos_est=sum( numpy.array(self.r_ins)*self.r_count/sum(self.r_count))
            self.r_edge_est=sum( numpy.array(self.r_edge)*self.r_count/sum(self.r_count))
        else:
            self.r_pos_est=numpy.nan
            self.r_edge_est=numpy.nan
        if self.l_ids!=[]:
            self.l_pos_est=sum( numpy.array(self.l_ins)*self.l_count/sum(self.l_count))
            self.l_edge_est=sum( numpy.array(self.l_edge)*self.l_count/sum(self.l_count))
        else:
            self.l_pos_est=numpy.nan
            self.l_edge_est=numpy.nan
        self.ids=list(self.l_ids) +list( self.r_ids)
        self.pos=numpy.nanmean([self.l_pos_est,self.r_pos_est])
    def rows(self,counts=False):
        self.IDs=','.join(self.ids)
        self.insertion_info= [self.IDs, self.pos, self.l_pos_est,  self.l_edge_est,self.r_pos_est, self.r_edge_est]
        if counts==True:
            return [self.insertion_info+[ sum( self.l_count),sum( self.r_count)]]
        if self.type=='single':
            self.dist_to_TSS=numpy.nanmin([ abs(self.l_pos_est-self.TSS_los), abs(self.r_pos_est-self.TSS_los)])
            return [self.insertion_info+[self.type, self.gene, self.dist_to_TSS]]
        if self.type=='tandem':
            row_list=[]
            for m in self.motif:
                row_dict=m.row()
                for key in row_dict.keys():
                    if len(set(self.ids)&set(key))==0: continue
                    motif_loc=int(row_dict[key][0])
                    self.dist_to_TSS= abs(motif_loc-self.TSS_los)
                    row_list.append(self.insertion_info+[self.type, self.gene, self.dist_to_TSS] +row_dict[key])
            return row_list

class ClusterInformation():
    """Class for organizing information about clusters of junctions and for
    combining clulsters."""
    def __init__(self):
        self.CN=0.
        self.pos_x=[]
        self.pos_y=[]
        self.IDs=[]
        self.mean_x=numpy.nan
        self.mean_y=numpy.nan
    def add_jxn(self, CN, pos_x, pos_y, ID):
        self.CN+=CN
        self.pos_x.append(str( pos_x))
        self.mean_x=numpy.mean(numpy.asarray(self.pos_x, float))
        self.pos_y.append(str( pos_y))
        self.mean_y=numpy.mean(numpy.asarray(self.pos_y, float))
        self.IDs.append(str( ID))
    def join_jxns(self, new_jxn):
        self.CN+=new_jxn.CN
        self.pos_x+=new_jxn.pos_x
        self.pos_y+=new_jxn.pos_y
        self.IDs+=new_jxn.IDs

#------------------------------------------------------------------------------#
#
#
#
#------------------------------------------------------------------------------#


def BuildMetaTable(indir, in_tail, outdir):
    MakeDir(outdir)
    MakeDir(outdir+'/cvg_dist')
    files=os.listdir(indir)
    image_directory={}
    for f in files:

        #Only process files with the appropriate tag
        filename=f.split('/')[-1]
        root=filename.split('.')[0]
        tail=root.split('_')[-1]
        name='_'.join(root.split('_')[:-1])
        if tail=='auto':
            src=indir+'/'+f
            dest=outdir+'/cvg_dist/'+f
            shutil.copy2(src, dest)
        if tail!=in_tail: continue
        print name

        inhandle=open(indir+'/'+f, 'r')
        intable=csv.reader(inhandle, delimiter= '\t')
        intable.next()
        for row in intable:
            line=cluster(row)
            image=line.Seq1
            if image_directory.has_key(image)==False:
                image_directory[image]=[]
            out_line=[name, line.ID, line.Seq1, line.Seq2, line.Quad1, line.Quad2, line.feature, line.count, line.j0_x, line.j0_y,line.err_1, line.err_2, line.count,row[-4], row[-3], row[-2], row[-1] ]
            image_directory[image].append(out_line)
            if line.Seq1!=line.Seq2:
                image=line.Seq2
                if image_directory.has_key(image)==False:
                    image_directory[image]=[]
                out_line=[name, line.ID, line.Seq2, line.Seq1, line.Quad2, line.Quad1, line.feature, line.count, line.j0_y, line.j0_x,line.err_1, line.err_2, line.count,row[-4], row[-3], row[-2], row[-1] ]
                image_directory[line.Seq2].append(out_line)
        inhandle.close()
    #Write the output files
    for key in image_directory.keys():
        outfile=outdir+'/'+key+'.tsv'
        outhandle=open(outfile, 'w')
        outtable=csv.writer(outhandle, delimiter='\t')
        for line in image_directory[key]:
            outtable.writerow(line)
        outhandle.close()


#------------------------------------------------------------------------------#
#
#Match junctions to identify TE insertions
#
#------------------------------------------------------------------------------#

def LoadMetaTables(indir, refnames):
    metatable_jxns={}
    file_list=os.listdir(indir)
    for f in sorted( file_list):
        file_parts=f.split('.')
        name,ext= '.'.join(file_parts[:-1]), file_parts[-1]
        if ext!='tsv': continue
        if file_parts[0]=='size_factors': continue
        inpath='{0}/{1}'.format(indir, f)
        print inpath
        inhandle=open(inpath, 'r')

        for textline in inhandle:
            row=textline.split('\t')
            try:

                line=MetaRow(row)
            except:
                print row
                print jabber
            if line.replicate=='2':continue
            if line.seq1!=line.seq2 and refnames.count(line.seq2)==0: continue
            if metatable_jxns.has_key(line.sample)==False:
                metatable_jxns[line.sample]={}
            if metatable_jxns[line.sample].has_key(line.seq1)==False:
                metatable_jxns[line.sample][line.seq1]={}
            if line.feature!='Tandem' and line.feature!='Junction':continue

            if metatable_jxns[line.sample][line.seq1].has_key(line.seq2)==False:
                metatable_jxns[line.sample][line.seq1][line.seq2]={}
                for quad in [('-','+'),('+','-'),('-','-'),('+','+')]:
                    metatable_jxns[line.sample][line.seq1][line.seq2][quad]=MetaJxn()

            quadrant=(line.quad1, line.quad2)
            if metatable_jxns[line.sample][line.seq1][line.seq2].has_key(quadrant)==False:
                metatable_jxns[line.sample][line.seq1][line.seq2][quadrant]=MetaJxn()
            metatable_jxns[line.sample][line.seq1][line.seq2][quadrant].add_jxn(line)
    return metatable_jxns

def FindUnmaskedPosInReference(ref_file):
    ref_seq=GetSeq(ref_file, 'True')
    array_dict={}
    for key in ref_seq:
        array_dict[key]=numpy.fromstring(ref_seq[key], '|S1')
        array_dict[key]=(array_dict[key]!='N')
    return array_dict

def FindInsertionsFromJunctions(tables, ref_file, distance_cutoff=300):
    corresponding_strands={'+': [('-','+'),('+','-')],'-': [('-','-'),('+','+')]}
    #Identify the unmasked positions in the referece
    unmasked_array=FindUnmaskedPosInReference(ref_file)

    insert_dict={}
    for sample in tables.keys():
        insert_dict[sample]={}
        for rpt in tables[sample].keys():
            insert_dict[sample][rpt]={}
            for target in tables[sample][rpt].keys():
                insert_dict[sample][rpt][target]={}
                if rpt==target: continue
                for strand in corresponding_strands.keys():

                    minus_ins=tables[sample][rpt][target][corresponding_strands[strand][0]]
                    minus_pos=minus_ins.Ys

                    plus_ins=tables[sample][rpt][target][corresponding_strands[strand][1]]
                    plus_pos=plus_ins.Ys

                    clusters_l, clusters_r=ClusterJunctions(minus_pos, plus_pos, unmasked_array[target])

                    insertion_list=[]
                    for insertions in clusters_l.keys():
                        insertion_list.append(Insertion())
                        if clusters_l[insertions]!=None:
                            for ind in clusters_l[insertions]:
                                insertion_list[-1].AddLeftJxn(minus_ins.IDs[ind], minus_ins.Ys[ind], minus_ins.Xs[ind], minus_ins.readcount[ind])
                        else:
                            r=12

                        if clusters_r[insertions]!=None:
                            for ind in clusters_r[insertions]:

                                insertion_list[-1].AddRightJxn(plus_ins.IDs[ind], plus_ins.Ys[ind], plus_ins.Xs[ind], plus_ins.readcount[ind])

                        insertion_list[-1].FinalizeInsertion()
                    insert_dict[sample][rpt][target][strand]=insertion_list
    return insert_dict


def ClusterJunctions(minus_pos, plus_pos, unmasked_array, distance_cutoff=300):

    #I want to match junctions on opposite strands of the reference that are near
    #each other.
    #However, the junctions of an insertion present in the reference will be separated by
    #a stretch of Ns, which will make them appear far apart.
    #To account for this I need to count the number of non-N nucleotides between
    #each pair of junctions.
    #I don't know how to do that sort of summation with efficient numpy operations
    #(if it is even possible). So, I first identify all pairs of junctions on opposite strands
    # within 50kb of each other, which requires N by M operations, but can be accomplished in numpy

    # Actually can be accomplished quicker with a Ball-Tree--not worth implementing now
    #however computing pw_dist is memory intensive, a

    #Join the plus and minus strands
    joined_array=numpy.hstack((minus_pos, plus_pos))
    original_array=numpy.array( [-1]*len(minus_pos)+[1]*len(plus_pos))
    original_ind=numpy.array( range(len(minus_pos))+range(len(plus_pos)))
    pw_dist=abs(numpy.subtract.outer(joined_array, joined_array))
    candidates_l, candidates_r= numpy.where(pw_dist<=50000)
    passed_threshold_l, passed_threshold_r=[],[]

    #I then count gthe number of unmasked nucleotides between each identified pair.
    #This is linear in the number of pairs with respect to non-numpy operations.

    for ind in range(len(candidates_l)):
        l,r=candidates_l[ind], candidates_r[ind]
        left_pos=int( min(joined_array[ l],joined_array[ r]))
        right_pos=int( max(joined_array[ l],joined_array[ r]))
        if unmasked_array[left_pos:right_pos].sum()<distance_cutoff:
            passed_threshold_l.append(l)
            passed_threshold_r.append(r)

    array_2=passed_threshold_r
    array_1=passed_threshold_l
    r_list=array_2
    l_list=array_1
    cluster_list=[]

    for i in range(len( r_list)):
        pair=set([r_list[i], l_list[i]])
        intersection=[]
        #Check whether either item in the pair already belongs to a cluster
        for c in range(len( cluster_list)):
            if len( pair&cluster_list[ c])!=0: intersection.append(c)

        if (len(intersection))!=0:

            new_cluster_list=[]
            new_cluster=pair
            for c in range(len(cluster_list)):
                if intersection.count(c)!=0:
                    new_cluster|=cluster_list[c]
                else:
                    new_cluster_list.append(cluster_list[c])


            new_cluster_list.append(new_cluster)
            cluster_list=new_cluster_list
        else:
            cluster_list.append(pair)

    #Separate out clusters into a paired lists
    final_clusters_l={}
    final_clusters_r={}
    for i in range(len(r_list)):

        for c in range(len( cluster_list)):
            if final_clusters_l.has_key(c)==False:
                final_clusters_l[c]=set()
                final_clusters_r[c]=set()

            if list(cluster_list[c]).count( r_list[i])==0: continue

            or_ind=original_ind[r_list[i]]
            or_arr=original_array[r_list[i]]
            if or_arr==1:
                final_clusters_r[c].add(or_ind)
            if or_arr==-1:
                final_clusters_l[c].add(or_ind)

    return final_clusters_l, final_clusters_r


def WriteInsertionTable(insertion_table, outfile):
    outhandle=open(outfile, 'w')
    outtable=csv.writer(outhandle, delimiter='\t')
    rpt_dict={}
    sample_list=['Repeat']

    for s in sorted( insertion_table.keys()):

        sample_list.append(s)

        for rpt in sorted(insertion_table[s].keys()):

            for target in insertion_table[s][rpt].keys():
                if target==rpt: continue

                for strand in ['-', '+']:
                    for insertion in insertion_table[s][rpt][target][strand]:
                        row=[s,rpt, target,strand]+insertion.rows(counts=True)[0]
                        outtable.writerow(row)

    outhandle.close()


#------------------------------------------------------------------------------
#
#Estimate the library-specific error
#
#------------------------------------------------------------------------------

def BuildErrorDict(infile, ref_file, cons_file):
    """This reads a """
    inhandle=open(infile, 'r')
    rpt_lengths= GetLengths(cons_file)
    ref_seq=FindUnmaskedPosInReference(ref_file)
    intable=csv.reader(inhandle, delimiter='\t')
    error_dict={}
    data_dict={}
    for row in intable:
        sample=row[0]
        if error_dict.has_key(sample)==False:
            error_dict[sample]={}
            data_dict[sample]=[]

        repeat=row[1]
        if rpt_lengths.has_key(repeat)==False: continue
        if error_dict[sample].has_key(repeat)==False:
            error_dict[sample][repeat]={'-':[], '+':[]}
        strand=row[3]
        chrom=row[2]
        count_1, count_2=float(row[10]), float(row[11])
        if count_1*count_2==0: continue
        x_1, x_2=float(row[6]), float(row[8])
        int_x1=int(x_1)
        int_x2=int(x_2)
        if ref_seq[chrom][int_x1-500:int_x1+500].sum()!=1000 or ref_seq[chrom][int_x2-500:int_x2+500].sum()!=1000: continue

        y_1, y_2=float(row[7]), float(row[9])
        X=x_1-x_2
        Y=y_1-y_2+rpt_lengths[repeat]

        #Change basis
        if strand=='+':
            new_X=X/2**.5-Y/2**.5
            new_Y=X/2**.5+Y/2**.5
        else:
            new_X=X/2**.5+Y/2**.5
            new_Y=X/2**.5-Y/2**.5

        error_dict[sample][repeat][strand].append([new_X,new_Y,count_1, count_2])
        if repeat.split('_')[-1]=='LTR':data_dict[sample].append([new_X,new_Y,count_1, count_2])
        elif repeat.find('LTR')!=-1:data_dict[sample].append([new_X,new_Y,count_1, count_2])

    for key in data_dict.keys():
        data_dict[key]=numpy.array(data_dict[key])
    return error_dict, data_dict

def RecomposeMatrix(eig_vec, eig_val_matrix):
    """Returns a covariance matrix from eigenvectors and eigenvalues"""
    return eig_vec*eig_val_matrix*numpy.linalg.inv(eig_vec)

def GetLibraryCovariances(infile):
    summary_dict={}
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t')
    #skip the first two rows; these are headers
    row=intable.next()
    row=intable.next()
    eig_vec=numpy.matrix([[ 0.70710678,  -0.70710678], [0.70710678,  0.70710678]])
    for row in intable:
        sample_name, sample_rep=row[0].split('_')
        if sample_rep=='2': continue
        scale_x, scale_y=float( row[7]), float(row[8])
        cov_opt=RecomposeMatrix(eig_vec, [[scale_x**2,0],[0,scale_y**2]])
        summary_dict[sample_name]=numpy.linalg.inv( cov_opt)
    return summary_dict

def EstimateCovariances( data_dict):
    cov_dict={}
    for key in sorted(data_dict.keys()):
        eig_1, sd_1=EstimateParam(data_dict[key][:,2:], data_dict[key][:,1]**2)
        eig_2,sd_2= EstimateParam(data_dict[key][:,2:], data_dict[key][:,0]**2)
        cov_dict[key]=[eig_1,eig_2, eig_1+sd_1, eig_2+sd_2]
    return cov_dict

def EstimateParam(x,y):
    """Fit the 2D power function to the data to find the best F: x-->y ."""
    x=numpy.array(x)
    y=numpy.array(y)
    k, cov = scipy.optimize.curve_fit(Power2D, x,y)

    return k**.5, cov**.25

def Power2D(x,sig):
    """"""
    return (sig/x).sum(1)

def PowerLaw(x,sig):
    return sig/x


#------------------------------------------------------------------------------#
#
#Clustering algorithms
#
#------------------------------------------------------------------------------#
def ClusterDirectory(indir,outdir,insertion_file, ref_file, rep_file ):
    err_dict, data_dict=BuildErrorDict(insertion_file, ref_file, rep_file)
    cov_dict=EstimateCovariances(data_dict)
    infile_list=os.listdir(indir)
    MakeDir(outdir)
    rep_len=GetLengths(rep_file)
    for infile in sorted( infile_list):
        if infile=='size_factors.tsv': continue
        ext=infile.split('.')[-1]
        name='.'.join( infile.split('.')[:-1])

        if ext!='tsv': continue
        if rep_len.has_key(name)==False: continue
        ClusterFile('{0}/{1}'.format(indir, infile), '{0}/{1}'.format(outdir, infile), name, cov_dict, rep_len[name]>1000)

def ClusterFile(infile, outfile, self_seq,  cov_dict, exclude_consensus=True ):

    outhandle=open(outfile, 'w')
    outtable=csv.writer(outhandle, delimiter='\t')
    full_file='.'.join( outfile.split('.')[:-1])+'_full.tsv'
    fullhandle=open(full_file, 'w')
    fulltable=csv.writer(fullhandle, delimiter='\t')
    quadrants=[('-','-'), ('-','+'),('+','-'), ('+','+')]
    #Identify target sequences

    data, errors, samples, ID, CN=ReadMetaTableAsDictionary(infile,cov_dict,emp= True)
    for seq in data.keys():
        for quad in data[seq].keys():
            #Read file

            if len(data[seq][quad])<2: continue
            #cluster data
            print "Clustering {0} {1}".format(seq, quad)

            cluster_output=CMeans(data[seq][quad], errors[seq][quad],'', 100, 200)
            cluster_assignments=cluster_output[1]
            #Organize the CMeans output
            jxn_tables, jxn_array, cluster_labels=BuildClusterTable(cluster_assignments, samples[seq][quad],data[seq][quad], ID[seq][quad], CN[seq][quad], cov_dict, exclude_consensus)

            #Fix mistakes
            jxn_tables, cluster_assignments, cluster_labels=CorrectErrors(jxn_tables, data[seq][quad], cluster_assignments, cluster_labels)


            for i in range(len(cluster_labels)):
                indices=cluster_assignments==cluster_labels[i]
                label=cluster_labels[i]
                pulled_data=data[seq][quad][indices, :]
##                try:
                if len(pulled_data)>1:
                    X,Y=numpy.mean(pulled_data, 0)
                else: X,Y=pulled_data[0]
##                except:
##                    print pulled_data
##                    print jabber
                try:
                    copy_number=[str( x.CN) for x in jxn_tables[label,:]]
                except:
                    print jxn_tables.shape
                    print label
                    print cluster_labels
                    print jabber
                pos_x=[','.join( x.pos_x) for x in jxn_tables[label,:]]
                pos_y=[','.join( x.pos_y) for x in jxn_tables[label,:]]
                ids=[','.join( x.IDs) for x in jxn_tables[label,:]]
                output_line=[';'.join((copy_number[index], pos_x[index], pos_y[index], ids[index])) for index in range (len(copy_number))]
                CN_array=numpy.asarray(copy_number, float)
                #Classify structure
                feature='Junction'
                if self_seq==seq:
                    if quad[0]!=quad[1]:
                        if X<Y: feature='Tandem'
                        else: feature='Deletion'
                    else:
                        feature='Inversion'
                additional_info=['','','']
                row=[self_seq, seq, quad[0], quad[1], X,Y,numpy.mean(CN_array>0), numpy.mean(CN_array), feature]+ additional_info +list(copy_number)
                outtable.writerow(row)
                row=[self_seq,seq, quad[0], quad[1], X,Y,numpy.mean(CN_array>0), numpy.mean(CN_array), feature]+ additional_info +output_line
                fulltable.writerow(row)
            outhandle.flush()
            fullhandle.flush()
    outhandle.close()
    fullhandle.close()

def ReadMetaTableAsDictionary (infile, err_dict, emp=True, exclude_cons=True):
    """Read a metatable and output list of junctions and corresponding uncertainty measures"""
    inhandle=open(infile, 'r')
    intable=csv.reader(inhandle, delimiter='\t',quoting=csv.QUOTE_NONE)
    row=intable.next()
    pos_list={}
    err_list={}
    sample_list={}
    id_list={}
    read_count_list={}

    for row in intable:
        line=MetaRow(row)
        quad=(line.quad1, line.quad2)

        if pos_list.has_key(line.seq2)==False:
            pos_list[line.seq2]={}
            err_list[line.seq2]={}
            sample_list[line.seq2]={}
            id_list[line.seq2]={}
            read_count_list[line.seq2]={}

        if pos_list[line.seq2].has_key(quad)==False:
            pos_list[line.seq2][quad]=[]
            err_list[line.seq2][quad]=[]
            sample_list[line.seq2][quad]=[]
            id_list[line.seq2][quad]=[]
            read_count_list[line.seq2][quad]=[]

        if numpy.isnan(line.x+line.y+line.err1+line.err2)==True: continue
        if line.feature=='Probable Artifact': continue

        read_count=min(10000, line.readcount)
        sample=line.sample

        if exclude_cons==True and line.seq1==line.seq2 and quad[0]!=quad[1] and (line.x-line.y>-400 and line.x -line.y<100) : continue
        if err_dict.has_key(sample)==False: continue

        if emp==True:
            err1=err_dict[sample][2] /read_count**.5
            err2=err_dict[sample][3]/read_count**.5

        else:
            err1=line.err1
            err2=line.err2

        if line.quad1==line.quad2:
                orientation=[[1,1],[1,1]]
        else:
                orientation=[[1,-1],[-1,1]]

        eig_vec=numpy.matrix([[ 0.70710678,  0.70710678], [-0.70710678,  0.70710678]])
        cov_mat=numpy.array( orientation* numpy.array( RecomposeMatrix(eig_vec, [[err1**2,0],[0,err2**2]])), float)
        err_list[line.seq2][quad].append(cov_mat)
        pos_list[line.seq2][quad].append([line.x, line.y])
        sample_list[line.seq2][quad].append(line.sample)
        id_list[line.seq2][quad].append(line.ID)
        read_count_list[line.seq2][quad].append(line.CN_map)

    for seq_key in pos_list.keys():
        for quad_key in pos_list[seq_key].keys():
            pos_list[seq_key][quad_key]=numpy.array(pos_list[seq_key][quad_key])
            err_list[seq_key][quad_key]=numpy.array( err_list[seq_key][quad_key])
            sample_list[seq_key][quad_key]=numpy.array(sample_list[seq_key][quad_key])
            id_list[seq_key][quad_key]=numpy.array(id_list[seq_key][quad_key])
            read_count_list[seq_key][quad_key]=numpy.array(read_count_list[seq_key][quad_key])

    return pos_list,err_list, sample_list, id_list, read_count_list

def BuildClusterTable(labels,samples,data, ids, CNs, cov_dict,exclude_consensus):
    """Organizes the C-means output in table."""
    sample_set=numpy.array( sorted(cov_dict.keys()))
    cluster_labels=numpy.array(list(set(list(labels))))
    count=0
    jxn_array=numpy.ndarray((len(cluster_labels), len(sample_set)))
    jxn_array_rich=numpy.ndarray((len(cluster_labels), len(sample_set)), object)
    jxn_array.fill(0)
    for i in range(len(cluster_labels)):
        for j in range( len(sample_set)):
            jxn_array_rich[i,j]=ClusterInformation()
    for l in cluster_labels:
##        print l
        cluster_ind=labels==l
        samples_with_jxn=samples[cluster_ind]
        jxn_CN=CNs[cluster_ind]
        jxn_data=data[cluster_ind,:]
        jxn_ids=ids[cluster_ind]
##        print numpy.mean(jxn_CN)
        for i in range(len(samples_with_jxn)):

            samp_ind=numpy.where(sample_set==samples_with_jxn[i])[0]
            for s in samp_ind:
##            jxn_array_rich[count,samp_ind]= [0,[],[], []]
                jxn_array[count,s]+=jxn_CN[i]
                jxn_array_rich[count, s].add_jxn(jxn_CN[i], jxn_data[i,0], jxn_data[i,1], str(jxn_ids[i]))

        count+=1
    return jxn_array_rich, jxn_array, cluster_labels

def CorrectErrors( data,data_locations, cluster_assignments, cluster_labels, threshold=100,cutoff=.01, FDR=.01):
    p_array=numpy.ndarray((len(data), len(data)))
    p_array.fill(1.)
    p_val_list=[]
    locations=[]
    loc_N=[]
    for i in range(len(data)):
        #Compute the average location of junctions in the cluster
        indices=cluster_assignments==cluster_labels[i]
        pulled_data=data_locations[indices, :]

        if len(pulled_data)>1:
            X,Y=numpy.mean(pulled_data, 0)
        else:

            X,Y=pulled_data[0]
        loc_N.append( len(pulled_data))
        locations.append(numpy.array( [X,Y]))

        for j in range(i):
            L2=((locations[i]-locations[j])**2).sum()**.5
            if L2>threshold: continue
            CN_i=numpy.array([d.CN for d in data[i,:]])
            CN_j=numpy.array([d.CN for d in data[j,:]])
            #Compute the p-value of association
            p=FisherExact(CN_i, CN_j)[1]
            p_array[i,j]=p
            p_val_list.append(p)

    #Hierarchical cluster-
    terminate=False
    count=0
    clusters,samples=data.shape

    #If an FDR corrected clustering is desired, the complete dendrogram could
    #be constructed here.
    #The complete dendrogram of the data contains all relevant p-values
    #which can be used to correct for multiple testing. That is: the dendrogram
    #reflects all decisions we might make. So first we cluster to completion
    #on a copy of the data and log the p-values of the Fisher exact tests

    data_copy=numpy.copy(data)
    p_copy=numpy.copy(p_array)
    p_vals=[]
    #Now construct the dendrogram again, but terminate once the p-value exceeds the
    #greatest p-value for which the null hypothesis should be rejected under
    #multiple test correction
    while terminate==False:
        #When we join two clusters, we add the second to the first, and set
        #the index flag of the second to false

        min_p=numpy.min(p_array)
        if min_p>cutoff or count>=len(data)-1:
            terminate==True
            break
        min_indices=numpy.where(p_array==min_p)
        print min_indices
        i,j=min_indices[0][0], min_indices[1][0]  #Indices with the minimum
        if i==j: break
        #Merge the clusters
        for s in range(samples):
                data[i,s].join_jxns( data[j,s])
                data[j,s]=[]
        #Update the cluster labels
        j_cluster_indices=cluster_assignments==cluster_labels[j]
        cluster_assignments[j_cluster_indices]=cluster_labels[i]
        cluster_labels[j]=cluster_labels[i]

        #Update locations
        locations[i]=(locations[i]*loc_N[i]+locations[j]*loc_N[j])/(loc_N[i]+loc_N[j])

        #Update the probailities
        for k in range(len(data)):
            if k==i: continue
            #If the cluster has been removed
            if data[k,0]==[]:
                p_array[k,:]=1.
                p_array[:,k]=1.
                continue
            if k<i:
                if p_array[i,k]==1:  continue
                #Only update if the average locations within the two clusters
                L2=((locations[i]-locations[k])**2).sum()**.5
                if L2>threshold:
                    p_array[i,k]=1
                    continue
                CN_i=numpy.array([d.CN for d in data[i,:]])
                CN_k=numpy.array([d.CN for d in data[k,:]])
                #Compute the p-value of association
                p=FisherExact(CN_i, CN_k)[1]
                p_array[i,k]=p
            if k>i:
                if p_array[k,i]==1: continue
                L2=((locations[i]-locations[k])**2).sum()**.5
                if L2>threshold:
                    p_array[k,i]=1
                    continue
                CN_i=numpy.array([d.CN for d in data[i,:]])
                CN_k=numpy.array([d.CN for d in data[k,:]])
                #Compute the p-value of association
                p=FisherExact(CN_i, CN_k)[1]
                p_array[k,i]=p
        count+=1
    print "Joined {0} cluster.".format(count)
    cluster_labels=list(set(list(cluster_assignments)))
    return data, cluster_assignments, numpy.array( cluster_labels)


def FisherExact(vec1, vec2):

    pres_1=vec1>0
    pres_2=vec2>0
    A_B=((pres_1==1)*(pres_2==1)).sum()
    A_notB=((pres_1==1)*(pres_2!=1)).sum()
    notA_B=((pres_1!=1)*(pres_2==1)).sum()
    notA_notB=((pres_1!=1)*(pres_2!=1)).sum()
    table=numpy.array( [[A_B, A_notB], [notA_B,notA_notB]])
    binom_p=scipy.stats.fisher_exact(table, 'less')
    ##    print expected_intersection
    ##    print observed_intersection

    return binom_p

#------------------------------------------------------------------------------#
#
#C-Means algorithm
#
#------------------------------------------------------------------------------#

def CMeans(data, errors,pdf=scipy.stats.norm, components=10, iterations=10):
    total_start=time.clock()
    time_dict={'Init':0,'Updates':0,'Total':0., 'Prune':0. }

    sample_size, dim=data.shape
    weight=numpy.array([1./components]*components)

    start=time.clock()
    means=InitializeModel(100, data,errors,1)
    time_dict['Init']+=time.clock()-start
    print len(means)

    components=len(data)
    weight=numpy.array([1./components]*components)

    weight/=sum(weight)
    mean_tracker=[]
    x_list=[]
    y_list=[]
    weight_tracker=[]
    it_tracker=[]
    theta_old=zip(weight, means)
    lk_list=[]
    dif_list=[]
    loss_list=[]
    lab_list=[]
    sd1=[]
    sd2=[]
    for e in errors:
        dec_vec, dec_val=DecomposeMatrix(e)
        err1,err2=dec_val[0][0]**.5, dec_val[1][1]**.5
        sd1.append(err1)
        sd2.append(err2)
    sd1=numpy.array(sd1)
    sd2=numpy.array(sd2)

    for i in range(iterations):

        wt, mean=zip(*theta_old)
        wt=numpy.array(wt)

        mean=numpy.array(mean)
        old_mean=mean

        x_list.append(mean[:,0])
        y_list.append(mean[:,1])

        #Update the model parameters
        start=time.clock()
        theta_old,theta_array, weight_array=UpdateModelParameters(data, errors, theta_old)
        time_dict["Updates"]+=time.clock()-start

        new_wt, new_mean=zip(*theta_old)
        new_wt=numpy.array(new_wt)
        new_mean=numpy.array(new_mean)

        labels=numpy.argmax(weight_array, 0)
        loss=(((data-new_mean[labels,:])**2).sum(1)**.5).sum()
        loss_list.append(loss)

        if i>0:
            label_change=numpy.mean( (labels==old_labels))
            lab_list.append(label_change)

        old_labels=labels
        if numpy.isnan(weight_array.sum())==True:

            print i
            weight_array=CalcBivariateNormal(data, numpy.array( mean), errors )
            weight_array=weight_array*numpy.array(wt)[:,None]

            weight_array/=weight_array.sum(0)
            labels=numpy.argmax( weight_array,0)
            print jabber
            return theta_old,labels, lk_list, weight_array, mean, numpy.array(weight_tracker), numpy.array( x_list),numpy.array( y_list),dif_list, loss_list, lab_list #,numpy.array( weight_tracker), mean_tracker, it_tracker

        #Prune the model to remove unnecessary clusters
        prune_start=time.clock()
        prune_list=[]

        for I in range(len(theta_old)):
            wt, mean=theta_old[I]
            if numpy.isnan(wt)+numpy.isnan(mean[0])+numpy.isnan(mean[1])==True:
                prune_list.append(I)
        for I in reversed(prune_list):
            theta_old.pop(I)

        wt, mean=zip(*theta_old)
        wt=numpy.array(wt)
        mean=numpy.array(mean)

        weight_array=CalcBivariateNormal(data, numpy.array( mean), errors )
        weight_array=weight_array*numpy.array(wt)[:,None]

        weight_array/=weight_array.sum(0)
        labels=numpy.argmax( weight_array,0)
        label_set= set(list(labels))

        #Check each component to see whether it explains at least one data point
        #If not, add it to the list of components to be removed
        for I in range(len(theta_old)):
            if I not in label_set:
                prune_list.append(I)
        #Beginning at the end of the prune list, remove unnecessary components
        #from the model
        for I in reversed(prune_list):
            theta_old.pop(I)

        wt, mean=zip(*theta_old)
        mean=numpy.array(mean)

        time_dict['Prune']+=time.clock()-prune_start

        #Check the termination criteria
        if len( lab_list)>10:
            if numpy.min( lab_list[-10:])>.999:
                print "Terminated after {0} iterations.".format( i)
                break


    wt, mean=zip(*theta_old)
    mean=numpy.array(mean)
    weight_array=CalcBivariateNormal(data, numpy.array( mean), errors )
    weight_array/=weight_array.sum(0)
    labels=numpy.argmax( weight_array,0)
    time_dict['Total']+=time.clock()-total_start
    print time_dict
    return theta_old,labels, lk_list, weight_array, mean, numpy.array(weight_tracker), numpy.array( x_list),numpy.array( y_list),dif_list, loss_list, lab_list #,numpy.array( weight_tracker), mean_tracker, it_tracker

def InitializeModel(components,data, errors, buffer_size=1 ):
    """Construct an initial model by centering components on randomly chosen
    data points, but with the additional requirement that no two components be within
    a set distance of each other. Place components until no eligible data points remain."""
    init_theta=[]
    dataLength=data.shape[0]

    random_seed=int(time.clock()*(10000))%(2**24)
    print "Random Seed: {0}".format(random_seed)
    numpy.random.seed(random_seed)

    prec=numpy.array( [numpy.linalg.inv(errors[i,:,:]) for i in range( errors.shape[0])])

    #Choose data points as the initial means
##    for i in range(components):
    while range(0, dataLength)!=[]:
        index=numpy.random.choice(range(0, dataLength))
        init_theta.append(data[index])

        pulledRead= data[index]
        distance_from_seed= MahalanobisDistance(pulledRead, data, prec)
        RemainingReads=numpy.where(distance_from_seed>buffer_size)

        data=data[RemainingReads]

        dataLength=data.shape[0]
##        print dataLength
##        if range(0, dataLength)==[]: break
    return init_theta

def UpdateModelParameters(data, errors, theta):
        #Unpack parameters
    wt, mean=zip(*theta)
    wt=numpy.array(wt)
    mean=numpy.array(mean)

    #Compute the fuzzy cluster assignments
    weight_array=CalcBivariateNormal(data, numpy.array( mean), errors )
    weight_array=weight_array*numpy.array(wt)[:,None]
    weight_array/=weight_array.sum(0)

    #Update the weights and means based on the cluster assignments
    new_weights=UpdateMixingProportions(weight_array)
##    new_weights=wt
    new_means=UpdateMeans(weight_array, data, errors, theta)

    #Repackage parameters
    theta=zip(new_weights, new_means)
    theta_array=numpy.hstack((new_weights[:,None],new_means))
    return theta, theta_array, weight_array

def CalcBivariateNormal(data,  means, cov ):

    """Compute the probability of each of N data points given each of the K components.
    Returns a NxK matrix of probabilities."""

    #Unpack the model parameters
    #Structure of the covariance matrix allows for some simplifications
    rho= cov[:,0,1]/( cov[:,0,0])
    sigma_x=abs(( cov[:,0,0]))

    alpha=(2*numpy.math.pi*sigma_x  *(1-(rho)**2)**.5)**-1  #Constant depending only on covariance
    beta=-1*(2*(1-(rho)**2)*sigma_x)**-1                 #Constant depending only on covariance

    #Compute distances from the mean
    x_dist = (means[:,0][:,None]- data[:,0]) #(x-E[X]) for all x in X and all k in K
    y_dist= (means[:,1][:,None]- data[:,1])  #(y-E[Y]) for all y in Y and all k in K

    # F(x,y)= alpha* exp( beta/sigma^2 *[(x-mean)^2 + (y-mean)^2 - 2 * correl * (x-mean)(y-mean)])
    #Compute the probability oft the bivariate normal
    PDF= numpy.exp (beta *( (x_dist**2) + (y_dist**2) - 2*rho* x_dist*y_dist))
    return  (PDF *alpha)

def UpdateMixingProportions(weight_array):
    N=weight_array.shape[1]
    weights=weight_array.sum(1)/N
    return weights

def UpdateMeans(weight_array, data, errors, theta):
    N=weight_array.shape[1]
    N_k=weight_array.sum(1)
    error_weights=errors/errors.sum()
    mean_x=(weight_array*( data[:,0]) ).sum(1)/(N_k)
    mean_y=(weight_array*( data[:,1]) ).sum(1)/(N_k)
    return numpy.vstack((mean_x, mean_y)).transpose()

def MahalanobisDistance(point, data, prec):
    mahab=[]
    for i in range(len(data)):
        mahab.append(scipy.spatial.distance.mahalanobis(point, data[i], prec[1]))
    return numpy.array( mahab)

def ClusterMetatables(indir, outdir, specification):

    MakeDir(outdir)

    spec_dict=ReadSpecificationFile(specification)
    metatable_dir=outdir+'/metatables'
    print "Building metatables."
    BuildMetaTable(indir,'outCN',metatable_dir)
    metatables=LoadMetaTables(metatable_dir, spec_dict['CV'])
    print ("Calling TE insertions in the reference chromosomes")
    insertion_file='{0}/TE_insertions.tsv'.format(outdir)
    insertions=FindInsertionsFromJunctions(metatables, spec_dict['Masked'])
    WriteInsertionTable(insertions, insertion_file)
##    BuildErrorDict(insertion_file, spec_dict['Masked'], spec_dict['Cons']  )
    cluster_dir=outdir+'/cluster_tables'
    ClusterDirectory(metatable_dir, cluster_dir, insertion_file, spec_dict['Masked'], spec_dict['Cons'] )


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
    indir=param['-i']
    outdir=param['-o']
    spec_file=param['-spec']
    ClusterMetatables(indir, outdir, spec_file)

    pass

if __name__ == '__main__':
    main(sys.argv)
