#!/usr/bin/env python
#a program to convert a bed file from one prokaryotic species to another based on a Mauve alignment
#Alexa McIntyre, 2017

from collections import defaultdict
import sys
import os
import re

def check_other_seq(seq_a,seq_b,int_b,factor_b):
    """updates the length of the gap depending on where the sequences next align without gaps"""
    if seq_a[0] == '-':
        max_gap = len(re.split('A|C|G|T',seq_a)[0])
        seq_a = seq_a[max_gap:]
        seq_b = seq_b[max_gap:]
        int_b += factor_b*max_gap
    return(seq_a,seq_b,int_b)

def min_recursive(first_seq,second_seq,first_ind,second_ind,pos_dict,factor1,factor2,factor_match,i=0):
    """adds indices of ungapped alignments to dictionary"""
    if len(first_seq) == 0 or len(second_seq) == 0:
        return(pos_dict)
    first_seq,second_seq,second_ind = check_other_seq(first_seq,second_seq,second_ind,factor2)
    second_seq,first_seq,first_ind = check_other_seq(second_seq,first_seq,first_ind,factor1)
    while first_seq[0] == '-' or second_seq[0] == '-':
        if first_seq[0] == '-':
            first_seq,second_seq,second_ind = check_other_seq(first_seq,second_seq,second_ind,factor2)
        if second_seq[0] == '-':
            second_seq,first_seq,first_ind = check_other_seq(second_seq,first_seq,first_ind,factor1)
    min_seq = min([len(first_seq.split('-')[0]),len(second_seq.split('-')[0])])
    
    #print first_seq[:10],second_seq[:10],first_ind,second_ind,min_seq
    
    range1 = [first_ind,first_ind+factor1*min_seq]
    range2 = [second_ind,second_ind+factor2*min_seq]
    seq1_inds = zip(first_seq[:min_seq],range(min(range1),max(range1)+1)[::factor1])
    seq2_inds = zip(second_seq[:min_seq],range(min(range2),max(range2)+1)[::factor2])

    pos_dict.update({x:(y,factor_match) for x,y in zip(zip(*seq1_inds)[1],zip(*seq2_inds)[1])})

    first_ind += factor1*min_seq
    second_ind += factor2*min_seq
    first_seq = first_seq[min_seq:]
    second_seq = second_seq[min_seq:]

    if len(first_seq) and len(second_seq) > 0:
        i+=1
        return(min_recursive(first_seq,second_seq,first_ind,second_ind,pos_dict,factor1,factor2,factor_match,i))
    else:
        #print i
        return(pos_dict)

def bed2bed(pos_dict,bedfi,outname):
    """converts original bed indices to aligned genome indices and writes to output file"""
    with open(bedfi,'r') as bed:
        with open(outname,'w') as outfi:
            for line in bed:
                cset = line.strip().split('\t')[:6]
                start = int(cset[1])
                end = int(cset[2])
                try:
                    strand = cset[5]
                except IndexError:
                    strand = '+'
                if start in pos_dict and end in pos_dict:
                    new_start = min(pos_dict[start][0],pos_dict[end][0])
                    new_end = max(pos_dict[start][0],pos_dict[end][0])
                    outfi.write('\t'.join([str(x) for x in [2,new_start,new_end,pos_dict[end][1],'.',strand,1,start,end]])+'\n')

def change_inds(bedfi,xmfafi,outname): #,partial_overlap,save):
    """main function - reads in xmfa and bed"""
    pos_dict = {}
    with open(xmfafi,'r') as xmfa:
        new_alignment = False
        first_seq = ''
        first_inds = [0,0]
        second_seq = ''
        second_inds = [0,0]
        for entry in xmfa.read().split('>')[1:]:
            if not new_alignment and entry[-2] != '=':
                new_alignment = True
                first_seq = ''.join(entry.strip().split('\n')[1:])
                first_inds = [int(x) for x in entry.split(':')[1].split()[0].split('-')]
                if entry.split()[1] == '-':
                    factor1 = -1
                else:
                    factor1 = 1
                #> 1:1-508662 + /Users/abm237/Desktop/enterobacter_sequences/IF2SWP3_hybrid_assembly.fasta
            elif new_alignment:
                if entry[-2] == '=':
                    second_seq = ''.join(entry.strip().split('\n')[1:-2])
                    if entry.split()[1] == '-':
                        factor2 = -1
                    else:
                        factor2 = 1
                    second_inds = [int(x) for x in entry.split(':')[1].split()[0].split('-')]
                    if factor1 < 0:
                        first_ind = first_inds[1]
                    else:
                        first_ind = first_inds[0]
                    if factor2 < 0:
                        second_ind = second_inds[1]
                    else:
                        second_ind = second_inds[0]
                    if factor1 == factor2:
                        factor_match = 1
                    else:
                        factor_match = 0
                    pos_dict = min_recursive(first_seq,second_seq,first_ind,second_ind,pos_dict,factor1,factor2,factor_match)
                    second_seq, second_inds = '',[]
                    first_seq, first_inds = '',[]
                    new_alignment = False
                else:
                    first_seq = ''.join(entry.strip().split('\n')[1:])
                    first_inds = [int(x) for x in entry.split(':')[1].split()[0].split('-')]
                    new_alignment = True
    #print len(pos_dict)
    bed2bed(pos_dict,bedfi,outname)

def main():
    """parses command line options"""
    #TODO: fix off-by-one errors in output! from conversion stage?
    #TODO: allow multiple contigs 
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Liftover bed positions')
    parser.add_argument('-x','--alignment',type=str,required=True,help='alignment file (xmfa format) from Mauve GUI for two species')
    parser.add_argument('-b','--bed',type=str,required=True,help='bed file with positions to exchange - if no strand column provided, chooses + by default')
    parser.add_argument('-o','--output',type=str,required=False,help='output file name')
    #parser.add_argument('-s','--save',action='store_true',required=False,help='TODO: save bed file of input positions with exact matches')
    #parser.add_argument('-p','--partial_overlap',action='store_true',required=False,help='TODO: allow partial overlaps, report separated')
    parser.add_argument('-v','--version',action='store_true',required=False,help='print version')
    args = parser.parse_args()

    if args.version:
        print 'proklift 0.2'
        sys.exit(0)

    if not args.output:
        output='liftover.bed'
    else:
        output=args.output

    change_inds(args.bed,args.alignment,output) #,args.partial_overlap,args.save)

if __name__ == "__main__":
    main()


