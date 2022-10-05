#!/usr/bin/env python3

#####################################################################################################
#Script Name    : ms_seqcon.py
#Description    : performs some useful MSA manipulations
#Author         : Michael A Sennett
#####################################################################################################

from Bio import SeqIO
from Bio import AlignIO
import numpy as np
import re
import argparse
import os

#################################################################################################
#
#different functions
#
#################################################################################################

#take an alignment and explode all it into individual sequence files
def msa_explode(direct, file, form):
    file = file[0]
    form = form[0]
    
    for record in SeqIO.parse(os.path.join(direct,file), form):
        seq_lst=list(record.seq)
        seq_text=[]
        for m in range(0, len(seq_lst), 70):
            n=m+70
            seq_text.append("".join(seq_lst[m:n])+"\n")
        newfile = str(record.id)+'.fst'
        with open(os.path.join(direct,newfile), 'w') as file:
            print(newfile)
            file.write(">"+str(record.id)+"\n")
            for text in seq_text:
                file.write(str(text))
            file.close()

#calculate the number of taxa and the length of the alignment
def align_dim(direct, file, form):
    file = file[0]
    form = form[0]

    c=0
    alignment = AlignIO.read(os.path.join(direct,file), form)
    for rec in alignment:
        c+=1
    print("taxa: ", c, "length: ", alignment.get_alignment_length())

#convert one alignment type into another
def change_format(direct, file, form, newform):
    file = file[0]
    form = form[0]
    newform = newform[0]

    alignment = AlignIO.read(os.path.join(direct,file),form)
    newfile=re.sub(r'([.].*)','.'+str(newform), str(file))
    AlignIO.write(alignment, os.path.join(direct,newfile), newform)

    print(newfile)

#any column that is all gaps, delete
def del_all_gaps(direct, file, form):
    
    file = file[0]
    form = form[0]

    alignment=AlignIO.read(os.path.join(direct,file),form)

    emp_arr=[]
    names=[]
    #record sequence    
    for rec in alignment:
        emp_arr.append(rec.seq)
        names.append(rec.id)
    seqs=np.asarray(emp_arr)
    seqs=np.transpose(seqs)
    columns=[]
    #determine columns that contain only gaps
    for i in range(len(seqs)):
        c=0
        l=len(seqs[i])
        for j in range(len(seqs[i])):
            if seqs[i][j] == '-':
                c+=1
        if c >= l:
            columns.append(i)
    
    seqs=np.transpose(seqs)
    #delete columns of only gaps
    seqs=np.delete(seqs,columns,1)
    new_app='_del.'+str(form)
    newfile=re.sub(r'([.].*)', new_app, str(file))

    #write new alignment file
    with open(os.path.join(direct,newfile), 'w+') as aln:
        for k in range(len(seqs)):
            aln.writelines(['>'+names[k]+'\n', "".join(list(seqs[k]))+'\n'])
        aln.close()

    print(newfile)

#applies gaps from MSA1 to MSA2, have to be in same order
def apply_gaps(direct, file, form):
    
    file1 = file[0]
    file2 = file[1]
    
    if len(form) >= 2:
        form1 = form[0]
        form2 = form[1]
    else:
        form1 = form[0]
        form2 = form[0]

    alignment1=AlignIO.read(os.path.join(direct, file1), form1)
    alignment2=AlignIO.read(os.path.join(direct, file2), form2)

    names=[]
    aln2=[]

    for i in range(len(alignment1)):
        tmp=str(alignment2[i].seq)
        for j in range(len(alignment1[i].seq)):
            if alignment1[i].seq[j] == '-' and tmp[j] != '-':
                print(alignment2[i].id, 'replace', tmp[j], 'with -')
                tmp=tmp.replace(tmp[j], '-')
            else:
                pass

        names.append(alignment2[i].id)
        aln2.append(tmp)
        
    newfile=re.sub(r'([.].*)','_gapped.'+str(form2), str(file2))

    with open(os.path.join(direct, newfile), 'w+') as msa:
        for k in range(len(names)):
            msa.writelines(['>'+names[k]+'\n', aln2[k]+'\n'])
        msa.close()

    print(newfile)

#get the number of identical residues for each sequence between two alignments (not including gaps)
def frac_corr(direct, file, form):

    file1 = file[0]
    file2 = file[1]
    if len(form) >=2:
        form1 = form[0]
        form2 = form[1]
    else:
        form1 = form[0]
        form2 = form[0]

    alignment1=AlignIO.read(os.path.join(direct, file1), form1)
    alignment2=AlignIO.read(os.path.join(direct, file2), form2)

    fc=[]

    for i in range(len(alignment1)):
        c=0
        C=0
        for j in range(len(alignment1[i].seq)):
            if alignment1[i].seq[j] != '-':
                C+=1
                if alignment1[i].seq[j] == alignment2[i].seq[j]:
                    c+=1
                else:
                    pass
            else:
                pass
        fc.append(round(c/C, 4))
    
    newfile='fraction_correct.txt'
    with open(os.path.join(direct, newfile), 'w+') as ff:
        for k in range(len(fc)):
            ff.write(str(fc[k])+'\n')
        ff.close()

    print(newfile)

#################################################################################################
#
#Argument Parser
#
#################################################################################################

parser = argparse.ArgumentParser()

parser.add_argument("-E", "--explode", action="store_true",
        help="take an alignment and explode all seqs into individual sequence files")
parser.add_argument("-M", "--Aln_Dim", action="store_true",
        help="calculate the number of taxa and the length of the alignment")
parser.add_argument("-R", "--Reformat", action="store_true",
        help="convert one alignment type into another")
parser.add_argument("-D", "--Delete", action="store_true", 
        help="delete any column that is all gaps")
parser.add_argument("-G", "--Gap", action="store_true",
        help="apply gaps from MSA1 to MSA2")
parser.add_argument("-C", "--Correct", action="store_true",
        help="get the fraction of identical residues between two alignments, MSA1 considered true")

parser.add_argument("--MSA", nargs='+', help="names of the MSA file(s)")
parser.add_argument("--form1", nargs='+', help="input formats of the MSA file(s)")
parser.add_argument("--form2", nargs='+', help="output format of the MSA file")
parser.add_argument("--direct", default=os.getcwd(), help="where to direct files, default cwd")

args = parser.parse_args()

#################################################################################################
#
#Running fxns
#
#################################################################################################

if args.explode:
    msa_explode(args.direct, args.MSA, args.form1)

elif args.Aln_Dim:
    align_dim(args.direct, args.MSA, args.form1)

elif args.Reformat:
    change_format(args.direct, args.MSA, args.form1, args.form2)

elif args.Delete:
    del_all_gaps(args.direct, args.MSA, args.form1)

elif args.Gap:
    apply_gaps(args.direct, args.MSA, args.form1)

elif args.Correct:
    frac_corr(args.direct, args.MSA, args.form1)

else:
    print("not an option, check usage: ./ms_seqcon.py -h")
