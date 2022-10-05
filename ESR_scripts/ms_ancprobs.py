#!/usr/bin/env python3

#####################################################################################################
#Script Name    : ms_ancprobs.py
#Description    : parse through iqtree statefile to generate ancestors
#Author         : Michael A Sennett
#####################################################################################################

from Bio import SeqIO
from Bio import AlignIO
import numpy as np
import re
import argparse
import os
import random as rd

#sample sequences randomly from prob dist
def sample_seqs(pp_arr, node, samples=0):
    
    rd.seed()
    seed=rd.random()
    
    AA=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    
    PP_vector_cum=np.empty((len(pp_arr),len(pp_arr[0])))
    for i in range(len(pp_arr)):
        PP_vector_cum[i]=np.cumsum(pp_arr[i])
        for j in range(len(PP_vector_cum[i])):
            PP_vector_cum[i][j] = np.round(PP_vector_cum[i][j], 4)
    
    seq_lst=[]
    seq_name=[]
    for k in range(samples):
        seq_name.append(">node_"+str(node)+"_s"+str(k+1))
        tmp_seq=[]
        for l in range(len(PP_vector_cum)):
            rd_num=rd.randint(1,10000)/10000
            for m in range(len(PP_vector_cum[l])):
                if rd_num <= PP_vector_cum[l,m]:
                    tmp_seq.append(AA[m])
                    break
                else:
                    continue
        seq_lst.append(tmp_seq)
    
    return seq_name, seq_lst, seed   

#get stats of prob dist; avg exp prob, std dev of exp prob, exp log-prob
#, avg site log-prob, and std dev of site log-prob
def get_arr_stats(pp_arr):
    
    pp_site=[]
    plnp_site=[]
    for i in range(len(pp_arr)):
        pp=[]
        plnp=[]
        for j in range(len(pp_arr[i])):
            prob=pp_arr[i][j]
            if prob == 0:
                plnp.append(0)
                pp.append(0)
            else:
                plnp.append(prob*np.log(prob))
                pp.append(prob*prob)

        pp_site.append(np.sum(pp))
        plnp_site.append(np.sum(plnp))

    avg_exp_prob=np.mean(pp_site)
    std_exp_prob=np.std(pp_site)
    exp_lnp=np.sum(plnp_site)
    avg_lnp=np.mean(plnp_site)
    std_lnp=np.std(plnp_site)
    
    return exp_lnp, avg_lnp, std_lnp, avg_exp_prob, std_exp_prob

#get stats of individ seq; avg of prob, geo avg of prob, log-prob
def get_seq_stats(seq, pp_arr):
    AA=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    res_dct = dict(zip(AA, range(len(AA))))
    
    probs=[]
    for i in range(len(seq)):
        res=seq[i]
        probs.append(pp_arr[i][res_dct[res]])

    avg=np.mean(probs)
    geo_avg=((np.prod(probs))/len(probs))**(1/len(probs))
    lnp=np.sum(np.log(probs))
    
    return avg, geo_avg, lnp

#goes through a statefile and pulls the relevant PP distribution for a node and outputs data in text file
def parse_state(direct, file, node, gap="None", samples=0):
    
    #create a list of all the lines printed?
    text=[]
    
    AA=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    res_dct = dict(zip(AA, range(len(AA))))
    
    ###############################################################################
    #create an array of post prob distribution, and array of gap and amino acid index
    ###############################################################################      
    with open(os.path.join(direct, file), "r") as state:
        state_probs = state.readlines()
        state.close()
    
    seq_state=[]
    
    for line in state_probs:
        if line.startswith("Node"+str(node)+'\t'):
            seq_state.append(line.replace('\n','').split('\t'))
        else:
            pass
    
    g_idx=[]
    if os.path.isfile(os.path.join(direct, gap)):
        path = os.path.join(direct, gap)
        alignment = list(SeqIO.parse(path, 'fasta'))
        for aa in range(len(alignment[0].seq)):
            if alignment[0].seq[aa] == '-':
                g_idx.append(aa)                
            else:
                pass
        smp_seq_un=['-']*(len(seq_state)-len(g_idx))
        smp_seq_gap=['-']*len(alignment[0].seq)
    else:
        smp_seq_un=['-']*len(seq_state)
        smp_seq_gap=smp_seq_un
        

    pp=[]
    ste=[]
    for site in range(len(seq_state)):
        if site not in g_idx:
            ste.append(int(site)+1)
            pp.append(seq_state[site][3:])
        else:
            pass
    
    pp = np.array(pp)    
    pp_arr = pp.astype(float)

    ###############################################################################
    #get the top three residue and prob values for each site
    ###############################################################################
    
    print("%12s %15s %6s %6s" % ("site","prob","prob2","prob3"))
    print("-"*42)
    
    text.append("%12s %15s %6s %6s\n" % ("site","prob","prob2","prob3"))
    text.append("-"*42+"\n")
    
    
    for i in range(len(pp_arr)):
        x=np.argmax(pp_arr[i])
        r1=AA[x]
        p1=pp_arr[i][x]
        y=pp_arr.argsort()[i,-2]
        r2=AA[y]
        p2=pp_arr[i][y]
        z=pp_arr.argsort()[i,-3]
        r3=AA[z]
        p3=pp_arr[i][z]
        j=ste[i]
        
        print("SITE: %6d  {0}  {1}  {2}  {3:.3f}  {4:.3f}  {5:.3f}".format(r1, r2, r3, p1, p2, p3) % (j))
        text.append("SITE: %6d  {0}  {1}  {2}  {3:.3f}  {4:.3f}  {5:.3f}\n".format(r1, r2, r3, p1, p2, p3) % (j))
    ent, site_lnp, std_lnp, avg_exp_prb, std_exp_prb=get_arr_stats(pp_arr)
    
    print("\nsite entropy:{0:13.3f}{1:8.3f}".format(-site_lnp, np.exp(-site_lnp)))
    print("std site entropy:{0:9.3f}\n".format(std_lnp))
    print("entropy:{0:18.3f} {1}".format(-ent, np.format_float_scientific(np.exp(ent), precision=3)))
    print("std entropy:{0:14.3f}\n".format(std_lnp))
    print("\nave exp prob: {0:.6f}".format(avg_exp_prb))
    print("std exp prob: {0:.6f}\n".format(std_exp_prb))

    text.append("\nsite entropy:{0:13.3f}{1:8.3f}\n".format(-site_lnp, np.exp(-site_lnp)))
    text.append("std site entropy:{0:9.3f}\n".format(std_lnp))
    text.append("entropy:{0:18.3f} {1}\n".format(-ent, np.format_float_scientific(np.exp(ent), precision=3)))
    text.append("std entropy:{0:14.3f}\n".format(std_lnp))
    text.append("\nave exp prob: {0:.6f}\n".format(avg_exp_prb))
    text.append("std exp prob: {0:.6f}\n".format(std_exp_prb))
    ###############################################################################
    #get the standard sequences, smp, altall, minKL
    ###############################################################################          
    
    #get the smp sequence gapped and ungapped
    print(">node_{0}_smp_gaps".format(node))
    text.append("\n>node_{0}_smp_gaps\n".format(node))
    smp_text=[]
    
    for j in range(len(pp_arr)):
        x=np.argmax(pp_arr[j])
        r1=AA[x]
        y=ste[j]-1
        smp_seq_gap[y] = r1
        smp_seq_un[j] = r1
    for k in range(0, len(smp_seq_gap), 70):
        l=k+70
        print("".join(smp_seq_gap[k:l]))
        text.append("".join(smp_seq_gap[k:l])+"\n")
    
    print(">node_{0}_smp".format(node))
    text.append("\n>node_{0}_smp\n".format(node))
    smp_text.append(">node_{0}_smp\n".format(node))

    for m in range(0, len(smp_seq_un), 70):
        n=m+70
        print("".join(smp_seq_un[m:n]))
        text.append("".join(smp_seq_un[m:n])+"\n")
        smp_text.append("".join(smp_seq_un[m:n])+"\n")

    smp_avg, smp_geo, smp_lnp=get_seq_stats(smp_seq_un, pp_arr)
    smp_prob=np.format_float_scientific(np.exp(smp_lnp), precision=3)
    
    print("\nlog prob: {0:.3f} {1}\nave prob: {2:.3f}\ngeo prob: {3:.3f}\n".format(smp_lnp, smp_prob, smp_avg, smp_geo))
    text.append("\nlog prob: {0:.3f} {1}\nave prob: {2:.3f}\ngeo prob: {3:.3f}\n".format(smp_lnp, smp_prob, smp_avg, smp_geo))

    smp_file="node_"+str(node)+"_smp.fst"
    with open(os.path.join(direct,smp_file),'w+') as smp_out:
        smp_out.writelines(smp_text)
    smp_out.close()

    #get the AltAll sequence
    AltA_text=[]
    AltA_seq_un=['-']*len(pp_arr)
    for o in range(len(pp_arr)):
        x=pp_arr[o][pp_arr.argsort()[o,-1]]
        y=pp_arr[o][pp_arr.argsort()[o,-2]]
        r=AA[pp_arr.argsort()[o,-1]]
        if x <= 0.8 and y >= 0.2:
            r=AA[pp_arr.argsort()[o,-2]]
            AltA_seq_un[o]=r
        else:
            r=AA[pp_arr.argsort()[o,-1]]
            AltA_seq_un[o]=r
    
    print(">node_{0}_AltAll".format(node))
    text.append("\n>node_{0}_AltAll\n".format(node))
    AltA_text.append(">node_{0}_AltAll\n".format(node))

    for k in range(0, len(AltA_seq_un), 70):
        l=k+70
        print("".join(AltA_seq_un[k:l]))
        text.append("".join(AltA_seq_un[k:l])+"\n")
        AltA_text.append("".join(AltA_seq_un[k:l])+"\n")

    AltA_avg, AltA_geo, AltA_lnp=get_seq_stats(AltA_seq_un, pp_arr)
    AltA_prob=np.format_float_scientific(np.exp(AltA_lnp), precision=3)
    
    print("\nlog prob: {0:.3f} {1}\nave prob: {2:.3f}\ngeo prob: {3}\n".format(AltA_lnp, AltA_prob, AltA_avg, AltA_geo))
    text.append("\nlog prob: {0:.3f} {1}\nave prob: {2:.3f}\ngeo prob: {3}\n".format(AltA_lnp, AltA_prob, AltA_avg, AltA_geo))
    AltA_file="node_"+str(node)+"_AltAll.fst"
    with open(os.path.join(direct,AltA_file),'w+') as AltA_out:
        AltA_out.writelines(AltA_text)
    AltA_out.close()
    ###############################################################################
    #get the sampled sequences and sats
    ############################################################################### 

    if samples > 0:
        samp_names=[]
        samp_seqs=[]
        samp_names, samp_seqs, seed_num=sample_seqs(pp_arr, node, samples)
        
        print('\nseed: {0}\n'.format(seed_num))
        text.append('\nseed: {0}\n'.format(seed_num))
        
        for p in range(len(samp_names)):
            samp_avg, samp_geo, samp_lnp = get_seq_stats(samp_seqs[p], pp_arr)
            samp_prob=np.format_float_scientific(np.exp(samp_lnp), precision=3)
            
            print(samp_names[p])
            text.append("\n"+samp_names[p]+"\n")
            
            for k in range(0, len(samp_seqs[p]), 70):
                l=k+70
                print("".join(samp_seqs[p][k:l]))
                text.append("".join(samp_seqs[p][k:l])+"\n")
            
            print("\nlog prob: {0:.3f} {1}\nave prob: {2:.3f}\ngeo prob: {3}\n".format(samp_lnp, samp_prob, samp_avg, samp_geo))
            text.append("\nlog prob: {0:.3f} {1}\nave prob: {2:.3f}\ngeo prob: {3}\n".format(samp_lnp, samp_prob, samp_avg, samp_geo))
    else:
        pass
    
    newfile="node_"+str(node)+"_ancprobs.out"
    with open(os.path.join(direct,newfile), 'w+') as outfile:
        for line in text:
            outfile.write(line)
        outfile.close()
#################################################################################################
#
#Argument Parser
#
#################################################################################################

parser = argparse.ArgumentParser()

parser.add_argument("-S", "--Seq_Sampled", type=int, default=0,
        help="take an alignment and explode all seqs into individual sequence files")
parser.add_argument("-G", "--Gap_file", default="None",
        help="apply gaps from provided seq in fasta format to ancestor")
parser.add_argument("-D", "--direct", default=os.getcwd(), 
        help="where to direct files, default cwd")

requiredNamed = parser.add_argument_group('required arguments')

requiredNamed.add_argument("-F", "--File", 
        help="iqtree .state file with prob distribution")
requiredNamed.add_argument("-N", "--Node", 
        help="ancestral node to reconstruct")


args = parser.parse_args()

#################################################################################################
#
#Running fxns
#
#################################################################################################

if args.Gap_file and args.Seq_Sampled:
    parse_state(args.direct, args.File, args.Node, args.Gap_file, args.Seq_Sampled)

elif args.Gap_file:
    parse_state(args.direct, args.File, args.Node, args.Gap_file)

elif args.Seq_Sampled:
    parse_state(args.direct, args.File, args.Node, args.Gap_file, args.Seq_Sampled)

elif args.File and args.Node:
    parse_state(args.direct, args.File, args.Node)
else:
    print("try again, ./ms_ancprobs.py -h")
