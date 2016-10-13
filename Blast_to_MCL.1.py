#!/usr/bin/python
#python2 script
#Blast_to_MCL.py <blastalnfile> 
import sys
import subprocess

Blastfile=sys.argv[1]


def make_start_end_dict(alnfile):
    allstartenddict={}
    q_allstartenddict={}
    q_len_dict={}
    s_len_dict={}
    bitscoredict={}
    with open(alnfile, "r") as aln:
	for line in aln:
	    linearray=line.strip().split('\t')
            sname=linearray[1]
	    qname=linearray[0]
	    if qname not in q_allstartenddict:
	        allstartenddict={}
                subbitscore={}
            subbitscore[sname]=0
	    allstartenddict.setdefault(sname, []).append([int(linearray[6]), int(linearray[7])])#should sort
            subbitscore[sname]=float(linearray[11])+subbitscore[sname]
	    q_len_dict[qname]=int(linearray[-1])#change depending on aln file
	    q_allstartenddict[qname]=allstartenddict
            bitscoredict[qname]=subbitscore
    return q_allstartenddict, q_len_dict, bitscoredict

def collapse1(allstartenddict):
    for genomes in allstartenddict:
	allstartenddict[genomes].sort(key=lambda x: x[0])
	x=0
        while x < len(allstartenddict[genomes])-1:
	    s1, s2, e1, e2 =allstartenddict[genomes][x][0], allstartenddict[genomes][x+1][0], allstartenddict[genomes][x][1], allstartenddict[genomes][x+1][1]
            if s2 >= s1 and e1>= e2:
                allstartenddict[genomes].remove(allstartenddict[genomes][x+1])
                x=x-1
            elif s2 <= e1 and e2 >= e1:
                allstartenddict[genomes][x][1]=e2
                allstartenddict[genomes].remove(allstartenddict[genomes][x+1])               
	        x=x-1
	    x=x+1
    return allstartenddict

def runcollapse(alnfile):
    """gets shared locations of each virus to reference phage"""
    allstartenddict, q_len_dict, bitscoredict=make_start_end_dict(alnfile)
    newdict={}
    for subject in allstartenddict:
	newdict[subject]=collapse1(allstartenddict[subject])
    #print newdict
    return newdict, q_len_dict,bitscoredict

def print_PID(collapsed_dict, q_len_dict, bitscoredict):
    for query in collapsed_dict:
	for subject in collapsed_dict[query]:
	    PID= sum(x[1]-x[0]+1 for x in collapsed_dict[query][subject])/float(q_len_dict[query])
            summation= sum(x[1]-x[0]+1 for x in collapsed_dict[query][subject])
            if subject not in query:
	        print '\t'.join([query, subject,str(PID), str(summation), str(bitscoredict[query][subject])])

try:
    open(sys.argv[1], "r")
    #print '\t'.join(['query', 'subject', 'totalaligned/lengthquery', 'totalaligned', 'lengthquery'])
    newdict, q_len_dict, bitscoredict = runcollapse(Blastfile)   
    print_PID(newdict, q_len_dict, bitscoredict)
except IOError:
    print "Blast_to_MCL.1.py <blastfile_with_qlen>\nCan't open blast file"

#print_PID(newdict, q_len_dict)



