#!/usr/bin/env python
import sys

blastfile=sys.argv[1]
percentID=int(sys.argv[2])
mindist=int(sys.argv[3])

query2line={}
querylist=[]
with open(blastfile) as blast:
    filteredarray=[line.strip() for line in blast if float(line.strip().split('\t')[2]) > percentID]
    queryarray=set([array.split('\t')[0] for array in filteredarray])
    for array in filteredarray:
        query2line.setdefault(array.split('\t')[0],[]).append(array)
        querylist.append(array.split('\t')[0])
        #print array

querycount=len(set(querylist))
# print querycount
# print "\n"
minimum=1
#print query2line
for query in query2line:
    subjectcoord={}
    totlocarray=[]
    dictcount={}
    for qline in query2line[query]: # run through subjects
        col=qline.split('\t')
        startendarray=[]
        if col[1] !=query:
            startendarray.append(int(col[6]))
            startendarray.append(int(col[7]))
            startendarray=range(int(col[6]), int(col[7])+1)
            #print startendarray
            subjectcoord.setdefault(col[1], []).append(startendarray)
            totlocarray=totlocarray+startendarray
            for item in startendarray:
                dictcount.setdefault(item,[]).append(1)
    subjecttotal=len(subjectcoord)
    totlocarray.sort()
    final=[]
    # print dictcount
    # print "\n"
    if len(dictcount)>0:
        for item in dictcount:
            if float(len(dictcount[item])+1)/querycount>=minimum:    
                final.append(item)
    final.sort()
    if len(final)>0:
        start=final[0]
#        print final
        for i in range(len(final)-1):
            if final[i+1] -final[i]> mindist: 
                end=final[i]
                print query+"\t"+str(start)+"-"+str(end)
                start=final[i+1]
            if i == len(final)-2:
                end=final[i+1]
                print query+"\t"+str(start)+"-"+str(end)
#    print query
                
