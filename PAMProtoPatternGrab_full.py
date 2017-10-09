#!/usr/bin/env python

from __future__ import division
import sys
#import glob
import os
import subprocess
import re

spacerfile=sys.argv[1]
dbfile=sys.argv[2]

alnfilename=os.path.basename(dbfile)[:]+'_vs_'+os.path.basename(spacerfile)

os.mkdir(alnfilename+'.dir')
alnfilename1=alnfilename+'.dir/'+alnfilename+'.aln'
subprocess.call(['blastn', '-task', 'blastn-short', '-query', spacerfile, '-db', dbfile, '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen btop', '-max_target_seqs', '1000000', '-evalue', '.1', '-out', alnfilename1])
print alnfilename1
filteredaln=alnfilename+'.dir/'+alnfilename+'.filt.aln'
extra=alnfilename+'.dir/'+alnfilename+'.extra.aln'
filteredalnfile=open(filteredaln, "w")
extraaln=open(extra, "w")
alnfile=open(alnfilename1, "r")

def checkgaps(BTOP):
	BTOParray=re.split('(\D+)', BTOP)
	x=1
	start=BTOParray[0]
	for characters in BTOParray[1].split():
		if characters=='-':
			if x%2==1:
				return BTOP
		else:
			if x%2==0:
				return BTOP
		x=x+1

spacerfiledict={}
spacers=open(spacerfile, "r")
for line in spacers:
	if line.startswith(">"):
		spacerID=line[1:].strip()
	else:
		spacerfiledict[spacerID]=line.strip()
#print spacerfiledict

def blastdbcm(linearray):
	if int(linearray[8])<int(linearray[9]): #positive strand
		range1=linearray[8]+'-'+linearray[9]
		entry=subprocess.Popen(['blastdbcmd', '-entry', linearray[1], '-db', dbfile,'-strand','plus','-range',range1], stdout=subprocess.PIPE)
        	fastaarray=entry.communicate()[0].split('\n')
        	protostring=fastaarray[1]
        if int(linearray[8])>int(linearray[9]):
                range1=linearray[9]+'-'+linearray[8] #negative strand
                entry=subprocess.Popen(['blastdbcmd', '-entry', linearray[1], '-db', dbfile,'-strand','minus','-range',range1], stdout=subprocess.PIPE)
                fastaarray=entry.communicate()[0].split('\n')
                protostring=fastaarray[1]

def check6thpos(protostring, spacerstring):
	proto=list(protostring)
	spacer=list(spacerstring)
	i=0
	distance=0
	if len(proto)==len(spacer):
		while len(proto)>i:
			if proto[i] !=spacer[i]:
				if i % 5 ==0:
					distance=distance
				else:
					distance=distance+1
			i=i+1
	else:
		distance=len(spacer)
	if distance > 0:
		return 1
	else:
		return 0

def gethamming(protostring, spacerstring):
	proto=list(protostring)
	spacer=list(spacerstring)
	i=0
	distance=0
	if len(proto)==len(spacer):
		while len(proto)>i:
			if proto[i]!=spacer[i]:
				distance=distance+1	
			i=i+1	
	else:	
		distance=len(spacer)
	return distance
def getPID(protostring, spacerstring):
        proto=list(protostring)
        spacer=list(spacerstring)
        i=0
        sim=0
        if len(proto)==len(spacer):
                while len(proto)>i:
                        if proto[i]==spacer[i]:
                                sim=sim+1
                        i=i+1
	return float(sim)/len(spacer)

CCcount=0
PIDdict={}
newalnfilearray=[]
CCcountEXT=0
countlines=0
for line in alnfile:
	countlines=countlines+1
	linearray=line.strip().split('\t')
	spacersize=len(spacerfiledict[linearray[0]])
	endextlen=spacersize-int(linearray[7]) #figure out
	startextlen=int(linearray[6])-1 #figure out
	BTOP=linearray[-1]
	BTOParray=re.split('(\D+)', BTOP)
	slen=int(linearray[-2])
	qlen=int(linearray[-3])
	spacerstring=spacerfiledict[linearray[0]]
	#print endextlen
	#print startextlen
        fastaarray=[]
        protostring="DOES NOT MEET CONDITION"
	if spacersize == int(BTOParray[0]): # fillupfile
		filteredalnfile.write(line)
                #extralnline='\t'.join([line.strip(),spacerfiledict[linearray[0]], spacerfiledict[linearray[0]], '1'])
                #extraaln.write(extralnline+'\n')
		if int(linearray[8])<int(linearray[9]) and int(linearray[8])-startextlen-3>0 and int(linearray[9])+endextlen+3<slen: #positive strand
                        range1=str(int(linearray[8])-startextlen-3)+'-'+str(int(linearray[9])+endextlen+3)
                        entry=subprocess.Popen(['blastdbcmd', '-entry', linearray[1], '-db', dbfile,'-strand','plus','-range',range1], stdout=subprocess.PIPE)
                        fastaarray=entry.communicate()[0].split('\n')
                        spacerstring=spacerfiledict[linearray[0]]
                        if len(fastaarray)<2:
                                protostring=fastaarray[0]
                        else:
                                protostring=fastaarray[1]
                elif int(linearray[8])>int(linearray[9]) and int(linearray[9])-endextlen-3>0 and int(linearray[8])+startextlen+3<slen: #negative strand
                        range1=str(int(linearray[9])-endextlen-3)+'-'+str(int(linearray[8])+startextlen+3) #negative strand
                        entry=subprocess.Popen(['blastdbcmd', '-entry', linearray[1], '-db', dbfile,'-strand','minus','-range',range1], stdout=subprocess.PIPE)
                        fastaarray=entry.communicate()[0].split('\n')
                        if len(fastaarray)<2:
				protostring=fastaarray[0]
			else:
				protostring=fastaarray[1]
			spacerstring=spacerfiledict[linearray[0]]
		if len(fastaarray)>1 and len(spacerstring)==len(protostring[3:-3]):# and check6thpos(protostring[3:-3], spacerstring)==0:# and protostring[1:3]=='CC':
			print linearray[1], '\t', protostring[3:-3], '\t', protostring[0:3], '\t', protostring[1:3],'\t', protostring[-3:-1], '\t', BTOP,'\t',0 
			#print linearray[0], '\t', spacerstring, '\t', [range]
			if "CC" in protostring[1:3]:
				CCcount=CCcount+1	
       #         extralnline='\t'.join([line.strip(),spacerfiledict[linearray[0]], protostring[3:-3], str(getPID(protostring[3:-3],spacerstring)), protostring[0:3], protostring[-3:-1]])
       #         extraaln.write(extralnline+'\n')
	#if they have mismatches	
	elif int(linearray[3])==spacersize and '-' not in BTOP and spacersize != int(BTOParray[0]):
	        if int(linearray[8])<int(linearray[9]) and int(linearray[8])-startextlen-3>0 and int(linearray[9])+endextlen<slen: #positive strand
                        range1=str(int(linearray[8])-startextlen-3)+'-'+str(int(linearray[9])+endextlen+3)
                        entry=subprocess.Popen(['blastdbcmd', '-entry', linearray[1], '-db', dbfile,'-strand','plus','-range',range1], stdout=subprocess.PIPE)
                        fastaarray=entry.communicate()[0].split('\n')
                        spacerstring=spacerfiledict[linearray[0]]
			if len(fastaarray)<2:
				protostring=fastaarray[0]
			else:
				protostring=fastaarray[1]
                elif int(linearray[8])>int(linearray[9]) and int(linearray[9])-endextlen-3>0 and int(linearray[8])+startextlen+3<slen: #negative strand
                        range1=str(int(linearray[9])-endextlen-3)+'-'+str(int(linearray[8])+startextlen+3) #negative strand
                        entry=subprocess.Popen(['blastdbcmd', '-entry', linearray[1], '-db', dbfile,'-strand','minus','-range',range1], stdout=subprocess.PIPE)
                        fastaarray=entry.communicate()[0].split('\n')
			if len(fastaarray)<2:
				protostring=fastaarray[0]
			else:
				protostring=fastaarray[1]
			spacerstring=spacerfiledict[linearray[0]]
		if len(fastaarray)>1 and len(spacerstring)==len(protostring[3:-3]):# and check6thpos(protostring[3:-3], spacerstring)==0:# and protostring[1:3]=='CC':
			print linearray[1], '\t', protostring[3:-3], '\t', protostring[0:3],'\t', protostring[1:3], '\t',protostring[-3:-1], '\t', BTOP, '\t', gethamming(protostring[3:], spacerstring) 
			#print linearray[0], '\t', spacerstring, '\t', [range]
			if "CC" in protostring[1:3]:
				CCcount=CCcount+1
		if gethamming(protostring, spacerstring)>4:
			filteredalnfile.write(line)
                extralnline='\t'.join([line.strip(),spacerfiledict[linearray[0]], protostring[3:-3], str(getPID(protostring[3:-3],spacerstring)), protostring[0:3], protostring[-3:-1]])
                extraaln.write(extralnline+'\n')
	elif int(linearray[3])<spacersize and '-' not in BTOP and float(linearray[-5])<float(.1): #partial matches and extension
		#print endextlen
		#print startextlen
		qstart=int(linearray[6])
		qend=int(linearray[7])
		if int(linearray[8])<int(linearray[9]) and int(linearray[8])-startextlen-3>0 and int(linearray[9])+endextlen+3<slen-1: #positive strand
                        range1=str(int(linearray[8])-startextlen-3)+'-'+str(int(linearray[9])+endextlen+3)
                        entry=subprocess.Popen(['blastdbcmd', '-entry', linearray[1], '-db', dbfile,'-strand','plus','-range',range1], stdout=subprocess.PIPE)
                        fastaarray=entry.communicate()[0].split('\n')
			if len(fastaarray)<2:
				protostring=fastaarray[0]
			if len(fastaarray)>1:
				protostring=fastaarray[1]
			    	if len(protostring[3:-3])==len(spacerstring):# and check6thpos(protostring[3:-3], spacerstring)==0:# and protostring[1:3]=='CC':
					spacerstring=spacerfiledict[linearray[0]]
					##print linearray[1], '\t', protostring[3:-3], ##'\t', protostring[0:3], '\t', protostring[1:3], '\t',protostring[-3:-1],'\t', BTOP, '\t', gethamming(protostring[3:],spacerstring)
					##print linearray[0], '\t', spacerfiledict[linearray[0]], '\t', [range]
					startprotoext=protostring[3:qstart+2]
					#print startprotoext
					startspacerext=spacerstring[0:qstart-1]
					#print startspacerext
					endprotoext=protostring[qend+3:-3]
					endspacerext=spacerstring[qend:qlen]
					if len(startprotoext)>0:# and "CC" in protostring[1:3]:
						#print len(startprotoext)
						startPID=getPID(startprotoext, startspacerext)
						PIDdict.setdefault(len(startprotoext), []).append(startPID)
						#print len(startprotoext),'\t',startPID
					if len(endprotoext)>0:# and "CC" in protostring[1:3]:
						#print len(endprotoext)
#						print protostring[1:3]
#						print protostring
						endPID=getPID(endprotoext, endspacerext)
						PIDdict.setdefault(len(endprotoext), []).append(endPID)
					#	print len(endprotoext),'\t',endPID
					if "CC" in protostring[1:3] and "CC" in protostring[1:3]:
						CCcountEXT=CCcountEXT+1
					#if len(startprotoext) ==11 and "CC" in protostring[1:3]: #or len(endprotoext)==5:
					#	print protostring[3:-3]
					#	print spacerstring
					#	print startprotoext
					#	#print endprotoext
					#if len(endprotoext)==11 and "CC" in protostring[1:3]:
					#	print protostring[3:-3]
					#	print spacerstring
					#	print endprotoext
			if gethamming(protostring, spacerstring)>4:
				filteredalnfile.write(line)
                elif int(linearray[8])>int(linearray[9]) and int(linearray[9])-endextlen-3>0 and int(linearray[8])+startextlen+3<slen-1: #negative strand
                        range1=str(int(linearray[9])-endextlen-3)+'-'+str(int(linearray[8])+startextlen+3) #negative strand
                        entry=subprocess.Popen(['blastdbcmd', '-entry', linearray[1], '-db', dbfile,'-strand','minus','-range',range1], stdout=subprocess.PIPE)
                        fastaarray=entry.communicate()[0].split('\n')
			if len(fastaarray)<2:	
				protostring=fastaarray[0]
			if len(fastaarray)>1:
				protostring=fastaarray[1]
			        if len(protostring[3:-3])==len(spacerstring):# and check6thpos(protostring[3:-3], spacerstring)==0:# and protostring[1:3]=='CC':
					startprotoext= protostring[3:qstart+2]
					startspacerext= spacerstring[0:qstart-1]
					endprotoext= protostring[qend+3:-3]
					endspacerext= spacerstring[qend:qlen]
					spacerstring=spacerfiledict[linearray[0]]
					#if len(startprotoext)==11 and "CC" in protostring[1:3]: #or len(endprotoext)==5:
					#	print protostring[3:-3]
					#	print spacerstring
					#	print startprotoext
					#	print endprotoext
					if len(startprotoext)>0:# and "CC" in protostring[1:3]:
						startPID=getPID(startprotoext, startspacerext)
						PIDdict.setdefault(len(startprotoext),[]).append(startPID)
					if len(endprotoext)>0:# and "CC" in protostring[1:3]:
						endPID=getPID(endprotoext, endspacerext)
						PIDdict.setdefault(len(endprotoext),[]).append(endPID)
					#print linearray[1], '\n', protostring[3:-3],'\n', spacerstring ##'\t', protostring[0:3], '\t', protostring[1:3],'\t',protostring[-3:-1], '\t', BTOP, '\t', gethamming(protostring[3:], spacerstring)
					##print linearray[0], '\t', spacerstring, '\t', [range]
					if "CC" in protostring[1:3]:
						CCcount=CCcount+1
			if gethamming(protostring, spacerstring)>4:
				filteredalnfile.write(line)
#                extralnline='\t'.join([line.strip(),spacerfiledict[linearray[0]], protostring[3:-3], str(getPID(protostring[3:-3],spacerstring)), protostring[0:3], protostring[-3:-1]])
#                extraaln.write(extralnline+'\n')
	elif '-' in BTOP:# check for gaps
                print 'gap found'
		x=1
		y=1
		querygaps=0
		subjectgaps=0
		BTOParray= re.split('(\D+)',BTOP)
		query=[]
		subject=[]
		if int(linearray[8])<int(linearray[9]) and int(linearray[8])-startextlen-3>0 and int(linearray[9])+endextlen<slen: #positive strand
                        range=str(int(linearray[8])-startextlen-3)+'-'+str(int(linearray[9])+endextlen)
                        entry=subprocess.Popen(['blastdbcmd', '-entry', linearray[1], '-db', dbfile,'-strand','plus','-range',range], stdout=subprocess.PIPE)
			fastaarray=entry.communicate()[0].split('\n')
		elif int(linearray[8])>int(linearray[9]) and int(linearray[9])-endextlen>0 and int(linearray[8])+startextlen+3<slen:
                        range=str(int(linearray[9])-endextlen)+'-'+str(int(linearray[8])+startextlen+3) #negative strand
                        entry=subprocess.Popen(['blastdbcmd', '-entry', linearray[1], '-db', dbfile,'-strand','minus','-range',range], stdout=subprocess.PIPE)
                        fastaarray=entry.communicate()[0].split('\n')
		charcount=0
		protospacerchars=list(str(fastaarray[1]))
		spacerchars=list(spacerfiledict[linearray[0]])
		#add gaps and get hamming distance
		for characters in list(BTOParray[1]):
			if characters=='-':
				if x%2==1: #gap in spacer
					insertloc=charcount+int(BTOParray[0])
					spacerchars.insert(insertloc, "-")
					charcount=charcount+1
				elif x%2==0: #gap in protospacer
					insertloc=charcount+int(BTOParray[0])
					protospacerchars.insert(insertloc, "-")
					charcount=charcount+1
			x=x+1
			#charcount=charcount+1
		print protostring
		print ''.join(protospacerchars)
		print spacerstring
		if len(fastaarray)>1 and len(spacerstring)==len(protostring[3:-3]):	
			protostring=''.join(protospacerchars)
			##print linearray[1], '\t', protostring[3:],'\t', protostring[0:3], '\t', protostring[1:3], '\t', BTOP, '\t', gethamming(protostring[3:],spacerstring)
			spacerstring=''.join(spacerchars)
			##print linearray[0], '\t', spacerstring, '\t', [range]
			if "CC" in protostring[1:3]:
				CCcount=CCcount+1
		distance=gethamming(protostring,spacerstring)
		if distance>4:
			filteredalnfile.write(line)
	extralnline='\t'.join([line.strip(),spacerfiledict[linearray[0]], protostring[3:-3], str(getPID(protostring[3:-3],spacerstring)), protostring[0:3], protostring[-3:-1]])
	extraaln.write(extralnline+'\n')
        protostring="DOES NOT MEET CONDITION"
###print "total PAM:", CCcount
###print "total alignments:",countlines
#print PIDdict

###totalext=0
######PIDsum=0
###for items in PIDdict:
###	totalext=len(PIDdict[items])+totalext
###	PIDsum=sum(PIDdict[items])+PIDsum
###	print items, '\t', sum(PIDdict[items])/len(PIDdict[items]),'\t', len(PIDdict[items])

###print PIDsum/totalext
