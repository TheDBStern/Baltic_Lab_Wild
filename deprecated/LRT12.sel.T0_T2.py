#-------------------------------------------------------------------------------
# Name: adapted from J.K. Kelly       
# Purpose:     Fits models and estimates LRT
#-------------------------------------------------------------------------------

import sys
import numpy
from math import log
from math import sqrt
from scipy.stats import norm
from scipy.stats import chi2

windowsize=50000

name = sys.argv[1] # "Scf_2L.varscan"
#chrom= sys.argv[2] # 2L
# Ne estimated from LLS P.planII
Ne = [1755,2011,1807,2861,2401,1963,2675,2597,1755,2011,1807,2861,2401,1963,2675,2597]

Method="Subtraction" # "Regression"

vn=map(lambda x: 1.0/(2*50)+1.0/(2*50)+10.0/(2*x), Ne)


#inZ=open("null.variance.estimates.txt", "rU")
#for line_idx, line in enumerate(inZ):
# B	2L	0.0170	0.0199
#        cols = line.replace('\n', '').split('\t') 
#	if line_idx>0:
#		if Method=="Subtraction" and cols[1]==chrom: 
#			vn.append(float(cols[2])) # v A, v B, v C
#		elif Method=="Regression" and cols[1]==chrom:
#			vn.append(float(cols[3])) # v A, v B, v C


src  =open(name+".z.stats.all.txt", "rU")
out2 =open(name+".LRT.txt", "w")
out1 =open(name+".sig.LRT.txt", "w")
out1a=open(name+".sig.LRT2.txt", "w")
out3 =open(name+".p0mean.txt", "w")
out4 =open(name+".pop.specific.responses.txt", "w")
pdist=[0 for j in range(101)]

m={"F0.BSE3":0,"F0.BSE4":0,"F0.BSE5":0,"F0.BSE6":0,"F0.BSE8":0,"F0.BSE9":0,"F0.BSE11":0,"F0.BSE12":0,"F10.BSE3":0,"F10.BSE4":0,"F10.BSE5":0,"F10.BSE6":0,"F10.BSE8":0,"F10.BSE9":0,"F10.BSE11":0,"F10.BSE12":0}
p_crit=1.0/1000.0
g_snps=0
LRT=[]
LRT2=[]
        
LRT1dist=[0 for j in range(11)]
LRT2dist=[0 for j in range(11)]

windowLRT1=[]
windowLRT2=[]
cwin=0

# evaluate contents of each line of input file
for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t') 

# Scf_2L	15997755_C_G	A0	110	0.9	B0	86	0.93023255814	C0	97	0.907216494845	A7	95	0.947368421053	B7	139	0.884892086331	C7	518	0.916988416988	2.4980915448	2.60697816774	2.5225476806	2.67863792565	2.44929928338	2.55707227019
	position=int( cols[1].split("_")[0] )
	g_snps+=1

	p0num=float(cols[3])*float(cols[4])+float(cols[6])*float(cols[7])+float(cols[9])*float(cols[10])+float(cols[12])*float(cols[13])+float(cols[15])*float(cols[16])+float(cols[18])*float(cols[19])+float(cols[21])*float(cols[22])+float(cols[24])*float(cols[25])
	p0=p0num/( float(cols[3])+float(cols[6])+float(cols[9])+float(cols[12])+float(cols[15])+float(cols[18])+float(cols[21])+float(cols[24]))
	pdist[int(100.0*p0)]+=1


	strx=cols[0]+'\t'+cols[1]
	for j in range(16):
		key = cols[j*3+2]
		m[key]=float(cols[j*3+3])		
		#strx+=('\t'+cols[j*3+2]+'\t'+cols[j*3+3]+'\t'+cols[j*3+4])


	dz=[float(cols[50])-float(cols[58]),float(cols[51])-float(cols[59]),float(cols[52])-float(cols[60]),float(cols[53])-float(cols[61]),float(cols[54])-float(cols[62]),float(cols[55])-float(cols[63]),float(cols[56])-float(cols[64]),float(cols[37])-float(cols[65])]
	vz=[vn[0]+1.0/m["F0.BSE3"]+1.0/m["F10.BSE3"],vn[1]+1.0/m["F0.BSE4"]+1.0/m["F10.BSE4"],vn[2]+1.0/m["F0.BSE5"]+1.0/m["F10.BSE5"],vn[3]+1.0/m["F0.BSE6"]+1.0/m["F10.BSE6"],vn[4]+1.0/m["F0.BSE8"]+1.0/m["F10.BSE8"],vn[5]+1.0/m["F0.BSE9"]+1.0/m["F10.BSE9"],vn[6]+1.0/m["F0.BSE11"]+1.0/m["F10.BSE11"],vn[7]+1.0/m["F0.BSE12"]+1.0/m["F10.BSE12"]]


	w=[1.0/vz[0],1.0/vz[1],1.0/vz[2],1.0/vz[3],1.0/vz[4],1.0/vz[5],1.0/vz[6],1.0/vz[7]]

	#null model (no change)
	LL0=0.0
	for k in range(8):
		LL0-= ( (dz[k]*dz[k])/(2.0*vz[k]) + log(vz[k])/2.0 )

	# alt model (delta Z allowed non-zero)
	dz0 = (w[0]*dz[0]+w[1]*dz[1]+w[2]*dz[2]+w[3]*dz[3]+w[4]*dz[4]+w[5]*dz[5]+w[6]*dz[6]+w[7]*dz[7])/(w[7]+w[6]+w[5]+w[4]+w[3]+w[1]+w[2]+w[0]) 
	LL1=0.0
	for k in range(8):
		LL1-= ( (dz[k]-dz0)*(dz[k]-dz0)/(2.0*vz[k]) + log(vz[k])/2.0 )
	LL3=0.0 # separate mu(dz) for each rep
	for k in range(8):
		LL3-= ( log(vz[k])/2.0 )

	lrt=2*(LL1-LL0)
	LRT.append(lrt)
	vax=int(lrt/10)
	if vax>10:
		vax=10
	LRT1dist[vax]+=1
	tail=1.0-chi2.cdf(lrt, 1) # correct
	out2.write(strx+'\t'+str(lrt)+'\t'+str(tail))
	lrt2=2*(LL3-LL1)
	vax=int(lrt2/10)
	if vax>10:
		vax=10
	LRT2dist[vax]+=1
	tail2=1.0-chi2.cdf(lrt2, 7) # correct
	out2.write('\t'+str(lrt2)+'\t'+str(tail2)+'\n')#+'\t'+str(dz0)+'\t'+str(dz[0])+'\t'+str(vz[0])+'\t'+str(dz[1])+'\t'+str(vz[1])+'\t'+str(dz[2])+'\t'+str(vz[2])+'\t'+str(dz[3])+'\t'+str(vz[3])+'\t'+str(dz[4])+'\t'+str(vz[4])+'\t'+str(dz[5])+'\t'+str(vz[5])+'\t'+str(dz[6])+'\t'+str(vz[6])+'\t'+str(dz[7])+'\t'+str(vz[7])+'\n')
	out4.write(strx+'\t'+str(lrt)+'\n')#+'\t'+str(dz[0]/sqrt(vz[0]))+'\t'+str(dz[1]/sqrt(vz[1]))+'\t'+str(dz[2]/sqrt(vz[2]))+'\t'+str(dz[3]/sqrt(vz[3]))+'\t'+str(dz[4]/sqrt(vz[4]))+'\t'+str(dz[5]/sqrt(vz[5]))+'\t'+str(dz[6]/sqrt(vz[6]))+'\t'+str(dz[7]/sqrt(vz[7]))+'\n')

	if tail < p_crit:
		out1.write(strx+'\t'+str(lrt)+'\t'+str(tail)+'\t'+str(dz0)+'\t'+str(dz[0])+'\t'+str(dz[1])+'\t'+str(dz[2])+'\t'+str(dz[3])+'\t'+str(dz[4])+'\t'+str(dz[5])+'\t'+str(dz[6])+'\t'+str(dz[7])+'\n')
	if tail2 < p_crit:
		out1a.write(strx+'\t'+str(lrt2)+'\t'+str(tail2)+'\t'+str(dz0)+'\t'+str(dz[0])+'\t'+str(dz[1])+'\t'+str(dz[2])+'\t'+str(dz[3])+'\t'+str(dz[4])+'\t'+str(dz[5])+'\t'+str(dz[6])+'\t'+str(dz[7])+'\n')



	# windows
	win = int(position/windowsize)
	if win != cwin:
		out5.write(cols[0]+'\t'+str(cwin*windowsize))
		nu=len(windowLRT1)
		if nu>0:
			out5.write('\t'+str(nu)+'\t'+str(sum(windowLRT1)/float(nu))+'\t'+str(max(windowLRT1)) )
		else:
			out5.write('\t'+str(nu)+'\t-99\t-99' )
		nu=len(windowLRT2)
		if nu>0:
			out5.write('\t'+str(nu)+'\t'+str(sum(windowLRT2)/float(nu))+'\t'+str(max(windowLRT2)) )
		else:
			out5.write('\t'+str(nu)+'\t-99\t-99' )
		out5.write('\n')
		windowLRT1=[]
		windowLRT2=[]
		cwin=win
	windowLRT1.append(lrt)
	windowLRT2.append(lrt2)
	

ax = numpy.array(LRT)
print "mean and var LRT ",numpy.mean(ax),numpy.var(ax)

for j in range(101):
	out3.write(str(j)+'\t'+str(pdist[j])+'\n')

print "total test ",g_snps

for j in range(51):
	print "LRT range ",10*j,LRT1dist[j],LRT2dist[j]


                        
