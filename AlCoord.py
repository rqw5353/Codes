#Removes all modifiers
#Loop over all oxygens
#For each oxygen, determine nearest neighbors within cutoff
#If only 1 nearbest neighbor, delete that oxygen from main data array
#Then for the remaining data, choose specific network former of interest
#Count nearest neighbors of that network former, this is your Q speciation
	#Coordination can be calculated following the same procedure but keeping all NBOs




import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.optimize
import os

def getData(file):
	data = np.loadtxt(file, delimiter = ' ', skiprows = 2)
	print ('Reading: ' + file)
	return data


cutoff = 2.1

def pbc(dist,L):

    dist -= np.around(dist/L) * L

    return dist




#1 = O
#2 = Si
#3 = Al
#4 = Ca

def process(data):
	data[:,1:] -= np.min(data[:,1:])

	Q_O = np.zeros(30)

	a = 0

#------- Removes modifiers (4) ------------------------------#
	locs = np.where(data[:,0] != 4)[0]
	data = data[locs,:]

#------- Look at specific network former and count nearest neighbors ---------------------------#
#------- Oxygen ---------------------#
	b = 0
	cutoff = 2.1
	while b < len(data[:,0]):
		intAtom = b
		if data[b,0] == 3:
			dist = (data[intAtom,:] - np.copy(data))[:,1:]
			dist = pbc(dist,np.max(data))
			d = np.linalg.norm(dist,axis=1)
			dist = dist[d < cutoff]
			Q_O[len(dist)-1] += 1.
			
		b = b + 1

	return (Q_O)





def read(dir,n):
	out = open('Al_Coord.'+dir+'.dat','w')
	f = glob.glob('./'+dir+'/dump.*.xyz')
	print(f)
	for a in range(len(f)):
		file = f[a]
		d = getData(file)
		Q = process(d)
		#Q = Q/ np.sum(Q)
		msg = str(f[a])

		for c in Q:
			msg += "\t" + str(c)

		print(msg)
		out.write(msg + '\n')
	out.close()

for i in glob.glob('3.5Ratio'):
	os.chdir(i)
	print(os.getcwd())
	for j in glob.glob('3500Temp'):
		read(j,j[0])
	os.chdir('..')
