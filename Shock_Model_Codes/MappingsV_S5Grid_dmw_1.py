# This program will serve to run Mappings V with some specified grid of values.

import subprocess
from os import linesep
from joblib import Parallel, delayed
import numpy as np

class MapRun:
	def __init__(self, OffSetsPath,deltaV,deltaN,v_init,v_final,density_init, density_final,IT = 1,cores = -1):
		self.V = np.arange(v_init,v_final + deltaV,deltaV)
		self.N = np.arange(density_init,density_final + deltaN,deltaN)
		self.num_cores = cores
		self.iterations = str(IT)
		# bk = 1.3807*(10**-16)
		# m = 0.7*2*1.6733*(10**-24) 
		# Tin = str(T_init)
		# ion = str(ion_frac)
		self.k = 1
		self.OffPath = str(OffSetsPath)

		self.abund = "/home/dmw/HM89.txt"

		subprocess.call(["mkdir","ShockFiles"])
 
	def run_mv(self, i,j):

		# if i == 0 and j == 0:
		# 	subprocess.call(["screen \n \n"])

		# else:
		# 	subprocess.call(["screen -d"])
		# 	subprocess.call(["screen \n \n"])
		# 	self.k += 1 


		self.k = i*len(self.N) + j + 1

		v = str(self.V[i])

		n0 = str(self.N[j])

		directory = "Run" + str(self.k) + "_" + v + "kms" + "_" + n0 + "gcm^3"
		dynofile = "ShockDynamics" + str(self.k) + ".txt"
		specfile = "ShockSpectrum" + str(self.k) + ".txt"
		subprocess.call(["mkdir",directory])

		Run = subprocess.Popen("map51",cwd=directory,stdin=subprocess.PIPE,universal_newlines=True)

		newline = linesep

		commands = ["y",self.abund,"n","y",self.OffPath,"n","n","n","S5","D","x","n","V","B","0","200",n0,v,"F","c","100",self.iterations,"Shock","C","A", directory,"E"]

		Run.communicate( newline.join( commands))

		# Run.wait()

		subprocess.call(["cat " + directory + "/" + "Shockdyn0001.sh5 > " + "ShockFiles" + "/" + dynofile], shell=True)
		subprocess.call(["cat " + directory + "/" + "Shockspec0001.csv > " + "ShockFiles" + "/" + specfile], shell=True)


		# subprocess.call(["echo " + v + ">>" + "ShockFiles" + "/" + file],shell=True)
		# subprocess.call(["echo " + n0 + ">>" + "ShockFiles" + "/" + file],shell=True)

	def run(self):
		Parallel(n_jobs = self.num_cores)(delayed(self.run_mv)(i, j) for i in range(len(self.V)) for j in range(len(self.N)) )

# mr = MapRun("/home/tpierrej/Data/OffSets/NoOffSets.txt", 1,100,30,75,100,6000)
# mr.run()

# def MapRun(OffSetsPath,deltaV,deltaN,v_init,v_final,density_init, density_final):

# 	V = np.arange(v_init,v_final + deltaV,deltaV)
# 	N = np.arange(density_init,density_final + deltaN,deltaN)
# 	# bk = 1.3807*(10**-16)
# 	# m = 0.7*2*1.6733*(10**-24) 
# 	# Tin = str(T_init)
# 	# ion = str(ion_frac)
# 	k = 1
# 	OffPath = str(OffSetsPath)

# 	Runnings = np.arange(0,2760,184)

# 	abund = "/home/tpierrej/Data/NewAbunds.txt"

# 	subprocess.call(["mkdir","ShockFiles"])
                     

# 	def RunMV(i,j):

# 		if i == 0 and j == 0:
# 			subprocess.call(["screen \n \n"])

# 		else:
# 			subprocess.call(["screen -d"])
# 			subprocess.call(["screen \n \n"])
# 			k += 1 



# 		v = str(V[i])
		
# 		n0 = str(N[j])

# 		directory = "Run" + str(k) + "_" + v + "kms" + "_" + n0 + "gcm^3"
# 		dynofile = "ShockDynamics" + str(k) + ".txt"
# 		specfile = "ShockSpectrum" + str(k) + ".txt"
# 		subprocess.call(["mkdir",directory])

# 		Run = subprocess.Popen("map51",cwd=directory,stdin=subprocess.PIPE,universal_newlines=True)

# 		newline = linesep

# 		commands = ["y",abund,"n","y",OffPath,"n","n","n","S5","D","x","n","V","B","0","200",n0,v,"F","c","100","1","Shock","C","A", directory,"E"]

# 		Run.communicate( newline.join( commands))

# 		# Run.wait()

# 		subprocess.call(["cat " + directory + "/" + "Shockdyn0001.sh5 > " + "ShockFiles" + "/" + dynofile], shell=True)
# 		subprocess.call(["cat " + directory + "/" + "Shockspec0001.csv > " + "ShockFiles" + "/" + specfile], shell=True)


# 		# subprocess.call(["echo " + v + ">>" + "ShockFiles" + "/" + file],shell=True)
# 		# subprocess.call(["echo " + n0 + ">>" + "ShockFiles" + "/" + file],shell=True)

		

# 	Parallel(n_jobs=-1)(delayed(MapRun.RunMV)(i, j) for i in range(len(V)) for j in range(len(N)) )
		

		

			



