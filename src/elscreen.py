from configuration import *
from couplib.constants import *
import math
import numpy as np

class ElScreen(object):
	"""Class for computing electrostatic screening and damping factors"""
#-------------------------------------------------------------------------------
	def __init__(self, cfg):
		"""Initialize with current configuration"""
		self.cfg = cfg
		#Electrostatic screening factor
		self.ElScrModel = int(self.cfg.MolSys[CFG_MSY_EL_SCR_MODEL])

		#No-electrostatic screening
		if (self.ElScrModel == EL_SCR_MODEL_NONE):
			self.ElScrFact = 1.0
		#Uniform screening
		elif (self.ElScrModel == EL_SCR_MODEL_UNIFORM):
			self.ElScrFact = float(self.cfg.MolSys[CFG_MSY_EL_SCR_FACTOR])
		#Exponential function
		elif (self.ElScrModel == EL_SCR_MODEL_EXPONENTIAL):
			self.ElScrFact = 0.0 #Not used
			(self.ElScrExp_A, self.ElScrExp_B,self.ElScrExp_S) = [float(x) for x in self.cfg.MolSys[CFG_MSY_EL_SCR_EXP_FUNC].split(LIST_SEPARATOR)]
		return
#-------------------------------------------------------------------------------
	def GetElScreenDamp(self, R):
		"""Calculates electrostatic screening and damping effects"""
		print("R (Bohr):", float(R))
		
		print("Electrostatic screening model:", end=' ')
		#No-electrostatic screening
		Factor = self.ElScrFact
		
		if (self.ElScrModel == EL_SCR_MODEL_NONE):
			print(" No screening")
		#Uniform screening
		elif (self.ElScrModel == EL_SCR_MODEL_UNIFORM):
			print("Uniform screening")
			if (Factor == EL_SCR_FACTOR_DEFAULT):
				#Refractive index
				n = float(self.cfg.MolSys[CFG_MSY_RFX])
				if (n != CFG_MSY_RFX_DEFAULT):
					print("Refractive index N=",n," will be used to calculate the uniform screening factor as 1/n^2")
					Factor = 1.0/(n*n)

		#Exponential function
		elif (self.ElScrModel == EL_SCR_MODEL_EXPONENTIAL):
			print("Exponential screening function:")
			print("Scr=A*exp(-Beta*R)+S where A=",self.ElScrExp_A,"Beta=",self.ElScrExp_B,"A-1, S=",self.ElScrExp_S)
			ElScrExp_B_Bohr_1 = self.ElScrExp_B/ATOB
			print("Beta=",ElScrExp_B_Bohr_1,"Bohr-1")
			print("exp(-Beta*R)=",math.exp(-ElScrExp_B_Bohr_1*R))
			Factor = self.ElScrExp_A*math.exp(-ElScrExp_B_Bohr_1*R)+self.ElScrExp_S
		print("Screening factor:",Factor)
		return Factor
		#return self.ElScreen*(1.0-math.exp(-1.0*self.ElExpDam*R))
#-------------------------------------------------------------------------------