from couplib.constants import *
from configuration import *
import math
import numpy as np
from interfaces import *
from forster import Forster

class Rate(object):
	"""Class for computing energy trasfer rates based on coupling and overlap"""
	def __init__(self, cfg):
		"""Initialize with current configuration"""
		self.cfg = cfg
		return
#-------------------------------------------------------------------------------
	def Boltzmann(self, Rate_GausOverlapFr, InvRate_GausOverlapFr, FragID1, ExStateID1, FragID2, ExStateID2):
		"""Compute Boltzmann factors
		FragID1, ExStateID1 - donor
		"""
		BlzComp = int(self.cfg.MolSys[CFG_MSY_BLZ_COMP])
		BlzTemp = float(self.cfg.MolSys[CFG_MSY_BLZ_TEMP])
		
		#BLZ_COMP_DEFAULT = 0 Do not compute Boltzmann factors
		#BLZ_COMP_EMS_INI = 1 Equliburum populations among emission levels (useful for calculation of initial rates only)
		#BLZ_COMP_ABS_INI = 2 Equliburum populations among absorption levels (useful for calculation of initial rates only)
		#BLZ_COMP_EMS_EQ  = 3 Equliburum populations are _MAINTAINED_ among emission levels (useful for kinetic equaitons)
		#BLZ_COMP_ABS_EQ  = 4 Equliburum populations are _MAINTAINED_ among absorption levels (useful for kinetic equaitons)
		
		if (BlzComp == BLZ_COMP_DEFAULT):
			return (1.0, 1.0)
		RT = GAS_CONST_CM1_K*BlzTemp
		#Boltzmann factor
		P = 0.0
		#Scan excited states
		P_NORM = 0.0 #Normalization factor
		
		SelectedExStateID = -1
		SelectedFragID = -1
		#Donor-based model
		if ((BlzComp == BLZ_COMP_EMS_INI) or (BlzComp == BLZ_COMP_EMS_EQ)):
			SelectedFragID = FragID1
			SelectedExStateID = ExStateID1
		#Acceptor-based model
		elif ((BlzComp == BLZ_COMP_ABS_INI) or (BlzComp == BLZ_COMP_ABS_EQ)):
			SelectedFragID = FragID2
			SelectedExStateID = ExStateID2
		else:
			print("Error: Unknown type of Bolzmann factor calculation ",CFG_MSY_BLZ_COMP," = ",BlzComp)

		#Scan all fragments
		for FragID, Fragment in self.cfg.Fragments.items():
			#Obtain transformed excited states
			TrnsExStates_FragID = self.cfg.Fragments[FragID][CFG_TRNS_EXS_EXSTATE]
			#Filter: states requested by the user in the input
			ExStIDListReq_FragID = ParseList(self.cfg.Fragments[FragID][CFG_FRG_EST])
			#if user requested all ex. states to be processed
			AllExStatesFrag = False
			if ( len(ExStIDListReq_FragID) == 1 ):
				if ( ExStIDListReq_FragID[0] == 0 ):
					AllExStatesFrag = True
			#Scan excited states
			for ExStateID, ExcitedState in TrnsExStates_FragID.items():
				#Emission-based donor factors
				if ((BlzComp == BLZ_COMP_EMS_INI) or (BlzComp == BLZ_COMP_EMS_EQ)):
					print(ExcitedState.Ems_cm1, end=' ')
					P_NORM += math.exp(-1.0*ExcitedState.Ems_cm1/RT)
				#Abosrption-based acceptor factors
				elif ((BlzComp == BLZ_COMP_ABS_INI) or (BlzComp == BLZ_COMP_ABS_EQ)):
					print(ExcitedState.Abs_cm1)
					P_NORM += math.exp(-1.0*ExcitedState.Abs_cm1/RT)
				#Selected state
				if ((SelectedExStateID == ExStateID) and (SelectedFragID == FragID) ):
					print("(selected state)")
					if ((BlzComp == BLZ_COMP_EMS_INI) or (BlzComp == BLZ_COMP_EMS_EQ)):
						P = math.exp(-1.0*ExcitedState.Ems_cm1/RT)
					#Abosrption-based factors
					elif ((BlzComp == BLZ_COMP_ABS_INI) or (BlzComp == BLZ_COMP_ABS_EQ)):
						P = math.exp(-1.0*ExcitedState.Abs_cm1/RT)
				else:
					print()
		if (P_NORM == 0.0):
			print("Bolzmann factor calculation error. Normalization is zero")
			exit(-1)

		print("Calculation of Boltzmann factors at T = ",BlzTemp,"K")
		if (BlzComp == BLZ_COMP_EMS_INI):
			print("Emission-based donor model of initial rates:")
		elif (BlzComp == BLZ_COMP_ABS_INI):
			print("Absorption-based acceptor model of initial rates:")
		elif (BlzComp == BLZ_COMP_EMS_EQ):
			print("Emission-based donor model of equilibirum maintaing rates:")
		elif (BlzComp == BLZ_COMP_ABS_EQ):
			print("Absorption-based acceptor model of equilibirum maintaing rates:")
		else:
			print("Error: Unknown type of Bolzmann factor calculation ",CFG_MSY_BLZ_COMP," = ",BlzComp)
		
		print("Unnormalized factor:",P)
		print("Normalization factor: ",P_NORM)
		P /= P_NORM
		print("Normalized Boltzmann factor:",P)
		
		#Rate correction factor alpha for equilibum maintaining rates only
		alpha = 1.0
		#Equilibum maintaining rates
		if ((BlzComp == BLZ_COMP_EMS_EQ) or (BlzComp == BLZ_COMP_ABS_EQ) ):
			K0 = Rate_GausOverlapFr/InvRate_GausOverlapFr
			print("r(A,B)/r(B,A) = ",K0)
			alpha = math.sqrt(1.0/K0)
			print("Rate correction fctor: alpha(A,B)= sqrt(r(B,A)/r(A,B)): ",alpha)
		
		return (P, alpha)
#-------------------------------------------------------------------------------
	def Rate(self, CouplingSummary, OverlapSummary, InvCouplingSummary, InvOverlapSummary, FragID1, ExStateID1, ExcitedState1, FragID2, ExStateID2, ExcitedState2):
		
		"""Calculates energy trasfer rate based on the Forster theory
		
		Input: 
		Coupling - electronic coupling in Hartree (used only in the legacy code)
		Overlap - spectral overlap 
		Inverted coupling and overlap are used to calculate Boltzmann factors that satisfy the detailed eqilibrium condition
		cfg - configuration with all excited states
		FragID1 - donor fragment
		ExStateID1 - donor excited state index
		ExcitedState1 - donor excited state 
		FragID2 - acceptor fragment
		ExStateID2 - acceptor excited state index
		ExcitedState2 - acceptor excited state
		
		Output:
		Rate in s-1 
		
		"""
		
		RateType= self.cfg.Methods[CFG_MET_RAT]
		#Only for  overlap based rates, otherwise return empty rates
		if ( RateType != CFG_MET_RAT_OVL):
			return RateInterface()
		
		Rate_s1 = 0.0
		Inv_Rate_s1 = 0.0
		
		BlzComp = int(self.cfg.MolSys[CFG_MSY_BLZ_COMP])
		
		if (LEGACY_LEVEL == LL_NONE): #Default	
			#Run rate calculation based on the Forster theory. Coupling is not used. Only Screening Information and donor excited state are needed
			#1st state is a donor
			Rate_s1 = Forster(self.cfg).Rate(CouplingSummary, OverlapSummary.Overlap_M1cm1nm4, ExcitedState1)
			#2nd state is a donor (inverted rate for Boltzmann equilibirum calculations)
			if (BlzComp != BLZ_COMP_DEFAULT):
				Inv_Rate_s1 = Forster(self.cfg).Rate(InvCouplingSummary, InvOverlapSummary.Overlap_M1cm1nm4, ExcitedState2)

		elif (LEGACY_LEVEL == LL_OVERLAP4): #Legacy code
			#Units of rate are s-1
			Rate_s1 = self.RateFromOverlap(CouplingSummary.ScreenedCoupl, OverlapSummary.Overlap_cm)
			#Inverted rate for calculation of Boltzmann factors that satisfy the detailed equilibirum condition
						#2nd state is a donor (inverted rate for Boltzmann equilibirum calculations)
			if (BlzComp != BLZ_COMP_DEFAULT):
				Inv_Rate_s1 = self.RateFromOverlap(InvCouplingSummary.ScreenedCoupl, InvOverlapSummary.Overlap_cm)
		else:
			print("Error: Unknown legacy level: LEGACY_LEVEL = ",LEGACY_LEVEL)
			exit(-1)

		(BoltzmannFactor, AlphaCorrection) = self.Boltzmann(Rate_s1, Inv_Rate_s1, FragID1, ExStateID1, FragID2, ExStateID2)

		return RateInterface(Rate_s1, BoltzmannFactor, AlphaCorrection)
#-------------------------------------------------------------------------------
	def RateFromOverlap(self, Coupling_Hartree, Overlap_cm):
		"""Calculates energy trasfer rate based on Forster theory  - old version
		
		Input: 
		Coupling_Hartree - electronic coupling in Hartree
		Overlap_cm - spectral overlap in inverse wavenumbers 1/(cm-1) = cm 
		
		Output:
		Rate in s-1 
		"""
		
		h_bar_au_sq = 1.0 #Squared Planck's constant divided by 2 pi in atomic units 
		#c_au = 137.035999139 - speed of light in atomic units
		
		#convert overlap from cm to atomic units of length (Bohrs)
		Overlap_Bohr = Overlap_cm/(1.0E2*Bohr_Radius_m)
		Rate_atom_units = (1.0/(h_bar_au_sq*c_au))*((Coupling_Hartree)**2)*Overlap_Bohr
		
		
		#Rate in inverse seconds
		Rate_s1 = Rate_atom_units/TimeAuToS

		"""
		#The code below is tested for correct units conversion (for legacy code with J in cm)
		print "RATE FORMULADIMENSION TEST:"
		overlap_cm = 2.83e-4
		print "overlap_cm:",overlap_cm
		coupling_cm = 80.0
		print "coupling_cm:",coupling_cm
		coupling_hartree = coupling_cm*4.55633E-6
		print "coupling_hartree:",coupling_hartree
		rate_s1= self.RateFromOverlap( coupling_hartree, overlap_cm)
		print "rate_s1:",rate_s1
		time_s = 1.0/rate_s1
		print "time_s:",time_s
		time_ps = time_s/1.0e-12
		print "time_ps:",time_ps
		rate_ps1 = 1.0/time_ps
		print "rate_ps1:",rate_ps1
		exit(-1)
		"""
		
		return (Rate_s1)
#-------------------------------------------------------------------------------