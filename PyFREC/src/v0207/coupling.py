import sys
from couplib.autoarr import AutoArr
from couplib.myreportservice import *
import numpy as np
from configuration import *
from resonances import Resonances
from overlap import Overlap
from rate import Rate
from couplib.constants import *
from math import *
from forster import Forster
from interfaces import *
#-------------------------------------------------------------------------------
class Coupling():
	"""Class computs couplings between fragments"""

	def __init__(self, cfg):
		"""Initialize with current configuration"""
		self.cfg = cfg
		self.CouplRep = AutoArr() #Summary of couplings calculations for 2 ex. states
		self.SurCouplRep = AutoArr() #Summary of couplings calculations for all requested states (survey calculation)
		self.ResonancesRep = AutoArr() #Summary of resonances among excited states
		self.RateRep = AutoArr() #Summary of Rates and Overlaps
		self.ps = MyReportService()
		return
	
	def ProcessExStPair(self, FragID1, ExStateID1, ExcitedState1, Center1, FragID2, ExStateID2, ExcitedState2, Center2):
		"""Process a pair of excited states (couplings, resonances, spectral overlaps, rates etc.)"""
		
		BlzComp = int(self.cfg.MolSys[CFG_MSY_BLZ_COMP])
		
		CouplingSummary = Forster(self.cfg).Coupling(FragID1, ExcitedState1, Center1, FragID2, ExcitedState2, Center2)
		
		#Inverted coupling for detailed equlibirum calculation ( State 1 - donor and State 2 - acceptor are switched)
		if (BlzComp != BLZ_COMP_DEFAULT):
			InvCouplingSummary = Forster(self.cfg).Coupling(FragID2, ExcitedState2, Center2, FragID1, ExcitedState1, Center1)
		
		#Overlaps are only calculated if tabulated spectra have been provided 
		OverlapSummary = Overlap(self.cfg).Overlap(ExcitedState1, ExcitedState2, self.cfg.MolSys[CFG_MSY_RTR])
		#Inverted overlap for detailed equlibirum calculation ( State 1 - donor and State 2 - acceptor are switched)
		if (BlzComp != BLZ_COMP_DEFAULT):
			InvOverlapSummary = Overlap(self.cfg).Overlap(ExcitedState2, ExcitedState1, self.cfg.MolSys[CFG_MSY_RTR])
		
		ResonanceSummary = Resonances().Resonance(ExcitedState1, ExcitedState2, self.cfg.MolSys[CFG_MSY_RES], self.cfg.MolSys[CFG_MSY_RTR], OverlapSummary)
		
		#Calculation of rates
		if (BlzComp != BLZ_COMP_DEFAULT):
			RateSummary = Rate(self.cfg).Rate(CouplingSummary, OverlapSummary, InvCouplingSummary, InvOverlapSummary, FragID1, ExStateID1, ExcitedState1, FragID2, ExStateID2, ExcitedState2)
		else:
			RateSummary = Rate(self.cfg).Rate(CouplingSummary, OverlapSummary, None, None, FragID1, ExStateID1, ExcitedState1, FragID2, ExStateID2, ExcitedState2)
		
		self.RateRep[FragID1][FragID2][ExStateID1][ExStateID2] = (CouplingSummary, OverlapSummary, RateSummary)
		self.SurCouplRep[FragID1][FragID2][ExStateID1][ExStateID2] = CouplingSummary
		self.ResonancesRep[FragID1][FragID2][ExStateID1][ExStateID2] = (ResonanceSummary, CouplingSummary, OverlapSummary)

	def IterateExStates(self, FragID1, Fragment1, FragID2, Fragment2):
		"""Iterate excited states for given pair of fragments"""

		#Obtain transformed excited states
		TrnsExStates_FragID1 = self.cfg.Fragments[FragID1][CFG_TRNS_EXS_EXSTATE]
		TrnsExStates_FragID2 = self.cfg.Fragments[FragID2][CFG_TRNS_EXS_EXSTATE]

		#Filter: states requested by the user in the input
		ExStIDListReq_FragID1 = ParseList(self.cfg.Fragments[FragID1][CFG_FRG_EST])
		ExStIDListReq_FragID2 = ParseList(self.cfg.Fragments[FragID2][CFG_FRG_EST])

		#Transformed origins
		#Orig1 = Fragment1[CFG_TRNS_ORIG]
		#Orig2 = Fragment2[CFG_TRNS_ORIG]
		#Transformed centers of fragments:
		Center1 = Fragment1[CFG_TRNS_CENTER]
		Center2 = Fragment2[CFG_TRNS_CENTER]
		

		#if user requested all ex. states to be processed
		AllExStatesFrag1 = False
		if ( len(ExStIDListReq_FragID1) == 1 ):
			if ( ExStIDListReq_FragID1[0] == 0 ):
				AllExStatesFrag1 = True

		AllExStatesFrag2 = False
		if ( len(ExStIDListReq_FragID2) == 1 ):
			if ( ExStIDListReq_FragID2[0] == 0 ):
				AllExStatesFrag2 = True

		self.SurCouplRep[FragID1][FragID2] = AutoArr()
		self.ResonancesRep[FragID1][FragID2] = AutoArr()
		
		FragTag1 = Fragment1[CFG_FRG_NAM]+str(Fragment1[CFG_FRG_ID])
		FragTag2 = Fragment2[CFG_FRG_NAM]+str(Fragment2[CFG_FRG_ID])
		
		self.ps.PrintDiv()
		#print("Fragments: ",str(FragID1+1),str(FragID2+1))
		print("Donor: ",FragTag1," Acceptor:",FragTag2)

		#Scan excited states
		for ExStateID1, ExcitedState1 in TrnsExStates_FragID1.items():
			self.SurCouplRep[FragID1][FragID2][ExStateID1] = AutoArr()
			self.ResonancesRep[FragID1][FragID2][ExStateID1] = AutoArr()
			for ExStateID2, ExcitedState2 in TrnsExStates_FragID2.items():
				#The list of the excited states accodring to the user's input (not the internal index)
				InputExStID1 = ExcitedState1.ExStID
				InputExStID2 = ExcitedState2.ExStID
				#Filter out states that were not requested by the user
				if    ( ( (int(InputExStID1) in ExStIDListReq_FragID1) or AllExStatesFrag1 ) and
					( (int(InputExStID2) in ExStIDListReq_FragID2) or AllExStatesFrag2 ) ):
					self.ProcessExStPair(FragID1, ExStateID1, ExcitedState1, Center1, FragID2, ExStateID2, ExcitedState2, Center2)

	def RunSurvey(self):
		"""Compute couplings for one pair for all ex. states"""
		
		#Scan fragments
		for FragID1, Fragment1 in self.cfg.Fragments.items():
			self.SurCouplRep[FragID1] = AutoArr()
			self.ResonancesRep[FragID1] = AutoArr()
			for FragID2, Fragment2 in self.cfg.Fragments.items():
				#off-diagonal elements of Hamiltonian
				#if ( FragID1 < FragID2 ): #only unique pairs of fragments are considered (triangular matrix)
				#All non-diagonal elements are needed rate D1->A2 is not equal to the rate D2->A1  (D= Donor, A=Acceptor)
				if ( FragID1 != FragID2 ): #only unique pairs of fragments are considered (triangular matrix)
					self.IterateExStates(FragID1, Fragment1, FragID2, Fragment2)
		return

	def PrintSurvey(self):
		"""Report couplings"""

		self.ps.PrintDiv()

		print("Transformed Excited States for Each Fragment:")
		print()
		self.PrintExStateHeader()
		for FragID, Fragment in self.cfg.Fragments.items():
			#Obtain transformed excited states
			TrnsExStates = self.cfg.Fragments[FragID][CFG_TRNS_EXS_EXSTATE]
			for ExStateID, ExcitedState in TrnsExStates.items():
				self.PrintExState(FragID, ExcitedState)
		print()
		self.ps.PrintDiv()

		res_cm1 = self.cfg.MolSys[CFG_MSY_RES]
		res_ev  = round(float(res_cm1)/eVToCM1,3)
		
		#Resonances based on overlaps
		res_thr = self.cfg.MolSys[CFG_MSY_RTR]

		print("Resonances satisfying condition", end=' ')
		
		if ( res_thr == 0.0 ):
			print(" |Ej-Ei| <= "+str(res_cm1)+" cm-1 or "+str(res_ev)+" eV")
		else:
			print(" overlap(Ej,Ei) >= "+str(res_thr))
		
		self.PrintResonancesHeader(res_thr)
		print()
		i = 0
		
		for FragID1, ResonancesReps in self.ResonancesRep.items():
			for FragID2, ResonancesReps1 in ResonancesReps.items():
				for ExStateID1, ResonancesReps2 in ResonancesReps1.items():
					for ExStateID2, ResonanceCouplingSummary in ResonancesReps2.items():
						i += self.PrintResonanceLine(FragID1, FragID2, ResonanceCouplingSummary, res_thr)

		self.ps.PrintDiv()
		print("Total resonances: ",i)
		print()
		print("Electronic couplings list:")
		print()
		print("Note that the orintation factor percentage is calculated as: Ori, % = 100 x Ori / 2 ")
		print()
		i = 0
		self.PrintCouplingHeader()
		for FragID1, CouplSurReps in self.SurCouplRep.items():
			for FragID2, CouplSurReps1 in CouplSurReps.items():
				for ExStateID1, CouplSurReps2 in CouplSurReps1.items():
					for ExStateID2, CouplingSummary in CouplSurReps2.items():
						i += self.PrintCouplingLine(FragID1, FragID2, CouplingSummary)

		self.ps.PrintDiv()
		print("Total couplings: ",i)
		print()
		print("Electronic couplings list (non-zero only) that satisfy the resonance condition. Pairs of states i,j and j,i are distinguished.", end=' ')
		
		if ( res_thr == 0.0 ):
			print(" |Ej-Ei| <= "+str(res_cm1)+" cm-1 or "+str(res_ev)+" eV")
		else:
			print(" overlap(Ej,Ei) >= "+str(res_thr))
		print()
		i = 0
		self.PrintCouplingHeader()
		for FragID1, ResonancesReps in self.ResonancesRep.items():
			for FragID2, ResonancesReps1 in ResonancesReps.items():
				for ExStateID1, ResonancesReps2 in ResonancesReps1.items():
					for ExStateID2, ResonanceCouplingSummary in ResonancesReps2.items():
						(ResonanceSummary, CouplingSummary, OverlapSummary) = ResonanceCouplingSummary
						if (ResonanceSummary.Flag):
							i += self.PrintCouplingLine(FragID1, FragID2, CouplingSummary)

		self.ps.PrintDiv()
		print("Total couplings that satisfy the resonance condition: ",i)
		print("Pairs of states i,j and j,i are distinguished.")
		
		
		if (self.cfg.Methods[CFG_MET_RAT] == CFG_MET_RAT_OVL ):
			print()
			print("Energy transfer rates:")
			print()
			self.PrintRateHeader()
			i = 0
			for FragID1, RateReps in self.RateRep.items():
				for FragID2, RateReps1 in RateReps.items():
					for ExStateID1, RateReps2 in RateReps1.items():
						for ExStateID2, RepSummary in RateReps2.items():
							(CouplingSummary, OverlapSummary, RateSummary) = RepSummary
							#if (ResonanceSummary.Flag):
							i += self.PrintRateLine(FragID1, FragID2, CouplingSummary, OverlapSummary, RateSummary)
			self.ps.PrintDiv()
			print("Total rates : ",i)
		return

	def Run(self):
		"""Compute couplings for each pair of excited states"""

		self.ps = MyReportService()
		#Number of fragments
		NFrags = len(self.cfg.Fragments)
		H = np.zeros(NFrags*NFrags,dtype=np.dtype('d')).reshape(NFrags,NFrags)
		#Scan fragments
		for FragID1, Fragment1 in self.cfg.Fragments.items():
			self.CouplRep[FragID1] = AutoArr()
			for FragID2, Fragment2 in self.cfg.Fragments.items():
				#off-diagonal elements of Hamiltonian
				if ( FragID1 < FragID2 ): #only unique pairs of fragments are considered (triangular matrix)
					#print "pairs:",FragID1," ",FragID2
					State1 = Fragment1[CFG_TRNS_EXS_EXSTATE][Fragment1[CFG_EXS_EXSTATE_ID]]
					State2 = Fragment2[CFG_TRNS_EXS_EXSTATE][Fragment2[CFG_EXS_EXSTATE_ID]]
					#Transformed origins
					#Orig1 = Fragment1[CFG_TRNS_ORIG]
					#Orig2 = Fragment2[CFG_TRNS_ORIG]
					#Transformed centers:
					Center1 = Fragment1[CFG_TRNS_CENTER]
					Center2 = Fragment2[CFG_TRNS_CENTER]
					#Couplings in Hartree, Distances in Bohrs
					CoupliSummary = Forster(self.cfg).Coupling(FragID1, State1, Center1, FragID2, State2, Center2)
					self.CouplRep[FragID1][FragID2] = CoupliSummary
					#R,A; Ori (-2...2); coupl (Hartree)
					#Ori_percent = (100.0*OriFact/2.0);
					#Compute Hamiltonian
					H[FragID1][FragID2] = CoupliSummary.ScreenedCoupl
					H[FragID2][FragID1] = H[FragID1][FragID2]
				#diagonal elements of Hamiltonian
				elif ( FragID1 == FragID2 ):
					#Target excited state
					#Eev = Fragment1[CFG_TRNS_EXS_EXSTATE][Fragment1[CFG_EXS_EXSTATE_ID]].Eev
					Abs_cm1 = Fragment1[CFG_TRNS_EXS_EXSTATE][Fragment1[CFG_EXS_EXSTATE_ID]].Abs_cm1
					#H[FragID1][FragID2] = Eev/HartreeToeV
					H[FragID1][FragID2] = Abs_cm1/HartreeToCM1
		return H
	
	def PrintRateHeader(self):
		""" Print rate header """

		if (LEGACY_LEVEL == LL_NONE):
			Overlap_Label = "J, M-1 cm-1 nm^4".rjust(STR_LEN_FLOAT)
		else:
			print("Error: Unknown legacy level: LEGACY_LEVEL = ",LEGACY_LEVEL)
			exit(-1)
			
		BlzComp = int(self.cfg.MolSys[CFG_MSY_BLZ_COMP])
		if (BlzComp != BLZ_COMP_DEFAULT):
			print("{} {} {} {} {} {} {} {} {} {} {} {}".format( "Frg i".rjust(STR_LEN_INT),
								"PDB TAG".rjust(STR_LEN_ALABEL),
								"Frg j".rjust(STR_LEN_INT),
								"PDB TAG".rjust(STR_LEN_ALABEL),
								"Ex.St i".rjust(STR_LEN_INT),
								"Ex.St j".rjust(STR_LEN_INT),
								"Scr.E,cm-1".rjust(STR_LEN_FLOAT),
								Overlap_Label,
								"Boltzmann ".rjust(STR_LEN_FLOAT),
								"Alpha".rjust(STR_LEN_FLOAT),
								"k, ps-1".rjust(STR_LEN_FLOAT),
								"t, ps".rjust(STR_LEN_FLOAT)))
		else:
			print("{} {} {} {} {} {} {} {} {} {}".format( "Frg i".rjust(STR_LEN_INT),
								"PDB TAG".rjust(STR_LEN_ALABEL),
								"Frg j".rjust(STR_LEN_INT),
								"PDB TAG".rjust(STR_LEN_ALABEL),
								"Ex.St i".rjust(STR_LEN_INT),
								"Ex.St j".rjust(STR_LEN_INT),
								"Scr.E,cm-1".rjust(STR_LEN_FLOAT),
								Overlap_Label,
								"k, ps-1".rjust(STR_LEN_FLOAT),
								"t, ps".rjust(STR_LEN_FLOAT)))

		return
	
	def PrintRateLine(self, FragID1, FragID2, CouplingSummary, OverlapSummary, RateSummary):
		"""Print rates for given pair of fragments"""
		ExStID1 = CouplingSummary.State1.ExStID
		St1Abs_cm1    = CouplingSummary.State1.Abs_cm1
		ExStID2 = CouplingSummary.State2.ExStID
		St2Abs_cm1    = CouplingSummary.State2.Abs_cm1
		
		i = FragID1+1
		j = FragID2+1
		tag1 = self.cfg.Fragments[FragID1][CFG_FRG_NAM]+str(self.cfg.Fragments[FragID1][CFG_FRG_ID])
		tag2 = self.cfg.Fragments[FragID2][CFG_FRG_NAM]+str(self.cfg.Fragments[FragID2][CFG_FRG_ID])
		
		#Screened coupling
		ScrCoupl = CouplingSummary.ScreenedCoupl
		if ( ScrCoupl == 0.0 ):
			return 0
		ScrCoupl_cm1 = ScrCoupl*HartreeToCM1
		int_round = 3
		
		Overlap = 0.0
		if (LEGACY_LEVEL == LL_NONE):
			Overlap = OverlapSummary.Overlap_M1cm1nm4	
		else:
			print("Error: Unknown legacy level: LEGACY_LEVEL = ",LEGACY_LEVEL)
			exit(-1)
		
		Rate_s1 = RateSummary.Rate_s1
		BoltzmannFactor = RateSummary.BoltzmannFactor
		AlphaCorrection = RateSummary.AlphaCorrection
		Time_s = np.inf
		Time_ps = np.inf
		Rate_ps1 = 0.0
		
		if (Rate_s1 != 0.0):
			Time_s = 1.0/Rate_s1
			Time_ps = Time_s/1.0e-12
			Rate_ps1 = 1.0/Time_ps
		
		l = str(STR_LEN_FLOAT).strip(); #length of a float number
		sf = str(4) #number of sig. figs. for scientific notation
		dm = str(3) #number of digits after the decimal point

		BlzComp = int(self.cfg.MolSys[CFG_MSY_BLZ_COMP])
		if (BlzComp != BLZ_COMP_DEFAULT):
			print(str("{} {} {} {} {} {} {:"+l+"."+sf+"g} {:"+l+"."+sf+"g} {:"+l+"."+sf+"e} {:"+l+"."+sf+"e} {:"+l+"."+sf+"g} {:"+l+"."+sf+"g}").format(
				str(i).rjust(STR_LEN_INT),
				tag1.rjust(STR_LEN_ALABEL),
				str(j).rjust(STR_LEN_INT),
				tag2.rjust(STR_LEN_ALABEL),
				str(ExStID1).rjust(STR_LEN_INT),
				str(ExStID2).rjust(STR_LEN_INT),
				ScrCoupl_cm1,
				Overlap,
				BoltzmannFactor,
				AlphaCorrection,
				Rate_ps1,
				Time_ps))
		else:
				print(str("{} {} {} {} {} {} {:"+l+"."+sf+"g} {:"+l+"."+sf+"g} {:"+l+"."+sf+"g} {:"+l+"."+sf+"g}").format(
				str(i).rjust(STR_LEN_INT),
				tag1.rjust(STR_LEN_ALABEL),
				str(j).rjust(STR_LEN_INT),
				tag2.rjust(STR_LEN_ALABEL),
				str(ExStID1).rjust(STR_LEN_INT),
				str(ExStID2).rjust(STR_LEN_INT),
				ScrCoupl_cm1,
				Overlap,
				Rate_ps1,
				Time_ps))
		return 1
	
	def PrintResonancesHeader(self,res_thr):
		""" Print resonance header"""
		if ( res_thr == 0.0 ):
			res_cond = "Ej-Ei".rjust(STR_LEN_FLOAT)
		else:
			res_cond = "Overlap".rjust(STR_LEN_FLOAT)
		
		print("{} {} {} {} {} {} {} {} {}".format("Frg i".rjust(STR_LEN_INT),
												  "PDB TAG".rjust(STR_LEN_ALABEL),
												  "Frg j".rjust(STR_LEN_INT),
												  "PDB TAG".rjust(STR_LEN_ALABEL),
												  "Ex.St i".rjust(STR_LEN_INT),
												  "Ex.St j".rjust(STR_LEN_INT),
												  "Ei cm-1".rjust(STR_LEN_FLOAT),
												  "Ej cm-1".rjust(STR_LEN_FLOAT),
												  res_cond))
		return

	def PrintResonanceLine(self, FragID1, FragID2, ResonanceCouplingSummary, res_thr):
		""" Print resonances summary for given pair of fragments"""

		(ResonanceSummary, CouplingSummary, OverlapSummary) = ResonanceCouplingSummary
		
		ExStID1 = CouplingSummary.State1.ExStID
		ExStID2 = CouplingSummary.State2.ExStID
		Ecm1_1  = ResonanceSummary.Ecm1_1
		Ecm1_2  = ResonanceSummary.Ecm1_2
		Diff = ResonanceSummary.Diff

		Overlap = 0.0
		if (LEGACY_LEVEL == LL_NONE):
			Overlap = OverlapSummary.Overlap_M1cm1nm4
		elif (LEGACY_LEVEL == LL_OVERLAP4):
			Overlap = OverlapSummary.Overlap_cm
		else:
			print("Error: Unknown legacy level: LEGACY_LEVEL = ",LEGACY_LEVEL)
			exit(-1)

		
		
		Flag = ResonanceSummary.Flag
		
		#No resonance
		if (not Flag):
			return 0

		i = FragID1+1
		j = FragID2+1

		tag1 = self.cfg.Fragments[FragID1][CFG_FRG_NAM]+str(self.cfg.Fragments[FragID1][CFG_FRG_ID])
		tag2 = self.cfg.Fragments[FragID2][CFG_FRG_NAM]+str(self.cfg.Fragments[FragID2][CFG_FRG_ID])

		l = str(STR_LEN_FLOAT).strip(); #length of a float number
		sf = str(4) #number of sig. figs. for scientific notation
		dm = str(0) #number of digits after the decimal point
		
		if (res_thr == 0.0):
			Diff_Overlap = Diff
			dm_res = str(0)
		else:
			Diff_Overlap = Overlap
			dm_res = str(3)
		
		print(str("{} {} {} {} {} {} {:"+l+"."+dm+"f} {:"+l+"."+dm+"f} {:"+l+"."+sf+"g}").format(
				str(i).rjust(STR_LEN_INT),
				tag1.rjust(STR_LEN_ALABEL),
				str(j).rjust(STR_LEN_INT),
				tag2.rjust(STR_LEN_ALABEL),
				str(ExStID1).rjust(STR_LEN_INT),
				str(ExStID2).rjust(STR_LEN_INT),
				Ecm1_1,
				Ecm1_2,
				Diff_Overlap))
		return 1

	def PrintCouplingHeader(self):
		""" Print coupling header"""
		print("{} {} {} {} {} {} {} {} {} {} {} {} {} {}".format("Frg i".rjust(STR_LEN_INT),
														"PDB TAG".rjust(STR_LEN_ALABEL),
														"Frg j".rjust(STR_LEN_INT),
														"PDB TAG".rjust(STR_LEN_ALABEL),
														"Ex.St i".rjust(STR_LEN_INT),
														"Ex.St j".rjust(STR_LEN_INT),
														"Abs. Ej-Ei,cm-1".rjust(STR_LEN_FLOAT),
														"R, A".rjust(STR_LEN_FLOAT),
														"AxB,(ea0)^2".rjust(STR_LEN_FLOAT),
														"Ori.".rjust(STR_LEN_FLOAT),
														"Ori.%".rjust(STR_LEN_FLOAT),
														"K',(ea0)^2".rjust(STR_LEN_FLOAT),
														"E,cm-1".rjust(STR_LEN_FLOAT),
														"Scr.E,cm-1".rjust(STR_LEN_FLOAT)))
		return

	def PrintCouplingLine(self, FragID1, FragID2, CouplingSummary):
		""" Print couplings for given pair of fragments"""
		ExStID1 = CouplingSummary.State1.ExStID
		St1Abs_cm1 = CouplingSummary.State1.Abs_cm1
		ExStID2 = CouplingSummary.State2.ExStID
		St2Abs_cm1 = CouplingSummary.State2.Abs_cm1
		
		EDiff = St2Abs_cm1-St1Abs_cm1
		
		R = CouplingSummary.R
		AMuAMuD = CouplingSummary.AMuAMuD
		OriFact = CouplingSummary.OriFact
		OriPercent = CouplingSummary.OriPercent
		K = CouplingSummary.K
		Coupl = CouplingSummary.Coupl
		ScrCoupl = CouplingSummary.ScreenedCoupl #Screeened coupling
		#Print non-zero couplings only
		if ( Coupl == 0.0 ):
			return 0
		
		i = FragID1+1
		j = FragID2+1

		tag1 = self.cfg.Fragments[FragID1][CFG_FRG_NAM]+str(self.cfg.Fragments[FragID1][CFG_FRG_ID])
		tag2 = self.cfg.Fragments[FragID2][CFG_FRG_NAM]+str(self.cfg.Fragments[FragID2][CFG_FRG_ID])

		Dist = np.linalg.norm(R)
		Coupl_cm1  = Coupl*HartreeToCM1
		ScrCoupl_cm1  = ScrCoupl*HartreeToCM1
		int_round  = 3 

		l = str(STR_LEN_FLOAT).strip(); #length of a float number
		sf = str(4) #number of sig. figs. for scientific notation
		dm = str(3) #number of digits after the decimal point
		print(str("{} {} {} {} {} {} {:"+l+"."+dm+"f} {:"+l+"."+dm+"f} {:"+l+"."+sf+"g} {:"+l+"."+sf+"g} {:"+l+"."+dm+"f} {:"+l+"."+sf+"g} {:"+l+"."+sf+"g} {:"+l+"."+sf+"g}").format(
				str(i).rjust(STR_LEN_INT),
				tag1.rjust(STR_LEN_ALABEL),
				str(j).rjust(STR_LEN_INT),
				tag2.rjust(STR_LEN_ALABEL),
				str(ExStID1).rjust(STR_LEN_INT),
				str(ExStID2).rjust(STR_LEN_INT),
				EDiff,
				Dist,
				AMuAMuD,
				OriFact,
				OriPercent,
				K,
				Coupl_cm1,
				ScrCoupl_cm1))
		return 1

	def ClosestNeighbours(self, FragID1, TryFragID2, CouplingSummary):
		""" Find the closest neighbour to fragment FragID1"""

		FragIDClosest = -1
		Min_Dist = sys.float_info.max

		for FragID2, CouplingSummary in self.CouplRep[FragID1].items():
			Dist = np.linalg.norm(CouplingSummary.R)
			if (Dist <= Min_Dist):
				Min_Dist = Dist
				FragIDClosest = FragID2

		if ( FragIDClosest == -1 ):
			print("Error: cannot find the closest neighbour to fragment "+str(FragID1+1))

		if (FragIDClosest == TryFragID2):
			return True
		else:
			return False

	def StrongestCouplings(self, FragID1, TryFragID2, CouplingSummary):
		""" Find the strongest coupling to fragment FragID1"""

		FragIDStrongest = -1
		Max_Coupling = 0.0

		for FragID2, CoupliSummary in self.CouplRep[FragID1].items():
			if ( abs(CouplingSummary.Coupl) >= Max_Coupling):
				Max_Coupling = abs(CouplingSummary.Coupl)
				FragIDStrongest = FragID2

		if ( FragIDStrongest == -1 ):
			print("Error: cannot find the strongest coupling to fragment "+str(FragID1+1))

		if (FragIDStrongest == TryFragID2):
			return True
		else:
			return False

	def PrintExStateHeader(self):
		print(str("{} {} {} {} {} {} {}").format("PDB TAG".rjust(STR_LEN_ALABEL), "Ex. St.".rjust(STR_LEN_INT),
										str("Tran. Dipole Moment (x,y,z) a.u.").center(3*STR_LEN_FLOAT+4),
										str("Abs. a.u.").rjust(STR_LEN_FLOAT_LARGE),
										str("Abs. Debye").rjust(STR_LEN_FLOAT_LARGE),
										str("Abs. En, nm").rjust(STR_LEN_FLOAT_LARGE),
										str("Ems. nm").rjust(STR_LEN_FLOAT_LARGE)))
		return

	def PrintExState(self, FragID, State):
		tag = self.cfg.Fragments[FragID][CFG_FRG_NAM]+str(self.cfg.Fragments[FragID][CFG_FRG_ID])
		l = str(STR_LEN_FLOAT_LARGE).strip(); #length of a float number (as a string)
		sf = str(4) #number of sig. figs. for scientific notation
		dm_nm = str(0)
		dm_au = str(2)
		ID = State.ExStID
		Abs_nm = State.Abs_nm
		x = State.x
		y = State.y
		z = State.z
		Norm = sqrt(x**2+y**2+z**2)
		NormD = Norm*TDM_auToDebye
		Ems_nm = State.Ems_nm
		print(str("{} {} ({} {} {}) {:"+l+"."+dm_au+"f} {:"+l+"."+dm_au+"f} {:"+l+"."+dm_nm+"f} {:"+l+"."+dm_nm+"f}").format(
																		tag.rjust(STR_LEN_ALABEL),
																		str(ID).rjust(STR_LEN_INT),
																		str(round(x,INT_ROUND)).ljust(STR_LEN_FLOAT),
																		str(round(y,INT_ROUND)).ljust(STR_LEN_FLOAT),
																		str(round(z,INT_ROUND)).ljust(STR_LEN_FLOAT),
																		Norm,
																		NormD,
																		Abs_nm,
																		Ems_nm))
		return

	def PrintCouplings(self):
		"""Report couplings"""

		self.ps.PrintDiv()
		print("Electronic couplings between fragments:")

		print()
		print("Transformed Excited States For Each Fragment:")

		self.PrintExStateHeader()
		for FragID, Fragment in self.cfg.Fragments.items():
			State= Fragment[CFG_TRNS_EXS_EXSTATE][Fragment[CFG_EXS_EXSTATE_ID]]
			self.PrintExState(FragID, State)
		print()
		self.ps.PrintDiv()
		print("Electronic couplings list (non-zero only):")
		print()
		self.PrintCouplingHeader()
		for FragID1, CouplReps in self.CouplRep.items():
			for FragID2, CouplingSummary in CouplReps.items():
				self.PrintCouplingLine(FragID1, FragID2, CouplingSummary)

		if (self.cfg.Methods[CFG_MET_VERBOSITY] < CFG_MET_VERB_TRJ):
			print()
			print("Electronic couplings (non-zero, closest neighbours only):")
			print("Note: only the last fragment among identical distances is printed")
			self.PrintCouplingHeader()
			for FragID1, CouplReps in self.CouplRep.items():
				for FragID2, CouplingSummary in CouplReps.items():
					if ( self.ClosestNeighbours(FragID1,FragID2,CouplingSummary) ):
						self.PrintCouplingLine(FragID1, FragID2, CouplingSummary)
	
			print()
			print("Electronic couplings (non-zero, strongest couplings for each fragment):")
			print("Note: only the last fragment among identical couplings is printed")
			self.PrintCouplingHeader()
			for FragID1, CouplReps in self.CouplRep.items():
				for FragID2, CouplingSummary in CouplReps.items():
					if ( self.StrongestCouplings(FragID1,FragID2,CouplingSummary) ):
						self.PrintCouplingLine(FragID1, FragID2, CouplingSummary)

		self.ps.PrintDiv()
		return
#-------------------------------------------------------------------------------
	def HomotransferProperties(self):
		"""Homotrasfer propoerties that depend only on propeties of one fragment (donor)"""
		self.ps.PrintDiv()
		print()
		print("Homotransfer Properties")
		for FragID, Fragment in list(self.cfg.Fragments.items()):

			print("Fragment:",str(Fragment[CFG_FRG_NAM]).strip()+str(Fragment[CFG_FRG_ID]).strip())

			for ExcitedStateID, ExcitedState in list(Fragment[CFG_EXS_EXSTATE].items()):
				print("State   :",ExcitedState.ExStID)
				#Compute overlap for a single fragment
				OverlapSummary =  None
				#Run "tabulated" overlap
				if ( (len(ExcitedState.Ems_Spec) > 0 ) and len(ExcitedState.Abs_Spec) > 0):
					OverlapSummary = Overlap(self.cfg).TabOverlap(ExcitedState, ExcitedState)
				
				#Empty as there is only one fragment, user provides kappa and n
				CouplingSummary = CouplingInterface()
				
				Rate_s1 = Forster(self.cfg).Rate(CouplingSummary, OverlapSummary.Overlap_M1cm1nm4, ExcitedState)
				
		self.ps.PrintDiv()
		return
#-------------------------------------------------------------------------------