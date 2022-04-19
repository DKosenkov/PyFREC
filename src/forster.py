from couplib.constants import *
from couplib.myreportservice import *
from configuration import *
from math import *
import numpy as np
from interfaces import *
from elscreen import ElScreen

class Forster(object):
	"""Class for Forster coupling and rates"""
#-------------------------------------------------------------------------------
	def __init__(self, cfg):
		"""Initialize with current configuration"""
		self.cfg = cfg
		self.ps = MyReportService()
		self.UserDefinedKappaSq = float(self.cfg.MolSys[CFG_MSY_KAPPA_SQ])
#-------------------------------------------------------------------------------
	def OscillatorsZones(self, State1, n):
		"""Calculates oscillator zones for given donor state"""
		if (State1.Ems_cm1 == 0.0):
			return
		DISTANCE_FORMAT = 2
		Lambda_nm = CM1_NM/State1.Ems_cm1 #Wavelenght of emission in nm
		b = Lambda_nm/(2*math.pi*n)
		print()
		print("Information on oscillator zones:") 
		print("b = Lambda/(2 PI n) =",round(b,DISTANCE_FORMAT),"nm")
		print()
		print("Dexter (contact zone)  0-0.01b (0-",round(b*0.01,DISTANCE_FORMAT),"nm)")
		print("Near field (Forster valid) 0.01b-0.1b (",round(b*0.01,DISTANCE_FORMAT),"-",round(b*0.1,DISTANCE_FORMAT),"nm)")
		print("Intermediate zone 0.1b-10 (",round(0.1*b,DISTANCE_FORMAT),"-",round(10*b,DISTANCE_FORMAT),"nm)")
		print("Radiation zone (far field) >10b (>",round(10*b,DISTANCE_FORMAT),"nm)")
		print()
		
		return
#-------------------------------------------------------------------------------
	def Rate(self, CouplingSummary, Overlap_M1cm1nm4, DonorExcitedState):
		"""Calculates Forster rates directly (without couplings)"""
		print()
		self.ps.PrintDiv()
		print("Direct calculation of Forster rates (electronic couplings are not explicitly used)")
		print()
		
		#format(,'.'+str(FD)+'E')
		FD = 2 #Format for double precision numbers display
		print("Spectral overlap:",format(Overlap_M1cm1nm4,'.'+str(FD)+'E'),"M-1 cm-1 nm^4")
		
		if (Overlap_M1cm1nm4 == 0.0):
			print("The Forster rate will be set to zero...")
			print()
			return 0.0
			
		#Use refractive index if user provided oterwise use  1/sqrt(ElScreening)
		n = float(self.cfg.MolSys[CFG_MSY_RFX])
		if (n != CFG_MSY_RFX_DEFAULT):
			print("User-specified refractive index:",n)
		else:
			#Electrostatic screening factor
			ElScreening = CouplingSummary.ElScreening
			print("Refractive index was not directly provided.")
			print("The electrostatic screening factor: c=", ElScreening," will be used to estimte the refraction index as 1/sqrt(c) =",n)
			n = 1.0/(pow(ElScreening, 0.5))
		
		over_n4 = 1.0/pow(n,4)
		
		
		
		kappa = self.GetOrientationFactor(CouplingSummary.OriFact)

		#Donor-acceptor distance in Angstrom
		#RAng = np.asscalar(np.linalg.norm(CouplingSummary.R))
		RAng = np.linalg.norm(CouplingSummary.R).item()
		R_nm = RAng/nm_to_ang
		
		Phi_D = DonorExcitedState.Phi_D
		
		#Quantum yield of the donor 
		if (Phi_D == 0.0):
			Phi_D = 1.0
			print("Warning: Quantum yield of the donor is not specified. Will be assumed:",Phi_D)
		else: 
			print("User-specified quantum yield of the donor:",Phi_D)
				
		
		#Fluorescence lifetime
		FlLifetime_s = DonorExcitedState.FlLifetime_s
		
		if (FlLifetime_s == 0.0):
			FlLifetime_s =  DonorExcitedState.FlLifetime_sb_s
			#Check if Strickler-Berg Lifetime is available:
			if (FlLifetime_s != 0.0):
				print("Estimated Strickler-Berg fluorescence lifetime of the donor:",format(FlLifetime_s,'.'+str(FD)+'E'),"s")
			else:
				print("Warning: Fluorescence lifetime of the donor is zero...")
				print("The Forster rate will be set to zero...")
				return 0.0
				#exit(-1)
		else:
			print("User-specified fluorescence lifetime of the donor:",format(FlLifetime_s,'.'+str(FD)+'E'),"s")
			
		#Calculate Forster radius
		#Pre-factor in calculation of R0 must be 0.02108
		#Factor of 1.0E17 is due to unit conversion of the spectral overlap from nm6 to mol-1 dm3 cm-1 nm4 ?
		PF = (9.0*np.log(10)*1.0E17)/(128*pow(math.pi,5.0)*AVOGADROS_Na)
		PF = pow(PF,1.0/6.0)
		#print("#Pre-factor in calculation of R0 [nm] and J [M-1 cm-1 nm^4] must be 0.02108 = ",PF)
		
		#kappa = math.sqrt(2.0/3.0)
		#Phi_D = 1.0
		#Overlap_M1cm1nm4 = 
		R0_nm = PF*pow(pow(kappa,2)*Phi_D*Overlap_M1cm1nm4*over_n4,1.0/6.0)
		
		#R0_6_nm = (9.0*np.log(10)*pow(kappa,2)*Phi_D*Overlap_M1cm1nm4*over_n4)/(128*pow(math.pi,5)*AVOGADROS_Na)
		#R0_nm = pow(R0_6_nm, 1.0/6.0)
		DISTANCE_FORMAT = 2 
		print("Forster radius:",round(R0_nm,DISTANCE_FORMAT),"nm")
		
		if ( R_nm > 0.0 ):
			print("Donor-acceptor distance:",round(RAng,DISTANCE_FORMAT),"Ang (",round(R_nm,DISTANCE_FORMAT)," nm)")
		else:
			R_nm = R0_nm
			print("Donor-acceptor distance is assumed to be the Forster radius:",round(R_nm,DISTANCE_FORMAT),"nm")

		#Zones estmates
		self.OscillatorsZones(DonorExcitedState, n)
		
		
		#Efficiency of the energy trasnfer (for D-A dimers only!)
		R_nm6 = pow(R_nm, 6.0)
		R0_6_nm = pow(R0_nm, 6.0)
		E = R0_6_nm/(R0_6_nm+R_nm6)
		print("Efficiency of the energy transfer (valid for dimers only) E =",E)
		
		#Rate of EET
		Rate_Fl_s1 = 1.0/FlLifetime_s
		Rate_s1 = Rate_Fl_s1*(R0_6_nm/R_nm6)
		print()
		print("Fluorescence rate:", format(Rate_Fl_s1,'.'+str(FD)+'E'),"s-1")
		print("Forster (non-radiative) rate:", format(Rate_s1,'.'+str(FD)+'E'),"s-1")
		print()
		#format(,'.'+str(FD)+'E')
	
		return Rate_s1
#-------------------------------------------------------------------------------
	def GetOrientationFactor(self, OriFact):
		"""Get orientation factor either calculated or user-specified"""
		kappa = 0.0

		if ( self.UserDefinedKappaSq != CFG_MSY_KAPPA_SQ_DEFAULT):
			kappa = math.sqrt(self.UserDefinedKappaSq)
			print("User-specified (NOT calculated) orientation factor kappa squared:",self.UserDefinedKappaSq,"( kappa =",kappa,")")
		else:
			kappa = OriFact
			print("Unnormalized orientation factor kappa:",kappa)

		return kappa
#-------------------------------------------------------------------------------
	def Coupling(self, FragID1, State1, Orig1, FragID2, State2, Orig2):
		"""Calculates dipole-dipole coupling usinCompute Forster couplings and form super-molecular Hamiltonian...g Forster theory"""
		#print "Test for checking forster:"
		
		tdm1 = np.array([State1.x, State1.y, State1.z])
		tdm2 = np.array([State2.x, State2.y, State2.z])
		#print "tdm1:",tdm1
		#print "tdm2:",tdm2

		o1 = np.array([Orig1.x, Orig1.y, Orig1.z])
		o2 = np.array([Orig2.x, Orig2.y, Orig2.z])
		#Angstroms
		RAng = np.subtract(o2, o1)
		
		print("R, Distance vector:")
		print(RAng)
		#print "R_ang:",RAng
		#Convert to Bohr
		R = RAng*ATOB
		
		#if (np.asscalar(np.linalg.norm(R))== 0.0):
		if (np.linalg.norm(R).item()== 0.0):
			print(("Error: Interfragment distance is zero for "+CFG_SEC_FRG+str(FragID1+1)+" and "+CFG_SEC_FRG+str(FragID2+1)))
			exit(-1)

		#Input:
		#tdm1 - transitoin dipole of donor (a.u.)
		#tdm2 - transitoin dipole of acceptor (a.u.)
		#R   - distance between donor and acceptor (Bohr)

		#Output:
		#AMuAMuD - product of absolute values of transtion dipole moments in a.u.^2
		#OriFact - orientation factor (-2...2) unitless
		#Coupl   - Forster coupling in Hartree

		#R1 = np.asscalar(np.linalg.norm(R))
		R1 = np.linalg.norm(R).item()

		if (R1 == 0.0):
			print("Error: interfragment distance is zero")
			exit(-1)
		#Distance factor
		DistFact = 1.0/R1**3

		#AMu1 = np.asscalar(np.linalg.norm(tdm1))
		#AMu2 = np.asscalar(np.linalg.norm(tdm2))
		AMu1 = np.linalg.norm(tdm1).item()
		AMu2 = np.linalg.norm(tdm2).item()

		AMuAMuD = AMu1*AMu2
		OriFact = 0.0
		Coupl = 0.0
		K = 0.0
		OriPercent = 0.0
		ScreenedCoupl = 0.0
		ElScreening = 1.0
		
		#Cannot detemine anything for zero transtion dipole moments (warnings are given before in ExStatesReader():)
		if (AMuAMuD != 0.0):
			#Normalized transition dipole moments and radius-vecotor of the interfragment distance:
			nMuA = tdm1/AMu1
			nMuD = tdm2/AMu2
			nR = R/R1
			#print "nMuA:",nMuA
			#print "nMuD:",nMuD
			#print "nR:",nR
			#Orientation factor
			OriFact = np.dot(nMuD,nMuA)-3.0*np.dot(nMuD,nR)*np.dot(nMuA,nR)
			
			OriFact = self.GetOrientationFactor(OriFact)
			
			#if (self.cfg.Methods[CFG_MET_VERBOSITY] < CFG_MET_VERB_TRJ):
			self.ps.PrintDiv()
			
			tag1 = self.cfg.Fragments[FragID1][CFG_FRG_NAM]+str(self.cfg.Fragments[FragID1][CFG_FRG_ID])
			tag2 = self.cfg.Fragments[FragID2][CFG_FRG_NAM]+str(self.cfg.Fragments[FragID2][CFG_FRG_ID])
			
			#Oriention factor decompositon works only if the factor is calculated (not user specified)
			if ( self.UserDefinedKappaSq == CFG_MSY_KAPPA_SQ_DEFAULT):
				#if (self.cfg.Methods[CFG_MET_VERBOSITY] < CFG_MET_VERB_TRJ):
				print(("Orientation Factor Decomposition for: "+str(tag1).strip()+" and "+str(tag2).strip()))
			
				STR_LEN = 16
			
				#if (self.cfg.Methods[CFG_MET_VERBOSITY] < CFG_MET_VERB_TRJ):
				print((str("{} {}").format(str("dot(nMuDnMuA) =").ljust(STR_LEN), np.dot(nMuD,nMuA))))
				print((str("{} {}").format(str("dot(nMuD,nR)  =").ljust(STR_LEN), np.dot(nMuD,nR))))
				print((str("{} {}").format(str("dot(nMuA,nR)  =").ljust(STR_LEN), np.dot(nMuA,nR))))
				print((str("{} {}").format(str("3.0*dot(nMuD,nR)*dot(nMuA,nR) ="), 3.0*np.dot(nMuD,nR)*np.dot(nMuA,nR))))
				print((str("{} {}").format(str("np.dot(nMuD,nMuA)-3.0*np.dot(nMuD,nR)*np.dot(nMuA,nR) ="), OriFact)))
				self.ps.PrintDiv()
			
			#Normalization of the orientation factor
			OriPercent = 100.0*OriFact/2.0
			#Distance-independent factor:
			K = AMuAMuD*OriFact
			#Forster coupling in Hartree
			Coupl = DistFact*AMuAMuD*OriFact
			#Electrostatic screening 
			ElScreening = ElScreen(self.cfg).GetElScreenDamp(R1)
			#Screened coupling in Hartree
			ScreenedCoupl = ElScreening*Coupl
		return CouplingInterface(State1, State2, RAng, AMuAMuD, OriFact, OriPercent, K, Coupl, ElScreening, ScreenedCoupl)
#-------------------------------------------------------------------------------