from couplib.myreportservice import *
from couplib.constants import *
from configuration import *
from math import *
import numpy as np

#-------------------------------------------------------------------------------
class AtomInterface():
	"""Interface class for the ionfromation about atoms (primarely read from PDB file)"""

	def __init__(self, x, y, z, ElSym, AtomPDBID = ATOM_PDB_ID_NONE, AtomName = ATOM_PDB_NAME_NONE, AltLoc = ATOM_PDB_ALT_LOC_NONE, ResID = ATOM_PDB_RES_ID_NONE, ResName=ATOM_PDB_RES_NAME_NONE, ChainID=ATOM_PDB_CHAIN_ID_NONE):
		"""Initialize with the PDB properties"""
		#Cartesian coordinates of the atom in Angstroms
		self.x = x
		self.y = y
		self.z = z
		self.ElSym = ElSym #Element symbol from the PDB file
		self.AtomPDBID = AtomPDBID #Atom PDB ID (unique)
		self.AtomName  = AtomName #Atom PDB name
		self.AltLoc    = AltLoc #Alternate location indicator
		self.ResID     = ResID # Atom PDB residue  ID (same for all atoms of given fragment)
		self.ResName   = ResName #Atom PDB residue  name (same for all atoms of given fragment)
		self.ChainID   = ChainID
		return

	def MyPrint(self, ID = -1, ALabelLen = STR_LEN_ALABEL, Round_XYZ = INT_ROUND, XYZ_Len = STR_LEN_FLOAT):
		"""Print atom"""
		if (ID == -1 ):
			ID = self.AtomPDBID
		print("{} {} {} {}".format( str(self.ElSym+str(ID)+str(self.AltLoc)).ljust(ALabelLen), str(round(self.x,Round_XYZ)).ljust(XYZ_Len),
									str(round(self.y,Round_XYZ)).ljust(XYZ_Len), str(round(self.z,Round_XYZ)).ljust(XYZ_Len)))
#-------------------------------------------------------------------------------
class OriginInterface():
	"""Interface class for the information about origins of fragments"""
	def __init__(self, x, y, z):
		
		#Cartesian coordinates of the atom in Angstroms
		self.x = x
		self.y = y
		self.z = z
		return

	def GetNPArray(self):
		"""Return NumPy Array"""
		return np.asarray([self.x, self.y, self.z])

	def MyPrint(self):
		"""Print origin"""
		print("({} {} {})".format(str(round(self.x,INT_ROUND)).ljust(STR_LEN_FLOAT),
					  str(round(self.y,INT_ROUND)).ljust(STR_LEN_FLOAT),
					  str(round(self.z,INT_ROUND)).ljust(STR_LEN_FLOAT)))
		return
#-------------------------------------------------------------------------------
class QDVibInterface():
	"""Interface class for the information about quamtum dynamics paramerters of fragments"""
	def __init__(self, ExStID, VibModeID, VibModeID_Internal, Vib_cm1, ElVibCoupl_cm1, Vib_Decay_ps1):
		
		self.ExStID = ExStID
		self.VibModeID = VibModeID
		self.VibModeID_Internal = VibModeID_Internal
		self.Vib_cm1 = Vib_cm1
		self.ElVibCoupl_cm1 = ElVibCoupl_cm1
		self.Vib_Decay_ps1 = Vib_Decay_ps1
		return

	def MyPrint(self):
		"""Print quantum dynamics parameters"""
		print("{} {} {} {} {}".format(str(self.ExStID).ljust(STR_LEN_FLOAT),
					      str(self.VibModeID).ljust(STR_LEN_FLOAT),
					  str(round(self.Vib_cm1,INT_ROUND)).ljust(STR_LEN_FLOAT),
					  str(round(self.ElVibCoupl_cm1,INT_ROUND)).ljust(STR_LEN_FLOAT),
					  str(round(self.Vib_Decay_ps1,INT_ROUND)).ljust(STR_LEN_FLOAT)))
		return
#-------------------------------------------------------------------------------
class ExciteStateInterface():
	"""Interface class for the excitate state information"""

	def __init__(self, ExStID, Abs_cm1, x, y, z, Ems_cm1, El_Deph_Rate_ps1, Epsilon_M1cm1,Phi_D,FlLifetime_s, FlLifetime_sb_s):
		"""Initialize with the excited state properties"""
		self.ExStID = ExStID #Excited state id (1-based index, read from an external file)
		self.Abs_cm1 = Abs_cm1 #Absorption maximum, cm-1
		
		if (Abs_cm1 != 0 ):
			self.Abs_nm  = CM1_NM/Abs_cm1 #Absorption maximum, nm
		else:
			self.Abs_nm  = 0.0
			
		#Transition dipole moment components (x,y,z) a.u.
		self.x = x
		self.y = y
		self.z = z
		self.Ems_cm1 = Ems_cm1 #Emission maximum, cm-1
		
		if (Ems_cm1 != 0 ):
			self.Ems_nm  = CM1_NM/Ems_cm1 #Emission maximum, nm
		else:
			self.Ems_nm  = 0.0
		self.El_Deph_Rate_ps1 = El_Deph_Rate_ps1
		self.Epsilon_M1cm1 = Epsilon_M1cm1 #Exctinction coefficient
		self.Phi_D = Phi_D #Fluorescence quantum yield
		self.FlLifetime_s = FlLifetime_s #fluorescnce lifetime from input in s 
		self.FlLifetime_sb_s = FlLifetime_sb_s #Strickler-Berg fluorescnce lifetime in s
		self.Abs_Spec = [] #Absorption spectrum from input
		self.Ems_Spec = [] #Emission spectrum from input
		self.Abs_Spec_nm = [] #Absorption spectrum from input  
		self.Ems_Spec_nm = [] #Emission spectrum from input
		self.Abs_Int_Lim_Low_nm = 0.0 #Lower integration limit of absorption spec. in nm
		self.Abs_Int_Lim_Up_nm  = 0.0 #Upper integration limit of absorption spec. in nm
		self.Ems_Int_Lim_Low_nm = 0.0 #Lower integration limit of emission spec. in nm
		self.Ems_Int_Lim_Up_nm  = 0.0 #Upper integration limit of emission spec. in nm
		self.warning = "Warning: zero transition dipole moment!"
		return
#-------------------------------------------------------------------------------
	def MyPrint(self,JobType=""):
		"""Print excited state properites (does not print quantum dynamics parameters e.g. el. dephasing rate)"""
		Norm = sqrt(self.x**2+self.y**2+self.z**2)
		NormD = Norm*TDM_auToDebye
		Format_CM1 = 2
		Format_nm = 1
		INTX_LEN = 5
		#Special print if calculation of lifetimes only
		if ( JobType == CFG_MET_ARX_LFT):
			print("{}\t{}\t{}\t{}\t{}".format(
															str(self.ExStID).ljust(INTX_LEN),
															str(round(self.Abs_nm,Format_nm)).ljust(STR_LEN_FLOAT),															
															str(round(self.Ems_nm,Format_nm)).ljust(STR_LEN_FLOAT),															
															str(round(self.Epsilon_M1cm1,Format_CM1)).ljust(STR_LEN_FLOAT),
															str(round(self.Phi_D,Format_CM1)).ljust(STR_LEN_FLOAT)))
		else:
			print("{} {} ({} {} {}) {} {} {}".format( str(self.ExStID).ljust(INTX_LEN),str(round(self.Abs_nm,Format_nm)).ljust(STR_LEN_FLOAT),
											str(round(self.x,INT_ROUND)).ljust(STR_LEN_FLOAT),
											str(round(self.y,INT_ROUND)).ljust(STR_LEN_FLOAT),
											str(round(self.z,INT_ROUND)).ljust(STR_LEN_FLOAT),
											str(round(Norm,INT_ROUND)).ljust(STR_LEN_FLOAT),
											str(round(NormD,INT_ROUND)).ljust(STR_LEN_FLOAT),
											str(round(self.Ems_nm,Format_nm)).ljust(STR_LEN_FLOAT)), end=' ')
			if ( Norm == 0.0):
				print(self.warning)
			else:
				print()
			return

	def MyTDMPrint(self):
		"""Print excited state transition dipole momment"""
		Norm = sqrt(self.x**2+self.y**2+self.z**2)
		print("{} ({} {} {}) {}".format( str(self.ExStID).ljust(3),
										str(round(self.x,INT_ROUND)).ljust(STR_LEN_FLOAT),
										str(round(self.y,INT_ROUND)).ljust(STR_LEN_FLOAT),
										str(round(self.z,INT_ROUND)).ljust(STR_LEN_FLOAT), round(Norm,INT_ROUND)), end=' ')
		return
#-------------------------------------------------------------------------------
class CouplingInterface():
	"""Interface class for the excitate state information"""

	def __init__(self, State1 = 0, State2 = 0, R = 0, AMuAMuD = 0, OriFact = 0, OriPercent = 0, K = 0,  Coupl = 0, ElScreening = 0, ScreenedCoupl = 0):
		"""Initialize with the coupling properties"""
		#Excited states under consideration
		self.State1 = State1 
		self.State2 = State2
		self.R = R #interfragment distnce angstrom, A
		self.OriFact = OriFact # Orientation factor (-2...2) unitless
		self.OriPercent = OriPercent # Normalized orientation percent
		self.AMuAMuD = AMuAMuD # Product of absolute values of transtion dipole moments in a.u.^2
		self.K = K # Distance independent factor (Product of absolute values of transtion dipole moments in a.u.^2 and Orientation factor)
		self.Coupl = Coupl #Forster coupling in Hartrees
		self.ElScreening = ElScreening #Electrostatic screening
		self.ScreenedCoupl = ScreenedCoupl #Screened coupilng in Hartrees
		return
#-------------------------------------------------------------------------------
class ResonanceInterface():
	"""Interface class for the information about detected resonances (matching excitation energies of fragments) """

	def __init__(self, Ecm1_1, Ecm1_2, Diff, Overlap, Flag):
		"""Initialize with the resonance properties"""
		#Excitation energies of states under consideration in cm-1
		self.Ecm1_1 = Ecm1_1
		self.Ecm1_2 = Ecm1_2
		self.Diff = Diff #Difference of excitation energies in cm-1
		self.Flag = Flag # true=resonance, false:=no resonance
		return
#-------------------------------------------------------------------------------
class OverlapInterface():
	"""Interface class for the information about overlap of spectra """
	def __init__(self, Overlap_M1cm1nm4 = 0.0, Flag = True):
		"""Initialize with overlaps"""
		self.Overlap_M1cm1nm4 = Overlap_M1cm1nm4 #Spectral overlap
		self.Flag = Flag #true = overlaps are available
		return
#-------------------------------------------------------------------------------
class RateInterface():
	"""Interface class for Forster rates """
	def __init__(self, Rate_s1 = 0.0, BoltzmannFactor = 0.0, AlphaCorrection=0.0):
		"""Initialize with rates"""
		
		self.Rate_s1 = Rate_s1 #Rate in s-1
		self.BoltzmannFactor = BoltzmannFactor
		self.AlphaCorrection = AlphaCorrection
		return
