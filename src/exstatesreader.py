from couplib.myfileservice import MyFileService
from couplib.myreportservice import MyReportService
from couplib.myreportservice import *
from couplib.constants import RE_FLOAT
from couplib.parselist import ParseList

from configuration import *
from interfaces import *

#-------------------------------------------------------------------------------
class ExStatesReader(MyFileService):
	"""Class reads standard orientation geometry, exciation energies, and transition dipole moments from EXSTATE_FILE"""
#-------------------------------------------------------------------------------
	def __init__(self, cfg):
		"""Initialize with current configuration"""
		self.cfg = cfg
		return
#-------------------------------------------------------------------------------
	def VerifyAtoms(self):
		"""Check number of atoms and their order that read from PDB file and excited state (CFG_FRG_GEF = "EXSTATE_FILE") file"""
		#Check number of atoms
		for FragID, Fragment in self.cfg.Fragments.items():

			FileName = Fragment[CFG_FRG_GEF]

			if ( Fragment[CFG_EXS_NATOMS] < 1 ):
				print("Error: There are no atoms in "+CFG_SEC_FRG+str(FragID+1)+" read from "+FileName)
				exit(-1)
			if ( Fragment[CFG_FRG_NATOMS] != Fragment[CFG_EXS_NATOMS] ):
				print("Error: Number of atoms in "+CFG_SEC_FRG+str(FragID+1)+": "+str(Fragment[CFG_FRG_NATOMS]), end=' ')
				print(" does not match the number of atoms: "+str(Fragment[CFG_EXS_NATOMS])+" read from "+FileName)

				if (Fragment[CFG_FRG_NATOMS] < Fragment[CFG_EXS_NATOMS]):
					print("(do not forget to delete any unused capping atoms in "+FileName+")")
				exit(-1)

			#Check atom order
			for AtomID, Atom in Fragment[CFG_FRG_ATOMS].items():
				ElSym = Atom.ElSym
				ElSymExs = self.cfg.Fragments[FragID][CFG_EXS_ATOMS][AtomID].ElSym
				if ( Atom.ElSym.strip().upper() != ElSymExs.strip().upper() ):
					print("Error: Atom "+ElSym+str(AtomID+1)+" in "+CFG_SEC_FRG+str(FragID+1)+" does not match atom ", end=' ')
					print(ElSymExs+str(AtomID+1)+" read from "+FileName)
					exit(-1)
		return
#-------------------------------------------------------------------------------
	def ReadExSGeo(self, FragID, Fragment, FileName, lines):
		"""Read EXSTATE_FILE file geometry section"""
		#Find GEOMETRY (CFG_EXS_SEC_GEO) section
		GeoAtomSectionStart = -1
		for LineCounter in range(0, len(lines)):
			line = lines[LineCounter]
			if ( re.match('^\s*'+CFG_EXS_SEC_GEO+'\s*$', line) ):
				GeoAtomSectionStart = LineCounter
				break

		if ( GeoAtomSectionStart < 0 ):
			print("Error: Cannot find "+CFG_EXS_SEC_GEO+" section in"+FileName+" for "+CFG_SEC_FRG+str(FragID+1))
			exit(-1)

		self.cfg.Fragments[FragID][CFG_EXS_NATOMS] = 0

		#Read GEOMETRY (CFG_EXS_SEC_GEO) section
		for line in lines[(GeoAtomSectionStart+1):]:
			#O   0.718236   0.718236   -4.951050
			matched = re.match('^\s*(\w+)\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s*$', line)
			if (matched):
				ElSymExs = str(matched.group(1)).strip().upper()
				x = float(matched.group(2))
				y = float(matched.group(3))
				z = float(matched.group(4))
				#print ElSymExs," ",x," ",y," ",z
				self.cfg.Fragments[FragID][CFG_EXS_ATOMS][self.cfg.Fragments[FragID][CFG_EXS_NATOMS]] = AtomInterface(x, y, z, ElSymExs)
				self.cfg.Fragments[FragID][CFG_EXS_NATOMS] += 1
			else:
				#End of section
				break
		return
#-------------------------------------------------------------------------------
	def ReadExS(self, FragID, Fragment, FileName, lines):
		""""Read excited states"""
		#Find EXCITED_STATES (CFG_EXS_SEC_ETD) section
		ExSectionStart = -1
		for LineCounter in range(0, len(lines)):
			line = lines[LineCounter]
			if ( re.match('^\s*'+CFG_EXS_SEC_ETD+'\s*$', line) ):
				ExSectionStart = LineCounter
				break
			
		if ( ExSectionStart < 0 ):
			print("Error: Cannot find "+CFG_EXS_SEC_ETD+" section in"+FileName+" for "+CFG_SEC_FRG+str(FragID+1))
			exit(-1)

		self.cfg.Fragments[FragID][CFG_EXS_NEXSTATE] = 0
		# FORMAT OF EXCITED STATES (values in the list):
		# 1.   Excited state number
		# 2.   Absorption maximum, nm*
		# 3,4,5 Transtion Dipole Moments a.u. x,y, and z components*
		# 6.   Emission maximum, nm*
		# 7.   Molar absorption (extinction) coefficient at absorption maximum, M-1 cm-1*
		# 8.  Quantum yield of the donor (unitless)*
		# 9.  Lifetime of the donor excited state in the abscense of the acceptor in seconds*
		#-------
		#*Put 0.0 if the value is unknown 
		
		#For quantum dyanamics simulations additional value(s) is (are) required:
		# 10.   Electronic depahsing rate, ps-1**

		for line in lines[(ExSectionStart+1):]:
			#No quantum dynamics, default behavior
			if (self.cfg.Methods[CFG_MET_QD] == CFG_MET_QD_NONE ):
				matched = re.match('^\s*(\d+)\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_DOUBLE+')\s+('+RE_FLOAT+')\s+('+RE_DOUBLE+')\s*$', line)
			else:
				matched = re.match('^\s*(\d+)\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_DOUBLE+')\s+('+RE_FLOAT+')\s+('+RE_DOUBLE+')\s+('+RE_FLOAT+')\s*$', line)
			#Additional parameters for quantum dynamics
			if (matched):
				ExStID = int(matched.group(1)) #Excited state ID
				Abs_nm = float(matched.group(2)) #Abosrption maximum, nm
				if (Abs_nm == 0.0):
					Abs_cm1	= 0.0
				else:
					Abs_cm1 = CM1_NM/Abs_nm

				#Trinsition dipole moment components
				x = float(matched.group(3))
				y = float(matched.group(4))
				z = float(matched.group(5))
				#print ExStID," ",Eev," ",x," ",y," ",z
				Ems_nm = float(matched.group(6)) #Emission maximum, nm=
				if ( Ems_nm == 0.0):
					Ems_cm1 = 0.0
				else:
					Ems_cm1 = CM1_NM/Ems_nm

				Alpha_M1cm1 = float(matched.group(7)) #Molar absorption (extinction) coefficient at absorption maximum, M-1 cm-1*
				#Quantum yield of the donor (unitless)
				Phi_D = float(matched.group(8)) #Molar absorption (extinction) coefficient at absorption maximum, M-1 cm-1* 
				FlLifetime_s = float(matched.group(9))#Lifetime of the donor excited state in the abscense of the acceptor in seconds
				#Strickler-Berg Lifetime of the donor excited state (will be calculated later)
				FlLifetime_sb_s = 0.0
				#Additional parameter for quamtum dynamics simulations
				El_Deph_Rate_ps1 = 0.0
				if (self.cfg.Methods[CFG_MET_QD] != CFG_MET_QD_NONE ):
					El_Deph_Rate_ps1 = float(matched.group(10)) #Electronic dephasing rate

				self.cfg.Fragments[FragID][CFG_EXS_EXSTATE][self.cfg.Fragments[FragID][CFG_EXS_NEXSTATE]] = ExciteStateInterface(ExStID, Abs_cm1, x, y, z, Ems_cm1, El_Deph_Rate_ps1,Alpha_M1cm1,Phi_D,FlLifetime_s,FlLifetime_sb_s)
				self.cfg.Fragments[FragID][CFG_EXS_NEXSTATE] += 1
				#save excited states
			else:
				#End of section
				break
		return
#-------------------------------------------------------------------------------
	def ReadCenter(self, FragID, Fragment, FileName, lines):
		""""Read coordinates of the fragment center"""
		#Find CENTER
		ExSectionStart = -1 
		for LineCounter in range(0, len(lines)):
			line = lines[LineCounter]
			if ( re.match('^\s*'+CFG_EXS_SEC_CNT+'\s*$', line) ):
				ExSectionStart = LineCounter
				break

		if ( ExSectionStart < 0 ):
			print("Error: Cannot find "+CFG_EXS_SEC_CNT+" section in"+FileName+" for "+CFG_SEC_FRG+str(FragID+1))
			exit(-1)

		#Read EXCITED_STATES (CFG_EXS_SEC_ETD) section
		for line in lines[(ExSectionStart+1):]:
			matched = re.match('^\s*('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s*$', line)
			if (matched):
				#Coordinates of the center:
				x = float(matched.group(1))
				y = float(matched.group(2))
				z = float(matched.group(3))
				self.cfg.Fragments[FragID][CFG_EXS_CENTER] = OriginInterface(x, y, z)
				break
			else:
				#End of section
				print("Error: Cannot find coordinates of origin in "+CFG_EXS_SEC_CNT+" section in"+FileName+" for "+CFG_SEC_FRG+str(FragID+1))
				exit(-1)
				break
		return
#-------------------------------------------------------------------------------
	def ReadSpectra(self, FragID, Fragment, FileName, lines, Section_Name):
		"""Reads tabulated absportion and emission spectra """

		FragTag = Fragment[CFG_FRG_NAM]+str(Fragment[CFG_FRG_ID])
		
		MyLines = lines
		State_ID = -1
		Int_Lim_Low_nm = 0.0
		Int_Lim_Up_nm  = 0.0
		Spec = []
		Spec_nm = []

		#Find start of the section
		ExSectionStart = -1 
		for LineCounter in range(0, len(lines)):
			line = lines[LineCounter]
			matched = re.match('^\s*'+Section_Name+'\s+(\d+)\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s*$', line)
			if ( matched ):
				State_ID = int(matched.group(1))
				Int_Lim_Low_nm = int(matched.group(2))
				Int_Lim_Up_nm  = int(matched.group(3))
				ExSectionStart = LineCounter
				break
		#No optional absorption spectrum was found
		if ( ExSectionStart < 0 ):
			return

		FirstLine = lines[(ExSectionStart+1)]
		matched = re.match('^\s*('+RE_DOUBLE+')\s+('+RE_DOUBLE+')\s*$', FirstLine)
		#if the first line does not contain data try to interpret it as a file name
		if (not matched):
			SpecFileName = str(FirstLine).strip()
			if not self.IsValidFile(SpecFileName,"Cannot open spectrum file: "+SpecFileName):
				exit(-1)
			#Open and read EXSTATE_FILE file content
			FILE = open(SpecFileName,'r')
			MyLines = FILE.readlines()
			FILE.close()
			ExSectionStart = 0

		for line in MyLines[(ExSectionStart+1):]:
			matched = re.match('^\s*('+RE_DOUBLE+')\s+('+RE_DOUBLE+')\s*$', line)
			if (matched):
				#Wavelengths  (nm) Absorbance (M-1 cm-1)
				Lambda_nm = float(matched.group(1))
				Wavenum_cm1 = CM1_NM/Lambda_nm;
				Value  = float(matched.group(2))
				#Read spectrum
				Spec.append((Wavenum_cm1, Value))
				Spec_nm.append((Lambda_nm, Value)) 
			else:
				#End of section
				break

		#Reverse the list of nm and cm-1 are inverse of each other
		Spec = Spec[::-1] 
		Spectrum_Recored_Flag = False
		#Find maximum
		(Max_cm1, Max_Val) = self.GetSpecMax(Spec)
		#Normalize spectrum
		Spec = self.SpecNorm(Spec, Max_Val)
		Format_CM1 = 2
		Format_nm  = 1
		NFrags = len(self.cfg.Fragments[FragID][CFG_EXS_EXSTATE])
		if (NFrags < 1):
			FragTag = Fragment[CFG_FRG_NAM]+str(Fragment[CFG_FRG_ID])
			print("Warning: Fragment ",FragTag," has no excited states. Check format of the ",CFG_EXS_SEC_ETD,"section")
			
		for ExStateID, ExcitedState, in self.cfg.Fragments[FragID][CFG_EXS_EXSTATE].items():
			#Record spectrum to the coresponding excited state
			if ( ExcitedState.ExStID == State_ID):
				if (Section_Name == CFG_EXS_SEC_ABS):
					ExcitedState.Abs_Int_Lim_Low_nm = Int_Lim_Low_nm
					ExcitedState.Abs_Int_Lim_Up_nm  = Int_Lim_Up_nm
					NormF = 1.0/ExcitedState.Epsilon_M1cm1
					Spec = self.SpecNorm(Spec, NormF)
					#Use exctintion for normalization
					(Max_cm1, Max_Val) = self.GetSpecMax(Spec)
					ExcitedState.Abs_Spec = Spec
					ExcitedState.Abs_Spec_nm = Spec_nm

					Max_cm1_str = round(Max_cm1,Format_CM1)
					Max_nm_str = round(CM1_NM/Max_cm1,Format_nm)

					print("Fragment:",FragTag,"tabulated ABSORPTION spectrum is available. Maximum:",Max_cm1_str," cm-1 (",Max_nm_str,"nm ) Eps=",Max_Val,"M-1 cm-1; Integration limits:", Int_Lim_Low_nm,"-",Int_Lim_Up_nm,"nm")
					if (ExcitedState.Abs_cm1 == 0.0):
						print("Maximum of the absorption is obtained from the tabulated spectrum.")
						ExcitedState.Abs_cm1 = Max_cm1
						ExcitedState.Abs_nm = CM1_NM/Max_cm1
					Spectrum_Recored_Flag = True

				elif (Section_Name == CFG_EXS_SEC_EMS): #Emission spectrum
					ExcitedState.Ems_Int_Lim_Low_nm = Int_Lim_Low_nm
					ExcitedState.Ems_Int_Lim_Up_nm  = Int_Lim_Up_nm
					#Check for normalization
					(Max_cm1, Max_Val)  = self.GetSpecMax(Spec)
					ExcitedState.Ems_Spec = Spec
					ExcitedState.Ems_Spec_nm = Spec_nm
					Max_cm1_str = round(Max_cm1,Format_CM1)
					Max_nm_str = round(CM1_NM/Max_cm1,Format_nm)

					print("Fragment:",FragTag,"tabulated normalized EMISSION spectrum is available. Maximum:",Max_cm1_str," cm-1 (",Max_nm_str,"nm ) I=",Max_Val," Arb. Units; Integration limits:", Int_Lim_Low_nm,"-",Int_Lim_Up_nm,"nm")
					if (ExcitedState.Ems_cm1 == 0.0):
						print("Maximum of the emission is obtained from the tabulated spectrum.")
						ExcitedState.Ems_cm1 = Max_cm1
						ExcitedState.Ems_nm = CM1_NM/Max_cm1

					Spectrum_Recored_Flag = True
				else:
					print("Error: unknown section name found while reading tabulated spectrum: ",Section_Name)
					exit(-1)
				break
		if (not Spectrum_Recored_Flag):
			print("Error: spectrum for the fragment: ",FragID," excited state: ",State_ID, " was not assigned to any excited state check sections: ",CFG_EXS_SEC_ETD,CFG_EXS_SEC_ABS,"or",CFG_EXS_SEC_EMS)
			exit(-1)
		return
#-------------------------------------------------------------------------------
	def SpecNorm(self, Spec, Max):
		"""" Normalize spectrum """
		SpecNorm = []
		for Wavenum_cm1, Value in Spec:
			SpecNorm.append((Wavenum_cm1, Value/Max))
		return SpecNorm
#-------------------------------------------------------------------------------
	def GetSpecMax(self, Spec):
		"""" Get maximum of the tabulated spectrum """

		Max_Val = -1.0
		Max_cm1 = -1.0
		
		for Wavenum_cm1, Value in Spec:
			if ( Value >= Max_Val):
				Max_Val = Value
				Max_cm1 = Wavenum_cm1

		return (Max_cm1, Max_Val)
#-------------------------------------------------------------------------------
	def ReadQDVib(self, FragID, Fragment, FileName, lines):
		""""Read molecular vibrations and el-vib. coupling for quantum dynamics"""
		#Find QDVIB section
		ExSectionStart = -1
		self.cfg.Fragments[FragID][CFG_EXS_NQDVID] = 0
		for LineCounter in range(0, len(lines)):
			line = lines[LineCounter]
			if ( re.match('^\s*'+CFG_EXS_SEC_QDV+'\s*$', line) ):
				ExSectionStart = LineCounter
				break

		if ( ExSectionStart < 0 ):
			print("Error: Cannot find "+CFG_EXS_SEC_QDV+" section in"+FileName+" for "+CFG_SEC_FRG+str(FragID+1))
			exit(-1)

		#Read QDVIB (CFG_EXS_SEC_QDV) section
		#For quantum dynamics simulations vibrational parameters are specified
		# 1. Moleuclar vibrational mode  number (starting from 1)
		# 2. Excited state nunmber (should match the first column in $EXCITED_STATES section) for which vibrational mode is specified
		# 2. Vibrational mode wavenumber, cm-1
		# 3. Electron-vibrational coupling (between this vib. mode and the excited state specified in the first column), cm-1
		# 4. Vibratioal decay rate, ps-1
		for line in lines[(ExSectionStart+1):]:
			matched = re.match('^\s*(\d+)\s+(\d+)\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s+('+RE_FLOAT+')\s*$', line)
			if (matched):
				#Coordinates of the center:
				ExStID = int(matched.group(1)) #Excited state ID
				VibModeID = int(matched.group(2)) #Vib. mode ID  - not used for indexing for information purposes only
				Vib_cm1 = float(matched.group(3))
				ElVibCoupl_cm1  = float(matched.group(4))
				Vib_Decay_ps1 = float(matched.group(5))
				self.cfg.Fragments[FragID][CFG_EXS_QDVID][self.cfg.Fragments[FragID][CFG_EXS_NQDVID]] = QDVibInterface(ExStID, VibModeID, self.cfg.Fragments[FragID][CFG_EXS_NQDVID], Vib_cm1, ElVibCoupl_cm1, Vib_Decay_ps1)
				self.cfg.Fragments[FragID][CFG_EXS_NQDVID] += 1
				break
			else:
				#End of section
				print("Error: Cannot find coordinates of origin in "+CFG_EXS_SEC_QDV+" section in"+FileName+" for "+CFG_SEC_FRG+str(FragID+1))
				exit(-1)
				break
		return
#-------------------------------------------------------------------------------
	def VerifyExcitedStates(self):
		"""Verify excited states (make sure that requested states were provided in EXSTATE_FILE) """

		for FragID, Fragment in self.cfg.Fragments.items():
			FileName = Fragment[CFG_FRG_GEF]
			#Check sanity of the excited state requested
			if ( int(Fragment[CFG_FRG_EST]) < 0 ):
				print("Error: An excited state number has cannot be negative. Instead given: "+str(Fragment[CFG_FRG_EST]), end=' ')
				print(" for "+CFG_SEC_FRG+str(FragID+1)+" in the configuration file (.ini)")
				exit(-1)

			#At least one excited state is required for each fragment
			if ( int(self.cfg.Fragments[FragID][CFG_EXS_NEXSTATE]) < 1 ):
				print("Error: at least one excited state is required. Instead found: "+str(self.cfg.Fragments[FragID][CFG_EXS_NEXSTATE]), end=' ')
				print(" for "+CFG_SEC_FRG+str(FragID+1)+" in "+FileName)
				exit(-1)

			#Make sure that ALL requested states are provided
			ExStIDListReq = ParseList(self.cfg.Fragments[FragID][CFG_FRG_EST])
			for ExStIDReq in ExStIDListReq:
				StateNotFound = True
				ExcitedStates = self.cfg.Fragments[FragID][CFG_EXS_EXSTATE]
				for ExcitedStateID, ExcitedState in ExcitedStates.items():
					# accept any/all excited state if user requested 0 (any excited state is acceptable)
					if ( (int(ExcitedState.ExStID) == int(ExStIDReq)) or (int(ExStIDReq) == 0) ):
						StateNotFound = False
				if ( StateNotFound ):
					print("Error: Requested in the configuration excited state"+CFG_FRG_EST+": "+str(ExStIDReq), end=' ')
					print(" for "+CFG_SEC_FRG+str(FragID+1)+" was not found in "+FileName)

			if (self.cfg.Methods[CFG_MET_QD] != CFG_MET_QD_NONE ):
				#Make sure that electroic excited state for which vibrational modes specified exist:
				#Scan vibrations
				QDVibs = self.cfg.Fragments[FragID][CFG_EXS_QDVID]
				for QDVibID, QDVib in QDVibs.items():
					#Check if mentioned excited state exist
					StateNotFound = True
					ExcitedStates = self.cfg.Fragments[FragID][CFG_EXS_EXSTATE]
					for ExcitedStateID, ExcitedState in ExcitedStates.items():
						if ( int(QDVib.ExStID) == int(ExcitedState.ExStID)):
							StateNotFound = False
					if ( StateNotFound ):
						print("Error: Quantum Dyanamics El./vib. modes requested in the configuration electronic excited state"+CFG_EXS_SEC_QDV+": "+str(QDVib.ExStID), end=' ')
						print(" for "+CFG_SEC_FRG+str(FragID+1)+" was not found in "+FileName)
						exit(-1)
		return
#-------------------------------------------------------------------------------
	def AssignTargetExcitedState(self, FragID, Fragment, FileName):
		"""Assign target excited state"""

		ExcitedStates = self.cfg.Fragments[FragID][CFG_EXS_EXSTATE]
		for ExcitedStateID, ExcitedState in ExcitedStates.items():
			if ( int(ExcitedState.ExStID) == int(Fragment[CFG_FRG_EST]) ):
				#Use our internal ID not number provided by the user
				self.cfg.Fragments[FragID][CFG_EXS_EXSTATE_ID] = ExcitedStateID
				break
		return
#-------------------------------------------------------------------------------
	def Read(self):
		"""Read EXSTATE_FILE file"""
		#Scan fragments
		for FragID, Fragment in self.cfg.Fragments.items():
			#Open EXSTATE_FILE file
			FileName = Fragment[CFG_FRG_GEF]
			if not self.IsValidFile(FileName,"Cannot find the EXSTATE_FILE file for "+CFG_SEC_FRG+str(FragID+1)):
				exit(-1)
			#Open and read EXSTATE_FILE file content
			FILE = open(FileName,'r')
			lines = FILE.readlines()
			FILE.close()

			#Do not read geometry and coordinates of the center if only lifetimes are requested
			if ( self.cfg.Methods[CFG_MET_APX] != CFG_MET_ARX_LFT):
				self.ReadExSGeo(FragID, Fragment, FileName, lines)
				self.ReadCenter(FragID, Fragment, FileName, lines)

			self.ReadExS(FragID, Fragment, FileName, lines)
			print()
			#Read tabulated absorption spectrum if provided (wavelengths in nm and Absorbance in M-1 cm-1 )
			self.ReadSpectra(FragID, Fragment, FileName, lines, CFG_EXS_SEC_ABS)
			#Read tabulated emission spectrum if provided (Wavelengths  in nm and Intensity in arbitrary units) with automatic normalization
			self.ReadSpectra(FragID, Fragment, FileName, lines, CFG_EXS_SEC_EMS)
			
			#Read this section only if quantum dynamics is requested
			if ( self.cfg.Methods[CFG_MET_QD] != CFG_MET_QD_NONE):
				self.ReadQDVib(FragID, Fragment, FileName, lines)
		
		#Do not check number of atoms if only lifetimes are calculated
		if ( self.cfg.Methods[CFG_MET_APX] != CFG_MET_ARX_LFT):
			self.VerifyAtoms()
		
		self.VerifyExcitedStates()
		return
#-------------------------------------------------------------------------------
	def AssignTargetExcitedStates(self):
		"""Assign which ex. state will be used for Forster calculation  (does not apply to survey calc.)"""
		for FragID, Fragment in self.cfg.Fragments.items():
			self.AssignTargetExcitedState(FragID, Fragment, Fragment[CFG_FRG_GEF])
#-------------------------------------------------------------------------------
	def PrintFragmens(self):
		"""Prints EXSTATE_FILE file data (geometries and excited states)"""

		self.ps = MyReportService()
		#Scan fragments
		for FragID, Fragment in self.cfg.Fragments.items():
			FileName = Fragment[CFG_FRG_GEF]
			self.ps.PrintSec(CFG_SEC_FRG+str(FragID+1)+" "+CFG_FRG_NAM+" = "+Fragment[CFG_FRG_NAM]+" "+CFG_FRG_ID+" = "+str(Fragment[CFG_FRG_ID]))
			print(CFG_FRG_GEF+" = "+FileName)

			#Do not print coordinates and origin only if lifetimes are calculated
			if ( self.cfg.Methods[CFG_MET_APX] != CFG_MET_ARX_LFT):
				for AtomID, Atom in Fragment[CFG_EXS_ATOMS].items():
					Atom.MyPrint(AtomID+1)
					#print str(ElSymExs)+str(AtomID+1)+" "+str(x)+" "+str(y)+" "+str(z)
				self.ps.PrintDiv()
				print("Coordinates of origin:")
				print(self.cfg.Fragments[FragID][CFG_EXS_CENTER].MyPrint())
				self.ps.PrintDiv()

			print("Excited states read from the input file:")
			INTX_LEN = 5
			#Print the information only relevant to the calculation of lifetimes
			if ( self.cfg.Methods[CFG_MET_APX] == CFG_MET_ARX_LFT):
				print("{}\t{}\t{}\t{}\t{}".format("State".ljust(INTX_LEN),
													  "Abs.Max. nm".ljust(STR_LEN_FLOAT),
													  "Ems.Max nm".ljust(STR_LEN_FLOAT),
													  "Eps M-1cm-1".ljust(STR_LEN_FLOAT),
													  "Q".ljust(STR_LEN_FLOAT)))
			else:
				print("{} {} ({} {} {}) {} {} {}".format("State".ljust(INTX_LEN),
											"Abs. Max, nm".ljust(STR_LEN_FLOAT),
											"x ea0".ljust(STR_LEN_FLOAT),
											"y ea0".ljust(STR_LEN_FLOAT),
											"z ea0".ljust(STR_LEN_FLOAT),
											"Norm ea0".ljust(STR_LEN_FLOAT),
											"Debye".ljust(STR_LEN_FLOAT),
											"Ems. Max, nm".ljust(STR_LEN_FLOAT)))
				#print("State, Abs. Max, nm, Transition Dipole Moment (x,y,z) ea0, abs. ea0, Debye, Ems. Max, nm", end=' ')

				#if (self.cfg.MolSys[CFG_MSY_RTR] != RESONANCE_THR_DEFAULT):
				#	print(", Ems. Max,nm,")
				#else:
				#	print()

			for ExcitedStateID, ExcitedState in Fragment[CFG_EXS_EXSTATE].items():
				ExcitedState.MyPrint(self.cfg.Methods[CFG_MET_APX])

		self.ps.PrintDiv()

		if (self.cfg.Methods[CFG_MET_QD] != CFG_MET_QD_NONE ):
			print("Quantum Dynamics Parameters of Fragments:")
			print("El. Ex. State, Vib. Mode, Vib. cm-1, El.-Vib. Coupl, cm-1, Vib. Decay ps-1")
			QDVibs = self.cfg.Fragments[FragID][CFG_EXS_QDVID]
			for QDVibID, QDVib in QDVibs.items():
				QDVib.MyPrint()
			self.ps.PrintDiv()
		return
#-------------------------------------------------------------------------------