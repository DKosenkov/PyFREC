#External system libraries
import os, sys, inspect
import re

#Printer for debugging purposes
#import pprint

#Command line parser
from   argparse import ArgumentParser
#Configuration (.ini file parser)
import configparser

#Location of our own shared libraries for couplib
STR_SHARED_LIBS = "."

#Add location of shared libraries  (STR_SHARED_LIB) to the system PATH
shared_lib_path = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe()))[0],STR_SHARED_LIBS)))
if shared_lib_path not in sys.path:
	sys.path.insert(0, shared_lib_path)

#Internal shared libraries
from couplib.autoarr import AutoArr
from couplib.myfileservice import MyFileService
from couplib.myreportservice import MyReportService
from couplib.parselist import ParseList

INT_CONFIG_TRJ_PDB_ON = False #Disable trajectory reader
INT_CONFIG_QD_ON  = False #Disable quantum dynamics

CFG_FREC_SAYS = r"FREC SAYS:" #Comment line in the make_fragment output which is easy to recognize later

#Internal keys
LL_NONE      =  1 #Non-comptatability with older versions
LEGACY_LEVEL = LL_NONE #Compatability with previos versions


#Sections of the configuration (ini) file
CFG_SEC_MET = r"METHODS"
CFG_SEC_MSY = r"MOLECULAR_SYSTEM"
CFG_SEC_QDY = r"QUANTUM_DYNAMICS"
CFG_SEC_FRG = r"FRAGMENT"
#Note that CFG_SEC_FRG (e.g. FRAGMENT1, FRAGMENT2, ... ) are EXCLUDED regular expressions are used instead 
CFG_SEC_LST_REQ = (CFG_SEC_MET, CFG_SEC_MSY, ) #required sections, no optional sections but FRAGMENTS that are checked elsewhere..
CFG_SEC_LST_ALW =  (CFG_SEC_MET, CFG_SEC_MSY, CFG_SEC_QDY, )  #allowed sections

#Keys in 'Methods' section
CFG_MET_VERBOSITY = r"VERBOSITY" #Controls amount of information printed in the output (may be critical for long MD trajectories)
CFG_MET_APX = r"JOB" #Calculation type
CFG_MET_RAT = r"RATES" #Calculation of exciton transfer rates
CFG_MET_QD = r"QD" #Quantum dynamics
CFG_MET_KEY_LST_REQ  = (CFG_MET_APX,) #required keys
CFG_MET_KEY_LST_ALW  = (CFG_MET_APX, CFG_MET_RAT, CFG_MET_QD, CFG_MET_VERBOSITY, ) #allowed keys

#Values of CFG_MET_VERBOSITY key
CFG_MET_VERB_MAX = 99 #Maximal output (default)
CFG_MET_VERB_TRJ = 0 #Succinct output for long MD trajectories
CFG_MET_VERB_VAL_LST_ALW = (CFG_MET_VERB_MAX, CFG_MET_VERB_TRJ, )

#Values of CFG_MET_APX key
CFG_MET_ARX_FRG  = r"MAKE_FRAGMENTS" #Initial fragmentation of the molecular system in the PDB file
CFG_MET_ARX_SUR  = r"SURVEY" #Forster coupling of ALL excited states for all fragments
CFG_MET_ARX_FRS  = r"FORSTER" #Forster coupling of SELECTED excited states (only one state per fragment is allowed) followed by variational calculation
CFG_MET_ARX_LFT  = r"LIFETIME" #Calculate lifetime(s) of the donor(s) only

CFG_MET_ARX_VAL_LST_ALW = (CFG_MET_ARX_FRG, CFG_MET_ARX_FRS, CFG_MET_ARX_SUR, CFG_MET_ARX_LFT, )

#Values of CFG_MET_RAT key
CFG_MET_RAT_NONE = r"NONE" #Do not calculate rates (default) CHECK IF THE VALUE ALLOWED?
CFG_MET_RAT_OVL  = r"OVERLAP" #Calculate rates based on spectral overlap CHECK IF VALUE ALLOWED?
CFG_MET_RAT_VAL_LST_ALW = (CFG_MET_RAT_NONE, CFG_MET_RAT_OVL, )

#Values of QD key
CFG_MET_QD_NONE = r"NONE"  #Do not run quantum dynamics (default)
CFG_MET_QD_MARKOVIAN = r"MARKOVIAN" #Markovian dynamics
CFG_MET_QD_NONMARKOVIAN = r"NONMARKOVIAN" #Non-Markovian dynamics (based on the trajectory)
CFG_MET_QD_VAL_LST_ALW = (CFG_MET_QD_NONE, CFG_MET_QD_MARKOVIAN, CFG_MET_QD_NONMARKOVIAN, ) #Allowed keywords 

#Keys in 'Quantum Dynamics' section
CFG_QDY_INT = r"INTEGRATION" #Integration parameters
CFG_QDY_LABELS = r"LABELS" #Labels for pyplot
CFG_QDY_INI_POPS = r"INI_POPS" ##Initial populations for each fragment 
CFG_QDY_KEY_LST_REQ  = () #required keys
CFG_QDY_KEY_LST_ALW  = (CFG_QDY_INT, CFG_QDY_LABELS, CFG_QDY_INI_POPS, ) #allowed keys

#Keys in 'MolecularSystem' section
CFG_MSY_GEO = r"GEOMETRY_FILE" #PDB file containing coordinates of the entire molecular system
CFG_MSY_TRJ = r"TRAJECTORY_FILE" #MD trajectory file (GROMACS .GRO formt) 

CFG_MSY_RES = r"RESONANCE_CM1" # resonance condition: two states are considered to be in resonace
				# if their excitation wavenumbers differ less or equal to the given value (cm-1)
				# this keywork is actually used only for SURVEY jobs only
CFG_MSY_RTR = r"RESONANCE_THR" # resonance condition: threshold value for overlaping Gaussian spectral bands

CFG_MSY_EL_SCR_MODEL = r"EL_SCR_MODEL" #Electrostatic screening model (
				# 0 = no screening (default)
				# 1 = uniform screening
				# 2 = exponential screening model 
CFG_MSY_EL_SCR_FACTOR = r"EL_SCR_FACTOR" #Uniform electrostatic screening factor (if EL_SCR_MODEL = 1, default EL_SCR_FACTOR = 0.8)
				#Exponential screening function (if EL_SCR_MODEL = 2) 
				#(defaults: pre-exponential  factor A=2.68, attenuation factor beta=0.27A-1, asymptotic value s0=0.54)
CFG_MSY_EL_SCR_EXP_FUNC = "EL_SCR_EXP_FUNC"

#Boltzmann factors
CFG_MSY_BLZ_COMP = r"BLZ_COMP" #Account for Boltzmann factors in calculation of rate
CFG_MSY_BLZ_TEMP = r"BLZ_TEMP" #Temperature for calculation of Boltzmann factors 

#Calculation of spectral overlaps lower and upper limits of integration in cm-1
CFG_MSY_OVL_LOW_LIM = r"OVRLP_INT_LOWER_LIM"
CFG_MSY_OVL_UPP_LIM = r"OVRLP_INT_UPPER_LIM"

#Calculation of spectral overlaps lower and upper limits of integration in nm for spectral overla in wavelengths domain
CFG_MSY_OVL_LOW_LIM_NM = r"OVRLP_INT_LOWER_LIM_NM"
CFG_MSY_OVL_UPP_LIM_NM = r"OVRLP_INT_UPPER_LIM_NM"

#Refraction index
CFG_MSY_RFX            = r"N"

#Orientation factor kappa squared
CFG_MSY_KAPPA_SQ = r"KAPPA_SQ"

#Strickler-Berg calculation of lifetime of donors
CFG_MSY_SB = r"SB"

CFG_MSY_KEY_LST_ALW = (CFG_MSY_GEO, CFG_MSY_TRJ, CFG_MSY_RES, CFG_MSY_RTR, CFG_MSY_EL_SCR_MODEL, CFG_MSY_EL_SCR_FACTOR, CFG_MSY_EL_SCR_EXP_FUNC, CFG_MSY_OVL_LOW_LIM, CFG_MSY_OVL_UPP_LIM, CFG_MSY_BLZ_COMP, CFG_MSY_BLZ_TEMP, CFG_MSY_OVL_LOW_LIM_NM, CFG_MSY_OVL_UPP_LIM_NM, CFG_MSY_RFX, CFG_MSY_SB, CFG_MSY_KAPPA_SQ, ) #Allowed keys
CFG_MSY_KEY_LST_REQ = ( ) #Required keys

#Keys in 'Fragment' sections (e.g. Fragment1)
CFG_FRG_NAM = r"NAME" #Fragment name (PDB residue name)
CFG_FRG_ID  = r"ID" #Fragment ID (PDB residue ID number)
CFG_FRG_CHAIN_ID  = r"CHAIN_ID" #Fragment ID (PDB chain ID letter code)
CFG_FRG_ALTLOC = r"ALTLOC" #Alternate location indicator. (the same atom may occupy serveral alternate locations see PDB file secs.)
CFG_FRG_AL  = r"ATOMLIST" #List of atoms that belong to the fragment
CFG_FRG_ATL = r"ATOMTRANS" #List of atoms that will be used for geometry matching
CFG_FRG_GEF = r"EXSTATE_FILE" #Geometry, Excitation energies, and transition dipole moments file
CFG_FRG_EST = r"EXSTATE" #Index of the excited state used for coupling calculations
CFG_FRG_KEY_LST_ALW = (CFG_FRG_NAM, CFG_FRG_ID, CFG_FRG_AL, CFG_FRG_ATL, CFG_FRG_GEF, CFG_FRG_EST, CFG_FRG_ALTLOC, CFG_FRG_CHAIN_ID, ) #Allowed keys

#Required keys depending on jobtype (CFG_MET_ARX_...)
CFG_FRG_KEY_LST_REQ_FRG = (CFG_FRG_NAM, CFG_FRG_ID, ) #Required keys for job type: make_fragments
CFG_FRG_KEY_LST_REQ_SUR = (CFG_FRG_NAM, CFG_FRG_ID, CFG_FRG_GEF, ) #Required keys for job type: survey
CFG_FRG_KEY_LST_REQ_FRS = (CFG_FRG_NAM, CFG_FRG_ID, CFG_FRG_GEF, CFG_FRG_EST, ) #Required keys for job type: Forster
CFG_FRG_KEY_LST_REQ_LFT = (CFG_FRG_NAM, CFG_FRG_ID, CFG_FRG_GEF, ) #Required keys for job type: lifetime calculation

#Default values for fragments
CFG_FRG_CHAIN_ID_DAFAULT = r""
CFG_FRG_ALTLOC_DAFAULT   = r""

#Sections in EXSTATE_FILE
#Dollar signs (S) are escaped for regular expression parsing
#see ExStatesReader() class for details
CFG_EXS_SEC_GEO = r"\$GEOMETRY"
CFG_EXS_SEC_ETD = r"\$EXCITED_STATES"
CFG_EXS_SEC_CNT = r"\$CENTER"
CFG_EXS_SEC_QDV = r"\$QDVIB"
CFG_EXS_SEC_ABS = r"\$ABS_SPEC"
CFG_EXS_SEC_EMS = r"\$EMS_SPEC" 
CFG_EXS_SEC_END = r"\$END"

#Empty PDB atom IDs and PDB atom and residue names for atoms read from  EXSTATE_FILE
ATOM_PDB_ID_NONE = -1
ATOM_PDB_NAME_NONE = ""
ATOM_PDB_RES_ID_NONE = -1
ATOM_PDB_RES_NAME_NONE = ""
ATOM_PDB_ALT_LOC_NONE = ""
ATOM_PDB_CHAIN_ID_NONE = ""

#Comment line symbol in input files (with exception of ccnfiguration .ini file)
COMM_SYMBOL = r"#"

#Separator in lists
LIST_SEPARATOR = r","

#Internal keys for each fragment (cannot be specified in the user's input):

CFG_FRG_NATOMS = r"NATOMS" #Number of atoms in a fragment (from PDB file)
CFG_FRG_ATOMS = r"ATOMS" #List of atoms in a fragment (from PDB file)

#Infromation obtained from quantum calulations and read from CFG_FRG_GEF = "EXSTATE_FILE"  file
CFG_EXS_NATOMS = r"EXS_NATOMS" #Number of atoms in a fragment (from CFG_FRG_GEF = "EXSTATE_FILE"  file)
CFG_EXS_ATOMS = r"EXS_ATOMS" #List of atoms in a fragment (from CFG_FRG_GEF = "EXSTATE_FILE"  file)

CFG_EXS_NEXSTATE = r"EXS_NEXSTATES" #Number of excited states in a fragment
CFG_EXS_EXSTATE = r"EXS_EXSTATES" #List of excited states in a fragment
CFG_EXS_CENTER = r"EXS_CENTER" #Center of the fragment

CFG_EXS_QDVID = r"EXS_QDVID" #Parameters of quatum dynamics for each fragment
CFG_EXS_NQDVID = r"EXS_NQDVID" #Number of vibrational mode per electronic excited state

CFG_EXS_EXSTATE_ID = r"EXS_EXSTATE_ID" #Target excited state (only one state for Forster and many states for survey)

CFG_TRNS = r"TRNS" #Tranformation (scaling, rotation, translation) parameters that match fragment coordinate system to the PDB file corrdinate system,
CFG_TRNS_EXS_ATOMS = r"TRNS_EXS_ATOMS" #List of atoms in a fragment (from CFG_FRG_GEF = "EXSTATE_FILE"  file) in TRANSFORMED coordinates (fragment matched to the PDB coordinate system)
CFG_TRNS_EXS_EXSTATE = r"TRNS_EXS_EXSTATES" #List of excited states in a fragment in TRANSFORMED coordinates (fragment matched to the PDB coordinate system)
CFG_TRNS_ORIG = r"TRNS_ORIG" #Coordinates of the fragment origin in TRANSFORMED coordinates (matched to PDB file coordinate system)
CFG_TRNS_CENTER = r"TRNS_CENTER" #Coordinates of the transformed center of the framnt used for coupling calculations

RESONANCE_CM1_DEFAULT = 1000 #Default value (cm-1) for the resonance condition
RESONANCE_THR_DEFAULT = 1.0  #Default value for the threshold condition (zero meams resonance condition is not used)

#Electrostatic screening model (no screening)
EL_SCR_MODEL_NONE = 0
EL_SCR_MODEL_UNIFORM = 1
EL_SCR_MODEL_EXPONENTIAL = 2

#Strickler-Berg calculation of lifetime of donors
SB_NONE = 0  #Do not run (Default)
SB_INTEGRATE = 1 #Evaluate integrals of the donor fluorescnce spectrum
SB_MAX_ONLY  = 2  #Estimate integrals based on its maximum of the donor fluorescnce spectrum 

#Default verbosity level
CFG_MET_VERB_DEFAULT = CFG_MET_VERB_MAX

#Default geometry file
CFG_MSY_GEO_NONE = ""

#Default trajectory file
CFG_MSY_TRJ_DAFAULT = ""

CFG_MSY_EL_SCR_MODEL_DEFAULT = EL_SCR_MODEL_NONE
#Electrostatic screening uniform factor:
EL_SCR_FACTOR_DEFAULT = 0.8
#Electrostatic screening expoential function (A,B,S):
#Pre-exponential factor A=2.68, attenuation factor beta=0.27A-1, asymptotic value S=0.54)
EL_SCR_EXP_FUNC_A_DEFAULT = 2.68
#Attenuation factor beta, A-1
EL_SCR_EXP_FUNC_B_DEFAULT = 0.27
#Asymptotic value, S
EL_SCR_EXP_FUNC_S_DEFAULT = 0.54

#Default lower and upper limits of numerical integration of
#spectral overlap integrals in cm-1
OVL_LOW_LIM_DEFAULT = 1.0e-6
OVL_UPP_LIM_DEFAULT = 1.0e5

#Default lower and upper limits of numerical integration of
#spectral overlap integrals in nm
OVL_LOW_LIM_NM_DEFAULT = 200.0
OVL_UPP_LIM_NM_DEFAULT = 2000.0

#number of points in the generated tabulated spectrum 
SPEC_TAB_POINTS = 200

#Default refraction index
CFG_MSY_RFX_DEFAULT = 1.0

#Default Strickler-Berg calculation
CFG_MSY_SB_DEFAULT = SB_NONE

#By default the kappa orinetation factor is not defined (value -1.0)
CFG_MSY_KAPPA_SQ_DEFAULT = -1.0

CFG_MET_QD_NONE = r"NONE" #No quantum dynamics by default

#Default values for quantum_dynamics section
QDY_INT_DEFAULT = r"0.0, 5.0, 1000" #Integration initial time 0.0, final time 5.0 ps, 1000 steps
QDY_LABELS_DEFAULT = r"" #No labels for pyplot by default


#Default calculation of Boltzmann factors
BLZ_COMP_DEFAULT = 0 #Do not calculate Boltzmann factors in calculation of rates
#Compute Boltzmann factors 
BLZ_COMP_EMS_INI = 1 #Equliburum populations among emission levels (useful for calculation of initial rates only)
BLZ_COMP_ABS_INI = 2 #Equliburum populations among absorption levels (useful for calculation of initial rates only)
BLZ_COMP_EMS_EQ  = 3 #Equliburum populations are _MAINTAINED_ among emission levels (useful for kinetic equaitons)
BLZ_COMP_ABS_EQ  = 4 #Equliburum populations are _MAINTAINED_ among absorption levels (useful for kinetic equaitons)
#
BLZ_TEMP_DEFAULT = 300.0 #Temperature, K for calculation of Boltzmann factors 

CFG_PRINT_TDM_VIS = False #Internal key to print coordinates for visualization of transition dipole moments
TDM_SCALE_VIS = 1.0 #Scaling factor for transition dipole moments used for visualization


#Visualization for quantum dynamics (charts) see qd.py for details imports matplotlib.pyplot 
#INT_CONFIG_WINVIS = False
INT_CONFIG_WINVIS = True

#-------------------------------------------------------------------------------
class Configuration(MyFileService):
	"""Class containing  configuration information read from ini file"""

	def __init__(self, IniFile=""):
		"""Setup empty parameters to be filled with info read from the file
		   Setup default values
		"""
		self.ps = MyReportService()
		#Configuration information:
		self.Methods = AutoArr() #Computational methods/approximations used
		self.QuantumDynamics  = AutoArr() #Quantum dynamics parameters
		self.Fragments = AutoArr() #Fragments
		self.MolSys = AutoArr() #Molecular system (e.g. dimer) data
		self.Trj = AutoArr() #MD trajectory
		
		#Set default verbosity (maximal)
		self.Methods[CFG_MET_VERBOSITY] = CFG_MET_VERB_DEFAULT
		
		#Do not calculate exciton transfer rates by default
		self.Methods[CFG_MET_RAT] = CFG_MET_RAT_NONE
		
		#Do not run quantum dynamics by default
		self.Methods[CFG_MET_QD] = CFG_MET_QD_NONE
		
		#Empty input geometry file by default 
		self.MolSys[CFG_MSY_GEO] = CFG_MSY_GEO_NONE

		#Default value for trajectory
		self.MolSys[CFG_MSY_TRJ] = CFG_MSY_TRJ_DAFAULT

		#Default value for the resonance condition
		self.MolSys[CFG_MSY_RES] = RESONANCE_CM1_DEFAULT
		
		#Default value for the resonance condition using spectral overlap (not used by default)
		self.MolSys[CFG_MSY_RTR] = RESONANCE_THR_DEFAULT
		
		self.MolSys[CFG_MSY_BLZ_TEMP] = BLZ_TEMP_DEFAULT

		#Default elctrostatic screening model 
		self.MolSys[CFG_MSY_EL_SCR_MODEL] = CFG_MSY_EL_SCR_MODEL_DEFAULT
		#Electrostatic screening uniform factor (not for default model though):
		self.MolSys[CFG_MSY_EL_SCR_FACTOR] = EL_SCR_FACTOR_DEFAULT
		#Exponential screening function (not for default model though):
		self.MolSys[CFG_MSY_EL_SCR_EXP_FUNC] = str(EL_SCR_EXP_FUNC_A_DEFAULT)+LIST_SEPARATOR+str(EL_SCR_EXP_FUNC_B_DEFAULT)+LIST_SEPARATOR+str(EL_SCR_EXP_FUNC_S_DEFAULT)

		#Default lower and upper limits of numerical integration of spectral overlap integrals in cm-1
		self.MolSys[CFG_MSY_OVL_LOW_LIM] = OVL_LOW_LIM_DEFAULT
		self.MolSys[CFG_MSY_OVL_UPP_LIM] = OVL_UPP_LIM_DEFAULT

		#Default lower and upper limits of numerical integration of spectral overlap integrals in nm
		self.MolSys[CFG_MSY_OVL_LOW_LIM_NM] = OVL_LOW_LIM_NM_DEFAULT
		self.MolSys[CFG_MSY_OVL_UPP_LIM_NM] = OVL_UPP_LIM_NM_DEFAULT
		
		#Default refraction index
		self.MolSys[CFG_MSY_RFX] = CFG_MSY_RFX_DEFAULT
		
		#Default Strickler-Berg calculation of lifetime of donors
		self.MolSys[CFG_MSY_SB] = CFG_MSY_SB_DEFAULT
		
		#Default Strickler-Berg calculation of lifetime of donors
		self.MolSys[CFG_MSY_SB] = CFG_MSY_SB_DEFAULT

		#The orientation factor is not defined by default
		self.MolSys[CFG_MSY_KAPPA_SQ] = CFG_MSY_KAPPA_SQ_DEFAULT
		
		#Default computation of Boltzmann factors
		self.QuantumDynamics[CFG_QDY_INT] = QDY_INT_DEFAULT
		self.QuantumDynamics[CFG_QDY_LABELS] = QDY_LABELS_DEFAULT
		
		#Configuration file name is passed as a command line argument
		self.ConfigFileName = self.ParseCommandLine(IniFile)
		#Read Configuration from the config file
		self.ReadConfig(self.ConfigFileName)
		
		self.MolSys[CFG_MSY_BLZ_COMP] = BLZ_COMP_DEFAULT
		
		return

	def ParseCommandLine(self,IniFile):
		"""Parse command line to get a name of the configuration file"""
		
		#Test mode - do not use the command line
		if ( len(IniFile.strip()) != 0 ):
			self.ConfigFileName = IniFile.strip()
			return self.ConfigFileName
		
		parser = ArgumentParser(description=r"")
		parser.add_argument(r"-f", dest=r"filename", required=True, help=r"provide a configuration (.ini) file", type=str )
		self.ConfigFileName = parser.parse_args().filename
		return self.ConfigFileName

	def IsInList(self, Option, List):
		"""Case insensitive check if given option (keyword or value) is in the list, no error messages"""
		if not (Option.upper() in (listed_option.upper() for listed_option in List)):
			return False
		return True

	def InList(self, Option, List, Comment = ''):
		"""Case insensitive check if given option (keyword or value) is in the listm with error message"""

		if not self.IsInList(Option, List):
			print(r"Error: item '"+Option+r"' is not in: ", end=' ')
			StrList = ', '.join([str(x) for x in List])
			print(StrList, end=' ')
			print(r" ",Comment)
			exit(-1)
		return Option

	def ListInList(self, List1, List2, Comment = ''):
		"""Each element from List1 is found in List2. Case insensitive, ignores order of items in both lists"""
		for Option in List1:
			self.InList(Option, List2, Comment)
		return

	def VerifyConfiguration(self, Configp):
		"""Verify validity of the user specified configuration in .ini file."""

		#Verify all requred sections but fragments "e.g. FRAGMENT1, FRAGMENT2, ... "
		#Each element from CFG_SEC_LST_REQ is found in list of sections
		self.ListInList(CFG_SEC_LST_REQ, Configp.sections(), "sections")

		#Verify required keys (CFG_MET_APX ) in METHODS section
		#Each element from CFG_MET_KEY_LST_REQ is found in METHODS section
		self.ListInList(CFG_MET_KEY_LST_REQ, Configp.options(CFG_SEC_MET), r"in section: "+CFG_SEC_MET)

		#Each element from METHODS section is found CFG_MET_KEY_LST_ALW
		self.ListInList(Configp.options(CFG_SEC_MET), CFG_MET_KEY_LST_ALW, r"in section: "+CFG_SEC_MET)

		#Check allowed jobtype CFG_MET_APX values
		JobType = self.Methods[CFG_MET_APX]
		self.InList(JobType, CFG_MET_ARX_VAL_LST_ALW, r"job types")
		
		#Check allowed rate calculation CFG_MET_RAT values
		RateType = self.Methods[CFG_MET_RAT]
		self.InList(RateType, CFG_MET_RAT_VAL_LST_ALW, r"rate calculation types")
		
		#Check allowed quantum dynamics calculation CFG_MET_QD values
		QuantumDynamicsType = self.Methods[CFG_MET_QD]
		self.InList(QuantumDynamicsType, CFG_MET_QD_VAL_LST_ALW, r"quantum dynamics types")

		if (self.Methods[CFG_MET_QD] != CFG_MET_QD_NONE ):
			#Check required keys in QUANTUM_DYNAMICS section
			self.ListInList(CFG_QDY_KEY_LST_REQ, Configp.options(CFG_SEC_QDY), r"in section: "+CFG_SEC_QDY)
			#Check allowed keys in QUANTUM_DYNAMICS section
			self.ListInList(Configp.options(CFG_SEC_QDY), CFG_QDY_KEY_LST_ALW, r"in section: "+CFG_SEC_QDY)

		#Check required keys in MOLECULAR_SYSTEM section
		self.ListInList(CFG_MSY_KEY_LST_REQ, Configp.options(CFG_SEC_MSY), r"in section: "+CFG_SEC_MET)
		#Check allowed keys in MOLECULAR_SYSTEM section
		self.ListInList(Configp.options(CFG_SEC_MSY), CFG_MSY_KEY_LST_ALW, r"in section: "+CFG_SEC_MET)
		
		#Check required keys for each jobtype
		ReqFrgKeys = None #list of required keywords depending on the job type
		if  ( JobType == CFG_MET_ARX_FRG ):
			ReqFrgKeys = CFG_FRG_KEY_LST_REQ_FRG
		elif( JobType == CFG_MET_ARX_SUR ):
			ReqFrgKeys = CFG_FRG_KEY_LST_REQ_SUR
		elif( JobType == CFG_MET_ARX_FRS ):
			ReqFrgKeys = CFG_FRG_KEY_LST_REQ_FRS
		elif( JobType == CFG_MET_ARX_LFT ):
			ReqFrgKeys = CFG_FRG_KEY_LST_REQ_LFT
		else:
			print(r"Error: unknown job type: "+JobType+r" in section: "+CFG_SEC_MET+r" (cannot verify fragments)")
			exit(-1)

		#In each fragment
		for FragID, Fragment in self.Fragments.items():
			#Check if ALL required keywords are present
			self.ListInList(ReqFrgKeys, Fragment, r"in "+CFG_SEC_FRG+str(FragID+1))

			if ( JobType != CFG_MET_ARX_FRG ): #Make fragments job does not require excited states
				#Verify excited states
				nFrgIDList = len(ParseList(Fragment[CFG_FRG_EST]))
				if ( JobType == CFG_MET_ARX_FRS ): #Forster calculation
					
					if ( nFrgIDList != 1 ):
						print(r"Error: exactly one excited state for each fragment is required for job type: "+JobType, end=' ')
						print(r"given "+nFrgIDList+r"in "+CFG_SEC_FRG+str(FragID+1))
						exit(-1)
				else: #Survey calculation
					if ( nFrgIDList < 1 ):
						print(r"Error: at least one excited state for each fragment is required for job type: "+JobType, end=' ')
						print(r"given "+nFrgIDList+"in "+CFG_SEC_FRG+str(FragID+1))
						exit(-1)
		return

	def ReadConfig(self, FileName):
		"""Read configuration from the config file."""
		if  not self.IsValidFile(FileName,r"Cannot find a configuration file"):
			exit(-1)

		self.ps.DateTimeStamp()
		self.ps.SysVersion()
		self.ps.PrintSec(r"CONFIGURATION", DivSym = "=")

		FILE = open(FileName,'r')
		self.ps.PrintKV(r"Configuration File",FileName)
		Configp = configparser.ConfigParser()
		#Configp.readfp(FILE)
		Configp.read_file(FILE)
		#Check sections
		SecFrgRE = r'^\s*'+CFG_SEC_FRG.strip().upper()+r'\d+\s*$'
		for Section in Configp.sections():
			if not re.match(SecFrgRE, Section.strip().upper()):
				self.InList(Section, CFG_SEC_LST_ALW, r'or '+CFG_SEC_FRG+r"1, "+CFG_SEC_FRG+r"2, etc.")
		#Read 'Methods' section
		self.ps.PrintSec(CFG_SEC_MET)
		for Option in Configp.options(CFG_SEC_MET):
			Option = Option.upper()
			self.Methods[Option] = Configp.get(CFG_SEC_MET, Option).strip()
			self.ps.PrintKV(Option,self.Methods[Option])

		if (self.Methods[CFG_MET_QD] != CFG_MET_QD_NONE ):
			#Read 'Quantum_Dynamics' section
			self.ps.PrintSec(CFG_SEC_QDY)
			for Option in Configp.options(CFG_SEC_QDY):
				Option = Option.upper()
				self.QuantumDynamics[Option] = Configp.get(CFG_SEC_QDY, Option).strip()
				self.ps.PrintKV(Option,self.QuantumDynamics[Option])

		#Read 'MolecularSystem' section
		self.ps.PrintSec(CFG_SEC_MSY)
		for Option in Configp.options(CFG_SEC_MSY):
			Option = Option.upper()
			#Do not change case (lower->upper) of values of the keywords
			#File names are case sensitive!
			self.MolSys[Option] = Configp.get(CFG_SEC_MSY, Option).strip()
			self.ps.PrintKV(Option,self.MolSys[Option])

		#Read and check 'FRAMGMENT' sections
		NFragments = 0
		for Section in Configp.sections():
			#Parse fragment sections
			if re.match(SecFrgRE, Section.strip().upper()):
				self.ps.PrintSec(Section)
				self.Fragments[NFragments] = AutoArr()
				for Option in Configp.options(Section):
					Option = Option.upper()
					#Do not change case (lower->upper) of values of the keywords
					self.Fragments[NFragments][Option] = Configp.get(Section, Option).strip()
					#Specify the number of atoms in the fragment
					self.ps.PrintKV(Option,self.Fragments[NFragments][Option])
					SpecifiedKeys = list(self.Fragments[NFragments].keys())
					#Set default options if not sepcified by user for fragments
					if ( not self.IsInList(CFG_FRG_AL, SpecifiedKeys) ):
						self.Fragments[NFragments][CFG_FRG_AL] = r"0"

					if ( not self.IsInList(CFG_FRG_ATL, SpecifiedKeys) ):
						self.Fragments[NFragments][CFG_FRG_ATL] = r"0"

					JobType = self.Methods[CFG_MET_APX]
					if ( (not self.IsInList(CFG_FRG_ATL, SpecifiedKeys)) and (JobType == CFG_MET_ARX_SUR) ):
						self.Fragments[NFragments][CFG_FRG_EST] = r"0"

					if ( not self.IsInList(CFG_FRG_ALTLOC, SpecifiedKeys) ):
						self.Fragments[NFragments][CFG_FRG_ALTLOC] = CFG_FRG_ALTLOC_DAFAULT
					
					if ( not self.IsInList(CFG_FRG_CHAIN_ID, SpecifiedKeys) ):
						self.Fragments[NFragments][CFG_FRG_CHAIN_ID] = CFG_FRG_CHAIN_ID_DAFAULT
						
				NFragments += 1

		if ( NFragments < 1):
			print(r"Error: no fragments were specified in the configuration file: "+FileName)
			exit(-1)

		self.ps.PrintDiv()
		self.ps.PrintKV(r"Number of Fragments",NFragments)
		self.ps.PrintDiv(DivSym = r"=")

		#Verify the configuration parameters
		self.VerifyConfiguration(Configp)
		FILE.close()

		#Prepare lists of atoms
		for FragID, Fragment in self.Fragments.items():
			Fragment[CFG_FRG_NATOMS] = 0 #The nuber of atoms in the fragment (will be updated duing actual PDB reading)
			Fragment[CFG_FRG_ATOMS]  = AutoArr() #Array of atoms read from the PDB file (will be updated during actual PDB file reading)
			Fragment[CFG_EXS_NATOMS] = 0 #Number of atoms in a fragment (from CFG_FRG_GEF = "EXSTATE_FILE"  file)
			Fragment[CFG_EXS_ATOMS] = AutoArr() #List of atoms in a fragment (from CFG_FRG_GEF = "EXSTATE_FILE"  file)
			Fragment[CFG_EXS_NEXSTATE] = 0 #Number of excited states in a fragment
			Fragment[CFG_EXS_EXSTATE] = AutoArr() #List of excited states in a fragment

		#a = self.Fragments[0]['NAMEz']
		#print a
		#pp = pprint.PrettyPrinter(indent=4)
		#pp.pprint(a)
		#print
		return
#-------------------------------------------------------------------------------