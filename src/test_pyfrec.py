from pytest import *
from pyfrec import *
import re
from couplib.constants import RE_FLOAT, RE_DOUBLE
from math import *
#from couplib.constants import RE_FLOAT

TEST_ROOR_DIR = r"./tests/"
TEST_PYFREC_LIFETIMES = TEST_ROOR_DIR+r"fluorescein.ini"
TEST_PYFREC_COUPLINGS = TEST_ROOR_DIR+r"bodipy_c_ac.ini"

FLOAT_MATCH_TOL = 1E-3 #Tolerance in reading float numbers  (exponents are accounted for separately)
#-------------------------------------------------------------------------------
def ParseCouplingOutput(out):
	"""Parse PyFREC output for accuracy: coupling calculation"""

	ReadingPDBStartedOkay = False
	ExcitedStatesPresentOkay = False
	ElectronicCouplingsOkay = False
	CalculationFinishedOkay = False

	for line in out.splitlines():
		#Reading PDB file...
		matched = re.match(r'^\s*Reading\s+PDB\s+file\.\.\.\s*$', line)
		if (matched):
			ReadingPDBStartedOkay = True
		
		#Last Atom of the 1st fragment
		#C20      -1.336       0.124        0.009
		#Last Atom of the 2nd fragment
		#C29      -9.676       -0.248       -0.303    
		matched = re.match(r'^\s*([\w\d]+)\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s*$', line)
		if (matched and ReadingPDBStartedOkay):
			FRG1_AtomName = "C20"
			FRG1_x = -1.336
			FRG1_y = 0.124
			FRG1_z = 0.009

			FRG2_AtomName = "C29"
			FRG2_x = -9.676
			FRG2_y = -0.248
			FRG2_z = -0.303    
			
			ErrorMSG = r"Coordiantes of the atom: "+FRG1_AtomName+" of the first fragment is not found: "
			if (str(matched.group(1)).strip() == FRG1_AtomName):
				if ( abs(float(matched.group(2)) - FRG1_x)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", x-coordiante problem"
				if ( abs(float(matched.group(3)) - FRG1_y)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", y-coordiante problem"
				if ( abs(float(matched.group(4)) - FRG1_z)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", z-coordiante problem"

			ErrorMSG = r"Coordiantes of the atom: "+FRG2_AtomName+" of the second fragment is not found: "
			if (str(matched.group(1)).strip() == FRG2_AtomName):
				if ( abs(float(matched.group(2)) - FRG2_x)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", x-coordiante problem"
				if ( abs(float(matched.group(3)) - FRG2_y)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", y-coordiante problem"
				if ( abs(float(matched.group(4)) - FRG2_z)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", z-coordiante problem"

		matched = re.match(r'^\s*Excited\s+states\s+read\s+from\s+the\s+input\s+file:\s*$', line)
		if (matched):
			ExcitedStatesPresentOkay = True

		#FRAGMENT2 NAME = TRZ ID = 2
		#matched = re.match(r'^\s*FRAGMENT(\d+)\s+NAME\s*=\s*(\w+)\s*ID\s*=\s*(\d+)\s*$', line)
		#if (matched):

		#1     0.0          (-3.3813      -0.0886      0.0094      ) 3.3825       8.5974       546.6
		#1     562.3        (0.002        0.0          0.3492      ) 0.3492       0.8876       0.0    
		matched = re.match(r'^\s*(\d+)\s+('+RE_FLOAT+r')\s+\(('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+\)\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s*$', line)
		if (ExcitedStatesPresentOkay and matched):
			ExcN1  = 1
			AbsNM1 = 0.0
			TDM_x1 = -3.3813
			TDM_y1 = -0.0886
			TDM_z1 = 0.0094
			TDMae1 = 3.3825
			TDMD1  = 8.5974
			EmsNM1 = 546.6
			ErrorMSG = " problem reading excited states first fragment, 1st excited state TMD"
			if ( (int(matched.group(1)) == ExcN1) and (abs(float(matched.group(2)) - AbsNM1)  < FLOAT_MATCH_TOL) ):
				if ( abs(float(matched.group(3)) - TDM_x1)*10 > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", x-componet"
				if ( abs(float(matched.group(4)) - TDM_y1)*10 > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", y-componet"
				if ( abs(float(matched.group(5)) - TDM_z1)*10 > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", z-componet"
				if ( abs(float(matched.group(6)) - TDMae1)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", magnitude a.u."
				if ( abs(float(matched.group(7)) - TDMD1)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", magnitude in Debye"
				if ( abs(float(matched.group(8)) - EmsNM1)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", emission maximum nm"

			ExcN2  = 1
			AbsNM2 = 562.3
			TDM_x2 = 0.002
			TDM_y2 = 0.0
			TDM_z2 = 0.3492
			TDMae2 = 0.3492
			TDMD2  = 0.8876 
			EmsNM2 = 0.0
			ErrorMSG = "  problem reading excited states second fragment, 1st excited state TMD"
			if ( (int(matched.group(1)) == ExcN2) and (abs(float(matched.group(2)) - AbsNM2)  < FLOAT_MATCH_TOL) ):
				if ( abs(float(matched.group(3)) - TDM_x2)*10 > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", x-componet of TDM"
				if ( abs(float(matched.group(4)) - TDM_y2)*10 > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", y-componet of TDM"
				if ( abs(float(matched.group(5)) - TDM_z2)*10 > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", z-componet of TDM"
				if ( abs(float(matched.group(6)) - TDMae2)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", magnitude in a.u."
				if ( abs(float(matched.group(7)) - TDMD2)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", magnitude in Debye"
				if ( abs(float(matched.group(8)) - EmsNM2)  > FLOAT_MATCH_TOL):
					assert False, ErrorMSG+r", emission maximum nm"

		#Electronic couplings list:
		matched = re.match(r'^\s*Electronic couplings list:\s*$', line)
		if (matched):
			ElectronicCouplingsOkay = True
		
		#1     BDP1            2     TRZ2            1            1    17783.468        9.308        1.181     -0.05109       -2.555     -0.06035       -2.434       -1.347
		matched = re.match(r'^\s*(\d+)\s+([\w\d]+)\s+(\d+)\s+([\w\d]+)\s+(\d+)\s+(\d+)\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s*$', line)
		if (ElectronicCouplingsOkay and matched):
			nFrg1 = 1
			Frg1Tag = "BDP1"
			nFrg1 = 2
			Frg1Tag = "TRZ2"
			nExc1 = 1
			nExc2 = 1
			dE = 17783.468
			R = 9.308
			mumu = 1.181
			Ori = -0.05109
			OriPrc = -2.555
			K = -0.06035
			E = -2.434
			EScr = -1.347
			
		
		#Done.
		matched = re.match(r'^\s*Done\.\s*$', line)
		if (matched):
			CalculationFinishedOkay = True
	
	assert ReadingPDBStartedOkay,r"PDB file reading problem..."
	assert ExcitedStatesPresentOkay,r"Excited states reading problem..."
	assert ElectronicCouplingsOkay,r"Electronic coupling calculation problem..."
	assert CalculationFinishedOkay,r"Calculation has NOT finished..."
	
#-------------------------------------------------------------------------------
def ParseLifetimeOutput(out):
	"""Parse PyFREC output for accuracy: fluorerscence lifetime calculation"""
	
	ConfigOkay = False
	AbsorptionTabSpecOkay = False
	FluorescenceTabSpecOkay = False
	QunatumYieldOKay = False
	StricklerBergOkay = False
	CalculationFinishedOkay = False
	
	nFrags = 1 #Number of fragments in the sample input
	FrgName = "FLR1" 

	for line in out.splitlines():
		#Number of Fragments     : 1
		matched = re.match(r'^\s*Number\s+of\s+Fragments\s+:\s+(\d+)\s*$', line)
		if (matched):
			ConfigOkay = True #Indicates the presence of the configuration information in the output
			nFragsOK = False
			if (int(matched.group(1).strip()) == nFrags):
				nFragsOK = True
			assert nFragsOK,r"Problem with reading the number of fragments, found:"+matched.group(1).strip()

		#Fragment: FLR1 tabulated ABSORPTION spectrum is available. Maximum: 19990.0  cm-1 ( 500.3 nm ) Eps= 92300.0 M-1 cm-1; Integration limits: 400 - 600 nm
		matched = re.match(r'^\s*Fragment:\s+(\w+)\s+tabulated\s+ABSORPTION\s+spectrum\s+is\s+available.\s+Maximum:\s+('+RE_FLOAT+r')\s+cm-1\s+\(\s*('+RE_FLOAT+r')\s*nm\s*\)\s+Eps\s*=\s*('+RE_FLOAT+r')\s*M-1\s*cm-1;\s*Integration\s+limits:\s*('+RE_FLOAT+r')\s*\-\s('+RE_FLOAT+r')\s*nm\s*$', line)
		if (matched):
			#Absorption spectrum parameters
			AbsMaximumCM1 = 19990
			AbsMaximumNM  = 500.3
			Eps = 92300.0
			AbsIntLim1 = 400
			AbsIntLim2 = 600
			AbsorptionTabSpecOkay = True #Indicates presence of the abs. spectrum in the output, following flags check the data 
			ErrorMSG = r"Problem reading tabulated ABSORPTION spectrum of "+FrgName
			if (str(matched.group(1)).strip() != FrgName):
				assert False, ErrorMSG+r", Fragment name error"
			if ( abs(float(matched.group(2)) - AbsMaximumCM1)  > FLOAT_MATCH_TOL):
				assert False, ErrorMSG+r", maximum in CM-1 error"
			if ( abs(float(matched.group(3)) - AbsMaximumNM) > FLOAT_MATCH_TOL):
				assert False, ErrorMSG+r", maximum in NM error"
			if ( abs(float(matched.group(4)) - Eps) > FLOAT_MATCH_TOL):
				assert False, ErrorMSG+r", Molar extinction (epsilon) error"
			if ( abs(float(matched.group(5)) - AbsIntLim1) > FLOAT_MATCH_TOL):
				assert False, ErrorMSG+r", Lower limit of integration error"
			if ( abs(float(matched.group(6)) - AbsIntLim2) > FLOAT_MATCH_TOL):
				assert False, ErrorMSG+r", Upper limit of integration error"
			
		#Fragment: FLR1 tabulated normalized EMISSION spectrum is available. Maximum: 18518.52  cm-1 ( 540.0 nm ) I= 1.0  Arb. Units; Integration limits: 500 - 700 nm
		matched = re.match(r'^\s*Fragment:\s+(\w+)\s+tabulated\s+normalized\s+EMISSION\s+spectrum\s+is\s+available.\s+Maximum:\s+('+RE_FLOAT+r')\s+cm-1\s+\(\s*('+RE_FLOAT+r')\s*nm\s*\)\s+I\s*=\s*('+RE_FLOAT+r')\s*Arb.\s*Units;\s*Integration\s+limits:\s*('+RE_FLOAT+r')\s*\-\s('+RE_FLOAT+r')\s*nm\s*$', line)
		if (matched):
			#Emission spectrum parameters
			EmsMaximumCM1 = 18518.52
			EmsMaximumNM  = 540.0
			EmsI = 1.0
			EmsIntLim1 = 500
			EmsIntLim2 = 700
			FluorescenceTabSpecOkay = True #Indicates presence of the ems. spectrum in the output, following flags check the data 
			ErrorMSG = r"Problem reading tabulated EMISSION spectrum of "+FrgName
			if (str(matched.group(1)).strip() != FrgName):
				assert False, ErrorMSG+r", Fragment name error"
			if ( abs(float(matched.group(2)) - EmsMaximumCM1)  > FLOAT_MATCH_TOL):
				assert False, ErrorMSG+r", maximum in CM-1 error"
			if ( abs(float(matched.group(3)) - EmsMaximumNM) > FLOAT_MATCH_TOL):
				assert False, ErrorMSG+r", maximum in NM error"
			if ( abs(float(matched.group(4)) - EmsI) > FLOAT_MATCH_TOL):
				assert False, ErrorMSG+r", maximal intensity error"
			if ( abs(float(matched.group(5)) - EmsIntLim1) > FLOAT_MATCH_TOL):
				assert False, ErrorMSG+r", Lower limit of integration error"
			if ( abs(float(matched.group(6)) - EmsIntLim2) > FLOAT_MATCH_TOL):
				assert False, ErrorMSG+r", Upper limit of integration error"
		#Quantum yield:  0.97
		matched = re.match(r'^\s*Quantum yield:\s*('+RE_FLOAT+r')\s*$', line)
		if (matched):
			QunatumYieldOKay = True #Indicates presence of the excited states information 
			QY = 0.97
			if ( abs(float(matched.group(1)) - QY)  > FLOAT_MATCH_TOL):
				assert False, r"Reading quantum yield problem"

		ErrorMSG = r"Strickler-Berg calculation, "
		#Accounts for the quantum  yield:  tau_1(obs) = 3.98E-09 s - this value will be used in further calculations...		
		matched = re.match(r'^\s*Accounts\s+for\s+the\s+quantum\s+yield:\s+tau_1\(obs\)\s*=\s*('+RE_DOUBLE+r')\s*s\s*\-\s*this\s+value\s+will\s+be\s+used\s+in\s+further\s+calculations\.\.\.\s*$', line)
		if (matched):
			StricklerBergOkay = True #Indicates presence of the Strickler-Berg calculation 
			tau1 = 3.98E-9
			if ( abs(float(matched.group(1)) - tau1)*1.0E9  > FLOAT_MATCH_TOL):
				assert False, ErrorMSGr+"problem with tau1(obs)"
		
		#Accounts for the quantum  yield:  tau_2(obs) = 4.37E-09 s
		matched = re.match(r'^\s*Accounts\s+for\s+the\s+quantum\s+yield:\s+tau_2\(obs\)\s*=s*('+RE_DOUBLE+r')s*s\s*$', line)
		if (matched):
			tau2 = 4.37E-9
			if ( abs(float(matched.group(1)) - tau2)*1.0E9  > FLOAT_MATCH_TOL):
				assert False, ErrorMSGr+"problem with tau2(obs)"
		
		ErrorMSG = r"Direct calculation of Forster rates, "
		#Spectral overlap: 1.29E+14 M-1 cm-1 nm^4
		matched = re.match(r'^\s*Spectral\s+overlap:\s*('+RE_DOUBLE+r')\s*M-1\s*cm-1\s*nm\^4\s*$', line)
		if (matched):
			Ovl = 1.29E14
			if ( abs(float(matched.group(1)) - Ovl)/1.0E14  > FLOAT_MATCH_TOL):
				assert False, ErrorMSGr+"problem with the spectral overlap"

		#Forster radius: 3.59 nm
		matched = re.match(r'^\s*Forster\s+radius:\s*('+RE_FLOAT+r')\s*nm\s*$', line)
		if (matched):
			FRadius = 3.59
			if ( abs(float(matched.group(1)) - FRadius )  > FLOAT_MATCH_TOL):
				assert False, ErrorMSGr+"problem with the Forster radius"

		#b = Lambda/(2 PI n) = 63.14 nm
		matched = re.match(r'^\s*b\s*=\s*Lambda\/\(2\s*PI\s*n\)\s*=\s*('+RE_FLOAT+r')\s*nm\s*$', line)
		if (matched):
			Radb = 63.14
			if ( abs(float(matched.group(1)) - Radb )  > FLOAT_MATCH_TOL):
				assert False, ErrorMSGr+"problem with the radiation zones calculation"

		#Forster (non-radiative) rate: 2.51E+08 s-1
		matched = re.match(r'^\s*Forster\s*\(non\-radiative\)\s*rate:\s*('+RE_FLOAT+r')\s*s\-1\s*$', line)
		if (matched):
			ForsterRate =  2.51E+08 
			if ( abs(float(matched.group(1)) - Radb )  > FLOAT_MATCH_TOL):
				assert False, ErrorMSGr+"problem with the radiation zones calculation"

		#Done.
		matched = re.match(r'^\s*Done\.\s*$', line)
		if (matched):
			CalculationFinishedOkay = True

	assert ConfigOkay,r"Fluorescence lifetime configuraton file/fragments reading problem..."
	assert AbsorptionTabSpecOkay,r"Fluorescence lifetime tabulated ABSORPTION spectrum reading problem..."
	assert QunatumYieldOKay,r"Reading excited states problem..."
	assert StricklerBergOkay,r"Strickler-Berg calculation problem..."
	assert CalculationFinishedOkay,r"Calculation has NOT finished..."

#-------------------------------------------------------------------------------
def test_pyfrec_lifetimes(capsys):
	"""Test PyFREC lifetimes calculations using sample inputs from test_pyfrec directory """

	try:
		#Print program header
		MyReportService().Header(PROGRAM_HEADER)
		#Load configuration file (data files, computational models, approximations, basis sets etc.) and call CalculationManager
		CalcManager(Configuration(TEST_PYFREC_LIFETIMES)).Run()
	except RuntimeError as e:
		num, message = e.args
		out, err = capsys.readouterr()
		ParseLifetimeOutput(out)
		assert False, "Error: "+str(num)+message
	except SystemExit as e:
		out, err = capsys.readouterr()
		ParseLifetimeOutput(out)
		assert e.code==0, "Fluorescence lifetime calculation "+TEST_PYFREC_LIFETIMES+" has failed with the code: "+str(e.code)

	out, err = capsys.readouterr()
	ParseLifetimeOutput(out)
#-------------------------------------------------------------------------------
def test_pyfrec_couplings(capsys):
	"""Test PyFREC lifetimes calculations using sample inputs from test_pyfrec directory """

	try:
		#Print program header
		MyReportService().Header(PROGRAM_HEADER)
		#Load configuration file (data files, computational models, approximations, basis sets etc.) and call CalculationManager
		CalcManager(Configuration(TEST_PYFREC_COUPLINGS)).Run()
	except RuntimeError as e:
		num, message = e.args
		out, err = capsys.readouterr()		
		ParseCouplingOutput(out)
		assert False, "Error: "+str(num)+message 
	except SystemExit as e:
		out, err = capsys.readouterr()
		ParseCouplingOutput(out)
		assert e.code==0, "Electronic coupling calculation "+TEST_PYFREC_COUPLINGS+" has failed with the code: "+str(e.code)
	
	out, err = capsys.readouterr()
	ParseCouplingOutput(out)
#-------------------------------------------------------------------------------
