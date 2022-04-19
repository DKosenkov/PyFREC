from couplib.constants import *
from couplib.myreportservice import *
from configuration import *
from math import *
import numpy as np
from interfaces import *
from overlap import Overlap

class SBCalculation(object):
	"""Strickler-Berg calculation of lifetime of donors"""
#-------------------------------------------------------------------------------
	def __init__(self, cfg):
		"""Initialize with current configuration"""
		self.cfg = cfg
		self.CalcType = int(self.cfg.MolSys[CFG_MSY_SB])
		self.ps = MyReportService()
#-------------------------------------------------------------------------------
	def CalcLifetimesOfFragments(self):
		"""Strickler-Berg calculation of the donor fluorescence  lifetime of all fragments"""

		#Do not run SB calculation
		if ( self.CalcType == SB_NONE):
			return

		print()
		print("Strickler-Berg (SB) calculation of the donor fluorescence  lifetime...")

		if ( self.CalcType == SB_INTEGRATE):
			print("Calculation method: Integration of the emission spectrum...")
		elif ( self.CalcType == SB_MAX_ONLY):
			print("Calculation method: Simplified, without integration of the emission spectrum...")
		else:
			print(("Error unbknown value of Strickler-Berg (SB) calculation type: ", self.CalcType))
			exit(-1)
		
		#Refractive index
		n = float(self.cfg.MolSys[CFG_MSY_RFX])
		
		#print("Refractive index:", n)
		#print("Refractive index squared:", n*n)

		#Scan fragments
		for FragID, Fragment in list(self.cfg.Fragments.items()):
			#FileName = Fragment[CFG_FRG_GEF]
			#print((CFG_SEC_FRG+str(FragID+1)+" "+CFG_FRG_NAM+" = "+Fragment[CFG_FRG_NAM]+" "+CFG_FRG_ID+" = "+str(Fragment[CFG_FRG_ID])))
			#print((CFG_FRG_GEF+" = "+FileName))
			print("Fragment:",str(Fragment[CFG_FRG_NAM]).strip()+str(Fragment[CFG_FRG_ID]).strip())
			
			for ExcitedStateID, ExcitedState in list(Fragment[CFG_EXS_EXSTATE].items()):
				self.ps.PrintDiv()
				print("State   :",ExcitedState.ExStID)
				ExcitedState.FlLifetime_sb_s = self.GetLifetimesOfFragments(FragID, ExcitedState, n)
				
		return 
#-------------------------------------------------------------------------------
	def GetLifetimesOfFragments(self, FragID, ExcitedState, n):
		"""Strickler-Berg (SB) calculation of the donor fluorescence lifetime of a given fragment"""

		Lifetime_s = 0.0
		#Run SB calculation with the integration

		Phi_D = ExcitedState.Phi_D #Fluorescence quantum yield
		print("Quantum yield: ", Phi_D)
		
		#Overlap class to run integrals of spectra
		Ovl = Overlap(self.cfg)
		if (len(ExcitedState.Ems_Spec) < 1):			
			ExcitedState.Ems_Spec = Ovl.FillSpecWN(ExcitedState, ExcitedState.Ems_nm, ExcitedState.Ems_FWHM_nm, ExcitedState.Epsilon_M1cm1, CFG_EXS_SEC_EMS)
		if (len(ExcitedState.Abs_Spec) < 1):
			ExcitedState.Abs_Spec = Ovl.FillSpecWN(ExcitedState, ExcitedState.Abs_nm, ExcitedState.Abs_FWHM_nm, ExcitedState.Epsilon_M1cm1, CFG_EXS_SEC_ABS) 
		
		#Ingegral of the absorption spectrum multiplied by 1/frec_cm1
		Abs_Tab =  Ovl.Int_Absorption_Over_Wavenumber_Tab(ExcitedState)
		self.ps.PrintDiv()
		
		#self.ps.PrintDiv()		
		Ems_Norm_Tab = 0.0
		#Integration of emission spectrum
		if ( self.CalcType == SB_INTEGRATE):
			#Ingegral of the emission spectrum multiplied by 1/frec_cm1^3						
			Ems_Norm_Tab = Ovl.Norm_Int_Emission_Over_Wavenumber_Tab(ExcitedState)
			self.ps.PrintDiv()

		#Simplified calculation is always performed
		#print("Simplified emission value 1/(nu_cm-1_max)^3:")
		Ems_Simple = pow(ExcitedState.Ems_cm1,3.0)
		#print("Emission cm3: ",ExcitedState.Ems_cm1)
		#print("Emission factor cm3: ",Ems_Simple)
		
		Stk_cm1 = ExcitedState.Abs_cm1 - ExcitedState.Ems_cm1
		
		#print("Stokes Shift:",Stk_cm1," cm-1 or ")
		#Pre-integral factor: 2.880E-9  (The factor of 100 at the end of the numerator is most likely due to conversion of the speed of light from m/s to cm/s)
		PIF = (8.0*np.log(10)*1000.0*math.pi*c*100)/AVOGADROS_Na
		#print("DEBUG: Pre-integral factor (must be 2.880E-9): ",PIF)
		print("Strickler-Berg Calculation summary:")
		print()
		Format_time = 2
		#Simplified Lifetime in s-1:
		if ( ExcitedState.Ems_cm1 != 0 ):		
			print("Simplified emission calculation based on the tabulated absorption lineshape:")
			#print("PIF*n*n*Ems_Simple*Abs_Tab:")
			#print(PIF,n,n,Ems_Simple,Abs_Tab)
			Simple_Tab_Rate_s1 = PIF*n*n*Ems_Simple*Abs_Tab
			Simple_Tab_tau_s = 1.0/Simple_Tab_Rate_s1
			Simple_Tab_tau_obs_s = Simple_Tab_tau_s*Phi_D
			print("tau_1 =",format(Simple_Tab_tau_s,'.'+str(Format_time)+'E'),"s (1/tau_1 =",format(Simple_Tab_Rate_s1,'.'+str(Format_time)+'E'),"1/s)")
			print("Accounts for the quantum  yield:  tau_1(obs) =",format(Simple_Tab_tau_obs_s,'.'+str(Format_time)+'E'),"s - this value will be used in further calculations...")
			#format(,'.'+str(Format_time)+'E')
			self.ps.PrintDiv()
		if (( self.CalcType == SB_INTEGRATE) and (ExcitedState.Ems_cm1 != 0) ):
			print("Integration based on tabulated spectra:")
			Rate_Tab_s1 = PIF*n*n*Abs_Tab/Ems_Norm_Tab
			tau_tab_s = 1.0/Rate_Tab_s1
			tau_tab_obs_s = tau_tab_s*Phi_D
			print("tau_2 =",format(tau_tab_s,'.'+str(Format_time)+'E'),"s (1/tau_2 =",format(Rate_Tab_s1,'.'+str(Format_time)+'E'),"1/s)")
			#print("Accounts for the quantum  yield:  tau_2(obs) =",format(tau_tab_obs_s,'.'+str(Format_time)+'E'),"s - this value will be used in further calculations...")
			print("Accounts for the quantum  yield:  tau_2(obs) =",format(tau_tab_obs_s,'.'+str(Format_time)+'E'),"s")
			return Simple_Tab_tau_obs_s
			#return tau_tab_obs_s  #this is tau_4(obs)
		elif( (self.CalcType == SB_MAX_ONLY) and (ExcitedState.Ems_cm1 != 0) ):
			print("Simplified lifetime will be used for calculation of Forster rates")
			print("if the emprical value is not specified in the input")
			return Simple_tau_obs_s
		else:
			print("The lifetime was not calculated (no emission informaiton was provided)...")
			return 0.0

#-------------------------------------------------------------------------------