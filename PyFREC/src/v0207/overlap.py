from couplib.constants import *
from configuration import *
from couplib.lookup import *
#from interfaces import OverlapInterface
from scipy.integrate import quad
from scipy.integrate import quad_explain
from interfaces import *
import numpy as np
import math
from exstatesreader import ExStatesReader

#Nummerical integration absolute and relative errors
INT_ERR_EPS_ABS = 1.0e-40 
INT_ERR_EPS_REL = 1.0e-3
INT_LIMIT = 250 #Number of subdivision in numerical integration

#-------------------------------------------------------------------------------
class Overlap(object):
	"""Class for calculation of spectral overlaps for resonance and rate calculations"""

	def __init__(self, cfg):
		"""Initialize with current configuration"""
		self.cfg = cfg
		#Overlap Integration limits (cm-1 scale):
		self.int_lim_low = float(self.cfg.MolSys[CFG_MSY_OVL_LOW_LIM])
		self.int_lim_upp = float(self.cfg.MolSys[CFG_MSY_OVL_UPP_LIM])
		
		#Overlap Integration limits (nm scale):
		self.int_lim_low_nm = float(self.cfg.MolSys[CFG_MSY_OVL_LOW_LIM_NM])
		self.int_lim_upp_nm = float(self.cfg.MolSys[CFG_MSY_OVL_UPP_LIM_NM])

		return
#-------------------------------------------------------------------------------
	def TabDiv(self, x, Wavenum_cm1_ta, Abs_M1cm1, pown):
		"""Value of the tabulated spectrum (Wavenum_cm1_ta, Abs_M1cm1) multiplied by 1/x^pown at the point x cm-1"""

		Tab_w = lookup(x, Wavenum_cm1_ta, Abs_M1cm1)

		dp = np.power(x,pown)
		if ( dp == 0.0 ): return np.inf;
		Tab_w = Tab_w/dp
		
		return Tab_w
#-------------------------------------------------------------------------------
	def TabMult(self, x, Ems_nm_tab, Ems_Val, Abs_nm_tab, Abs_M1cm1, pown):
		"""Value of the product of emission and absorption spectra multiplied by the power of the wavelenght sampled at the point x"""
		
		w = lookup(x, Ems_nm_tab, Ems_Val)*lookup(x, Abs_nm_tab, Abs_M1cm1)
		dp = np.power(x,pown)
		w = w*dp
		return w
#-------------------------------------------------------------------------------
	def Overlap(self, ExcitedState1, ExcitedState2, Resonace_Threshold_Gaussian):
		"""Class for calculation of spectral overlaps for resonance and rate calculations"""
		if (Resonace_Threshold_Gaussian == 0.0):
			#Empty overlap
			return OverlapInterface(False)
		
		#ExcitedState1 - Donor (emission spectrum is used)
		Ems_cm1 = ExcitedState1.Ems_cm1 #Emission maximum 
		Ems_a1 = 0.0 #Intensities of peaks are zeroes (normalization condition is used)
		
		#ExcitedState2 - Acceptor (absorption spectrum is used)
		Abs_cm1 = ExcitedState2.Abs_cm1 #Absorption maximum 

		#Molar exctinction coefficient of acceptor M-1 cm-1
		Epsilon_M1cm1 = ExcitedState2.Epsilon_M1cm1

		#print("Ex. States Donor(Emission):",str(ExcitedState1.ExStID),"Acceptor(Absorption):",str(ExcitedState2.ExStID))
		#print("Ex. Energies (cm-1):", Ems_cm1, Abs_cm1)
		#print("FWHM (cm-1):", Ems_FWHM_cm1, Abs_FWHM_cm1)

		#First state is assumed to be a donor (emission spectrum) and
		#Second state is assumed to be an acceptor (absorption spectrum)

		Overlap_M1cm1nm4 = 0.0

		if (LEGACY_LEVEL == LL_NONE): #Default
			#Run "tabulated" overlap
			if ( (len(ExcitedState1.Ems_Spec) > 0 ) and len(ExcitedState2.Abs_Spec) > 0):
				print("Empirical tabulated spectra for the donor and acceptor are provided numerical overlap will be calculated...")

				Overlap = self.TabOverlap(ExcitedState1, ExcitedState2)
				Overlap_M1cm1nm4 = Overlap.Overlap_M1cm1nm4
		else:
			print("Error: Unknown legacy level: LEGACY_LEVEL = ",LEGACY_LEVEL)
			exit(-1)
		
		return OverlapInterface(Overlap_M1cm1nm4)
#-------------------------------------------------------------------------------
	def Int_Absorption_Over_Wavenumber_Tab(self, ExcitedState):
		"""Integration of the molar extinction coefficient (absorption) in [M-1 cm-1] over cm-1 wavenubmer domain"""
		
		#print("Tabulated interpolation!")
		#print(ExcitedState.Abs_Spec)
		#Work only if the tabulated spectrum is available
		if (len(ExcitedState.Abs_Spec) < 1):
			return 0.0
		
		print()
		print("Absorption (tabulated spectrum):")
		
		#Wavenumber values
		Wavenum_cm1_ta = [i[0] for i in ExcitedState.Abs_Spec]
		#print(Wavenum_cm1_ta)

		#Correponding absorbance values
		Abs_M1cm1 = [i[1] for i in ExcitedState.Abs_Spec]
		#print(Abs_M1cm1)
		
		print("Numerical integration tolerances relative:",INT_ERR_EPS_REL)
		self.int_lim_low = CM1_NM/ExcitedState.Abs_Int_Lim_Up_nm
		self.int_lim_upp = CM1_NM/ExcitedState.Abs_Int_Lim_Low_nm

		print("Integration limits:",self.int_lim_low," ",self.int_lim_upp, "cm-1")
		pow_x = 1 #Integral 1/frec_cm1
		
		#print("Interpolation:")
		#for w in np.linspace(self.int_lim_low,self.int_lim_upp,400):
			#Abs_w = lookup(w, Wavenum_cm1_ta, Abs_M1cm1)
		#	Abs_w =self.TabDiv(w, Wavenum_cm1_ta, Abs_M1cm1, pow_x)
		#	print(w,Abs_w)
		#print("-------------")
		#quad_explain()
		#epsabs=INT_ERR_EPS_ABS
		#epsrel=INT_ERR_EPS_REL
		(Abs, Err) = quad(self.TabDiv, self.int_lim_low, self.int_lim_upp, args=(Wavenum_cm1_ta, Abs_M1cm1, pow_x),epsrel=INT_ERR_EPS_REL,limit=INT_LIMIT)
		
		print("Integral int[a(nu)/nu]d[nu] = ",Abs)
		Abs_Ext = Abs
		print("Epsilon x Integral int[a(nu)/nu]d[nu] = ",Abs_Ext)
		print("Actual integration error: ",Err,"(",100*(Err/Abs),"%)")
		
		return Abs_Ext
#-------------------------------------------------------------------------------
	def Norm_Int_Emission_Over_Wavenumber_Tab(self, ExcitedState):
		"""Integration of the emission in arb. units over cm-1 wavenubmer domain"""
		
		#print("Tabulated interpolation!")
		#print(ExcitedState.Ems_Spec)
		#Work only if the tabulated spectrum is available
		if (len(ExcitedState.Ems_Spec) < 1):
			return 0.0
		
		print()
		print("Emission (tabulated spectrum):")
		
		#Wavenumber values
		Wavenum_cm1_ta = [i[0] for i in ExcitedState.Ems_Spec]
		#print(Wavenum_cm1_ta)

		#Correponding emission values
		Ems_Val = [i[1] for i in ExcitedState.Ems_Spec]
		#print(Ems_Val)

		#print("-------------")
		print("Numerical integration tolerances relative:",INT_ERR_EPS_REL)
		self.int_lim_low = CM1_NM/ExcitedState.Ems_Int_Lim_Up_nm
		self.int_lim_upp = CM1_NM/ExcitedState.Ems_Int_Lim_Low_nm

		print("Integration limits cm-1",self.int_lim_low," ",self.int_lim_upp)
		#-----
		#print("Interpolation:")
		#for w in np.linspace(self.int_lim_low,self.int_lim_upp,200):
		#	Ems_w = lookup(w, Wavenum_cm1_ta, Ems_Val)
		#	print(w,Ems_w)
		
		(Norm, Err) = quad(lookup, self.int_lim_low, self.int_lim_upp, args=(Wavenum_cm1_ta, Ems_Val), epsrel=INT_ERR_EPS_REL,limit=INT_LIMIT)
		print("Integral int[e(nu)d[nu] = ",Norm)
		print("Actual integration error:",Err,"(",100*(Err/Norm),"%)")

		pow_x = 3 #Integral 1/frec_cm1^3
		(Ems, Err) = quad(self.TabDiv, self.int_lim_low, self.int_lim_upp, args=(Wavenum_cm1_ta, Ems_Val, pow_x),epsrel=INT_ERR_EPS_REL,limit=INT_LIMIT)

		print("Integral int[e(nu)/nu^3]d[nu] = ",Ems)
		print("Actual integration error:",Err,"(",100*(Err/Ems),"%)")
		
		#Normalization of the emission integral:
		Ems_Norm = Ems/Norm
		#print("Noramalized integral of emission:",Ems_Norm)
		
		print("Integral (int[e(nu)/nu^3]d[nu])/(int[e(nu)]dnu) = ",Ems_Norm)
		#Inv_Ems_Norm = 1.0/Ems_Norm
		#print("Inverted integral of emission: ",Inv_Ems_Norm)
		return Ems_Norm
#-------------------------------------------------------------------------------
	def TabOverlap(self, ExcitedState_Donor, ExcitedState_Acceptor):
		"""Overlap of the tabulated abrorption and emission spectra: normalized donor emission and molar absorption in M-1cm-1 of the acceptro with Lambda^4 factor (numerical)"""
		#Convert spectra from the wavenumber to wavelenght domain
		Ems_Spec = [ (CM1_NM/WL, I) for WL, I in ExcitedState_Donor.Ems_Spec]
		Ems_Spec.reverse()

		Abs_Spec = [ (CM1_NM/WL, A) for WL, A in ExcitedState_Acceptor.Abs_Spec]
		Abs_Spec.reverse()

		#Waveleghts of the emission spectrum
		Ems_nm_tab = [i[0] for i in Ems_Spec]

		#Emission integration limits
		Ems_Int_Lim_Low_nm = ExcitedState_Donor.Ems_Int_Lim_Low_nm
		Ems_Int_Lim_Up_nm = ExcitedState_Donor.Ems_Int_Lim_Up_nm
		#Ems_Int_Lim_Low_nm = min(Ems_nm_tab)
		#Ems_Int_Lim_Up_nm = max(Ems_nm_tab)

		#Intensity of the emission spectrum
		Ems_Val = [i[1] for i in Ems_Spec]
		
		WAVELENROUND = 1
		print("Numerical integration relative tolerance:",INT_ERR_EPS_REL)
		print()
		print("Integration limits:")
		print("Tabulated emission spectrum:                    ",round(Ems_Int_Lim_Low_nm, WAVELENROUND),"-",round(Ems_Int_Lim_Up_nm, WAVELENROUND),"nm")
		#-----
		#print("Interpolation EMS:")
		#for w in np.linspace(self.int_lim_low,self.int_lim_upp,200):
		#	Ems_w = lookup(w, Wavenum_cm1_ta, Ems_Val)
		#	print(w,Ems_w)
		
		#Absorbance
		Abs_nm_tab = [i[0] for i in Abs_Spec]
		Abs_M1cm1  = [i[1] for i in Abs_Spec]

		#Absorption integration limits
		#Abs_Int_Lim_Low_nm = min(Abs_nm_tab)
		#Abs_Int_Lim_Up_nm = max(Abs_nm_tab)
		Abs_Int_Lim_Low_nm = ExcitedState_Acceptor.Abs_Int_Lim_Low_nm
		Abs_Int_Lim_Up_nm = ExcitedState_Acceptor.Abs_Int_Lim_Up_nm
		
		print("Tabulated absorption spectrum:                  ",round(Abs_Int_Lim_Low_nm, WAVELENROUND),"-",round(Abs_Int_Lim_Up_nm, WAVELENROUND),"nm")

		#Spectral overlap integration limits
		Int_Lim_Low_nm = max(Abs_Int_Lim_Low_nm,Ems_Int_Lim_Low_nm)
		Int_Lim_Up_nm  = min(Abs_Int_Lim_Up_nm,Ems_Int_Lim_Up_nm)
		print("Overlap:                                        ",round(Int_Lim_Low_nm, WAVELENROUND),"-",round(Int_Lim_Up_nm, WAVELENROUND),"nm")

		(EmsNorm, Err) = quad(lookup, Ems_Int_Lim_Low_nm, Ems_Int_Lim_Up_nm, args=(Ems_nm_tab, Ems_Val), epsrel=INT_ERR_EPS_REL,limit=INT_LIMIT)
		
		print()
		FD = 2
		print("EMS = int[F_D(Lambda)]d(Lambda) =",format(EmsNorm,'.'+str(FD)+'E'))
		print("Actual integration error:",format(Err,'.'+str(FD)+'E'),"(",format(100*(Err/EmsNorm),'.'+str(FD)+'E'),"%)")

		pow_n = 4 #Integral Lambda^4
		print
		
		#Step_nm = 1
		#print("Tabulated integration: Lambda, Lambda^pown, Ems, Abs, Ems*Abs*Lambda^pown")
		#for Lambda_i in range(int(Int_Lim_Low_nm),int(Int_Lim_Up_nm), Step_nm):
		#for Lambda_i in range(251,847, Step_nm):
			#Lambda = float(Lambda_i)
			#dp = np.power(Lambda,pow_n)
			#Ems = lookup(Lambda, Ems_nm_tab, Ems_Val)
			#Abs = lookup(Lambda, Abs_nm_tab, Abs_M1cm1)
			#Mult = self.TabMult(Lambda, Ems_nm_tab, Ems_Val, Abs_nm_tab, Abs_M1cm1, pow_n)
			#print(Lambda, dp, Ems, Abs, Mult)
			#print(Lambda, Ems)
			#print(Lambda, Abs)
		
		(EmsAbsLambda4, Err) = quad(self.TabMult, Int_Lim_Low_nm, Int_Lim_Up_nm, args=(Ems_nm_tab, Ems_Val, Abs_nm_tab, Abs_M1cm1, pow_n),epsrel=INT_ERR_EPS_REL,limit=INT_LIMIT)
		print("ABSEMS = int[epsilon(Lambda)xF_D(Lambda)xLambda^4]d(Lambda) =",format(EmsAbsLambda4,'.'+str(FD)+'E'))
		if (EmsAbsLambda4 != 0):
			print("Actual integration error:",Err,"(",format(100*(Err/EmsAbsLambda4),'.'+str(FD)+'E'),"%)")

		#Normalization of the emission integral:
		EmsAbsLambda4_Norm = EmsAbsLambda4/EmsNorm
		#format(,'.'+str(FD)+'E')
		print("Resulting overlap: ABSEMS/EMS =",format(EmsAbsLambda4_Norm,'.'+str(FD)+'E'),"M-1 cm-1 nm^4")
		print()
		Overlap_M1cm1nm4 = EmsAbsLambda4_Norm

		return  OverlapInterface(Overlap_M1cm1nm4)
#-------------------------------------------------------------------------------