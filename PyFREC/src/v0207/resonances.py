from couplib.constants import *
from configuration import *
from interfaces import ResonanceInterface
from overlap import Overlap
import math

#-------------------------------------------------------------------------------
class Resonances(object):
	"""Class for identification of resonances of excited states of fragments"""

	def __init__(self):
		"""Initialize with current configuration"""
		return
	
	def Resonance(self, ExcitedState1, ExcitedState2, ResonanceCondition_cm1, Resonace_Threshold, OverapSummary):
		"""Resonance condition based on either frequency difference or Gaussian overlap"""
		
		Diff = ExcitedState2.Abs_cm1-ExcitedState1.Abs_cm1

		if ( Resonace_Threshold == 0.0):
			print("Resonances are computed based on a difference of absorption maxima")
			Overlap = 0.0
			Flag = True if (abs(Diff)<= float(ResonanceCondition_cm1)) else False
			return ResonanceInterface(ExcitedState1.Abs_cm1, ExcitedState2.Abs_cm1, Diff, Overlap, Flag)
		else:
			print("Resonances are computed based on spectral overlaps")
		
		Overlap = 0.0
		
		if (LEGACY_LEVEL == LL_NONE): #Default 
			#Overlap based on the normalizad emission of the donor and molar absorption (M-1 cm-1) of the acceptor with Lambda^4 factor
			Overlap = OverapSummary.Overlap_M1cm1nm4
		else:
			print(("Error: Unknown legacy level: LEGACY_LEVEL = ",LEGACY_LEVEL))
			exit(-1)

		Flag = True if (abs(Overlap)>= float(Resonace_Threshold)) else False
		
		return ResonanceInterface(ExcitedState1.Abs_cm1, ExcitedState2.Abs_cm1, Diff, Overlap, Flag)


