from sbcalc import SBCalculation
from pdbreader import PDBReader

from exstatesreader import ExStatesReader
from transformfrags import TransformFrags
from coupling import Coupling
from varmethod import VariationMethod

from configuration import *

if (INT_CONFIG_QD_ON):
	from qd import QuantumDynamics

if (INT_CONFIG_TRJ_PDB_ON):
	from trjreader import TRJReader
#-------------------------------------------------------------------------------
class CalcManager():
	"""Calculation manager which decides how to perform requested calculation"""

	def __init__(self,cfg):
		"""Initialization"""
		self.cfg = cfg #Data from the configuration file (.ini)
		return

	def Run(self):
		# Determine job type
		JobType = self.cfg.Methods[CFG_MET_APX]

		#Couplings
		cou = None
		#Excited state energies and couplings
		Hamiltonian = None

		#Make fragments only if the input file was provided by the user
		if (self.cfg.MolSys[CFG_MSY_GEO] != CFG_MSY_GEO_NONE):
			pdb = PDBReader(self.cfg)
			pdb.Read()
			pdb.PrintFragmens()
		
		if (INT_CONFIG_TRJ_PDB_ON):
			#Read trajectory file (if any)
			trj =  TRJReader(self.cfg)
			NMDFrames = trj.Read()
		else:
			NMDFrames = 0
		
		#Process trajectory only for Forster jobs for now
		if ( (NMDFrames > 0) and (JobType == CFG_MET_ARX_FRS) ):
			#for FrameID, Frame in trj.Frames:
			for FrameID in range(NMDFrames):
				print(("Process MD Frame = ",FrameID))
				print("Reading excited states file...")
				#Process current frame
				exs = ExStatesReader(self.cfg)
				exs.Read()
				#Assign target excited state for the Forster calculation
				exs.AssignTargetExcitedStates()
				#Calculate Strickler-Berg lifetimes
				SBCalculation(self.cfg).CalcLifetimesOfFragments()

				if (self.cfg.Methods[CFG_MET_VERBOSITY] < CFG_MET_VERB_TRJ):
					exs.PrintFragmens()
				print("Transforming fragments...")
				TransformFrags(self.cfg).Run()
				#Prepare to compute Forster couplings
				cou = Coupling(self.cfg)
				print("Compute Forster couplings and form super-molecular Hamiltonian...")
				#Couplings
				Hamiltonian = cou.Run()
				cou.PrintCouplings()
				trj.LoadFrame(FrameID)
		else:
			if(JobType == CFG_MET_ARX_LFT):
				#Lifetime calculations
				print("Calculation of lifetimes...")
				print("Reading excited states file...")
				exs = ExStatesReader(self.cfg)
				exs.Read()
				exs.PrintFragmens()
				#Calculate Strickler-Berg lifetimes
				SBCalculation(self.cfg).CalcLifetimesOfFragments()
				
				#If user specified kappa squared run homotransfer properties
				UserDefinedKappaSq = float(self.cfg.MolSys[CFG_MSY_KAPPA_SQ])
				if ( UserDefinedKappaSq != CFG_MSY_KAPPA_SQ_DEFAULT):
					#Calculate Forster theory-related properties of fragments for homotransfer
					cou = Coupling(self.cfg)
					cou.HomotransferProperties()

			#Forster coupling or survey calculation
			if( (JobType == CFG_MET_ARX_FRS) or (JobType == CFG_MET_ARX_SUR) ):
				#Read standard orientation geometry, exciation energies, and transition dipole moments from EXSTATE_FILE
				print("Reading excited states file...")
				exs = ExStatesReader(self.cfg)
				exs.Read()

				#Assign target excited state for Forster calculation
				if (JobType == CFG_MET_ARX_FRS):
					exs.AssignTargetExcitedStates()

				#Calculate Strickler-Berg lifetimes
				SBCalculation(self.cfg).CalcLifetimesOfFragments()
				exs.PrintFragmens()
				#Transform coordinates ess->pdb
				print("Transforming fragments...")
				TransformFrags(self.cfg).Run()
				#Prepare to compute Forster couplings
				cou = Coupling(self.cfg)

				if( JobType == CFG_MET_ARX_SUR ):
					print("Survey Forster couplings and rates...")
					#Survey Couplings
					cou.RunSurvey()
					cou.PrintSurvey()
				#Forster coupling
				elif( JobType == CFG_MET_ARX_FRS ):
					print("Compute Forster couplings and form super-molecular Hamiltonian...")
					#Couplings
					Hamiltonian = cou.Run()
					cou.PrintCouplings()
					print("Variational calculation...")
					#Calculates orbitals (C) and their energies (E) using variational principle solving an eigenvalue problem in the Huckel theory
					vm = VariationMethod(Hamiltonian)
					(Oribtals, Energies) = vm.Calc()

					QuantumDynamicsType = self.cfg.Methods[CFG_MET_QD]
					
					if (QuantumDynamicsType == CFG_MET_QD_MARKOVIAN):
						print("Markovian Quantum Dynamics...")
						qd = QuantumDynamics(self.cfg)
						qd.Calc(Hamiltonian)
					
		print("Done.")
		return
#-------------------------------------------------------------------------------