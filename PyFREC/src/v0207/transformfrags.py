from couplib.autoarr import AutoArr
from couplib.myreportservice import *
from couplib.vrrotmat2vec import Vrrotmat2Vec

from couplib.parselist import ParseList
import numpy as np
import math
import copy

from configuration import *
from interfaces import AtomInterface
from interfaces import ExciteStateInterface
from interfaces import OriginInterface

#-------------------------------------------------------------------------------
class TransformFrags():
	"""Class transforms coordinates of fragments"""

	def __init__(self, cfg):
		"""Initialize with current configuration"""
		self.cfg = cfg
		self.ps = MyReportService()
		self.TrnsFragments = AutoArr() #Fragments used for matching coordinates (may have less atoms than fragemtns in cfg.Fragemtns)
		return

	def GeometryDiff(self, FragmentPDB, FragmentEXS, Filter = 0):
		"""Computes RMS difference between two sets of atoms FragmentPDB (must have AtomPDBID) and FragmentB selecting only allowed by filter"""

		if (len(FragmentPDB) == 0):
			print("Error: cannot calculated RMS for fragments the number of atoms is zero!")
			exit(-1)
		natoms = 0
		sumsq = 0.0
		CfgAtomIDList = ParseList(Filter)
		for AtomID, Atom in FragmentPDB.items():
			AtomEx = FragmentEXS[AtomID]
			if ( (Atom.AtomPDBID in CfgAtomIDList) or (CfgAtomIDList == [0]) ):
				sumsq += ((Atom.x-AtomEx.x)**2+(Atom.y-AtomEx.y)**2+(Atom.z-AtomEx.z)**2)
				natoms += 1
		rms = math.sqrt(sumsq/natoms)
		return (sumsq, rms, natoms)

#-------------------------------------------------------------------------------
	def Procrustes(self, X, Y):
		"""
		X - fragment to be transformed
		Y - coordiantes from PDB (target)
		"""
	
		muX = X.mean(0)
		muY = Y.mean(0)
	
		X0 = X - muX
		Y0 = Y - muY
	
		ssX = (X0**2.).sum()
		ssY = (Y0**2.).sum()
	
		normX = np.sqrt(ssX)
		normY = np.sqrt(ssY)
	
		X0 /= normX
		Y0 /= normY
	
		# optimum rotation matrix of Y
		A = np.dot(X0.T, Y0)
		U,s,Vt = np.linalg.svd(A,full_matrices=False)
		V = Vt.T
		T = np.dot(V, U.T)
	
		traceTA = s.sum()
	
		b = traceTA * normX / normY
		c = muX - np.dot(muY, T)
		#d = 1 - traceTA**2
	
		return (T, b, c)
#-------------------------------------------------------------------------------

	def TransformFrag(self, FragID, Fragment):
		"""Transform given fragment"""

		NCRD = 3 #Number of coordinates in 3D space (x,y,z)
		self.ps = MyReportService()

		if (Fragment[CFG_FRG_NATOMS] < NCRD):
			print("Error: in fragments should be at least 3 atoms (not linear) given "+str(Fragment[CFG_FRG_NATOMS])+" in "+CFG_SEC_FRG+str(FragID+1))
			print("(use additional dummy atoms 'X' for diatomic and linear molecules)")
			exit(-1)

		#Determine atoms that will be used for the transformation (CFG_FRG_ATL = "ATOMTRANS")
		TrnsFragment = AutoArr()
		TrnsFragment[CFG_FRG_NATOMS] = 0 #The number of atoms in the fragment (will be updated duing actual PDB reading)
		TrnsFragment[CFG_FRG_ATOMS]  = AutoArr()
		TrnsFragment[CFG_EXS_NATOMS] = 0 #The number of atoms in the fragment (CFG_FRG_GEF = "EXSTATE_FILE" )
		TrnsFragment[CFG_EXS_ATOMS]  = AutoArr()


		CfgAtomIDList = ParseList(Fragment[CFG_FRG_ATL])
		for AtomID, Atom in Fragment[CFG_FRG_ATOMS].items():
			AtomPDBID = Atom.AtomPDBID

			for CfgAtomID in CfgAtomIDList:
				if ( (CfgAtomID == 0) or (CfgAtomID==AtomPDBID) ):
					#Add atom from the pdb file
					TrnsFragment[CFG_FRG_ATOMS][TrnsFragment[CFG_FRG_NATOMS]] = Fragment[CFG_FRG_ATOMS][AtomID]
					TrnsFragment[CFG_FRG_NATOMS] += 1
					#Add corresponding atom from the CFG_FRG_GEF = "EXSTATE_FILE" to the fragment.
					#The order of atoms in EXSTATE_FILE has to match exactly the order of atoms in pdb file!
					TrnsFragment[CFG_EXS_ATOMS][TrnsFragment[CFG_EXS_NATOMS]] = Fragment[CFG_EXS_ATOMS][AtomID]
					TrnsFragment[CFG_EXS_NATOMS] += 1

		#Should be at least 3 atoms to run fragment matching
		if ( (TrnsFragment[CFG_FRG_NATOMS] < NCRD) ):
			print("Error: At least 3 atoms are required to match fragments. Given "+str(TrnsFragment[CFG_FRG_NATOMS])+" in "+CFG_SEC_FRG+str(FragID+1))
			exit(-1)

		#Transform coordinates
		self.ps.PrintDiv()
		print("Original coordinates used for matching (independent atom index):")
		#TrnsFragment from CFG_FRG_GEF = "EXSTATE_FILE"  file

		npExsAtoms = np.zeros(NCRD*TrnsFragment[CFG_EXS_NATOMS],dtype=np.dtype('d')).reshape(TrnsFragment[CFG_EXS_NATOMS],NCRD)

		for AtomID, Atom in TrnsFragment[CFG_EXS_ATOMS].items():
			Atom.MyPrint(AtomID+1)
			npExsAtoms[AtomID,0:NCRD] = np.array([Atom.x,Atom.y,Atom.z])
		print()
		print("Mol. Sys. (target) coordinates used for matching:")
		#Fragment from PDB file
		npFrgAtoms = np.zeros(NCRD*TrnsFragment[CFG_FRG_NATOMS],dtype=np.dtype('d')).reshape(TrnsFragment[CFG_FRG_NATOMS],NCRD)
		for AtomID, Atom in TrnsFragment[CFG_FRG_ATOMS].items():
			Atom.MyPrint()
			npFrgAtoms[AtomID,0:NCRD] = np.array([Atom.x,Atom.y,Atom.z])

		self.ps.PrintDiv()

		#npExsAtoms -> npFrgAtoms
		(T, b, c) = self.Procrustes(npFrgAtoms, npExsAtoms)
		
		#Save transformation globally
		self.cfg.Fragments[FragID][CFG_TRNS] = (T, b, c)

		#Transform coordinates globally (ignoring scaling unlike Z)
		for AtomID, Atom in Fragment[CFG_EXS_ATOMS].items():
			Pnt = np.array((Atom.x, Atom.y, Atom.z),dtype=np.dtype('d'))
			TrnsPnt = c + np.dot(Pnt , T)
			#with scaling: TrnsPnt = c + b*np.dot(Pnt, T)
			self.cfg.Fragments[FragID][CFG_TRNS_EXS_ATOMS][AtomID] = AtomInterface(TrnsPnt[0], TrnsPnt[1], TrnsPnt[2], Atom.ElSym)
		#Transform and save origin globally
		Orig = np.zeros(NCRD,dtype=np.dtype('d'))
		#Fragment center
		Center = self.cfg.Fragments[FragID][CFG_EXS_CENTER].GetNPArray()
		print("User defined center within coordinate system of the fragment:", Center)
		
		#TrnsOrig = c + b*np.dot(Orig, T)
		#Scaling is not used
		TrnsOrig = c + np.dot(Orig, T)
		self.cfg.Fragments[FragID][CFG_TRNS_ORIG] = OriginInterface(TrnsOrig[0], TrnsOrig[1], TrnsOrig[2])
		#Save center of the fragment (not origin)
		TrnsCenter = c + np.dot(Center, T)
		self.cfg.Fragments[FragID][CFG_TRNS_CENTER] = OriginInterface(TrnsCenter[0], TrnsCenter[1], TrnsCenter[2])

		self.ps.PrintKV("Transformed origin", TrnsOrig, LineKeyLen=19)
		self.ps.PrintKV("Norm", np.linalg.norm(TrnsOrig), LineKeyLen=19)
		if (CFG_PRINT_TDM_VIS):
			self.ps.PrintDiv()
			print("Transition Dipole Moments Visualization: State ID, TDM(X,Y,Z)  TransTDM(X,Y,Z)")
		#Transform and save transition dipole moments globally
		for ExcitedStateID, ExcitedState in Fragment[CFG_EXS_EXSTATE].items():
			Pnt = np.array((ExcitedState.x, ExcitedState.y, ExcitedState.z),dtype=np.dtype('d'))
			print(Pnt)
			TrnsPnt = c + np.dot(Pnt, T) - TrnsOrig
			TrnsPntVis = TrnsPnt+TrnsOrig
			if (CFG_PRINT_TDM_VIS):
				print(ExcitedStateID,Pnt[0]*TDM_SCALE_VIS, Pnt[1]*TDM_SCALE_VIS, Pnt[2]*TDM_SCALE_VIS, end=' ')
				print(TrnsPntVis[0]*TDM_SCALE_VIS, TrnsPntVis[1]*TDM_SCALE_VIS, TrnsPntVis[2]*TDM_SCALE_VIS)
				print("Norm:",np.linalg.norm(Pnt),np.linalg.norm(TrnsPnt))
			#with scaling: TrnsPnt = c + b*np.dot(Pnt, T) - TrnsOrig
			ExcitedState.x = TrnsPnt[0]
			ExcitedState.y = TrnsPnt[1]
			ExcitedState.z = TrnsPnt[2]
			self.cfg.Fragments[FragID][CFG_TRNS_EXS_EXSTATE][ExcitedStateID] = copy.deepcopy(ExcitedState)
			#self.cfg.Fragments[FragID][CFG_TRNS_EXS_EXSTATE][ExcitedStateID] = ExciteStateInterface(ExcitedState.ExStID, ExcitedState.Abs_cm1, ExcitedState.Abs_FWHM_cm1, TrnsPnt[0], TrnsPnt[1], TrnsPnt[2], ExcitedState.Ems_cm1, ExcitedState.Ems_FWHM_cm1, ExcitedState.El_Deph_Rate_ps1, ExcitedState.Alpha_M1cm1, ExcitedState.Phi_D, ExcitedState.FlLifetime_s)
			#ExcitedState.MyTDMPrint()
			#print
			#self.cfg.Fragments[FragID][CFG_TRNS_EXS_EXSTATE][ExcitedStateID].MyTDMPrint()
			#print

		if (self.cfg.Methods[CFG_MET_VERBOSITY] > CFG_MET_VERB_TRJ):
			self.ps.PrintDiv()
			print("Mol. Sys. (target) coordinates all atoms:")
			for AtomID, Atom in Fragment[CFG_FRG_ATOMS].items():
				Atom.MyPrint(AtomID+1)
			self.ps.PrintDiv()
			print("Transformed coordinates (independent index) compare against target...")
			for AtomID, Atom in Fragment[CFG_TRNS_EXS_ATOMS].items():
				Atom.MyPrint(AtomID+1)

		self.ps.PrintDiv()
		print("TRANSFOMRATION")
		print("{}: {} (not used info. only)".format("Scaling".ljust(STR_LEN_FLOAT),b))
		print("{}: {}".format("Translation".ljust(STR_LEN_FLOAT),c))
		print("{}: {}".format("Norm".ljust(STR_LEN_FLOAT),np.linalg.norm(c)))
		print("{}:".format("Rotation".ljust(STR_LEN_FLOAT)))
		print(T)
		#(w, x, y, z) = GetQuaternion(T)
		ax_ang = Vrrotmat2Vec(T)
		sqax_ang =  np.squeeze(ax_ang)
		ang_deg = sqax_ang[3]*180.0/math.pi
		print("{}: ({} {} {}) {} deg".format("Axis-Angl".ljust(STR_LEN_FLOAT),sqax_ang[0],sqax_ang[1], sqax_ang[2],ang_deg))
		#print sqax_ang 		
		#print "{}: {} {} {} {}".format("Quaternion".ljust(STR_LEN_FLOAT),w,x,y,z)
		print()

		self.ps.PrintSec("METRICS REPORT")

		STR_LEN_LBL = 34
		NumLbl = 14
		print("{}  {} {} {}".format("Transformation errors, Angstrom".ljust(STR_LEN_LBL),"SUM".ljust(NumLbl),"RMS".ljust(NumLbl),"N".ljust(NumLbl)))
		#(s, d, n) = self.GeometryDiff(Fragment[CFG_FRG_ATOMS], Fragment[CFG_EXS_ATOMS])
		#print "{}: {} {} {}".format("All atoms PDB vs. initial EXS".ljust(STR_LEN_LBL),str(round(s,INT_ROUND)).ljust(STR_LEN_FLOAT),str(round(d,INT_ROUND)).ljust(STR_LEN_FLOAT),str(n).ljust(STR_LEN_FLOAT))
		#(s, d, n) = self.GeometryDiff(Fragment[CFG_FRG_ATOMS], Fragment[CFG_EXS_ATOMS], Fragment[CFG_FRG_ATL])
		#print "{}: {} {} {}".format("Only matched atoms PDB vs. initial EXS".ljust(STR_LEN_LBL),str(round(s,INT_ROUND)).ljust(STR_LEN_FLOAT),str(round(d,INT_ROUND)).ljust(STR_LEN_FLOAT),str(n).ljust(STR_LEN_FLOAT))
		(s, d, n) = self.GeometryDiff(Fragment[CFG_FRG_ATOMS], Fragment[CFG_TRNS_EXS_ATOMS])
		print("{}: {} {} {}".format("All atoms:".ljust(STR_LEN_LBL),str(s).ljust(NumLbl),str(d).ljust(NumLbl), str(n).ljust(NumLbl)))
		(s, d, n) = self.GeometryDiff(Fragment[CFG_FRG_ATOMS], Fragment[CFG_TRNS_EXS_ATOMS], Fragment[CFG_FRG_ATL])
		print("{}: {} {} {}".format("Only matched atoms:".ljust(STR_LEN_LBL),str(s).ljust(NumLbl),str(d).ljust(NumLbl),str(n).ljust(NumLbl)))

		#pp = pprint.PrettyPrinter(indent=4)
		#pp.pprint(npFrgAtoms)
		#Transform transition dipole moment
		#Transform origin
		return

	def Run(self):
		"""Transform fragments"""

		#Scan fragments
		for FragID, Fragment in self.cfg.Fragments.items():
			self.ps.PrintSec(CFG_SEC_FRG+str(FragID+1)+" "+CFG_FRG_NAM+" = "+Fragment[CFG_FRG_NAM]+" "+CFG_FRG_ID+" = "+str(Fragment[CFG_FRG_ID]))
			self.TransformFrag(FragID, Fragment)
		
		if (self.cfg.Methods[CFG_MET_VERBOSITY] > CFG_MET_VERB_TRJ):
			self.ps.PrintDiv()
			print("All fragments PDB coordinates:")
			for FragID, Fragment in self.cfg.Fragments.items():
				for AtomID, Atom in Fragment[CFG_FRG_ATOMS].items():
					Atom.MyPrint()
	
			self.ps.PrintDiv()
			print("All fragments transformed coordinates (independent index):")
			for FragID, Fragment in self.cfg.Fragments.items():
				for AtomID, Atom in Fragment[CFG_TRNS_EXS_ATOMS].items():
					Atom.MyPrint(AtomID+1)
			self.ps.PrintDiv()
		return
#-------------------------------------------------------------------------------
