"""
################################################################################
# Python FRagment Electronic Coupling (PyFREC) v1.4                            #
# Program for calculation of electronic couplings                              #
################################################################################
"""
from couplib.myfileservice import MyFileService
from couplib.myreportservice import MyReportService
from couplib.parselist import ParseList

from configuration import *
from interfaces import *

CIF_FORMAT = ".CIF"

#-------------------------------------------------------------------------------
class PDBReader(MyFileService):
	"""Class reads a complete molecular system from the PDB file according to the user's input in the configuration file"""

	def __init__(self, cfg):
		"""Initialize with current configuration"""
		self.cfg = cfg
		return

	def ParseAtom(self,Atom):
		"""Parse atom records forming fragments and eliminating unused atoms"""
		#Scan fragments to match ResID
		for FragID, Fragment in list(self.cfg.Fragments.items()):
			#check chain IDs
			FragChainID = str(Fragment[CFG_FRG_CHAIN_ID]).strip().upper()
			AtomChainID = str(Atom.ChainID).strip().upper()
			#Found residue with matching ID
			if ( (int(Fragment[CFG_FRG_ID]) == Atom.ResID) and ( (FragChainID == AtomChainID) or (FragChainID == '') ) ):
				#Check if residue names match
				if str(Fragment[CFG_FRG_NAM]).strip().upper() == Atom.ResName:
					#Check if atom IDs match or all atoms from the residue are included (CFG_FRG_AL == 0) in the configuration
					CfgAtomIDList = ParseList(Fragment[CFG_FRG_AL])
					#Alternate location indicator
					CfgAtomAltLocList = Fragment[CFG_FRG_ALTLOC].split(LIST_SEPARATOR)
					for CfgAtomID in CfgAtomIDList:
						if ( (CfgAtomID == 0) or (CfgAtomID==Atom.AtomPDBID) ):
							AltLoc = str(Atom.AltLoc).strip().upper()
							IsAtomAcceptable = True
							if ( AltLoc != ''):
								#Check Alternate location indicators
								IsAtomAcceptable = False
								for CfgAtomAltLoc in CfgAtomAltLocList:
									CfgAltLoc = str(CfgAtomAltLoc).strip().upper()
									if (CfgAltLoc == AltLoc):
										IsAtomAcceptable = True
								if ( not IsAtomAcceptable ):
									print(("WARNING: In "+CFG_SEC_FRG+str(FragID+1)+" atom "+str(Atom.ElSym)+str(Atom.AtomPDBID)+" in "+str(Atom.ResName).strip()+" "+str(Atom.ResID)+" with alternate location id: '"+str(AltLoc)+"' will be ignored. (use keyword '"+str(CFG_FRG_ALTLOC)+"' to add the alternate location)"))
							#Add atom to the fragment
							if ( IsAtomAcceptable ):
								self.cfg.Fragments[FragID][CFG_FRG_ATOMS][self.cfg.Fragments[FragID][CFG_FRG_NATOMS]] = Atom
								self.cfg.Fragments[FragID][CFG_FRG_NATOMS] += 1
								break
		return

	def ReadPDB(self,FileName):
		"""Read PDB input file"""
		#Open and read PDB file content
		'''
		01 - 06       Record name    "HETATM"
		07 - 11       Integer        serial        Atom serial number.
		13 - 16       Atom           name          Atom name.
		17            Character      altLoc        Alternate location indicator.
		18 - 20       Residue name   resName       Residue name.
		22            Character      chainID       Chain identifier.
		23 - 26       Integer        resSeq        Residue sequence number.
		27            AChar          iCode         Code for insertion of residues.
		31 - 38       Real(8.3)      x             Orthogonal coordinates for X.
		39 - 46       Real(8.3)      y             Orthogonal coordinates for Y.
		47 - 54       Real(8.3)      z             Orthogonal coordinates for Z.
		55 - 60       Real(6.2)      occupancy     Occupancy.
		61 - 66       Real(6.2)      tempFactor    Temperature factor.
		77 - 78       LString(2)     element       Element symbol; right-justified.
		79 - 80       LString(2)     charge        Charge on the atom.
		'''
		FILE = open(FileName,'r')
		for line in FILE.readlines():
			SecFrgRE = r'^\s*'+CFG_SEC_FRG.strip().upper()+r'\d+\s*$'
			if ( re.match(r'^ATOM\s*', line) or re.match(r'^HETATM\s*', line) ):
				#The following numerical values are taken from PDB format specification
				#http://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html
				AtomPDBID = int(line[6:11])  # Atom serial number
				AtomName  = str(line[12:16]).strip().upper() # Atom name
				AltLoc    = str(line[16:17]).strip().upper() # Alternate location indicator (the same atom at alternate locations)
				ResName   = str(line[17:20]).strip().upper() # Residue name
				ChainID   = str(line[21:22]).strip().upper() # Chain id name
				ResID     = int(line[22:26]) # Residue sequence number
				x         = float(line[30:38]) # Orthogonal coordinates for X in Angstroms
				y         = float(line[38:46]) # Orthogonal coordinates for Y in Angstroms
				z         = float(line[46:54]) # Orthogonal coordinates for Y in Angstroms
				ElSym     = str(line[76:78]).strip() # Atom name
				#print "{},{},{},{},{},{},{},{}".format(AtomPDBID,AtomName,AltLoc,ResName,ResID,x,y,z,ElSym)
				self.ParseAtom(AtomInterface(x, y, z, ElSym, AtomPDBID, AtomName, AltLoc, ResID, ResName, ChainID))
		FILE.close()
		return

	def ReadCIF(self,FileName):
		"""Read CIF input file"""
		FILE = open(FileName,'r')
		for line in FILE.readlines():
			SecFrgRE = r'^\s*'+CFG_SEC_FRG.strip().upper()+r'\d+\s*$'
			RegExline = r'^(ATOM|HETATM)\s+(\d+)\s+([a-zA-Z0-9_]+)\s+([a-zA-Z0-9_]+)\s+.\s+([a-zA-Z0-9_]+)\s+([a-zA-Z0-9_]+)\s+(\d+)'
			matched = re.match(RegExline, line)
			if(matched):
				AtomPDBID = int(matched.group(2)) # Atom ID
				ElSym = str(matched.group(3)) # Element
				AtomName = str(matched.group(4)) # Atom type
				ResName = str(matched.group(5)) # Residue name
				ChainID = str(matched.group(6)) # Chain name?
				#_ChainCode = int(matched.group(7)) # Chain code?
				RegExline2 = r'\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+('+RE_FLOAT+r')\s+.\s+(\d+)'
				matched2 = re.search(RegExline2, line)
				if(matched2):
					x = float(matched2.group(1)) # x, A
					y = float(matched2.group(2)) # y, A
					z = float(matched2.group(3)) # z, A
					#n = float(matched2.group(4)) # Occupation ?
					#b = float(matched2.group(5)) # Betha?
					ResID = int(matched2.group(6)) # Residue number
					AltLoc = ""
					#print(line)
					self.ParseAtom(AtomInterface(x, y, z, ElSym, AtomPDBID, AtomName, AltLoc, ResID, ResName, ChainID))
					 
		FILE.close()
		return

	def Read(self):
		"""Read PDB file"""
		#Open PDB file
		FileName = self.cfg.MolSys[CFG_MSY_GEO]
		if not self.IsValidFile(FileName,"Error: Cannot find the PDB file"):
			exit(-1)

		filename, file_extension = os.path.splitext(FileName)

		if ( file_extension.strip().upper() == CIF_FORMAT ): #CIF file
			print("Reading CIF file...")
			self.ReadCIF(FileName)
		else: #PDB file
			print("Reading PDB file...")
			self.ReadPDB(FileName)
		return

	def ReadingReport(self, ReadingReport):
		"""Reports summary of PDB reading (warnings etc.)"""
		return

	def PrintFragmens(self):
		"""Print fragments read from PDB file"""

		self.ps = MyReportService()
		#Scan fragments
		for FragID, Fragment in list(self.cfg.Fragments.items()):
			self.ps.PrintSec(CFG_FREC_SAYS+" "+CFG_SEC_FRG+str(FragID+1)+" "+CFG_FRG_NAM+" = "+Fragment[CFG_FRG_NAM]+" "+CFG_FRG_ID+" = "+str(Fragment[CFG_FRG_ID])+" "+CFG_FRG_AL+" = "+Fragment[CFG_FRG_AL]+" "+CFG_FRG_NATOMS+" = "+str(Fragment[CFG_FRG_NATOMS]) )
			Atoms = self.cfg.Fragments[FragID][CFG_FRG_ATOMS]
			for AtomID, Atom in list(Atoms.items()):
				Atom.MyPrint()

		self.ps.PrintDiv()
		return
#-------------------------------------------------------------------------------