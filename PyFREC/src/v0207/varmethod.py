from configuration import *
from couplib.myreportservice import *
import numpy as np
from couplib.constants import *

#-------------------------------------------------------------------------------
class VariationMethod(object):
	"""Calculation manager that runs variational calculations"""

	def __init__(self, Hamiltonian):
		"""Initialization"""
		self.ps = MyReportService()
		self.Hamiltonian = Hamiltonian
		return

	def Calc(self):
		"""
		Variational method calculations:
		eigenvalue problem HC = EC
		E - daigonal matrix of energies
		C - eigenvectors (orbitals)
		H - Hamiltonian in a matrix form
		"""

		#print "Hamiltonian:"
		#pp = pprint.PrettyPrinter(indent=4)
		#pp.pprint(self.Hamiltonian)

		Energies, Orbitals = np.linalg.eig(self.Hamiltonian);

		#Sort eigenvalues and vectors ascending
		EnergiesPerm = np.argsort(Energies)
		Energies = Energies[EnergiesPerm]
		Orbitals = Orbitals[:,EnergiesPerm]

		self.ps.PrintDiv()
		print("Energies of coupled states in the basis of molecular fragment excitations:")
		print(("{} {} {}".format("i".ljust(STR_LEN_INT),"E_coup., cm-1".ljust(STR_LEN_FLOAT),"Eigen Vectors")))
		i = 0
		for E in Energies:
			print(("{} {} {}".format(str(i+1).ljust(STR_LEN_INT),str(round(E*HartreeToCM1,INT_ROUND)).ljust(STR_LEN_FLOAT),np.array_str(Orbitals[:,i],max_line_width=MAX_LINE))))
			i += 1

		self.ps.PrintDiv()
		print("Super-molecular excitations in the basis of molecular fragment exciations in columns:")
		print((np.array_str(Orbitals,max_line_width=MAX_LINE)))
		self.ps.PrintDiv(DivSym = "=")

		return (Energies, Orbitals)
#-------------------------------------------------------------------------------