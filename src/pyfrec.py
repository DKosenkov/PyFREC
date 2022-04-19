PROGRAM_HEADER = """
################################################################################
#                                                                              #
# Python FRagment Electronic Coupling (PyFREC) v2.7                            #
# Program for calculation of electronic couplings                              #
#                                                                              #
# Dmitri Kosenkov,  2016-2022                                                  #
# Monmouth University                                                          #
#                                                                              #
#########################################################################DVK####
"""

from configuration import Configuration
from couplib.myreportservice import MyReportService
from calcmanager import CalcManager
#-------------------------------------------------------------------------------
def main():
	#Report builder

	#Print program header
	MyReportService().Header(PROGRAM_HEADER)

	#Load configuration file (data files, computational models, approximations, basis sets etc.) and call CalculationManager
	CalcManager(Configuration()).Run()

	return 0
#-------------------------------------------------------------------------------
#Call the main program only if the script is called directly. Not through the 'import' statement
if __name__ == "__main__":
	main()
