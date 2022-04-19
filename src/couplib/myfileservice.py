import sys
import os.path
#-------------------------------------------------------------------------------
class MyFileService(object):
	"""Class containing service routines to work with files"""

	def __init__(self):
		return

	def IsValidFile(self,FileName="", ErrMsg=""):
		"""Check if the file exists and prints error message if the file doesn't"""
		if not os.path.exists(FileName):
			print("Error: ")
			print(ErrMsg)
			print("The file %s does not exist!"%FileName)
			return False
		return True

	def GetLines(self,FileName="", ErrMsg=""):
		"""Read content of the text file"""
		if not self.IsValidFile(FileName,"Cannot find file with molecular geometry"):
			exit(-1)

		FILE = open(FileName,'r')
		Lines = FILE .readlines()
		FILE.close()
		return Lines
#-------------------------------------------------------------------------------