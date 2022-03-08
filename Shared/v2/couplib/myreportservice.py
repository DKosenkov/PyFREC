import time
import datetime
import sys

#Number of symbols in a terminal window
REP_LINE_LEN = 80
#Max line for priojnting long arrays
MAX_LINE = 500
#Mximal length of keys  to display
REP_LINE_KEY_SIZE = 24
#Default symbol to break sections
REP_DEF_DIV = "-"
#Default separator for key-value pairs
REP_DEF_SEP = ":"
#Length of line to represent float number
STR_LEN_FLOAT = 12
STR_LEN_FLOAT_LARGE = 16
#Length of line to represent integer number
STR_LEN_INT = 12
#Length of the line to repretent atomic label (element symblol + indexes)
STR_LEN_ALABEL  = 8
#Number of digits after the decimal point
INT_ROUND = 4
#Space for Boolean
STR_LEN_BOOL = 3
#-------------------------------------------------------------------------------
class MyReportService(object):
	"""Service class for nice terminal output"""

	@staticmethod
	def Header(ProgramHeader = ""):
		"""Print title of the program"""
		print(ProgramHeader)
		return

	def __init__(self, LineLen=REP_LINE_LEN, LineKeyLen=REP_LINE_KEY_SIZE):
		"""Setup line of the terminal line"""
		self.LineLen = LineLen
		self.LineKeyLen = LineKeyLen
		return

	def DateTimeStamp(self):
		"""Print date-time stamp"""
		ts = time.time()
		st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
		print(st)
		return

	def SysVersion(self):
		"""Print Python version"""
		print (sys.version)
		return


	def PrintSec(self, SectionTitle, LineLen = REP_LINE_LEN, DivSym = REP_DEF_DIV):
		"""Print titles of sections"""
		print(DivSym*LineLen)
		print(SectionTitle.upper())
		return

	def PrintKV(self, Key, Value, LineKeyLen = REP_LINE_KEY_SIZE ,Separator = REP_DEF_SEP):
		"""Print key-value pairs"""
		print("{}{} {}".format( Key.ljust(LineKeyLen),Separator,Value))
		return

	def PrintDiv(self, LineLen = REP_LINE_LEN, DivSym = REP_DEF_DIV):
		"""Print line divider"""
		print(DivSym*LineLen)
		return
#-------------------------------------------------------------------------------