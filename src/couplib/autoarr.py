from collections import OrderedDict
#-------------------------------------------------------------------------------
class AutoArr(OrderedDict):
	"""Impementation of automatic associative arrays (auto vivification)"""
	def __getitem__(self, element):
		try:
			return OrderedDict.__getitem__(self, element)
		except KeyError:
			self[element] = type(self)()
			return self[element]
#-------------------------------------------------------------------------------