#-------------------------------------------------------------------------------
def ParseList(List):
	"""Parse lists (e.g. 0,1,3-8,9)"""

	Res = list()

	if ( str(List).isdigit() ):
		Res.append(int(List))
	else:
		for Element in List.split(','):
			if '-' in Element:
				x1,x2 = Element.split('-')
				Res.extend(list(range(int(x1), int(x2)+1)))
			else:
				Res.append(int(Element))
	return Res
#-------------------------------------------------------------------------------