rem Test 
rem
set "PYTHON_CMD=C:\Anaconda3\python -u"
set "PYFREC_CMD=.\..\..\PyFREC\src\v0207\pyfrec.py -f"
rem 
%PYTHON_CMD% %PYFREC_CMD% .\BODIPY_OPEN.ini 1> .\BODIPY_OPEN.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\Fluorescein.ini 1> .\Fluorescein.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\RhodamineB.ini 1> .\RhodamineB.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\POPOP.ini 1> .\POPOP.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\pTerphenyl.ini 1> .\pTerphenyl.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\R6G.ini 1> .\R6G.log 2>&1


