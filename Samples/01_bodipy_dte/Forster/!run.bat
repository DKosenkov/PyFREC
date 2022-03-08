rem
set "PYTHON_CMD=C:\Anaconda3\python -u"
set "PYFREC_CMD=.\..\..\..\PyFREC\src\v0207\pyfrec.py -f"
rem 
%PYTHON_CMD% %PYFREC_CMD% .\BODIPY_DTE_closed.ini 1> .\BODIPY_DTE_closed.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\BODIPY_DTE_open.ini 1> .\BODIPY_DTE_open.log 2>&1



