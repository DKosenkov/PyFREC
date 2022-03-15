rem
set "PYTHON_CMD=C:\Anaconda3\python -u"
set "PYFREC_CMD=.\..\..\..\PyFREC\src\v0207\pyfrec.py -f"
rem 
%PYTHON_CMD% %PYFREC_CMD% .\bodipy_dte_open.ini 1> .\bodipy_dte_open.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\bodipy_dte_closed.ini 1> .\bodipy_dte_closed.log 2>&1
