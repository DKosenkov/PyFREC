rem BODIPY-TRZ
set "PYTHON_CMD=C:\Anaconda3\python -u"
set "PYFREC_CMD=.\..\..\..\PyFREC\src\v0207\pyfrec.py -f"

%PYTHON_CMD% %PYFREC_CMD% .\bodipy_trz_a.ini 1> .\bodipy_trz_a.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\bodipy_trz_b.ini 1> .\bodipy_trz_b.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\bodipy_trz_c.ini 1> .\bodipy_trz_c.log 2>&1
