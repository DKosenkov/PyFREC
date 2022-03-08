rem BODIPY-TRZ 
set "PYTHON_CMD=C:\Anaconda3\python -u"
set "PYFREC_CMD=.\..\..\..\PyFREC\src\v0207\pyfrec.py -f"

%PYTHON_CMD% %PYFREC_CMD% .\BODIPY_a_ac.ini 1> .\BODIPY_a_ac.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\BODIPY_b_ac.ini 1> .\BODIPY_b_ac.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\BODIPY_c_ac.ini 1> .\BODIPY_c_ac.log 2>&1


