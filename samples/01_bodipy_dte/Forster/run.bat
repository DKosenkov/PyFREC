rem Sample BODIPY-DTE calculations: electronic coupings and Forster energy transfer modeling
set "PYTHON_CMD=python -u"
set "PYFREC_CMD=.\..\..\..\src\pyfrec.py -f"
%PYTHON_CMD% %PYFREC_CMD% .\bodipy_dte_closed.ini 1> .\bodipy_dte_closed.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\bodipy_dte_open.ini 1> .\bodipy_dte_open.log 2>&1



