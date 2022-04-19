rem Sample BODIPY-TRZ calculations: make molecular fragments
set "PYTHON_CMD=python -u"
set "PYFREC_CMD=.\..\..\..\src\pyfrec.py -f"
conda activate pyfrec
%PYTHON_CMD% %PYFREC_CMD% .\bodipy_trz_a.ini 1> .\bodipy_trz_a.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\bodipy_trz_b.ini 1> .\bodipy_trz_b.log 2>&1
%PYTHON_CMD% %PYFREC_CMD% .\bodipy_trz_c.ini 1> .\bodipy_trz_c.log 2>&1
conda deactivate
