rem conda install -c anaconda pytest
set "PYTHON_CMD=python -u"
set "PYFREC_CMD=.\..\src\pyfrec.py -f"
rem location of the tested codes (NOT data files ...)
cd ..
rem save results of tests to the "tests" directory
conda run -n pyfrec pytest test_pyfrec.py 1> .\tests\tests.log 2>&1
cd tests


