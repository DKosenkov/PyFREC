#!/bin/bash 
# Run the script with: '/bin/bash -x' for debugging
#
cd .\..\src
#rem save results of tests to the "tests" directory
#%PYTHON_CMD% %PYTEST_MY_CMD% 1> .\..\tests\tests.log 2>&1
#conda run -n pyfrec pytest tests.py --junitxml result_Tests.xml
conda run -n pyfrec pytest tests.py 1> .\..\tests\tests.log 2>&1
#call .\..\bin\pyfrec-dea.bat
#rem conda install -c anaconda pytest
cd .\tests


