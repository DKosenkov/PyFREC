#!/bin/bash 
#conda install -c anaconda pytest
#location of the tested codes (NOT data files ...)
cd ..
# save results of tests to the "tests" directory
conda run -n pyfrec pytest test_pyfrec.py 1> ./tests/tests.log 2>&1
cd tests

#!/bin/bash
# Run the script with: '/bin/bash -x' for debugging
#
cd ..
#rem save results of tests to the "tests" directory
conda run -n pyfrec pytest test_pyfrec.py 1> ./tests/tests.log 2>&1
cd ./tests

