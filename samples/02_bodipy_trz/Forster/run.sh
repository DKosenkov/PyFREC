#!/bin/bash
# Run the script with: '/bin/bash -x' for debugging
#
#Sample BODIPY-TRZ calculations: electronic coupings and Forster energy transfer modeling
#
PYTHON_CMD='python -u'
PYFREC_CMD='./../../../src/pyfrec.py -f'
#
# Sample BODIPY-DTE calculations: electronic coupings and Forster energy transfer modeling
$PYTHON_CMD $PYFREC_CMD ./bodipy_a_ac.ini 1> ./bodipy_a_ac.log 2>&1
$PYTHON_CMD $PYFREC_CMD ./bodipy_b_ac.ini 1> ./bodipy_b_ac.log 2>&1
$PYTHON_CMD $PYFREC_CMD ./bodipy_c_ac.ini 1> ./bodipy_c_ac.log 2>&1
