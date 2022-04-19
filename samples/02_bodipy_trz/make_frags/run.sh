#!/bin/bash
# Run the script with: '/bin/bash -x' for debugging
#
#Sample BODIPY-TRZ calculations: make molecular fragments
#
PYTHON_CMD='python -u'
PYFREC_CMD='./../../../src/pyfrec.py -f'
#
# Sample BODIPY-DTE calculations: electronic coupings and Forster energy transfer modeling
$PYTHON_CMD $PYFREC_CMD ./bodipy_trz_a.ini 1> ./bodipy_trz_a.log 2>&1
$PYTHON_CMD $PYFREC_CMD ./bodipy_trz_b.ini 1> ./bodipy_trz_b.log 2>&1
$PYTHON_CMD $PYFREC_CMD ./bodipy_trz_c.ini 1> ./bodipy_trz_c.log 2>&1