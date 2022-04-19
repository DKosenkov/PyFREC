#!/bin/bash
# Run the script with: '/bin/bash -x' for debugging
#
#Sample BODIPY-DTE calculations: electronic coupings and Forster energy transfer modeling
#
PYTHON_CMD='python -u'
PYFREC_CMD='./../../../src/pyfrec.py -f'
#
# Sample BODIPY-DTE calculations: electronic coupings and Forster energy transfer modeling
$PYTHON_CMD $PYFREC_CMD ./bodipy_dte_closed.ini 1> ./bodipy_dte_closed.log 2>&1
$PYTHON_CMD $PYFREC_CMD ./bodipy_dte_open.ini 1> ./bodipy_dte_open.log 2>&1
