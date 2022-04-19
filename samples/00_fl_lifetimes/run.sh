#!/bin/bash 
# Run the script with: '/bin/bash -x' for debugging
# Sample fluorescence lifetime calculations
#
PYTHON_CMD='python -u'
PYFREC_CMD='./../../src/pyfrec.py -f'
#
$PYTHON_CMD $PYFREC_CMD ./bodipy_open.ini 1> ./bodipy_open.log 2>&1
$PYTHON_CMD $PYFREC_CMD ./fluorescein.ini 1> ./fluorescein.log 2>&1
$PYTHON_CMD $PYFREC_CMD ./rhodamineb.ini 1> ./rhodamineb.log 2>&1
$PYTHON_CMD $PYFREC_CMD ./popop.ini 1> ./popop.log 2>&1
$PYTHON_CMD $PYFREC_CMD ./pterphenyl.ini 1> ./pterphenyl.log 2>&1
$PYTHON_CMD $PYFREC_CMD ./r6g.ini 1> ./r6g.log 2>&1
