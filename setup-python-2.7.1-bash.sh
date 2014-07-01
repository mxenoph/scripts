#!/bin/bash

# Gets the name of the operating system i.e Linux_x86_64
export OS=`uname`_`uname -m`
export STOW_HOME=/ebi/research/software/${OS}/opt/stow

export PYTHON_ROOT=$STOW_HOME/python-2.7.1-shared
export PYTHON_VERSION=2.7.1
export PYTHON=$PYTHON_ROOT/bin/python

#export PYTHONPATH_DEF=$PYTHONPATH
#export LD_LIBRARY_PATH_DEF=$LD_LIBRARY_PATH

export PYTHONPATH=/homes/kostadim/myrto/lib/python2.7/site-packages

#unset PYTHONPATH

export PATH=$PYTHON_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PYTHON_ROOT/lib:$STOW_HOME/gcc-4.5.2/lib64:$STOW_HOME/lib:$STOW_HOME/centos5-libraries/lib

#export PYTHONPATH=${PYTHON_ROOT}/lib/python2.7/site-packages:$STOW_HOME/libsbml-4.3.0-xerces/lib/python2.6/site-packages:$STOW_HOME/python-2.6.1/lib/python2.6/site-packages:$STOW_HOME/ecell-3.1.106/lib/python2.5/site-packages:${PYTHON_ROOT}/lib/python2.7/site-packages/copasi_python
#export PYTHONPATH=${PYTHON_ROOT}/lib/python2.7/site-packages:$STOW_HOME/libsbml-4.3.0-xerces/lib/python2.6/site-packages:$STOW_HOME/ecell-3.1.106/lib/python2.5/site-packages:${PYTHON_ROOT}/lib/python2.7/site-packages/copasi_python
# export PYTHONPATH=$PYTHON_ROOT/lib/python2.5/site-packages/gtk-2.0:$PYTHON_ROOT/lib/python2.5/site-packages/cairo
