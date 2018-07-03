#!/usr/bin/env bash

cd ~/source/pyatactools-dev
python setup.py install --prefix=~/local/
mv $HOME/local/lib/python2.7/site-packages/pyatactools-0.0.1-py2.7.egg $HOME/local/lib/python2.7/site-packages/pyatactools-0.0.1-py2.7.egg.zip
unzip $HOME/local/lib/python2.7/site-packages/pyatactools-0.0.1-py2.7.egg.zip
touch $HOME/local/lib/python2.7/site-packages/pyatactools/__init__.py
mv $HOME/local/lib/python2.7/site-packages/pyatactools-0.0.1-py2.7.egg.zip $HOME/local/lib/python2-7/site-packages/pyatactools-0.0.1-py2.7.egg
#export PYTHONPATH=$HOME/local/lib/python2.7/site-packages/pyatactools-0.0.1-py2.7.egg:$PYTHONPATH
