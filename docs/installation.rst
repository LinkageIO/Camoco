
.. _installation:


Installation
############

Camoco can be installed using the `Python Package Index (pypi) <https://pypi.org/>`__. This
process is typical for python packages and uses a shell tool called :code:`pip` to manage 
software dependencies. Python requires that :code:`numpy` already be installed.::

  $ pip install numpy
  $ pip install camoco

Camoco does have a Python version requirement as well as a few dependencies in order to 
install successfully.

.. note:: 
  Code blocks starting with $ designate shell commands.


Operating System Compatibilty
=============================
Camoco was built on and tested using Linux (Ubuntu 18.04). It is also suggested that you
have a machine with at least 8Gb of usable RAM.

Python Version Requirement
==========================
Camoco requires python version 3.6+. If this differs from your system version of python,
it is recommended that you install Camoco within a python virtual environment. We have 
had success using miniconda, which can be installed `here <https://conda.io/docs/user-guide/install/index.html>`__.

Once miniconda is installed, create a virtual environment: ::

  $ conda create -n camoco_env python=3.6

Follow the prompts on the screen. Next, activate the python virtual environment: ::
  
  $ source activate camoco_env

Python should now be 3.6:

.. command-output:: python --version

Install Camoco as described above: ::
  
  $ pip install numpy
  $ pip install camoco



