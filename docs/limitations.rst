
.. _limitations:
.. currentmodule:: camoco

Strengths and Limitations of Camoco 
###################################
Camoco excels at building and integrating different genomic data types due to
its internal architecture and design. This *build once, use many times* 
philosophy allows you to effectively explore your datasets and to design
experiments that utilize and compare different input parameters and options.

Camoco also straddles the divide between being high performance and high utility.
Being written in Python, Camoco is very flexible in how it can be used. Critical
parts of the code leverage optimized numpy/scipy/pandas code and several crucial 
parts of Camoco are written in Cython to optimize performance. 

The CLI provides basic functionality, but using Camoco interactively through IPython
or by importing camoco as a module allows for more sophisticated and specialized
scripts to be written. Scripts are typically straight forward to write as Camoco
shares a python API similar to that of numpy/scipy/pandas.

Camoco was built and designed with a modern Python workflow in mind. These
Due to some design choices and dependencies, Camoco has some limitations that 
make it difficult to run on some systems / platforms.

*Older python releases (<3.6)*
  Several newer features released in newer versions of Python are heavily 
  used internally within Camoco. Requiring Python versions 3.6+ may conflict
  the system installed version of Python. Its recommended that Camoco is 
  installed using a python virtual environment. We personally have had success
  with miniconda which can be installed `here <https://conda.io/docs/user-guide/install/index.html>`_.
  See the installation documentation for more details.

*Network file systems*
  Camoco heavily utilizes SQLite for data persistence. As SQLite is file based
  and does not support concurrent writing. Some network based file systems that
  are commonly used within academic departments and super computer institutes 
  (e.g. NFS) do not play well with SQLite and can `corrupt the underlying
  database <https://www.sqlite.org/howtocorrupt.html#_filesystems_with_broken_or_missing_lock_implementations>`_.
  By default, Camoco stores its databases in `~/.camoco`. If you are running
  Camoco from a computer that utilizes a network file system, it is recommended
  that you create a symbolic link from `~/.camoco` to a local (i.e. non-network)
  directory that is stored somewhere else on the machine.

*Parallelizing Camoco*
  Being written in Python means that Camoco is limited by the GIL. Long story
  short, Camoco is inherently single-threaded. This issue is slowly being 
  addressed by converting some code to utilize the newer python async 
  functionality. Until that is fully supported, the easiest way to parallelize
  your Camoco workflow is by running many instances of Camoco from the command
  line. As Camoco leverages databases (which are not as limited by parallelization
  as python), running many instances of Camoco in parallel is feasible. For more
  a short tutorial showcasing this, please read this `blog post <http://blog.linkage.io/speeding-things-up-with-gnu-parallel/>`_.
 
