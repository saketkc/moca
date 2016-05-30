==========================================
MoCA: Tool for MOtif Conservation Analysis
==========================================

.. image:: https://img.shields.io/pypi/v/moca.svg
        :target: https://testpypi.python.org/pypi/moca/0.1.0

.. image:: https://img.shields.io/travis/saketkc/moca.svg
        :target: https://travis-ci.org/saketkc/moca

.. image:: https://coveralls.io/repos/github/saketkc/moca/badge.svg?branch=master
        :target: https://coveralls.io/github/saketkc/moca?branch=master

.. image:: https://landscape.io/github/saketkc/moca/master/landscape.svg?style=flat
        :target: https://landscape.io/github/saketkc/moca/master

.. image:: https://requires.io/github/saketkc/moca/requirements.svg?branch=master
        :target: https://requires.io/github/saketkc/moca/requirements/?branch=master


Python rewrite of `MoCA0.1.0`_

LICENSE
-------
ISC


API Documentation
-----------------

http://saketkc.github.io/moca/


Installation
------------
``moca`` is most compatible with the `conda`_ environment.

::

    $ git clone https://github.com:saketkc/moca.git
    $ cd moca
    $ conda create env -f environment.yml python=2.7
    $ source activate mocatest
    $ pip install .


A sample configuration file is available: `tests/data/application.cfg`

Workflow
--------

MoCA makes use of PhyloP/PhastCons/GERP scores to assess the quality of a
motif, the hypothesis being a 'true motif' would evolve slower as compared
to its surrounding(flanking sequences).

.. image:: https://raw.githubusercontent.com/saketkc/moca_web/master/docs/abstract/workflow.png

Usage
-----

::

    $ mocacli --help
    Usage: mocacli [OPTIONS]

    Run moca

    Options:
      -i, --bedfile TEXT            Bed file input  [required]
      -o, --oc TEXT                 Output Directory
      -c, --configuration TEXT      Configuration file  [required]
      --flank-seq INTEGER           Flanking sequence length  [required]
      --flank-motif INTEGER         Length of sequence flanking motif  [required]
      -g, -gb, --genome-build TEXT  Key denoting genome build to use in
                                configuration file  [required]
      --help                        Show this message and exit.


Example
-------

::

    $ mocacli -i tests/data/ENCFF002CDP.ctcf.bed\
        -g hg19
        -c tests/data/application.cfg\
        -o output_dir

.. image:: http://www.saket-choudhary.me/moca/_static/img/ENCFF002CEL.png

Tests
-----
``moca`` is mostly extensively tested. See `code-coverage`_. 

Run tests locally

::

    $ ./runtests.sh
      

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _`MoCA0.1.0`: https://github.com/saketkc/moca_web
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _`conda`: http://conda.pydata.org/docs/using/using.html
.. _`code-coverage`: https://coveralls.io/github/saketkc/moca?branch=master
