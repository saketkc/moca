===============================
MoCA
===============================

.. image:: https://img.shields.io/pypi/v/moca.svg
        :target: https://testpypi.python.org/pypi

.. image:: https://img.shields.io/travis/saketkc/moca.svg
        :target: https://travis-ci.org/saketkc/moca


.. image:: https://codecov.io/github/saketkc/moca/coverage.svg?branch=master
        :target: https://codecov.io/github/saketkc/moca?branch=master

.. image:: https://coveralls.io/repos/github/saketkc/moca/badge.svg?branch=master
        :target: https://coveralls.io/github/saketkc/moca?branch=master


Tool for motif conservation analysis
Python rewrite of `MoCA0.1.0`

* Free software: BSD license

Documentation
-------------

http://saketkc.github.io/moca/


Installation
------------
``moca`` is most compatible with the `conda`_ environment.

::

    $ git clone https://github.com:saketkc/moca.git
    $ cd moca
    $ conda create env -f environment.yml python=2.7
    $ source activate mocatest
    $ mocacli --bedfile tests/data/ENCFF002CDP.ctcf.bed --phylop /media/data1/genomes/hg19/phylop/hg19.100way.phyloP100way.bw --gerp /media/data1/genomes/hg19/gerp/All_hg19_RS.bw -gt /media/data1/genomes/hg19/fasta/hg19.sizes -gf /media/data1/genomes/hg19/fasta/hg19.fa --configuration tests/data/application.cfg


Tests
-----
``moca`` is mostly extensively tested. See `code-coverage`_. 

Run tests locally

::

    $ nosetests -v
     nosetests -v
     Test load broadPeak ... ok
     Test generate fasta ... ok
     Test load macsPeak ... ok
     Test load narrowPeak ... ok
     test_scorefile (tests.test_bedoperations.TestBedoperations) ... ok
     Test configuration genomes ... ok
     Test configuration sections ... ok
     bits            2.3                    * ... ok
     Test fimo runner ... ok
     Test fimo_to_sites ... ok
     Test meme runner ... ok
     Test load wig ... ok
     Test if query is out of bounds ... ok
     Test wig query ... ok

     ----------------------------------------------------------------------
     Ran 14 tests in 2.506s




Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _`MoCA0.1.0`: https://github.com/saketkc/moca_web
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _`conda`: http://conda.pydata.org/docs/using/using.html
.. _`code-coverage`: https://coveralls.io/github/saketkc/moca?branch=master
