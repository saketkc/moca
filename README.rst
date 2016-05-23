===============================
MoCA
===============================

.. image:: https://img.shields.io/pypi/v/moca.svg
        :target: https://testpypi.python.org/pypi/moca/0.1.0

.. image:: https://img.shields.io/travis/saketkc/moca.svg
        :target: https://travis-ci.org/saketkc/moca

.. image:: https://codecov.io/github/saketkc/moca/coverage.svg?branch=master
        :target: https://codecov.io/github/saketkc/moca?branch=master

.. image:: https://coveralls.io/repos/github/saketkc/moca/badge.svg?branch=master
        :target: https://coveralls.io/github/saketkc/moca?branch=master

.. image:: https://landscape.io/github/saketkc/moca/master/landscape.svg?style=flat
        :target: https://landscape.io/github/saketkc/moca/master

.. image:: https://requires.io/github/saketkc/moca/requirements.svg?branch=master
        :target: https://requires.io/github/saketkc/moca/requirements/?branch=master

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
    $ export MEME_ETC_DIR="/home/user/anaconda2/envs/mocatest/etc"
    $ export MEME_BIN_DIR="/home/user/anaconda2/envs/mocatest/bin"
    $ pip install .
    $ mocacli --bedfile tests/data/ENCFF002CDP.ctcf.bed --phylop /media/data1/genomes/hg19/phylop/hg19.100way.phyloP100way.bw --gerp /media/data1/genomes/hg19/gerp/All_hg19_RS.bw -gt /media/data1/genomes/hg19/fasta/hg19.sizes -gf /media/data1/genomes/hg19/fasta/hg19.fa --configuration tests/data/application.cfg


Usage
-----

::

    $ mocacli [OPTIONS]

      Run moca

      Options:
        -i, --bedfile TEXT            Bed file input  [required]
        -o, --oc TEXT                 Output Directory
        -c, --configuration TEXT      Configuration file  [required]
        --flank-seq INTEGER           Flanking sequence length  [required]
        --flank-motif INTEGER         Length of sequence flanking motif  [required]
        -g, -gb, --genome-build TEXT  Key denoting genome build to use in
                                      configuration file  [required]
Example
-------

::

    mocacli -i tests/data/ENCFF002CDP.ctcf.bed -c tests/data/application.cfg -g hg19

.. image:: http://www.saket-choudhary.me/moca/_static/img/ENCFF002CEL.png



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
