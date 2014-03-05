.. cyclopy documentation master file, created by
   sphinx-quickstart on Wed Mar  5 13:06:15 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to cyclopy's documentation!
===================================
.. automodule:: cyclopy
    :members:


cyclopy is a python module thought to be used from ipython_

.. _ipython: http://ipython.org/


cyclopy can be downlaoded and installed in your linux system following these steps:

.. sourcecode:: bash

    cd place/where/you/keep/code
    git clone git@github.com:luca-penasa/cyclopy.git

    cd cyclopy
    python setup.py instal [-u]

    # use -u flag if you whant to istall cyclopy for your user only


Then you can try it from within an ipython session (code linked), create a random series and plot its 4-pi multitaper sectrum.


.. plot::

    import numpy as np

    from matplotlib.pyplot import subplot, figure, title, subplots_adjust, ylabel

    import cyclopy
    from cyclopy import Elements


    series = cyclopy.Elements.TimeSeriesEven(np.random.randn(100))
    figure()
    subplot(121)
    series.plot()
    subplot(122)
    title("Spectra")
    ylabel("Spectral Power")
    series.plotMTMSpectrum(pi=4)
    subplots_adjust()


Refer to the tutorial to see what cyclopy can do: :ref:`tutorial`


.. _tutorial:

Tutorial
========
    .. toctree::
        tutorial.rst

Api doc
=======
    .. toctree::
        :maxdepth: 1

        apidoc/modules.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`