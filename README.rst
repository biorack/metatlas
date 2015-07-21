Metabolite Atlas
================

Metabolomics is the comprehensive profiling of the small molecule composition of a biological sample. This approach is being used to provide new insights into a variety of biological systems (clinical, bioenergy, etc.). A grand challenge for metabolomics is the complexity of the data, which often include many experimental artifacts. This is compounded by the tremendous chemical diversity of metabolites where identification of each uncharacterized metabolite is in many ways its own puzzle. The Metabolite Atlas project will provide easy to use tools that enable a scientist to quickly capture knowledge about what compounds they have observed and to propagate that knowledge to future experiments. Instead of having to sift through billions of data points, a scientist will use pre-existing method specific metabolite atlas results to suggest candidate identifications as soon as the data is generated.


Features
--------
- Parse LCMS mzML files to Pytables (HDF)
- Query and plot LCMS data
    - Extracted-ion chromatogram (XIC) data
    - Spectrogram data
    - HeatMap of Retention Time vs (m/z)
    - Generic query for (rt, mz, i) data
- Custom plotters for XIC, Chromatogram, and Heatmap


Installation
------------

.. code-block:: bash

    $ pip install metatlas


Documentation
-------------

Documentation is available online_.

See example notebook_.

For version information, see the Revision History_.


.. _online: http://metabolite-atlas.github.io/metatlas/

.. _notebook: http://nbviewer.ipython.org/github/metabolite-atlas/metatlas/blob/master/docs/example_notebooks/data_access_examples.ipynb

.. _History: https://github.com/metabolite-atlas/metatlas/blob/master/HISTORY.rst
