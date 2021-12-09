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


Documentation
-------------

`Targeted analysis workflow <https://github.com/biorack/metatlas/blob/main/docs/Targeted_Analysis.md>`_
