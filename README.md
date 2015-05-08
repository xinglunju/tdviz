#A simple tool of visualizing 3-D datacubes in FITS format

The visualization is implemented by MayaVI and the interface is made with TraitsUI.
I removed those 'enthought' prefixes of modules to make it compatible with other installations (e.g. Anaconda).

#Known issues

1. Cannot reload a different fits file if current one includes path: have to delete current 'fitsfile' selection first...
