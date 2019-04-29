# A simple tool of visualizing 3-D datacubes in FITS format

The visualization is implemented by MayaVI and the interface is made with TraitsUI.
It stems from my homework for the [ay201b](https://ay201b.wordpress.com/topical-modules/2013-topical-modules/) course.

![interface](https://github.com/xinglunju/tdviz/blob/master/screenshot.png)

## Denpendencies

* Mayavi,
* pyface,
* Traits,
* TraitsUI,
* astropy,
* numpy,
* scipy,
* (optional) ImageMagick to make GIF movies.

## Versions

Now it is full compatible with Python 3 (using TDViz3.py). The old file (TDViz.py) may still work with Python 2 but is not guarenteed.

I tested it with the following environment and everything works well:

* Python==3.7.3
* mayavi==4.6.2
* pyface==6.0.0
* traits==5.1.1
* traitsui==6.0.0
* astropy==3.1.2
* numpy==1.16.3
* scipy==1.2.1

## Render P-P-V cubes with three options

1. Iso-surfaces colored by velocities.
2. Iso-surfaces colored by intensities.
3. Scattered dots colored by intensities.

The first two can be saved in a 'mesh' file and be imported to softwares such as Blender or Meshlab, and be uploaded to Sketchfab.

## Known issues

1. Cannot load a different fits file if the current file name includes path: have to delete the current 'fitsfile' name, then nagivate to the new fits file.
2. To be found...
