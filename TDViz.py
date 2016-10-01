from traits.api import HasTraits, Button, Instance, List, Str, Enum, Float, File, Int
from traitsui.api import View, Item, VGroup, HSplit, CheckListEditor, HGroup, Group
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene
from mayavi import mlab
from astropy import coordinates as coord
from astropy import units as u
import numpy as np
import os
import glob
try:
	import astropy.io.fits as pyfits
except ImportError:
	import pyfits

class TDViz(HasTraits):
	fitsfile    = File(filter=[u"*.fits"])
	plotbutton1 = Button(u"Plot")
	plotbutton2 = Button(u"Plot")
	plotbutton3 = Button(u"Plot")
	clearbutton = Button(u"Clear")
	scene = Instance(MlabSceneModel, ())
	rendering = Enum("Surface-Spectrum", "Surface-Intensity", "Volume-Intensity")
	save_the_scene = Button(u"Save")
	save_in_file = Str("test.wrl")
	add_cut = Button(u"Cutthrough")
	remove_cut = Button(u"Remove the Last Cut")
	movie = Button(u"Movie")
	iteration	= Int(0)
	quality		= Int(8)
	delay		= Int(0)
	angle		= Int(360)
	spin = Button(u"Spin")
	zscale = Float(1.0)
	xstart = Float(0.0)
	xend   = Float(1.0)
	ystart = Float(0.0)
	yend   = Float(1.0)
	zstart = Float(0.0)
	zend   = Float(1.0)
	datamin= Float(0.0)
	datamax= Float(1.0)
	opacity= Float(0.3)

	view = View(
      HSplit(
		VGroup(
			Item("fitsfile", label=u"Select a FITS datacube", show_label=True),
			Item("rendering", tooltip=u"Choose the rendering type you like", show_label=True),
			Item('plotbutton1', tooltip=u"Plot the 3D scene with surface rendering, colored by spectrum", visible_when="rendering=='Surface-Spectrum'"),
			Item('plotbutton2', tooltip=u"Plot the 3D scene with surface rendering, colored by intensity", visible_when="rendering=='Surface-Intensity'"),
			Item('plotbutton3', tooltip=u"Plot the 3D scene with volume rendering, colored by intensity", visible_when="rendering=='Volume-Intensity'"),
			HGroup(Item('xstart', tooltip=u"starting pixel in X axis", show_label=True, springy=True),
				   Item('xend', tooltip=u"ending pixel in X axis", show_label=True, springy=True)),
			HGroup(Item('ystart', tooltip=u"starting pixel in Y axis", show_label=True, springy=True),
				   Item('yend', tooltip=u"ending pixel in Y axis", show_label=True, springy=True)),
			HGroup(Item('zstart', tooltip=u"starting pixel in Z axis", show_label=True, springy=True),
				   Item('zend', tooltip=u"ending pixel in Z axis", show_label=True, springy=True)),
			HGroup(Item('datamax', tooltip=u"Maximum datapoint shown", show_label=True, springy=True),
				   Item('datamin', tooltip=u"Minimum datapoint shown", show_label=True, springy=True)),
			Item('zscale', tooltip=u"Stretch the datacube in Z axis", show_label=True),
			Item('opacity', tooltip=u"Opacity of the scene", show_label=True),
			Item("add_cut", tooltip="Add a cutthrough view"),
			Item("remove_cut", tooltip="Remove all cutthroughs"),
			Item("spin", tooltip=u"Spin 360 degrees"),
			"clearbutton",
			Item('_'),
			Item("movie", tooltip="Make a GIF movie", show_label=False),
			HGroup(Item('iteration', tooltip=u"number of iterations, 0 means inf.", show_label=True),
				   Item('quality', tooltip=u"quality of plots, 0 is worst, 8 is good.", show_label=True)),
			HGroup(Item('delay', tooltip=u"time delay between frames, in millisecond.", show_label=True),
				   Item('angle', tooltip=u"angle the cube spins", show_label=True)),
			Item('_'),
			Item("save_the_scene", tooltip=u"Save current scene in a .wrl file", visible_when="rendering=='Surface-Spectrum' or rendering=='Surface-Intensity'"),
			Item("save_in_file", tooltip=u"3D model file name", show_label=True),
			show_labels=False
		),
		VGroup(
		    Item(name='scene',
                editor=SceneEditor(scene_class=MayaviScene),
                resizable=True,
                height=600,
                width=600
            ), show_labels=False
		)
      ),
	  resizable=True,
	  title=u"TDViz"
	)

	def _fitsfile_changed(self):
		img = pyfits.open(self.fitsfile)     # Read the fits data
		dat = img[0].data
		self.hdr  = img[0].header
		
		naxis = self.hdr['NAXIS']
		## The three axes loaded by pyfits are: velo, dec, ra
		## Swap the axes, RA<->velo
		if naxis == 4:
			self.data = np.swapaxes(dat[0],0,2)*1000.0
		elif naxis == 3:
			self.data = np.swapaxes(dat,0,2)*1000.0
		#onevpix = self.hdr['CDELT3']
		self.data[np.isnan(self.data)] = 0.0
		self.data[np.isinf(self.data)] = 0.0

		self.datamax = np.asscalar(np.max(self.data))
		self.datamin = np.asscalar(np.min(self.data))
		self.xend    = self.data.shape[0] - 1 
		self.yend    = self.data.shape[1] - 1
		self.zend    = self.data.shape[2] - 1

		self.data[self.data<self.datamin] = self.datamin

	def loaddata(self):
		channel = self.data
		## Reset the range if it is beyond the cube:
		if self.xstart < 0:
			print 'Wrong number!'
			self.xstart = 0
		if self.xend > channel.shape[0]-1:
			print 'Wrong number!'
			self.xend = channel.shape[0]-1
		if self.ystart < 0:
			print 'Wrong number!'
			self.ystart = 0
		if self.yend > channel.shape[1]-1:
			print 'Wrong number!'
			self.yend = channel.shape[1]-1
		if self.zstart < 0:
			print 'Wrong number!'
			self.zstart = 0
		if self.zend > channel.shape[2]-1:
			print 'Wrong number!'
			self.zend = channel.shape[2]-1
		## Select a region, use mJy unit
		region=channel[self.xstart:self.xend,self.ystart:self.yend,self.zstart:self.zend]

		## Stretch the cube in V axis
		from scipy.interpolate import splrep
		from scipy.interpolate import splev
		vol=region.shape
		stretch=self.zscale
		## Stretch parameter: how many times longer the V axis will be
		sregion=np.empty((vol[0],vol[1],vol[2]*stretch))
		chanindex=np.linspace(0,vol[2]-1,vol[2])
		for j in range(0,vol[0]-1):
		    for k in range(0,vol[1]-1):
		        spec=region[j,k,:]
		        tck=splrep(chanindex,spec,k=1)
		        chanindex2=np.linspace(0,vol[2]-1,vol[2]*stretch)
		        sregion[j,k,:]=splev(chanindex2,tck)
		self.sregion = sregion
		# Reset the max/min values
		if self.datamin < np.asscalar(np.min(self.sregion)):
			print 'Wrong number!'
			self.datamin = np.asscalar(np.min(self.sregion))
		if self.datamax > np.asscalar(np.max(self.sregion)):
			print 'Wrong number!'
			self.datamax = np.asscalar(np.max(self.sregion))
		self.xrang = abs(self.xstart - self.xend)
		self.yrang = abs(self.ystart - self.yend)
		self.zrang = abs(self.zstart - self.zend)*stretch

		## Keep a record of the coordinates:
		crval1 = self.hdr['crval1']
		cdelt1 = self.hdr['cdelt1']
		crpix1 = self.hdr['crpix1']
		crval2 = self.hdr['crval2']
		cdelt2 = self.hdr['cdelt2']
		crpix2 = self.hdr['crpix2']
		crval3 = self.hdr['crval3']
		cdelt3 = self.hdr['cdelt3']
		crpix3 = self.hdr['crpix3']

		ra_start = (self.xstart + 1 - crpix1) * cdelt1 + crval1
		ra_end = (self.xend + 1 - crpix1) * cdelt1 + crval1
		#if ra_start < ra_end:
		#	ra_start, ra_end = ra_end, ra_start
		dec_start = (self.ystart + 1 - crpix2) * cdelt2 + crval2
		dec_end = (self.yend + 1 - crpix2) * cdelt2 + crval2
		#if dec_start > dec_end:
		#	dec_start, dec_end = dec_end, dec_start
		vel_start = (self.zstart +1 - crpix3) * cdelt3 + crval3
		vel_end = (self.zend + 1 - crpix3) * cdelt3 + crval3
		#if vel_start < vel_end:
		#	vel_start, vel_end = vel_end, vel_start
		vel_start /= 1e3
		vel_end /= 1e3

		## Flip the V axis
		if cdelt3 > 0:
			self.sregion = self.sregion[:,:,::-1]
			vel_start, vel_end = vel_end, vel_start

		self.extent =[ra_start, ra_end, dec_start, dec_end, vel_start, vel_end]
	
	def labels(self):
		'''
		Add 3d text to show the axes.
		'''
		fontsize = max(self.xrang, self.yrang)/40.
		tcolor = (1,1,1)
		mlab.text3d(self.xrang/2,-40,self.zrang+40,'R.A.',scale=fontsize,orient_to_camera=True,color=tcolor)
		mlab.text3d(-40,self.yrang/2,self.zrang+40,'Decl.',scale=fontsize,orient_to_camera=True,color=tcolor)
		mlab.text3d(-40,-40,self.zrang/2-10,'V (km/s)',scale=fontsize,orient_to_camera=True,color=tcolor)
		# Label the coordinates of the corners
		# Lower left corner
		ra0 = self.extent[0]; dec0 = self.extent[2]
		c = coord.ICRS(ra=ra0, dec=dec0, unit=(u.degree, u.degree))
		RA_ll = str(int(c.ra.hms.h))+'h'+str(int(c.ra.hms.m))+'m'+str(round(c.ra.hms.s,1))+'s'
		mlab.text3d(0,-20,self.zrang+20,RA_ll,scale=fontsize,orient_to_camera=True,color=tcolor)
		DEC_ll = str(int(c.dec.dms.d))+'d'+str(int(abs(c.dec.dms.m)))+'m'+str(round(abs(c.dec.dms.s),1))+'s'
		mlab.text3d(-80,0,self.zrang+20,DEC_ll,scale=fontsize,orient_to_camera=True,color=tcolor)
		# Upper right corner
		ra0 = self.extent[1]; dec0 = self.extent[3]
		c = coord.ICRS(ra=ra0, dec=dec0, unit=(u.degree, u.degree))
		RA_ll = str(int(c.ra.hms.h))+'h'+str(int(c.ra.hms.m))+'m'+str(round(c.ra.hms.s,1))+'s'
		mlab.text3d(self.xrang,-20,self.zrang+20,RA_ll,scale=fontsize,orient_to_camera=True,color=tcolor)
		DEC_ll = str(int(c.dec.dms.d))+'d'+str(int(abs(c.dec.dms.m)))+'m'+str(round(abs(c.dec.dms.s),1))+'s'
		mlab.text3d(-80,self.yrang,self.zrang+20,DEC_ll,scale=fontsize,orient_to_camera=True,color=tcolor)
		# V axis
		if self.extent[5] > self.extent[4]:
			v0 = self.extent[4]; v1 = self.extent[5]
		else:
			v0 = self.extent[5]; v1 = self.extent[4]
		mlab.text3d(-20,-20,self.zrang,str(round(v0,1)),scale=fontsize,orient_to_camera=True,color=tcolor)
		mlab.text3d(-20,-20,0,str(round(v1,1)),scale=fontsize,orient_to_camera=True,color=tcolor)

		mlab.axes(self.field, ranges=self.extent, x_axis_visibility=False, y_axis_visibility=False, z_axis_visibility=False)
		mlab.outline()

	def _plotbutton1_fired(self):
		mlab.clf()
		self.loaddata()
		self.sregion[np.where(self.sregion<self.datamin)] = self.datamin
		self.sregion[np.where(self.sregion>self.datamax)] = self.datamax

		# The following codes from: http://docs.enthought.com/mayavi/mayavi/auto/example_atomic_orbital.html#example-atomic-orbital
		field = mlab.pipeline.scalar_field(self.sregion)     # Generate a scalar field
		colored = self.sregion
		vol=self.sregion.shape
		for v in range(0,vol[2]-1):
			colored[:,:,v] = self.extent[4] + v*(-1)*abs(self.hdr['cdelt3'])
		new = field.image_data.point_data.add_array(colored.T.ravel())
		field.image_data.point_data.get_array(new).name = 'color'
		field.image_data.point_data.update()

		field2 = mlab.pipeline.set_active_attribute(field, point_scalars='scalar')
		contour = mlab.pipeline.contour(field2)
		contour2 = mlab.pipeline.set_active_attribute(contour, point_scalars='color')

		mlab.pipeline.surface(contour2, colormap='jet', opacity=self.opacity)
		
		## Insert a continuum plot
		##im = pyfits.open('g28_SMA1.cont.image.fits')
		##dat = im[0].data
		##dat0 = dat[0]
		##channel = dat0[0]
		##region = np.swapaxes(channel[self.xstart:self.xend,self.ystart:self.yend]*1000.,0,1)
		##field = mlab.contour3d(region, colormap='gist_ncar')
		##field.contour.minimum_contour = 5

		self.field = field
		self.labels()
		mlab.view(azimuth=0, elevation=0, distance='auto')
		mlab.show()

	def _plotbutton2_fired(self):
		mlab.clf()
		self.loaddata()
		#field=mlab.contour3d(self.sregion,colormap='gist_ncar')     # Generate a scalar field
		field=mlab.contour3d(self.sregion)     # Generate a scalar field
		field.contour.maximum_contour = self.datamax
		field.contour.minimum_contour = self.datamin
		field.actor.property.opacity = self.opacity

		self.field = field
		self.labels()
		mlab.view(azimuth=0, elevation=0, distance='auto')
		mlab.show()
	
	def _plotbutton3_fired(self):
		mlab.clf()
		self.loaddata()
		field = mlab.pipeline.scalar_field(self.sregion) # Generate a scalar field
		mlab.pipeline.volume(field,vmax=self.datamax,vmin=self.datamin)
		
		self.field = field
		self.labels()
		mlab.view(azimuth=0, elevation=0, distance='auto')
		mlab.show()

#	def _datamax_changed(self):
#		if hasattr(self, "field"):
#			self.field.contour.maximum_contour = self.datamax

	def _add_cut_fired(self):
		self.cut=mlab.pipeline.scalar_cut_plane(self.field,plane_orientation="x_axes")
		self.cut.enable_contours=True
		self.cut.contour.number_of_contours=5

	def _remove_cut_fired(self):
		self.cut.stop()

	def _save_the_scene_fired(self):
		mlab.savefig(self.save_in_file)

	def _movie_fired(self):
		if os.path.exists("./tenpfigz"):
			print "The chance of you using this name is really small..."
		else:
			os.system("mkdir tenpfigz")

		if filter(os.path.isfile, glob.glob("./tenpfigz/*.jpg")) != []:
			os.system("rm -rf ./tenpfigz/*.jpg")

		i = 0
		## Quality of the movie: 0 is the worst, 8 is ok.
		self.field.scene.anti_aliasing_frames = self.quality
		self.field.scene.disable_render = True
		mlab.savefig('./tenpfigz/screenshot0'+str(i)+'.jpg')
		while i<(self.angle/5):
			self.field.scene.camera.azimuth(5)
			self.field.scene.render()
			i += 1
			if i<10:
				mlab.savefig('./tenpfigz/screenshot0'+str(i)+'.jpg')
			elif 9<i<100:
				mlab.savefig('./tenpfigz/screenshot'+str(i)+'.jpg')
		self.field.scene.disable_render = False

		os.system("convert -delay "+str(self.delay)+" -loop "+str(self.iteration)+" ./tenpfigz/*.jpg ./tenpfigz/animation.gif")
	
	def _spin_fired(self):
		i = 0
		self.field.scene.disable_render = True
		@mlab.animate
		def anim():
			while i<72:
				self.field.scene.camera.azimuth(5)
				self.field.scene.render()
				yield
		a = anim()
		#while i<72:
		#	self.field.scene.camera.azimuth(5)
		#	self.field.scene.render()
		#	i += 1
		#	#mlab.savefig('./'+str(i)+'.png')
		self.field.scene.disable_render = False

	def _clearbutton_fired(self):
		mlab.clf()

app = TDViz()
app.configure_traits()   
