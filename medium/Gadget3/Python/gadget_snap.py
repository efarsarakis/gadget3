import numpy as np
import calcGrid
import os.path as path
import pylab
import struct
import time
from const import *
from math import *
from utilities import *

class gadget_snapshot:
	def __init__( self, filename, verbose=False, loadonly=False, onlyHeader=False, nommap=False ):
		if type( loadonly ) == np.ndarray:
			loadlist = loadonly
			loadonly = True
		else:
			loadlist = False

		self.filecount = 1
		if path.exists( filename ):
			self.files = [filename]
		elif path.exists( filename + '.0' ):
			self.files = [filename + '.0']
			while path.exists( filename + ".%d" % self.filecount ):
				self.files += [filename + ".%d" % self.filecount]
				self.filecount += 1

			if not nommap:
				print "Multiple files detected, thus mmap is deactivated."
				nommap = True
		else:
			raise Exception( "Neither %s nor %s.0 exists." % (filename, filename) )

		self.center = pylab.array( [0,0,0], dtype="float64" )

		self.datablocks_skip = ["HEAD"]
		self.datablocks_int32 = ["ID  "]

		self.get_blocks( 0, verbose=verbose )
		self.load_header( 0, verbose=verbose )

		if onlyHeader:
			return

		if self.filecount == 1 and not nommap:
			self.load_data_mmap( verbose=verbose )
		else:
			self.load_data( verbose=verbose, loadonly=loadonly, loadlist=loadlist )

		self.convenience()
		return

	def get_blocks( self, fileid, verbose=False ):
		swap, endian = endianness_check( self.files[fileid] )
		
		self.blocks = {}
		f = open( self.files[fileid], 'r' )
		
		fpos = f.tell()
		s = f.read(16)
		while len(s) > 0:
			fheader, name, length, ffooter = struct.unpack( endian + "i4sii", s )
			self.blocks[ name ] = fpos

			print "Block %s, length %d, offset %d." %(name, length, fpos)

			f.seek( length, 1 )
			fpos = f.tell()
			s = f.read(16)

		if verbose:
			print "%d blocks detected." % len(self.blocks.keys())
		return

	def load_header( self, fileid, verbose=False ):
		swap, endian = endianness_check( self.files[fileid] )
		
		f = open( self.files[fileid], 'r' )
		s = f.read(16)
		while len(s) > 0:
			fheader, name, length, ffooter = struct.unpack( endian + "i4sii", s )

			if name != "HEAD":
				f.seek( length, 1 )
				s = f.read( 16 )
			else:
				f.seek( 4, 1 ) # skip fortran header of data block
				s = f.read(24)
				self.nparticles = struct.unpack( endian + "6i", s )
				s = f.read(48)
				self.masses = struct.unpack( endian + "6d", s )
				s = f.read(24)
				self.time, self.redshift, self.flag_sfr, self.flag_feedback = struct.unpack( endian + "ddii", s )
				s = f.read(24)
				self.nparticlesall = struct.unpack( endian + "6i", s )
				s = f.read(16)
				self.flag_cooling, num_files, self.boxsize = struct.unpack( endian + "iid", s )
				s = f.read(120)

				self.nparticles = pylab.array( self.nparticles )
				self.nparticlesall = pylab.array( self.nparticlesall )

				sort = self.nparticlesall.argsort()
				if self.nparticlesall[ sort[::-1] ][1] == 0:
					self.singleparticlespecies = True
				else:
					self.singleparticlespecies = False

				if self.nparticlesall[1] > 0:
					if self.singleparticlespecies:
						print "Pure DARKNESS detected."
					else:
						print "DARKNESS detected."

				f.close()
				s = "" # stop reading

				self.npart = pylab.array( self.nparticles ).sum()
				self.npartall = pylab.array( self.nparticlesall ).sum()

				if fileid == 0:
					if (num_files != self.filecount) and not (num_files == 0 and self.filecount == 1):
						raise Exception( "Number of files detected (%d) and num_files in the header (%d) are inconsistent." % (self.filecount, num_files) )
		print "Snapshot contains %d particles." % self.npartall
		return

	def load_data_mmap( self, verbose=False ):
		# memory mapping works only if there is only one file
		swap, endian = endianness_check( self.files[0] )

		self.data = {}

		self.get_blocks( 0, verbose=verbose )
		f = open( self.files[0], 'r' )
		for block in self.blocks.keys():
			# skip some blocks
			if block in self.datablocks_skip:
				continue

			f.seek( self.blocks[block], 0 )
			fheader, name, length, ffooter = struct.unpack( endian + "i4sii", f.read(16) )

			# we assume that the output contains only float32 / int32 fields
			elements = (length-8)/4

			if self.singleparticlespecies:
				npart = self.npart
			else:
				npart, npartall = get_block_size_from_table( name )
			dim = elements / npart

			if verbose:
				print "Loading block %s, offset %d, length %d, elements %d, dimension %d, particles %d." % (block, self.blocks[block], length, elements, dim, npart)

			offset = f.tell() + 4

			if block in self.datablocks_int32:
				blocktype = "i4"
			else:
				blocktype = "f4"

			if dim == 1:
				self.data[ block.strip().lower() ] = np.memmap( self.files[0], offset=offset, mode='c', dtype=endian+blocktype, shape=(npart) )
			else:
				self.data[ block.strip().lower() ] = np.memmap( self.files[0], offset=offset, mode='c', dtype=endian+blocktype, shape=(npart, dim) )
		return

	def load_data( self, verbose=False, nommap=False, loadonly=False, loadlist=[] ):
		swap, endian = endianness_check( self.files[0] )
		self.data = {}

		lb = 0
		for fileid in range( self.filecount ):
			self.load_header( fileid, verbose=verbose )
			self.get_blocks( fileid, verbose=verbose )

			ub = lb + self.npart
			
			f = open( self.files[fileid], 'r' )
			for block in self.blocks.keys():
				# skip some blocks
				if block in self.datablocks_skip:
					continue

				if loadonly and not block in loadlist:
					continue

				if verbose:
					print "Loading block %s of file %s." % (block, fileid)

				f.seek( self.blocks[block], 0 )
				fheader, name, length, ffooter = struct.unpack( endian + "i4sii", f.read(16) )

				# we assume that the output contains only float32 / int32 fields
				elements = (length-8)/4
				if self.singleparticlespecies:
					npart = self.npart
					npartall = self.npartall
				else:
					npart, npartall = get_block_size_from_table( name )
				dim = elements / npart

				if block in self.datablocks_int32:
					blocktype = 'i4'
				else:
					blocktype = 'f4'

				blockname = block.strip().lower()
				if not self.data.has_key( block ):
					if dim == 1:
						self.data[ blockname ] = pylab.zeros( npartall, dtype=blocktype )
					else:
						self.data[ blockname ] = pylab.zeros( (npartall, dim), dtype=blocktype )

					if blocktype == 'f4':
						self.data[ blockname ] = self.data[ blockname ].astype( 'float64' ) # change array type to float64

				if verbose:
					print "Block contains %d elements." %elements

				f.seek( 4, 1 ) # skip fortran header of data field
				if dim == 1:
					self.data[ blockname ][lb:ub] = np.fromfile( f, dtype=endian+blocktype, count=elements )
				else:
					self.data[ blockname ][lb:ub,:] = np.fromfile( f, dtype=endian+blocktype, count=elements ).reshape( npart, dim )

			lb = ub
		
		if self.filecount > 1:
			self.npart = self.npartall
		print '%d particles loaded.' % self.npartall
		return

	def get_block_size_from_table( self, blockname ):
		npart = 0
		npartall = 0
		if block in ["POS ", "VEL ", "ID  "]:
			# present for all particles
			return self.nparticles.sum(), self.nparticlesall.sum()
		if block in ["U   ","RHO ", "NE  ", "NH  ", "HSML", "SFR ", "Z   ", "XNUC", "PRES"]:
			# present for hydro particles
			npart += self.nparticles[0]
			npartall += self.nparticlesall[0]
		if block in ["AGE ", "Z   "]:
			# present for star particles
			npart += self.nparticles[4]
			npartall += self.nparticlesall[4]
		if block in ["BHMA", "BHMD"]:
			#present for black hole particles
			npart += self.nparticles[5]
			npartall += self.nparticlesall[5]
		return npart, npartall

	def convenience( self ):
		if self.data.has_key( "pos" ):
			self.pos = self.data["pos"]
		
		if self.data.has_key( "rho" ):
			self.rho = self.data["rho"]
		
		if self.data.has_key( "vel" ):
			self.vel = self.data["vel"]
		
		if self.data.has_key( "id" ):
			self.id = self.data["id"]
		
		if self.data.has_key( "mass" ):
			self.mass = self.data["mass"]
		return

	def set_center( self, center ):
		if type( center ) == list:
			self.center = pylab.array( center )
		elif type( center ) != np.ndarray:
			self.center = center
		else:
			raise Exception( "center has to be of type list or numpy.ndarray" )
		return
		
	def r( self, center=False ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center

		radius = pylab.sqrt( (self.data[ "pos" ][:,0]-center[0])**2+(self.data[ "pos" ][:,1]-center[1])**2+(self.data[ "pos" ][:,2]-center[2])**2 )
		return radius

	def get_slice( self, value, box=[0,0], nx=200, ny=200, center=False, axes=[0,1] ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center
		
		dim0 = axes[0]
		dim1 = axes[1]
		
		if (box[0] == 0 and box[1] == 0):
			box[0] = max( abs( self.data[ "pos" ][:,dim0] ) ) * 2
			box[1] = max( abs( self.data[ "pos" ][:,dim1] ) ) * 2

		if (value == "mass"):
			return calcGrid.calcDensSlice( self.data["pos"].astype('float64'), self.data["hsml"].astype('float64'), self.data[value].astype('float64'), nx, ny, box[0], box[1], center[0], center[1], center[2], dim0, dim1 )
		else:
			return calcGrid.calcSlice( self.data["pos"].astype('float64'), self.data["hsml"].astype('float64'), self.data["mass"].astype('float64'), self.data["rho"].astype('float64'), self.data[value].astype('float64'), nx, ny, box[0], box[1], center[0], center[1], center[2], dim0, dim1 )

	def get_raddens( self, nshells=200, dr=0, center=False ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center
		return calcGrid.calcRadialProfile( self.pos.astype('float64'), self.data["mass"].astype('float64'), 1, nshells, dr, center[0], center[1], center[2] )

	def get_radprof( self, value, nshells=200, dr=0, center=False ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center
		return calcGrid.calcRadialProfile( self.pos.astype('float64'), self.data[value].astype('float64'), 2, nshells, dr, center[0], center[1], center[2] )
	
	def get_radmassprof( self, nshells, dr=0, center=False, solarmass=False ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center		

		p = calcGrid.calcRadialProfile( self.pos.astype('float64'), self.data["mass"].astype('float64'), 0, nshells, dr, center[0], center[1], center[2] )
		for i in range( 1, nshells ):
			p[0,i] += p[0,i-1]
		if solarmass:
			p[0,:] /= 1.989e33
		return p

	def plot_raddens( self, log=False, nshells=200, dr=0, center=False, color='k' ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center	
		
		p = self.get_raddens( nshells, dr, center )
		if log:
			pylab.semilogy( p[1,:], p[0,:], color )
		else:
			pylab.plot( p[1,:], p[0,:], color )

	def plot_radprof( self, value, log=False, nshells=200, dr=0, center=False, color='k' ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center	
		
		p = self.get_radprof( value, nshells, dr, center )
		if log:
			pylab.semilogy( p[1,:], p[0,:], color )
		else:
			pylab.plot( p[1,:], p[0,:], color )

	def plot_radmassprof( self, log=False, nshells=200, dr=0, center=False, color='k', solarmass=False ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center
		
		p = self.get_radmassprof( nshells, dr, center, solarmass )
		if log:
			pylab.semilogy( p[1,:], p[0,:], color )
		else:
			pylab.plot( p[1,:], p[0,:], color )
	
	def plot_radvecprof( self, value, log=False, nshells=200, dr=0, center=False, color='k' ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center
		
		vec = (self.data[value]*self.data['pos']).sum(axis=1)/self.r()
		p = calcGrid.calcRadialProfile( self.pos.astype('float64'), vec.astype('float64'), 2, nshells, dr, center[0], center[1], center[2] )
		if log:
			pylab.semilogy( p[1,:], p[0,:], color )
		else:
			pylab.plot( p[1,:], p[0,:], color )

	def plot_pos( self, axes=[0,1] ):
		pylab.plot( self.pos[:,axes[0]], self.pos[:,axes[1]], ',' )
		pylab.axis( "scaled" )

	def plot_slice( self, value, logplot=True, colorbar=False, box=[0,0], nx=200, ny=200, center=False, axes=[0,1], minimum=1e-8, newfig=True ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center
		
		dim0 = axes[0]
		dim1 = axes[1]
		
		if (box[0] == 0 and box[1] == 0):
			box[0] = max( abs( self.data[ "pos" ][:,dim0] ) ) * 2
			box[1] = max( abs( self.data[ "pos" ][:,dim1] ) ) * 2

		slice = self.get_slice( value, box, nx, ny, center, axes )
		x = (pylab.array( range( nx+1 ) ) - nx/2.) / nx * box[0]
		y = (pylab.array( range( ny+1 ) ) - ny/2.) / ny * box[1]

		if (newfig):
			fig = pylab.figure( figsize = ( 13, int(12*box[1]/box[0] + 0.5) ) )
			pylab.spectral()
		
		if logplot:
			pc = pylab.pcolor( x, y, pylab.transpose( pylab.log10( pylab.maximum( slice, minimum ) ) ), shading='flat' )
		else:
			pc = pylab.pcolor( x, y, pylab.transpose( slice ), shading='flat' )
		if colorbar:
			cb = pylab.colorbar()
		pylab.axis( "image" )

		xticklabels = []
		for tick in pc.axes.get_xticks():
			if (tick == 0):
				xticklabels += [ r'$0.0$' ]
			else:
				xticklabels += [ r'$%.2f \cdot 10^{%d}$' % (tick/10**(ceil(log10(abs(tick)))), ceil(log10(abs(tick)))) ]
		pc.axes.set_xticklabels( xticklabels, size=16, y=-0.1, va='baseline' )

		yticklabels = []
		for tick in pc.axes.get_yticks():
			if (tick == 0):
				yticklabels += [ r'$0.0$' ]
			else:
				yticklabels += [ r'$%.2f \cdot 10^{%d}$' % (tick/10**(ceil(log10(abs(tick)))), ceil(log10(abs(tick)))) ]
		pc.axes.set_yticklabels( yticklabels, size=16, ha='right' )
		return pc

	def plot_cylav( self, value, logplot=True, box=[0,0], nx=512, ny=512, center=False, minimum=1e-8 ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center
		
		if (box[0] == 0 and box[1] == 0):
			box[0] = max( abs( self.data[ "pos" ][:,0] ) ) * 2
			box[1] = max( abs( self.data[ "pos" ][:,1:] ) ) * 2

		grid = calcGrid.calcGrid( self.pos.astype('float64'), self.data["hsml"].astype('float64'), self.data["mass"].astype('float64'), self.data["rho"].astype('float64'), self.data[value].astype('float64').astype('float64'), nx, ny, ny, box[0], box[1], box[1], 0, 0, 0 )
		cylav = calcGrid.calcCylinderAverage( grid )
		x = (pylab.array( range( nx+1 ) ) - nx/2.) / nx * box[0]
		y = (pylab.array( range( ny+1 ) ) - ny/2.) / ny * box[1]

		fig = pylab.figure( figsize = ( 13, int(12*box[1]/box[0] + 0.5) ) )
		pylab.spectral()
		
		if logplot:
			pc = pylab.pcolor( x, y, pylab.transpose( pylab.log10( pylab.maximum( cylav, minimum ) ) ), shading='flat' )
		else:
			pc = pylab.pcolor( x, y, pylab.transpose( slice ), shading='flat' )

		pylab.axis( "image" )
		xticklabels = []
		for tick in pc.axes.get_xticks():
			if (tick == 0):
				xticklabels += [ r'$0.0$' ]
			else:
				xticklabels += [ r'$%.2f \cdot 10^{%d}$' % (tick/10**(ceil(log10(abs(tick)))), ceil(log10(abs(tick)))) ]
		pc.axes.set_xticklabels( xticklabels, size=16, y=-0.1, va='baseline' )

		yticklabels = []
		for tick in pc.axes.get_yticks():
			if (tick == 0):
				yticklabels += [ r'$0.0$' ]
			else:
				yticklabels += [ r'$%.2f \cdot 10^{%d}$' % (tick/10**(ceil(log10(abs(tick)))), ceil(log10(abs(tick)))) ]
		pc.axes.set_yticklabels( yticklabels, size=16, ha='right' )
		return pc

	def getbound( self, center=False, vel=[0,0,0] ):
		if type( center ) == list:
			center = pylab.array( center )
		elif type( center ) != np.ndarray:
			center = self.center
		
		start = time.time()
		radius = pylab.zeros( self.npartall )
		for i in range( 3 ):
			radius += (self.data[ "pos" ][:,i] - center[i])**2
		radius = pylab.sqrt( radius )
		rs = radius.argsort()
		
		mass = 0.
		bcount = 0.
		bmass = 0.
		bcenter = [0., 0., 0.]
		bparticles = []
		for part in range( self.npart ):
			if (part == 0) or (( ( self.vel[rs[part],:] - vel )**2. ).sum() < 2.*G*mass/radius[rs[part]]):
				bcount += 1.
				bmass += self.data['mass'][rs[part]]
				bcenter += self.pos[rs[part],:]
				bparticles += [self.id[rs[part]]]
			mass += self.data['mass'][rs[part]]
		
		bobject = {}
		bobject['mass'] = bmass
		bobject['center'] = bcenter / bcount
		bobject['count'] = bcount
		bobject['particles'] = bparticles
		
		print "Calculation took %gs." % (time.time()-start)
		return bobject
