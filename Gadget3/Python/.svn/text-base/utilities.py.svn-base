import struct
import os.path as path
import pylab
from pylab import rcParams

class anyobject:
	def __init__( self ):
		return

def dict2obj( dict ):
	o = anyobject()
	for key in dict.keys():
		o.__dict__[ key ] = dict[ key ]
	
	return o

def obj2dict( obj ):
	d = {}
	for key in obj.__dict__.keys():
		d[ key ] = obj.__dict__[ key ]
	
	return d

def endianness_local():
	s = struct.pack( "bb", 1, 0 )
	r, = struct.unpack( "h", s )

	if r == 1:
		return '<'
	else:
		return '>'

def endianness_check( filename ):
	if not path.exists( filename ):
		print "File %s does not exist." % filename
		return False

	filesize = path.getsize( filename )
	
	endian_local = endianness_local()
	endian_data = endian_local

	f = open( filename )
	s = f.read(4)
	if len(s) > 0:
		size, = struct.unpack( "<i", s )
		size = abs( size )

		if size < filesize:
			f.seek( size, 1 )
			s = f.read(4)
		else:
			s = ""
		
		if (len(s) > 0) and (struct.unpack( "<i", s )[0] == size):
			f.close()
			endian_data = "<"
		else:
			f.seek( 0, 0 )
			size, = struct.unpack( ">i", f.read(4) )
			size = abs( size )
			
			if size < filesize:
				f.seek( size, 1 )
				s = f.read(4)
			else:
				s = ""
			
			if (len(s) > 0) and (struct.unpack( ">i", s )[0] == size):
				f.close()
				endian_data = ">"
			else:
				f.close()
				print "File %s is corrupt." % filename
				return False
	else:
		f.close()
		print "File %s is empty." % filename

	return (endian_local != endian_data, endian_data)

def imsave( filename, data, **kwargs ):
	figsize = (pylab.array(data.shape)/100.0)[::-1]
	rcParams.update( {'figure.figsize':figsize} )
	fig = pylab.figure( figsize=figsize )
	pylab.axes( [0,0,1,1] )
	pylab.axis( 'off' )
	fig.set_size_inches( figsize )
	pylab.imshow( data, origin='lower', **kwargs )
	pylab.savefig( filename, facecolor='black', edgecolor='black', dpi=100 )
	pylab.close( fig )


