import gadget_snap
import pylab
import struct
import os
import numpy
import scipy.io
from utilities import *
	
def gadget_readsnap( snapshot, snappath=".", snapbase="snapshot_", verbose=False, onlyHeader=False, loadonly=False, nommap=False ):
	return gadget_snap.gadget_snapshot( "%s/%s%03d" % (snappath, snapbase, snapshot), verbose=verbose, loadonly=loadonly, nommap=nommap, onlyHeader=onlyHeader )

def gadget_readsnapname( snapshot, verbose=False, loadonly=False, nommap=False, onlyHeader=False ):
	return gadget_snap.gadget_snapshot( snapshot, verbose=verbose, loadonly=loadonly, nommap=nommap, onlyHeader=onlyHeader )
	
def gadget_snapexists( snapshot, snappath=".", snapbase="snapshot_" ):
	return os.path.exists( "%s/%s%03d" % (snappath, snapbase, snapshot) ) or os.path.exists( "%s/%s%03d.0" % (snappath, snapbase, snapshot) )

def gadget_snapcount( snappath=".", snapbase="snapshot_" ):
	i = 0
	while gadget_snapexists( i, snappath=snappath, snapbase=snapbase ):
		i += 1
	return i

def gadget_merge_ics( outfile, filename1, filename2, offset1, offset2, voffset1=[0.,0.,0.], voffset2=[0.,0.,0.] ):
	snap1 = gadget_readsnapname( filename1 )
	snap2 = gadget_readsnapname( filename2 )

	for i in range(3):
		snap1.pos[:,i] += offset1[i]
		snap2.pos[:,i] += offset2[i]
		
	for i in range(3):
		snap1.vel[:,i] += voffset1[i]
		snap2.vel[:,i] += voffset2[i]

	npart = snap1.npart + snap2.npart
	data = {}
	data[ 'count' ] = npart
	data[ 'pos' ] = pylab.zeros( [npart, 3] )
	data[ 'pos' ][ 0:snap1.npart, : ] = snap1.pos
	data[ 'pos' ][ snap1.npart:npart, : ] = snap2.pos
	data[ 'vel' ] = pylab.zeros( [npart, 3] )
	data[ 'vel' ][ 0:snap1.npart, : ] = snap1.vel
	data[ 'vel' ][ snap1.npart:npart, : ] = snap2.vel
	data[ 'mass' ] = pylab.zeros( npart )
	data[ 'mass' ][ 0:snap1.npart ] = snap1.data["mass"]
	data[ 'mass' ][ snap1.npart:npart ] = snap2.data["mass"]
	data[ 'u' ] = pylab.zeros( npart )
	data[ 'u' ][ 0:snap1.npart ] = snap1.data["u"]
	data[ 'u' ][ snap1.npart:npart ] = snap2.data["u"]
	nxnuc = pylab.shape( snap1.data["xnuc"] )[1]
	data[ 'xnuc' ] = pylab.zeros( [npart, nxnuc] )
	data[ 'xnuc' ][ 0:snap1.npart, : ] = snap1.data["xnuc"]
	data[ 'xnuc' ][ snap1.npart:npart, : ] = snap2.data["xnuc"]

	gadget_write_ics( outfile, data, transpose=False )
	return
	
def gadget_write_ics( filename, data, transpose=True, time=0., skipxnuc=False ):
	f = open( filename, 'wb' )
	
	label = 'HEAD'
	blksize = 256
	npart = pylab.zeros( 6, dtype=pylab.int32 )
	npart[0] = data['count']
	massarr = pylab.zeros( 6, dtype=pylab.float64 )
	redshift = 0
	flag_sfr = 0
	flag_feedback = 0
	la = pylab.zeros( 160, dtype=pylab.int8 )
	
	header = struct.pack( "i4sii", 8, label, blksize+8, 8 )
	f.write( header )
	
	f.write( struct.pack( "i", blksize ) )
	npart.tofile( f )
	massarr.tofile( f )
	f.write( struct.pack( "ddii", time, redshift, flag_sfr, flag_feedback ) )
	la.tofile( f )
	f.write( struct.pack( "i", blksize ) )
	
	for field in ['POS ','VEL ']:
		d = data[ field.strip().lower() ]
		header = struct.pack( "i4sii", 8, field, d.size*4+8, 8 )
		f.write( header )
		f.write( struct.pack( "i", d.size*4 ) )
		if d.shape[0] == data['count']:
			d.astype( pylab.float32 ).tofile( f )
		elif d.shape[1] == data['count']:
			d.transpose().astype( pylab.float32 ).tofile( f )
		else:
			print 'Error: ' + field + ' has the wrong dimensions.'
		f.write( struct.pack( "i", d.size*4 ) )
	
	header = struct.pack( "i4sii", 8, 'ID  ', data['count']*4+8, 8 )
	f.write( header )
	f.write( struct.pack( "i", data['count']*4 ) )
	if data.has_key( "id" ):
		data["id"].astype( pylab.uint32 ).tofile( f )
	else:
		pylab.arange( 1, data['count']+1 ).astype( pylab.uint32 ).tofile( f )
	f.write( struct.pack( "i", data['count']*4 ) )
	
	if skipxnuc:
		fields = ['MASS', 'U   ']
	else:
		fields = ['MASS', 'U   ', 'XNUC']
	
	for field in fields:
		d = data[ field.strip().lower() ]
		header = struct.pack( "i4sii", 8, field, d.size*4+8, 8 )
		f.write( header )
		f.write( struct.pack( "i", d.size*4 ) )
		if d.shape[0] == data['count']:
			d.astype( pylab.float32 ).tofile( f )
		elif d.ndim > 1 and d.shape[1] == data['count']:
			d.transpose().astype( pylab.float32 ).tofile( f )
		else:
			d.astype( pylab.float32 ).tofile( f )
			print 'Error: ' + field + ' has the wrong dimensions.'
		f.write( struct.pack( "i", d.size*4 ) )
	
	f.close()
	return

def gadget_readenergy( snappath="./" ):
	if not snappath.endswith( "/" ):
		snappath += "/"

	f = open( snappath + "energy.txt", "r" )
	e = scipy.io.read_array( f )
	f.close()

	data = {}
	data['time'] = e[:,0]
	data['ein' ] = e[:,1]
	data['epot'] = e[:,2]
	data['ekin'] = e[:,3]
	data['etot'] = e[:,1:4].sum(axis=1)
	data['eadd'] = e[:,pylab.shape(e)[1]-1]

	return dict2obj( data )

