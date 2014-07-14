# all constants in cgs

# gravitational constant
G = 6.674e-8		# [ cm^3 g^-1 s^-1 ]

# avogadro constant
NA = 6.0221418e23	# [ ]

# boltzmann constant
KB = 1.3806504e-16	# [ erg K^-1 ]

# planck constant
H = 6.62606896e-27	# [ erg s ]

# speed of light in vacuum
c = 2.99792458e10	# [ cm s^-1 ]

# solar mass
msol = 1.989e33		# [ g ]

# solar radius
rsol = 6.955e10		# [ cm ]

# solar luminosity
lsol = 3.839e33		# [ erg s^-1 ]

# electron charge
qe = 4.80320427e-10 # [ esu ]

# ev2erg
ev2erg = 1.602177e-12 # [ erg eV^-1 ]

# colors
colorset = [ 'k', 'b', 'r', 'g', 'c', 'm', 'y' ]

# atomic symbols
asymb = [ 'neutron', 'h', 'he',
          'li', 'be',  'b',  'c',  'n',  'o',  'f', 'ne',
          'na', 'mg', 'al', 'si',  'p',  's', 'cl', 'ar',
           'k', 'ca', 'sc', 'ti',  'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
          'rb', 'sr',  'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te',  'i', 'xe',
          'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf',
          'ta',  'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn',
          'fr', 'ra', 'ac', 'th', 'pa',  'u', 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr', 'rf',
          'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg','uub','uut','uuq','uup','uuh','uus','uuo' ]

# atomic numbers to symbols
asymb2num = {}
for i in range( len( asymb ) ):
	if (asymb2num.has_key( asymb[i] )):
		print "asymb2num key error: key `%s` already exists." % asymb[i]
	else:
		asymb2num[ asymb[i] ] = i

class nuc:
	def __init__( self, name ):
		self.name = name

		nucname = ""
		nuca = ""

		for j in range( len( self.name ) ):
			if self.name[j].isalpha():
				nucname += self.name[j]
			else:
				nuca += self.name[j]
		
		self.symb = nucname

		self.nz = asymb2num[ nucname ]
		if len( nuca ) > 0:
			self.na = int( nuca )
		else:
			if ( nucname == "n" ):
				self.na = 1
				self.nz = 0
			elif ( nucname == "p" ):
				self.na = 1
				self.nz = 1
			elif ( nucname == "d" ):
				self.na = 2
				self.nz = 1
			elif ( nucname == "t" ):
				self.na = 3
				self.nz = 1
			else:
				print "unknown nucleus:", name
				self.na = 0
				self.nz = 0
		
		self.nn = self.na - self.nz
		return

          
