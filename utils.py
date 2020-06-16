import pickle
import numpy as np
from scipy.stats import binned_statistic as binStat

# Convenient named constants for indices for various parameters pertaining to an individual in the VIPER model
IDX_SUS = 0 # Susceptibility
IDX_VLD = 1 # Viral load
IDX_RLD = 2 # Recovery load
IDX_RST = 3 # Resistance (never changed - intrinsic to the individual)
IDX_X = 4   # X coordinate (only changed if the individual is travelling)
IDX_Y = 5   # Y coordinate (only changed if the individual is travelling)
IDX_STA = 6 # State
IDX_QRN = 7 # Is this individual quarantined or not
IDX_DPR = 8 # Disease progression rate
IDX_TI = 9  # Timestamp of infection (only changed at time of infection)
IDX_TQ = 10 # Timestamp of quarantine
IDX_TR = 11 # Timestamp of recovery/expiry (only changed at time of recovery/expiry)

nParams = 12

sqrt2 = np.sqrt(2)

# Convenient named constants for various states of an individual
# Being quarantined is not a medical state and so it is tracked
# separately (see IDX_QRN above) -- this allows VIPER to quarantine
ST_S = 0	# Susceptible
ST_E = 1	# Exposed
ST_I = 2	# Infectious
ST_R = 3	# Recovered
ST_X = 4	# Expired

nStates = 5

# Convenient named constants useful for reporting stats
SIDX_S = 0	# Number of susceptible but non-recovered individuals
SIDX_E = 1	# Number of exposed individuals
SIDX_I = 2	# Number of infectious individuals
SIDX_Q = 3	# Number of quarantined individuals
SIDX_R = 4	# Number of recovered individuals
SIDX_X = 5	# Number of expired individuals
SIDX_V = 6	# Average virulence of the viral strains in E and I populations
SIDX_EI = 7	# Number of infected individuals
SIDX_D = 8	# Number of individuals infected each day

nStats = 9

# This class provides an encapsulation for the VIPER model as well as methods
# to perform epidemiological simulations under lock-down, quarantining etc
class Population:
	# The base constructor -- almost never used directly
	def __init__( self, N, BIR, BCR, INC, BVL, VMR, QTH, BQP, XTH, BXP, BTR, BTD, BDPR, alphaInit, popIdx, ind, SUS, RST, X, Y ):
		self.N = N
		self.BIR = BIR
		self.BCR = BCR
		self.INC = INC
		self.BVL = BVL
		self.VMR = VMR
		self.QTH = QTH
		self.BQP = BQP
		self.XTH = XTH
		self.BXP = BXP
		self.BTR = BTR
		self.BTD = BTD
		self.BDPR = BDPR
		self.alphaInit = alphaInit
		
		self.popIdx = popIdx
		self.ind = ind
		self.SUS = SUS
		self.RST = RST
		self.X = X
		self.Y = Y
	
	# Constructor that initializes using scalar values alone
	# This provides a generic population with location, SUS, RST values
	# set uniformly randomly. To use India-specific initialization, use
	# the patch provided by the class Demographics (see below and setup.py)
	
	# N: the size of the population before the pandemic began
	# BIR: the (base) probability that a contact event will lead to a successful infection
	# BCR: the (base) contact radius
	# INC: the incubation period for the virus
	# BVL: the base viral load presented in individuals at the end of the exposed period
	# VMR: the rate at which the DPR mutates upon contact
	# QTH: the viral load over which an individual's chances of getting quarantined increase linearly
	# BQP: the (base) quarantining probability for a person with viral load at the QTH threshold
	# XTH: the viral load over which an individual's chances of getting expired increase linearly
	# BXP: the (base) expiry probability for a person with viral load at the XTH threshold
	# BTR: the (base) fraction of population that travels at any time instant
	# BTD: the (base) distance to which people travel
	# BDPR: the (base) disease progression rate for the initial strain of the virus
	# alphaInit: the initial fraction of population that is infected with the virus
	@classmethod
	def valInit( cls, N, BIR = 0.5, BCR = 0.25, INC = 3, BVL = 0.05, VMR = 0.0, QTH = 0.3, BQP = 0.0, XTH = 0.7, BXP = 0.0, BTR = 0.01, BTD = 1.0, BDPR = 0.1, alphaInit = 0.01 ):
		# Give the individuals unique identifiers
		popIdx = np.arange(N)
		
		# Initialize the individuals
		ind = np.zeros( (N, nParams) )

		# Initially no one is infected or quarantined
		ind[ :, IDX_VLD ] = 0
		ind[ :, IDX_RLD ] = 0
		ind[ :, IDX_STA ] = ST_S
		ind[ :, IDX_DPR ] = 0
		ind[ :, IDX_QRN ] = 0
		
		# People have susceptibilities uniformly in a range
		# The patch offered by the Demographics class allows this to be changed
		SUS = np.random.uniform( 0.01, 0.99, (N,) )
		ind[ :, IDX_SUS ] = SUS
		
		# People have resistances uniformly in a range
		# The patch offered by the Demographics class allows this to be changed
		RST = 0.1 * np.random.uniform( 0.01, 0.99, (N,) )
		ind[ :, IDX_RST ] = RST

		# Set their locations uniformly at random within the [0,1] x [0,1] square
		# The patch offered by the Demographics class allows this to be changed
		X = np.random.uniform( 0, 1, (N,) )
		ind[ :, IDX_X ] = X
		Y = np.random.uniform( 0, 1, (N,) )
		ind[ :, IDX_Y ] = Y
		
		# Set the timers to invalid values since nothing has happened yet
		ind[ :, IDX_TI ] = 0
		ind[ :, IDX_TQ ] = 0
		ind[ :, IDX_TR ] = 0
		
		return cls( N, BIR, BCR, INC, BVL, VMR, QTH, BQP, XTH, BXP, BTR, BTD, BDPR, alphaInit, popIdx, ind, SUS, RST, X, Y )
	
	# Constructor that initializes using data loaded from file
	@classmethod
	def fileInit( cls, filename ):
		# Get data from the file
		with open( filename, 'rb' ) as file:
			data = pickle.load( file )
		
		# Create an object out of this data and return it
		return cls( *data )
		
	# Save a copy of this object to a file
	def dump( self, filename ):
		data = ( self.N, self.BIR, self.BCR, self.INC, self.BVL, self.VMR, self.QTH, self.BQP, self.XTH, self.BXP, self.BTR, self.BTD, self.BDPR, self.alphaInit, self.popIdx, self.ind, self.SUS, self.RST, self.X, self.Y )
		with open( filename, 'wb' ) as file:
			pickle.dump( data, file )
			
	# Turn the clock back to before the pandemic started
	def reset( self ):
		# Reinitialize the individuals
		self.ind.fill(0)
		
		# Reset their susceptibilities to their original values (they do change for recovered and expired individuals)
		self.ind[ :, IDX_SUS ] = self.SUS
							 
		# Initially no one is infected or quarantined
		self.ind[ :, IDX_VLD ] = 0
		self.ind[ :, IDX_RLD ] = 0
		self.ind[ :, IDX_STA ] = ST_S
		self.ind[ :, IDX_DPR ] = 0
		self.ind[ :, IDX_QRN ] = 0
							 
		# Reset their resistances to their original value (although they should not have changed)
		self.ind[ :, IDX_RST ] = self.RST
				  
		# Reset their locations
		self.ind[ :, IDX_X ] = self.X
		self.ind[ :, IDX_Y ] = self.Y
				  
		# Reset the timers to a happier time before the pandemic began
		self.ind[ :, IDX_TI ] = 0
		self.ind[ :, IDX_TQ ] = 0
		self.ind[ :, IDX_TR ] = 0
		
	# T: the number of time steps for which to run the simulation
	# LKP: the lock-down policy
	def simulate( self, T, LKP, minimal = False, checkSanity = False ):
		numInfInit = int( self.alphaInit * self.N )
		
		# Randomly select the unfortunate souls who are going to get infected initially
		idxInit = np.random.permutation( self.N )[ : numInfInit ]
		
		# Infect them!
		self.ind[ idxInit, IDX_STA ] = ST_E			# These individuals are incubating right now
		self.ind[ idxInit, IDX_DPR ] = self.BDPR	# The virus has not mutated yet
		self.ind[ idxInit, IDX_TI ] = 1				# They got infected at t = 1
		
		# Store the intial stats
		stats = np.zeros( (nStats, T+1) )
		stats[ :, 0 ] = np.array( [ self.N - numInfInit, numInfInit, 0, 0, 0, 0, self.BDPR, numInfInit, numInfInit ] )
		
		for t in range(T):
			if checkSanity:
				self.doSanityChecks()
			# Track the progress of disease in infected + exposed individuals
			self.letDiseaseProgress( t+2 )
			# As much as permitted by the lock-down level at this time instant, allow individuals to travel
			self.letIndividualsTravel( LKP[t] )
			# Allow interactions and fresh infections and get a daily count
			daily = self.letIndividualsInfect( LKP[t], t+2 )
			# As dictated by the quarantine policy, catch and quarantine individuals
			self.applyQuarantinePolicy( t+2 )
			
			# Compute essential statistics
			idxEorI = (self.ind[ :, IDX_STA ].astype( int ) == ST_E) | (self.ind[ :, IDX_STA ].astype( int ) == ST_I)
			stats[ SIDX_EI, t+1 ] = np.count_nonzero( idxEorI )
			
			# Compute optional statistics
			if not minimal:
				stats[ SIDX_S, t+1 ] = np.count_nonzero( self.ind[ :, IDX_STA ].astype( int ) == ST_S )
				stats[ SIDX_E, t+1 ] = np.count_nonzero( self.ind[ :, IDX_STA ].astype( int ) == ST_E )
				stats[ SIDX_I, t+1 ] = np.count_nonzero( self.ind[ :, IDX_STA ].astype( int ) == ST_I )
				stats[ SIDX_Q, t+1 ] = np.count_nonzero( self.ind[ :, IDX_QRN ].astype( int ) == 1 )
				stats[ SIDX_R, t+1 ] = np.count_nonzero( self.ind[ :, IDX_STA ].astype( int ) == ST_R )
				stats[ SIDX_X, t+1 ] = np.count_nonzero( self.ind[ :, IDX_STA ].astype( int ) == ST_X )
				
				stats[ SIDX_D, t+1 ] = daily
				
				if len(idxEorI) == 0:
					stats[ SIDX_V, t+1 ] = 0
				else:
					stats[ SIDX_V, t+1 ] = np.mean( self.ind[ idxEorI, IDX_DPR ] )
		
		# Simulation over!!
		
		# Compute disease progression statistics and return final results of the simulation
		if minimal:
			return stats
		else:
			idxEverInfected = self.ind[ :, IDX_TI ].astype( int ) > 0
			tInfect = self.ind[ idxEverInfected, IDX_TI ]
			
			idxQuarantined = self.ind[ :, IDX_TQ ].astype( int ) > 0
			tQuarantine = self.ind[ idxQuarantined, IDX_TQ ] - self.ind[ idxQuarantined, IDX_TI ]
			
			idxRecovered = self.ind[ :, IDX_STA ].astype( int ) == ST_R
			tRecovery = self.ind[ idxRecovered, IDX_TR ] - self.ind[ idxRecovered, IDX_TI ]
			
			idxExpired = self.ind[ :, IDX_STA ].astype( int ) == ST_X
			tExpiry = self.ind[ idxExpired, IDX_TR ] - self.ind[ idxExpired, IDX_TI ]
			
			return ( stats, tInfect, tQuarantine, tRecovery, tExpiry )
			
	
	# Simulate incubation period, update viral and recovery loads, and process expiries and recoveries
	def letDiseaseProgress( self, t ):
		# Process recoveries of infectious individuals
		# TODO: recovery does not necessarily grant immunity: high levels of immunity may
		# only be available to those in whom the disease did progress sufficiently
		idxRecovered = ((self.ind[ :, IDX_STA ].astype( int ) == ST_I) & (self.ind[ :, IDX_VLD ] < self.BVL))
		if len(idxRecovered) > 0:
			self.ind[ idxRecovered, IDX_TR ] = t		# These people recovered at time t
			self.ind[ idxRecovered, IDX_SUS ] = 0		# They are now immune from the disease
			self.ind[ idxRecovered, IDX_VLD ] = 0		# Their remaining viral load vanishes in an instant
			self.ind[ idxRecovered, IDX_RLD ] = 0		# Their recovery load is made good instantly
			self.ind[ idxRecovered, IDX_STA ] = ST_R	# They are now recovered
			self.ind[ idxRecovered, IDX_DPR ] = 0		# They are free of the virus
			self.ind[ idxRecovered, IDX_QRN ] = 0		# They are released from any quarantine in which they may have been
		
		# Process expiries of infectious individuals
		# If expiry threshold is too high, assume no one is going to get expired!
		if self.XTH < 1 - 1e-6:
			idxExpiryRisk = self.popIdx[ (self.ind[ :, IDX_STA ].astype( int ) == ST_I) & (self.ind[ :, IDX_VLD ] > self.XTH) ]
			if len(idxExpiryRisk) > 0:
				expiryProb = self.BXP + ( 1 - self.BXP ) * (self.ind[ idxExpiryRisk, IDX_VLD ] - self.XTH)/(1 - self.XTH)
				tosses = np.random.uniform( 0, 1, ( len(expiryProb), ) )
				idxExpired = idxExpiryRisk[ tosses < expiryProb ]
				if len(idxExpired) > 0:
					self.ind[ idxExpired, IDX_TR ] = t		# These people expired at time t
					self.ind[ idxExpired, IDX_SUS ] = 0		# Either way, they are now immune from the disease
					self.ind[ idxExpired, IDX_VLD ] = 0		# Their remaining viral load vanishes in an instant
					self.ind[ idxExpired, IDX_RLD ] = 0		# Their recovery load does not matter any more
					self.ind[ idxExpired, IDX_STA ] = ST_X	# They are removed from the population
					self.ind[ idxExpired, IDX_DPR ] = 0		# Either way, they are now free of the virus
					self.ind[ idxExpired, IDX_QRN ] = 0		# They are released from any quarantine in which they may have been
		
		# Update viral and recovery loads
		idxInfected = self.popIdx[ self.ind[ :, IDX_STA ].astype( int ) == ST_I ]
		if len(idxInfected) > 0:
			# Some of the viral load is shed into recovery load, depending on the resistance of the person
			VLD_DeltaNeg = self.ind[ idxInfected, IDX_RST ] * self.ind[ idxInfected, IDX_VLD ]
			# The viral load also goes up, depending on the disease progression rate
			NormalLoad = 1 - self.ind[ idxInfected, IDX_VLD ] - self.ind[ idxInfected, IDX_RLD ]
			VLD_DeltaPos = NormalLoad * self.ind[ idxInfected, IDX_DPR ]
			
			self.ind[ idxInfected, IDX_RLD ] += VLD_DeltaNeg
			self.ind[ idxInfected, IDX_VLD ] += (VLD_DeltaPos - VLD_DeltaNeg)
		
		# Process incubation maturity cases
		idxExposed = self.popIdx[ self.ind[ :, IDX_STA ].astype( int ) == ST_E ]
		if len(idxExposed) > 0:
			tExposed = self.ind[ idxExposed, IDX_TI ].astype( int )
			idxMatured = idxExposed[ tExposed == t - self.INC ]
			if len(idxMatured) > 0:
				self.ind[ idxMatured, IDX_VLD ] = self.BVL	# These people now have a base viral load
				self.ind[ idxMatured, IDX_STA ] = ST_I	# They are now infectious
				
	# Let individuals travel according to the current lock-down level
	def letIndividualsTravel( self, level ):
		# Only non-quarantined and non-expired individuals can travel
		idxTravelAllowed = self.popIdx[ (self.ind[ :, IDX_QRN ].astype( int ) == 0) & (self.ind[ :, IDX_STA ].astype( int ) != ST_X) ]
		if len(idxTravelAllowed) > 0:
			nTravelAllowed = len(idxTravelAllowed)
			nTravellers = int( nTravelAllowed * self.BTR )
			if nTravellers > 0:
				# Allow a random set of people who are allowed to travel, to do so
				idxTravellers = idxTravelAllowed[ np.random.permutation( nTravelAllowed )[ :nTravellers ] ]
				
				# The travel distance is limited by the lock-down level
				tDist = self.BTD * np.exp( -level )
				
				# Dividing by sqrt2 since perturbations are made to two dimensions
				travelX = tDist / sqrt2 * np.random.uniform( -1, 1, (nTravellers,) )
				travelY = tDist / sqrt2 * np.random.uniform( -1, 1, (nTravellers,) )
				
				# Make sure individuals do not travel to another planet or something
				self.ind[ idxTravellers, IDX_X ] = ( self.ind[ idxTravellers, IDX_X ] + travelX ) % 1
				self.ind[ idxTravellers, IDX_Y ] = ( self.ind[ idxTravellers, IDX_Y ] + travelY ) % 1
		
	# Let individuals infect each other based on the current lock-down level
	def letIndividualsInfect( self, level, t ):
		# All non-quarantined non-expired individuals can participate in interactions
		# Individuals going through the incubation period do participate but neither infect nor can be infected
		idxParticipants = self.popIdx[ (self.ind[ :, IDX_QRN ].astype( int ) == 0) & (self.ind[ :, IDX_STA ].astype( int ) != ST_X) ]
		# Find a mask for the (non-quarantined) infectious participants
		maskInfectious = self.ind[ idxParticipants, IDX_STA ].astype( int ) == ST_I
		# Find a mask for non-immune susceptible/recovered individuals
		maskSusceptible = (self.ind[ idxParticipants, IDX_SUS ] > 0) & ( (self.ind[ idxParticipants, IDX_STA ].astype( int ) == ST_S) | (self.ind[ idxParticipants, IDX_STA ].astype( int ) == ST_R) )
		
		# If there is either no one infecting or else no one to infect, nothing left to do
		if ( np.count_nonzero( maskInfectious ) == 0 ) | ( np.count_nonzero( maskSusceptible ) == 0 ):
			return 0
		
		# Jiggle the coordinates pf the participants a bit to simulate random interactions
		# The jiggle dies down with lock-down level
		XParticipants = self.ind[ idxParticipants, IDX_X:IDX_Y + 1 ]
		jitter = self.BCR * np.exp( -level ) / sqrt2 * np.random.uniform( 0, 1, (2,) )
		XParticipants += jitter
		
		# Create interaction bins out of these jittered coordinates
		# The jitter added earlier allows individuals to fall into potentially different bins each time
		# TODO: this is slow at the moment, especially at high lock-down levels when number of bins
		# skyrockets since everyone is trapped inside their virtual houses. Speed this up by using a
		# reverse hash to iterate only over non-empty bins (there can be at most N such bins)
		rContact = self.BCR * np.exp( -level ) / sqrt2
		XParticipants = np.floor( XParticipants / rContact )
		nBinsPerCoord = np.max( XParticipants ) + 1
		BParticipants = ( XParticipants[:, 0] * nBinsPerCoord + XParticipants[:, 1] ).astype( int )
		
		# Find the bins for infectious and susceptible individuals
		BInfectious = BParticipants[ maskInfectious ]
		BSusceptible = BParticipants[ maskSusceptible ]
		
		# Find the infection probabilities in each bin
		bins = np.arange( np.max(BParticipants) + 2 )
		# print( "At time %d, we have jitter (%f, %f) and %d bins" % ( t, jitter[0], jitter[1], len(bins) ) )
		binParticipantCount = np.histogram( BParticipants, bins )[0]
		binParticipantCount[ binParticipantCount < 1 ] = 1 # Avoid a divide-by-zero error
		binInfectiousCount = np.histogram( BInfectious, bins )[0]
		binInfectionProb = binInfectiousCount / binParticipantCount * self.BIR
		
		# Find the mean DPR within each bin -- this will be used to decide the DPR for newly infected people
		idxInfectious = idxParticipants[maskInfectious]
		binDPR = binStat( BInfectious, self.ind[ idxInfectious, IDX_DPR ], statistic = "mean", bins = bins )[0]
		binDPR[ np.isnan( binDPR ) ] = 0 # Fix values for bins where there were no infectious people
		
		# Find the infection probabilities for each susceptible individual
		idxSusceptible = idxParticipants[maskSusceptible]
		infectionProb = binInfectionProb[BSusceptible] * self.ind[ idxSusceptible, IDX_SUS ]
		
		# Find new DPRs for the mutated viruses for the susceptible individuals
		# Currently this bit of code has been inactivated since we have set VMR = 0.0
		newDPR = binDPR[BSusceptible] * ( 1 + self.VMR * ( 2 * np.random.randint( 0, 2, ( len(BSusceptible), ) ) - 1 ) )
		newDPR[ newDPR < 0.001 ] = 0.001 # Make sure DPR never dips too low
		newDPR[ newDPR > 0.999 ] = 0.999 # Make sure DPR never goes too high
		
		# Get some random bits to decide who gets infected among the susceptible
		tosses = np.random.uniform( 0, 1, ( len(infectionProb), ) )
		maskInfected = tosses < infectionProb
		idxInfected = idxSusceptible[maskInfected]
		if len(idxInfected) > 0:
			self.ind[ idxInfected, IDX_STA ] = ST_E						# These individuals are incubating right now
			self.ind[ idxInfected, IDX_DPR ] = newDPR[maskInfected]		# These individuals get the mutated version
			self.ind[ idxInfected, IDX_TI ] = t							# They got infected at t = t
			
		return len(idxInfected)
		
	# Quarantine individuals based on quarantine policy
	def applyQuarantinePolicy( self, t ):
		# If quarantine threshold is too high, assume no one is going to get quarantined!
		if self.QTH > 1 - 1e-6:
			return
			
		# Only non-quarantined and non-expired individuals with viral loads above the quarantine threshold can be quarantined
		idxPotential = self.popIdx[ ( self.ind[ :, IDX_QRN ].astype( int ) == 0 ) & ( self.ind[ :, IDX_STA ].astype( int ) != ST_X ) & ( self.ind[ :, IDX_VLD ] > self.QTH ) ]
		if len(idxPotential) > 0:
			quarantineProb = self.BQP + ( 1 - self.BQP ) * ( self.ind[ idxPotential, IDX_VLD ] - self.QTH ) / ( 1 - self.QTH )
			tosses = np.random.uniform( 0, 1, ( len(quarantineProb), ) )
			idxQuarantined = idxPotential[ tosses < quarantineProb ]
			if len(idxQuarantined) > 0:
				self.ind[ idxQuarantined, IDX_QRN ] = 1		# Put them in quarantine
				self.ind[ idxQuarantined, IDX_TQ ] = t		# Note down the time of quarantining
				
	# Do sanity checks to make sure there are no anomalies
	def doSanityChecks( self ):
		# Make sure there are no negative values anywhere
		minVal = np.min( self.ind )
		assert minVal >= 0, "Negative values detected in the parameter values: %f" % minVal
	
		# Make sure everyone is in one of the five states S, E, I, R, X
		states = set( np.unique( self.ind[ :, IDX_STA ] ).astype( int ) )
		gold = set( np.arange( nStates ) )
		assert len( states - gold ) == 0, "Invalid states detected %" % "".join([str(i)+", " for i in states])
		
		# Make sure no one except those in I state have a non-zero viral and recovery load
		idxNonInfected = (self.ind[ :, IDX_STA ].astype( int ) != ST_I)
		assert max( self.ind[ idxNonInfected, IDX_VLD ] ) < 1e-6, "Non-infectious individuals with non-zero viral load detected"
		assert max( self.ind[ idxNonInfected, IDX_RLD ] ) < 1e-6, "Non-infectious individuals with non-zero recovery load detected"
		
		# Make sure no one other than those in state I are in quarantine
		# This may have to be modified if VIPER is extended to enable false positives in quarantining
		idxInvalidQuarantine = self.popIdx[ (self.ind[ :, IDX_QRN ].astype( int ) == 1) & (self.ind[ :, IDX_STA ].astype( int ) != ST_I) ]
		assert len(idxInvalidQuarantine) == 0, "Non-infectious individuals quarantined"
		
		# Make sure all individual locations are within bounds
		minX = np.min( self.ind[ :, IDX_X ] )
		assert minX >= 0, "Negative values detected in X coordinates %f" % minX
		minY = np.min( self.ind[ :, IDX_Y ] )
		assert minY >= 0, "Negative values detected in Y coordinates %f" % minY
		maxX = np.max( self.ind[ :, IDX_X ] )
		assert maxX <= 1, "Large values detected in X coordinates %f" % maxX
		maxY = np.max( self.ind[ :, IDX_Y ] )
		assert maxY <= 1, "Large values detected in Y coordinates %f" % maxY
		
# This class provides patches to the SUS, RST and location parameters for an object of the
# Population class. RST and SUS values are fitted to Indian statistics and locations are
# clustered into cities for a certain fraction of the population.
class Demographics:
	# The base constructor -- almost never used directly
	def __init__( self, N, urbanRatio, numCities, cityRadius, loc, isUrban, cityCenters, cityIdx, demoData ):
		self.N = N
		self.urbanRatio = urbanRatio
		self.numCities = numCities
		self.cityRadius = cityRadius
		self.loc = loc
		self.isUrban = isUrban
		self.cityCenters = cityCenters
		self.cityIdx = cityIdx
		self.demoData = demoData
	
	# Constructor that initializes using scalar values for various parameters
	# This provides a patch to the generic population generated using the class
	# Population (see above and setup.py) with a certain fraction living in cities
	# and RST and SUS values distributed to fit India statistics
	
	# N: the size of the population before the pandemic began
	# urbanRatio: what fraction of this population lives in cities?
	# numCities: how many cities do we have
	# cityRadius: the soft geographical extent of each city. ~99.7 of the population
	# of any city should live within 3 x cityRadius distance of the city center
	@classmethod
	def valInit( cls, N, urbanRatio = 0.34, numCities = 4, cityRadius = 0.1 ):
		( loc, isUrban, cityCenters, cityIdx ) = cls.getInitLocations( N, urbanRatio, numCities, cityRadius )
		demoData = cls.getDemographics(N)
		return cls( N, urbanRatio, numCities, cityRadius, loc, isUrban, cityCenters, cityIdx, demoData )
		
	# Constructor that initializes using data loaded from file
	@classmethod
	def fileInit( cls, filename ):
		# Get data from the file
		with open( filename, 'rb' ) as file:
			data = pickle.load( file )
		
		# Create an object out of this data and return it
		return cls( *data )
		
	# Save a copy of this object to a file
	def dump( self, filename ):
		data = ( self.N, self.urbanRatio, self.numCities, self.cityRadius, self.loc, self.isUrban, self.cityCenters, self.cityIdx, self.demoData )
		with open( filename, 'wb' ) as file:
			pickle.dump( data, file )
			
	# Distribute people into cities and non-urban areas randomly as requested
	@classmethod
	def getInitLocations( cls, N, urbanRatio, numCities, cityRadius ):
		loc = np.zeros( (N, 2) )
		isUrban = np.zeros( (N,) )
		
		# Find out who all is living in the cities
		numUrban = int( N * urbanRatio )
		permLoc = np.random.permutation( N )
		isUrban[ permLoc[ : numUrban ] ] = 1
		
		# Non-urban populations are distributed uniformly randomly
		loc[ permLoc[ numUrban : ], : ] = np.random.uniform( 0, 1, ( N - numUrban, 2 ) )
		
		# Urban populations are distributed in clusters around city centers
		# First, let us find the city centers
		cityCenters = np.random.uniform( 0, 1, ( numCities, 2 ) )
		# Next, assign each city dweller, a city uniformly randomly
		cityIdx = np.random.randint( 0, numCities, (numUrban,) )
		# Finally, set the location of each city dweller as a jitter around their respective city centers
		loc[ permLoc[ : numUrban ], : ] = cityCenters[ cityIdx, : ]
		loc[ permLoc[ : numUrban ], : ] += np.random.normal( 0, cityRadius / sqrt2, ( numUrban, 2 ) )
		
		# Make sure no one leaves the planet!
		loc[ loc < 0 ] = 0
		loc[ loc > 1 ] = 1
		
		return ( loc, isUrban, cityCenters, cityIdx )
	
	# Get RST and SUS values for the population that are derived from India statistics. For references
	# mentioned below e.g. (Ramakrishnan et al. 2019), please refer to the ESOP paper
	# https://arxiv.org/abs/2005.11257
	
	# Age groups are 0-20, 20-30, 30-40, 40-45, 45-50, 50-55, 55-60, 60-65, 65-70, 70-75, 75-80, 80+
	# There are 12 age groups that are indexed 0 through 11
	# This funny splitting of groups is due to the way in which various papers report their statistics
	# Whereas (Group 2003, Table 2) uses age intervals such as 20-30 years, 30-40 years, on the other
	# hand, (Ramakrishnan et al. 2019, Table 1) uses 45-55 years, 55-65 years as intervals instead
	
	# Genders are Female, Male -- these are indexed as Female = 0, Male = 1
	# The above convention is chosen to simplfy the susceptibility calculations
	# In all tables below, the first column is for the female gender and the second column for males
	@classmethod
	def getDemographics( cls, N ):
		# Taken from (Ministry Home Affairs India 2016, Detailed Tables)
		ageSplitbyGender = np.array([
			[ 0.363, 0.380 ],
			[ 0.205, 0.197 ],
			[ 0.152, 0.151 ],
			[ 0.061, 0.061 ],
			[ 0.054, 0.053 ],
			[ 0.043, 0.044 ],
			[ 0.037, 0.035 ],
			[ 0.031, 0.030 ],
			[ 0.022, 0.021 ],
			[ 0.015, 0.014 ],
			[ 0.009, 0.008 ], 
			[ 0.008, 0.006 ]
		])
		
		cumSplitbyGender = np.cumsum( ageSplitbyGender, axis = 0 )
		
		# According to the 2011 census, India's gender ratio is 51.5% males to 48.5% females
		genderSplit = 0.515
		
		# Taken from (Group 2003, Table 2)
		diabetesSplit = np.array([
			[ 0.000, 0.000 ],
			[ 0.000, 0.000 ],
			[ 0.125, 0.114 ],
			[ 0.227, 0.218 ],
			[ 0.227, 0.218 ],
			[ 0.326, 0.330 ],
			[ 0.326, 0.330 ],
			[ 0.346, 0.410 ],
			[ 0.346, 0.410 ],
			[ 0.330, 0.326 ],
			[ 0.330, 0.326 ],
			[ 0.177, 0.242 ]
		])
		
		# Taken from (Ramakrishnan et al. 2019, Table 1)
		hypertensionSplit = np.array([
			[ 0.062, 0.161 ],
			[ 0.140, 0.267 ],
			[ 0.140, 0.267 ],
			[ 0.140, 0.267 ],
			[ 0.346, 0.424 ],
			[ 0.346, 0.424 ],
			[ 0.454, 0.490 ],
			[ 0.454, 0.490 ],
			[ 0.514, 0.515 ],
			[ 0.514, 0.515 ],
			[ 0.513, 0.522 ],
			[ 0.513, 0.522 ],
		])
		
		demoData = np.zeros( ( N, 5 ) )
		# Sample random numbers in bulk for speed
		randomCoins = np.random.uniform( 0, 1, ( N, 4 ) )
		
		IDX_AGE = 0		# The age group of the person
		IDX_SEN = 1		# The seniority of the person
		IDX_GEN = 2		# The gender of the person
		IDX_DBT = 3		# The diabetic status of the person
		IDX_HTN = 4		# The hypertensive status of the person
		
		# TODO: The below for loop can surely be sped-up using vectorization/broadcasting techniques
		for i in range( N ):
			# Since our primary data only gives us conditional statistics e.g. gender ratio, age split by gender
			# and not joint statistics, we should not sample individual attributes jointly but rather sequentially
			
			thisGender = 0
			# First sample the gender of the person
			if randomCoins[ i, 0 ] <= genderSplit:
				thisGender = 1							# Its a boy!
			demoData[ i, IDX_GEN ] = thisGender
			
			# Next, sample the age-group of the person - use fresh randomness
			thisAgeGroup = np.digitize( randomCoins[ i, 1 ], cumSplitbyGender[ :, thisGender ], right = True )
			demoData[ i, IDX_AGE ] = thisAgeGroup
			
			# Next, decide the seniority of the person
			thisSeniority = 0
			# Persons above 55 years of age are classified as senior according to (Joshi 2020, Table 1)
			if thisAgeGroup > 5:
				thisSeniority = 1
			demoData[ i, IDX_SEN ] = thisSeniority
			
			# Next, decide the diabetic status of the person - use fresh randomness
			thisDiabetic = 0
			if randomCoins[ i, 2 ] <= diabetesSplit[ thisAgeGroup, thisGender ]:
				thisDiabetic = 1
			demoData[ i, IDX_DBT ] = thisDiabetic
			
			# Finally, decide the hypertensive status of the person - use fresh randomness
			thisHypertensive = 0
			if randomCoins[ i, 3 ] <= hypertensionSplit[ thisAgeGroup, thisGender ]:
				thisHypertensive = 1
			demoData[ i, IDX_HTN ] = thisHypertensive
			
		return demoData