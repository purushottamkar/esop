from utils import (Population, Demographics)

# The size of the population
N = 20000

# Instantiate a population with defaults for all VIPER parameters
# This population has naive distributions for initial locations of
# individuals as well as demographic distributions e.g. RST, SUS.
# All VIPER paramaters are customizable e.g. BIR, BCR, INC etc.
# Please see the valInit method in the class Population in utils.py
# in case you wish to modify these parameters
pop = Population.valInit( N = N )
pop.dump( "pop_generic" )

# To be able to run experiments with more "realistic" populations, create
# a patch that assigns population demographics that fit to Indian statistics
# as well as assigns initial locations for individuals that are clustered
# into cities. Urban population ratio, num cities, city radius are customizable.
# Please see the valInit method in the class Demographics in utils.py
# in case you wish to modify these parameters
demloc = Demographics.valInit( N = N )
demloc.dump( "patch_India_fit" )