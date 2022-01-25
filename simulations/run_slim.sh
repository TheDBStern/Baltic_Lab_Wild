## Set parameter values ##
fitfunc=1 #1 = expotential (population genetic epistasis unless alpha=0),
          #3 = QT directional epistasis, 4 = truncating, 5 = shifted optimum

#multiplicative / exponential model parameters
a=8
s=0.05
b=-0.395
#quantitative fitness function parameters
mu=0.435
std=0.0175

slim -d fmin=0 -d fmax=1 -d npops=10 -d nloci=121 -d popsize=1750 -d scaleT0=0 -d scales=0 -d seed=$RANDOM \
             -d fitnessFunction=1 \
             -d a=$a -d s=$s -d b=$b -d mu=$mu -d std=$std epistasis_simulations.slim | \
             tail -n +14 > allele_frequencies.csv
