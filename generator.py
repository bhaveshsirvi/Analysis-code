#This file does the generator level generation according to the decay.dec file and saves the n-tuples.

#!/usr/bin/env python

import sys

import basf2 as b2

import modularAnalysis as ma

from variables import variables as vm  # shorthand for VariableManager

import variables.collections as vc

import variables.utils as vu

import vertex

import stdPi0s

import stdCharged as stdc

import variables

import generators as ge

main = b2.Path()

# Define number of events and experiment number
main.add_module('EventInfoSetter', evtNumList=[20000], expList=[0])

# Generate B0B0bar events
ge.add_evtgen_generator(
    path=main,
    finalstate='signal',
    signaldecfile=b2.find_file('decay.dec')
)

ma.fillParticleListFromMC("pi+", "", path=main)
ma.fillParticleListFromMC("pi-", "", path=main)
ma.fillParticleListFromMC("K+", "", path=main)
ma.fillParticleListFromMC("K-", "", path=main)


## D mesons
ma.reconstructMCDecay(decayString="D+ -> K- pi+ pi+" ,cut="", path=main)

ma.matchMCTruth("D+", path=main)


## B meson decay channels 
ma.reconstructMCDecay(decayString = "B-:ch1 -> D+ pi- pi-"  , cut = "", path=main)

ma.matchMCTruth("B-:ch1", path=main)

b_minus_had_vars1 =  vu.create_aliases_for_selected(
    list_of_variables=['M','deltaE','genMotherPDG']+vc.mc_kinematics + vc.mc_truth,
    decay_string='^B-:ch1 -> ^D+ ^pi- ^pi-',
    prefix=['B1','Dplus','pim11','pim12'])


ma.variablesToNtuple(
    decayString="B-:ch1",
    variables=b_minus_had_vars1 + ['M','Mbc','deltaE','isSignal'],
    filename="generated-output.root",
    treename="tree",
    path=main,
)

b2.process(main)

print(b2.statistics)