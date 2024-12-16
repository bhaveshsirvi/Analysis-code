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

main = b2.Path()

#Replace the collection name if you want to run for different collections.
ma.inputMdstList(
    filelist=["/belle/collection/MC/MC15rd_charged_exp20-26_4S_v2"],
    path=main
)

ma.fillParticleList("pi+", "thetaInCDCAcceptance and nCDCHits > 20 and pionID > 0.3 and abs(dr) < 2 and abs(dz) < 4", path=main)
ma.fillParticleList("pi-", "thetaInCDCAcceptance and nCDCHits > 20  and pionID > 0.3 and abs(dr) < 2 and abs(dz) < 4", path=main)
ma.fillParticleList("K+", "thetaInCDCAcceptance and nCDCHits > 20  and kaonID > 0.6 and abs(dr) < 2 and abs(dz) < 4", path=main)
ma.fillParticleList("K-", "thetaInCDCAcceptance and nCDCHits > 20 and kaonID > 0.6 and abs(dr) < 2 and abs(dz) < 4", path=main)


## D mesons
ma.reconstructDecay(decayString="D+ -> K- pi+ pi+" ,cut="1.83965 < M < 1.89965", path=main)

ma.matchMCTruth("D+", path=main)


## B meson decay channels 
ma.reconstructDecay(decayString = "B-:ch1 -> D+ pi- pi-"  , cut = "Mbc > 5.24 and abs(deltaE) < 0.25", path=main)

ma.matchMCTruth("B-:ch1", path=main)

vertex.treeFit("B-:ch1",conf_level=0.0,updateAllDaughters=True,path=main)

#The following is required in case you want to suppress the continuum background. In my case this was not required as we used sweights to subtract the background.

ma.buildRestOfEvent(target_list_name="B-:ch1", path=main) 
cleanMask = (
    "cleanMask",
    "nCDCHits > 0 and useCMSFrame(p)<=3.2",
    "p >= 0.05 and useCMSFrame(p)<=3.2",
)
ma.appendROEMasks(list_name="B-:ch1", mask_tuples=[cleanMask], path=main)

ma.buildContinuumSuppression(list_name="B-:ch1", roe_mask="cleanMask", path=main)  

simpleCSVariables = [
    "R2",
    "thrustBm",
    "thrustOm",
    "cosTBTO",
    "cosTBz",
    "KSFWVariables(et)",
    "KSFWVariables(mm2)",
    "KSFWVariables(hso00)",
    "KSFWVariables(hso01)",
    "KSFWVariables(hso02)",
    "KSFWVariables(hso03)",
    "KSFWVariables(hso04)",
    "KSFWVariables(hso10)",
    "KSFWVariables(hso12)",
    "KSFWVariables(hso14)",
    "KSFWVariables(hso20)",
    "KSFWVariables(hso22)",
    "KSFWVariables(hso24)",
    "KSFWVariables(hoo0)",
    "KSFWVariables(hoo1)",
    "KSFWVariables(hoo2)",
    "KSFWVariables(hoo3)",
    "KSFWVariables(hoo4)",
    "CleoConeCS(1)",
    "CleoConeCS(2)",
    "CleoConeCS(3)",
    "CleoConeCS(4)",
    "CleoConeCS(5)",
    "CleoConeCS(6)",
    "CleoConeCS(7)",
    "CleoConeCS(8)",
    "CleoConeCS(9)",
]



b_minus_had_vars1 =  vu.create_aliases_for_selected(
    list_of_variables=vc.kinematics  + ['M' ,'genMotherPDG']+vc.pid ,
    decay_string='^B-:ch1 -> ^D+ ^pi- ^pi-',
    prefix=['B1','Dplus','pim11','pim12'])


ma.variablesToNtuple(
    decayString="B-:ch1",
    variables=["M","isContinuumEvent" , 'isSignalAcceptMissing','isSignal'] + b_minus_had_vars1+vc.kinematics +vc.mc_kinematics + vc.mc_truth + vc.deltae_mbc +  simpleCSVariables,
    filename="B1.root",
    treename="tree",
    path=main,
)

b2.process(main)

print(b2.statistics)