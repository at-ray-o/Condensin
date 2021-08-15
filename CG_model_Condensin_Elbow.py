#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 17:10:38 2019

@author: Ryota
"""




from simtk import openmm
from simtk import unit
from simtk.openmm import app
from simtk.openmm import vec3
import numpy as np
# import matplotlib.pyplot as plt
from mdtraj.reporters import DCDReporter
import mdtraj as md
import os
import time

class ForceReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, True, False)

    def report(self, simulation, state):
        forces = state.getForces().value_in_unit(unit.kilojoules / unit.mole / unit.nanometer)
        i = 0.0
        for f in forces:
            # self._out.write('%g %g %g\n' % (f[0], f[1], f[2]))
            print('%g %g %g %g' % (i, f[0], f[1], f[2]))  # Gives value in kcal/mol/angstrom
            i += 1
        print('finished printing forces')

start = time.time()
np.random.seed()


input_file="coordRestart.pdb"
input_coords = md.load(input_file)
nparticles = input_coords.n_atoms
#Creat topology
topology=input_coords.topology.to_openmm()

indexH1=0
indexH2=82
indexHi=41
indexEl1=20
indexEl2=62


# No need to write angstrom as mdtraj import already converted the units
positions = unit.Quantity(np.array(input_coords.xyz.tolist()[0][0:nparticles]), unit.nanometer)

# Input parameters in openmm  units#

kBoltzmann = 1.380648 * 10 ** (-23) * 6.023 * 10 ** (23) / 1000 * unit.kilojoule_per_mole / unit.kelvin
viscosity = 8.90 * 10 ** -4 * 1000 * 6.023 * 10 ** 23 / (10 ** 9 * 10 ** 12) * \
            unit.amu / unit.nanometer / unit.picoseconds

kBT = 1 * unit.kilojoule_per_mole


# Basic units
sigmaBase = 1 * unit.nanometer  # 1 Angstrom

cutoffBase = (4*2.5) * sigmaBase
epsilonBase = 2.8 * kBT
temperature = kBT/unit.MOLAR_GAS_CONSTANT_R


# Set box size
edge = 10000000 * sigmaBase * 10

# Create Forces
# Load bonds
fBonds = open("bondList.dat", "r")
bondsDat = fBonds.read()
bondsDat = bondsDat.split('\n')
bondList = [[bondsDat[i].split()[0], bondsDat[i].split()[1], bondsDat[i].split()[2]] for i in range(len(bondsDat) - 1)]
fBonds.close()

fAnglesCondensin = open("CondensinAngles.dat", "r")
condAnglesDat = fAnglesCondensin.read()
condAnglesDat = condAnglesDat.split('\n')
condAnglesList = [[condAnglesDat[i].split()[0], condAnglesDat[i].split()[1], condAnglesDat[i].split()[2]] for i in range(len(condAnglesDat) - 1)]
fAnglesCondensin.close()


fAnglesDNA = open("DNAangles.dat", "r")
dnaAnglesDat = fAnglesDNA.read()
dnaAnglesDat = dnaAnglesDat.split('\n')
dnaAnglesList = [[dnaAnglesDat[i].split()[0], dnaAnglesDat[i].split()[1], dnaAnglesDat[i].split()[2]] for i in range(len(dnaAnglesDat) - 1)]
fAnglesDNA.close()

fparticles = open("particleList.dat", "r") ## Elbow
particlesDat = fparticles.read()
particlesDat = particlesDat.split('\n')
particleList = [particlesDat[i].split() for i in range(len(particlesDat)-1)]

for i in range(len(particleList)):
    particleList[i][2] = float(particleList[i][2])
fparticles.close()


#fEpsNonNat = open("epsilonMatrixCCMepigen.dat", "r")
#matrix = np.loadtxt(fEpsNonNat)
#epsilonMatrixNonNativeSopsc = unit.Quantity(matrix, unit.kilojoules_per_mole)


#fEpsNonNat.close()

massArray = [1] * nparticles

print("Loaded everything")

# Using data from Liu Reddy Thirumalai 2012
unitMass = 1 * unit.amu
massArray = [2*unitMass] * nparticles


unitLen = sigmaBase
unitEn = 1 * kBT
unitTimeHighFric = 1 * unit.picosecond
diffCoeff =  1 * unitLen**2/unitTimeHighFric
#high_fric = 0.10052277275 *kBT*unitTimeHighFric/unitLen**2
#high_fric_per_m = high_fric/unitMass

unitTimeLowFric = 1 * unit.picosecond * 0.4

tau_L = unitTimeLowFric.value_in_unit(unit.picosecond)
low_fric = 0.01 / tau_L
timeStepLD = 0.01 * tau_L

tau_H = 148 * unit.picoseconds

timeStepBD = 0.0001 * unitTimeHighFric



# Create an empty system object
system = openmm.System()

system.setDefaultPeriodicBoxVectors([edge, 0, 0], [0, edge, 0], [0, 0, edge])

nCondensin = 83
nDNA = nparticles - nCondensin


# =============================================================================
R0_fene = 1.5 * sigmaBase
#R0_fene_other = 3 * sigmaBase
k_fene = 20 * kBT/sigmaBase**2
feneBonds = openmm.CustomBondForce("-k/2.0*(R0)^2*log(1-((r-r0)/(R0))^2)")
feneBonds.addGlobalParameter('k', k_fene)
feneBonds.addGlobalParameter('R0', R0_fene)
feneBonds.addPerBondParameter('r0')

for i in range(len(bondList)):
    num = feneBonds.addBond(int(bondList[i][0]), int(bondList[i][1]),[float(bondList[i][2])])

system.addForce(feneBonds)


# =============================================================================

# Create an epsilon and sigma matrix

sigmaMat = unit.Quantity(np.zeros((nparticles, nparticles)), unit.nanometer)
epsilonAttMat = unit.Quantity(np.zeros((nparticles, nparticles)), unit.kilojoules_per_mole)
epsilonRepMat = unit.Quantity(np.zeros((nparticles, nparticles)), unit.kilojoules_per_mole)

epsilon_repel = 2 * unit.kilojoules_per_mole
epsilon_attract = 2 * unit.kilojoules_per_mole

for i in range(nparticles-1):
    for j in range(i+1,nparticles):
        sigmaMat[i][j]=0.5*(particleList[i][2]+particleList[j][2])*unit.nanometer
        sigmaMat[j][i]=sigmaMat[i][j]

for i in range(nparticles-1):
    for j in range(i+1,nparticles):
        #If either i or j is not hinge and DNA
        if i==41 and j>=nCondensin:
            epsilonAttMat[i][j]=epsilon_attract
            epsilonAttMat[j][i]=epsilon_attract
        else:
            epsilonRepMat[i][j]=epsilon_repel
            epsilonRepMat[j][i]=epsilon_repel


# =============================================================================


sigmaMatrixFlat = [item for sublist in sigmaMat for item in sublist]
epsilonMatrixAttFlat = [item for sublist in epsilonAttMat for item in sublist]
epsilonMatrixRepFlat = [item for sublist in epsilonRepMat for item in sublist]

repelBonds = openmm.CustomNonbondedForce("4*epsRep*(sigma/r)^12+4*epsAtt*((sigma/r)^12-(sigma/r)^6);\
                                         sigma=fsigma(type1,type2);\
                                         epsAtt=fepsilonAtt(type1,type2);\
                                         epsRep=fepsilonRep(type1,type2)")
repelBonds.addTabulatedFunction('fsigma',
                                         openmm.Discrete2DFunction(nparticles, nparticles, sigmaMatrixFlat))
repelBonds.addTabulatedFunction('fepsilonAtt', openmm.Discrete2DFunction(nparticles, nparticles,
                                                                                     epsilonMatrixAttFlat))
repelBonds.addTabulatedFunction('fepsilonRep', openmm.Discrete2DFunction(nparticles, nparticles,
                                                                                     epsilonMatrixRepFlat))
repelBonds.addPerParticleParameter('type')
repelBonds.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)
repelBonds.setCutoffDistance(cutoffBase)

for i in range(len(particleList)):
    repelBonds.addParticle([i])

system.addForce(repelBonds)
# =============================================================================


# =============================================================================

# =============================================================================




lp=150
lpEl=4
k_bend = lp * kBT/sigmaBase**2 # 150 is the literature value
k_bend_stiff = 100 * kBT/sigmaBase**2
k_bend_weak = lpEl * kBT/sigmaBase**2

angleForceCondensin = openmm.CustomAngleForce("kbend*(1+(cos(theta)))")
angleForceCondensin.addPerAngleParameter('kbend')

for i in range(len(condAnglesList)):
    if (int(condAnglesList[i][1])==28 or int(condAnglesList[i][1])==53):
        angleForceCondensin.addAngle(int(condAnglesList[i][0]), int(condAnglesList[i][1]), int(condAnglesList[i][2]),[k_bend_weak])
    else:
        angleForceCondensin.addAngle(int(condAnglesList[i][0]), int(condAnglesList[i][1]), int(condAnglesList[i][2]),[k_bend])
system.addForce(angleForceCondensin)

# =============================================================================

angleForceDNA = openmm.CustomAngleForce("kDNA*(theta-theta0)^2")
angleForceDNA.addGlobalParameter('kDNA',7.775*unitEn)
angleForceDNA.addGlobalParameter('theta0',np.pi)
for i in range(len(dnaAnglesList)):
    angleForceDNA.addAngle(int(dnaAnglesList[i][0]), int(dnaAnglesList[i][1]), int(dnaAnglesList[i][2]))

system.addForce(angleForceDNA)

# =============================================================================

# =============================================================================
theta_conf=0.833070358+np.pi/2      #O state
theta_conf2=4        #B state

#confAngle = openmm.CustomCompoundBondForce(3, "kangle*(cos(theta)-cos(angle(p1,p2,p3)))^2")
confAngle = openmm.CustomCompoundBondForce(3, "kangle*(thetaBend-angle(p1,p2,p3))^2")
#angle.addGlobalParameter('pi',np.pi)
confAngle.addGlobalParameter('kangle',k_bend_stiff)
confAngle.addGlobalParameter('thetaBend',theta_conf)
#angle.addPerBondParameter('d')
confAngle.addBond([82,0,1])
confAngle.addBond([0,82,81])
system.addForce(confAngle)
# =============================================================================


cm = openmm.CMMotionRemover()
system.addForce(cm)
print("Setup forces")

Marray = [1] * nparticles
Marray[indexH1]=4
Marray[indexH2]=4
Marray[indexHi]=4

## Elbow
Marray[indexEl1]=0.4
Marray[indexEl1-1]=0.4
Marray[indexEl1+1]=0.4
Marray[indexEl2]=0.4
Marray[indexEl2-1]=0.4
Marray[indexEl2+1]=0.4

## For the DNA part
for i in range(nCondensin,nparticles):
    Marray[i] = 3.4

for index in range(nparticles):
    #nonBondedTotalForce.addParticle([index])
    system.addParticle(Marray[index])


# Let's try to simulate it

# =============================================================================
from simtk.openmm import LangevinIntegrator, VerletIntegrator, BrownianIntegrator
from simtk.openmm.app import Simulation, PDBReporter, StateDataReporter, PDBFile
from sys import stdout

platform = openmm.Platform.getPlatformByName('CPU')

#integratorBD = BrownianIntegrator(temperature, high_fric, timeStepBD)
integratorLD = LangevinIntegrator(temperature, low_fric, timeStepLD)

 #import openmmtools as mmtools
 #Create custom integrator




integratorVariableD = openmm.CustomIntegrator(timeStepBD);
integratorVariableD.addUpdateContextState();
integratorVariableD.addPerDofVariable("gaussPrev",0)
integratorVariableD.addPerDofVariable("gaussNow",0)


integratorVariableD.addGlobalVariable("boltzmann",kBT)
integratorVariableD.addGlobalVariable("dt",timeStepBD)
integratorVariableD.addComputePerDof("gaussNow","gaussian")
integratorVariableD.addComputePerDof("x", "x+preFac1*f+preFac2*(gaussPrev+gaussNow); preFac1= dt/m; preFac2=sqrt(boltzmann*dt/(2*m))")
integratorVariableD.addComputePerDof("gaussPrev","gaussNow")



simulation = Simulation(topology, system, integratorLD)
#simulation = Simulation(topology, system, integratorBD)
#print("check")
#input('a')


simulation.context.setPositions(positions)
simulation.minimizeEnergy()


#simulation.minimizeEnergy(tolerance=0.001*kBT)
nsteps = 150000000
freq = 10000

simulation.reporters.append(DCDReporter('output.dcd', freq))
#simulation.reporters.append(DCDReporter('output' + '_lp=' + str(lp) + '.dcd', freq))
simulation.reporters.append(StateDataReporter(stdout,freq,potentialEnergy=True,\
                                              totalEnergy=True,step=True,time=True,separator='   '))
#simulation.reporters.append(ForceReporter('forces.txt', 1))
with open('topology.pdb', 'w') as f:
    PDBFile.writeFile(topology, positions, f)



#simulation.step(10000)
#simulation.context.setParameter("theta",theta_conf)
simulation.context.setParameter("thetaBend",theta_conf)
simulation.step(10*freq)
#simulation.step(confSteps)
simulation.context.setParameter("thetaBend",theta_conf2)
simulation.step(10*freq)
simulation.context.setParameter("thetaBend",theta_conf)
simulation.step(10*freq)
#simulation.step(nsteps)

state = simulation.context.getState(getPositions=True)
finalpos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)

end = time.time()
print("Time Elapsed : ", end - start)
# =============================================================================
