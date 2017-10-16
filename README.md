# cell-economy-model
Cellular economy models based on global optimization using GAMS

## overview

version: 1.0
date: 07-08-2017
author: Michael Jahn
affiliation: Science for Life Laboratory (KTH), Stockholm, Sweden
based on: R. Burnap, 2015, Molenaar et al., 2009
characteristics: protein economy model of a photoautotrophic cell
status: beta

## changelog
### 28-08-2017
- removed STA, substrate transport and assimilation and merged it with CBM, carbon metabolism. CBM includes CO2 fixation and other pathways that generate precursors 'pre'. 
- implemented light inhibition term for PSET rate v('PSET'). The term includes a substrate inhibition constant Ki similar to Km. v=c_enz*kcat*[S]/(Km+[S]*(1+[S]/Ki))

### 07-08-2017
- simulation of photosystem inhibition by adding DCMU. Implemented simply as reduction of kcat('PSET') by a certain factor, e.g. from 10 to 2 = -80%. DCMU is an inhibitor of PS2, the first of the two PS complexes, driving electron transfer to PQ and oxidation of water (2H2O -> 4H+ + 4e- + O2). It inhibits the reduction of PQ between the Q_A and Q_B reaction centers, preventing the primary step of photosynthesis. However, cyclic electron flow via PS1 may continue in real cells.

### 13-07-2017
- simulation of mixotrophic behavior. Added additional enzyme GLM representing glc uptake and catabolism. Added glc as a substrate. This model shows identical behavior to the autotrophic model when glc=0

### 07-07-2017
- changed cell volume balance: the cell surface area is now entirely based on the cytoplasmic membrane. The thylakoid membrane and its proteins were eliminated from the equation to prevent infinitey low beta (volume to surface ratio). Very low beta values would mean a highly convoluted cell shape. 
- thylakoid proteins (thyP) are now part of the intracellular protein  balance which is more realistic, i.e. they take up intracellular space.
- a perfectly spherical cell with V=1 has radius r=0.620 and surface area A=4.836. Ideally, the volume to surface ratio beta should therefore not exceed beta=V/A=1/4.836=0.206

### 02-03-2017
- changed membrane structure: Now split into thylakoid and cytoplasmic membrane. Thylakoid membrane houses the LHC and PSET, the photosynthetic machinery, while cytoplasmic membrane houses GLM
- membrane integrity conditions state that membrane protein concentration cannot exceed membrane lipid  concentration c(cpmP) <= c(cpm)

### 17-02-2017
- added shadow prices to the report. A shadow price is the increase in the objective function when a constraint is relaxed by 1 unit, in other words: shadow prices show which variable has the biggest effect on mu and represents the bottleneck

### 30-01-2017
- changed constraints from experimental proteomics data from ranges to a linear model, for each protein x, a(x)=m*mu + n. slope and intercept of this equation can be constrained by fixed values or optimized by the solver.

### 26-01-2017
- added constraints from experimental proteomics data, a(x) ranges for each protein group x

### 23-01-2017
- removed ATP/NADPH balance, implemented as metabolites
- substrate imports costs NADPH
- precursor generating enzyme changed to CBM, carbon metabolism,
- optionally included AAB, amino acid biosynthesis, niche proteins changed to MAI, regulation and maintenance
