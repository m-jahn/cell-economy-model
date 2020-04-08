# cell-economy-model
Cellular economy models based on global optimization using GAMS

-----

**Note: This version is deprecated. please check out the [new implementation based on python](https://github.com/m-jahn/cell-economy-models)**

-----

This cellular economy model is a 'coarse-grain' model originally conceived by Molenaar et al., 2009, and R. Burnap, 2015. The purpose of the model is not to reflect cellular behavior in its entire complexity, but rather to reduce complexity to an amount that still allows drawing significant conclusions while keeping the number of components and reactions as small as possible. Following this assumption, a cellular economy model may not contain all known metabolic pathways, enzymes or known kinetic parameters thereof. It rather bundles fundamental cellular processes in 'super-enzymes' (sectors). These are single catalytic units that serve as proxies for several similar or related pathways. However, the output from a coarse-grain model is very useful as it illustrates metabolic tradeoffs under different conditions, without getting lost in details. The model is a non-linear mixed-integer optimization problem formulated in GAMS, a programming language for optimization problems. Although different aspects of cellular behaviour can be probed using this model, it is focused on growth as the objective function.

## overview

- type: cellular economy model of a photoautotrophic cell
- version: 1.0
- date: 2018-04-06
- author: Michael Jahn
- affiliation: Science for Life Laboratory (KTH), Stockholm, Sweden
- original concept and models by:
-- Burnap, R.L. (2015). Systems and Photosystems: Cellular Limits of Autotrophic Productivity in Cyanobacteria. Front. Bioeng. Biotechnol. 3, 1.
-- Molenaar, D., Van Berlo, R., De Ridder, D., and Teusink, B. (2009). Shifts in growth strategies reflect tradeoffs in cellular economics. Mol. Syst. Biol. 5, 323.


## changelog

### 2018-05-20
- simulation of photoinhibition by including the according term in the rate equation for the light harvesting sector LHC.
- new branch of the model including the possibility to fix dynamic fractions of each protein as un-utilized, that can be seen as a protein 'reserve'. Thsi behavior is known for ribosomes, whose mass expands with growth rate but is not directly proportional to it. Even at µ=0, cells keep a large fraction of the proteome as inactive/un-utilized ribosomes.

### 2017-08-28
- removed STA, substrate transport and assimilation and merged it with CBM, carbon metabolism. CBM includes CO2 fixation and other pathways that generate precursors 'pre'. 
- implemented light inhibition term for PSET rate v('PSET'). The term includes a substrate inhibition constant Ki similar to Km. v=c_enz*kcat*[S]/(Km+[S]*(1+[S]/Ki))

### 2017-08-07
- simulation of photosystem inhibition by adding DCMU. Implemented simply as reduction of kcat('PSET') by a certain factor, e.g. from 10 to 2 = -80%. DCMU is an inhibitor of PS2, the first of the two PS complexes, driving electron transfer to PQ and oxidation of water (2H2O -> 4H+ + 4e- + O2). It inhibits the reduction of PQ between the Q_A and Q_B reaction centers, preventing the primary step of photosynthesis. However, cyclic electron flow via PS1 may continue in real cells.

### 2017-07-13
- simulation of mixotrophic behavior. Added additional enzyme GLM representing glc uptake and catabolism. Added glc as a substrate. This model shows identical behavior to the autotrophic model when glc=0

### 2017-07-07
- changed cell volume balance: the cell surface area is now entirely based on the cytoplasmic membrane. The thylakoid membrane and its proteins were eliminated from the equation to prevent infinitey low beta (volume to surface ratio). Very low beta values would mean a highly convoluted cell shape. 
- thylakoid proteins (thyP) are now part of the intracellular protein  balance which is more realistic, i.e. they take up intracellular space.

### 2017-03-02
- changed membrane structure: Now split into thylakoid and cytoplasmic membrane. Thylakoid membrane houses the LHC and PSET, the photosynthetic machinery, while cytoplasmic membrane houses GLM
- membrane integrity conditions state that membrane protein concentration cannot exceed membrane lipid  concentration c(cpmP) <= c(cpm)

### 2017-02-17
- added shadow prices to the report. A shadow price is the increase in the objective function when a constraint is relaxed by 1 unit, in other words: shadow prices show which variable has the biggest effect on mu and represents the bottleneck

### 2017-01-30
- changed constraints from experimental proteomics data from ranges to a linear model, for each protein x, a(x)=m*mu + n. slope and intercept of this equation can be constrained by fixed values or optimized by the solver.

### 2017-01-26
- added constraints from experimental proteomics data, a(x) ranges for each protein group x

### 2017-01-23
- removed ATP/NADPH balance, implemented as metabolites
- substrate imports costs NADPH
- precursor generating enzyme changed to CBM, carbon metabolism,
- optionally included AAB, amino acid biosynthesis, niche proteins changed to MAI, regulation and maintenance
