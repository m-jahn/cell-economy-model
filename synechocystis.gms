$ontext
PHOTOAUTOTROPHIC CELL MODEL
version: 1.0
subversion: 'light and CO2 limitation'
date: 2018-04-06
author: Michael Jahn
affiliation: Science for Life Laboratory (KTH), Stockholm, Sweden
based on: R. Burnap, 2015, Molenaar et al., 2009
characteristics: protein economy model of a photoautotrophic cell
$offtext
$offlisting
$offsymxref offsymlist

OPTION
    limrow = 0,
    limcol = 0,
    solprint = off,
    sysout = off,
    decimals = 4,
    resLim=1000;

SETS
cmp "all cell components" / LHC    "light harvesting complex, photosystems"
                            PSET   "photosynthesis electron transport enzyme"
                            CBM    "Carbon metabolism"
                            LPB    "lipid biosynthesis enzyme"
                            RIB    "ribosomal proteins"
                            MAI    "regulation/maintenance/environmental homeostasis, other proteins"
                            GLM    "glucose import and catabolism"
                            hvi    "absorbed light (hv), intracellular"
                            atp    "ATP"
                            nadph  "NADPH, a universal reductant"
                            pre    "universal metabolic precursor"
                            lip    "lipids and other cell membrane components"
                            cpm    "cytoplasmic membrane"
                            thy    "thylakoid membrane"/


  pro(cmp)  "all proteins" / LHC, PSET, CBM, LPB, RIB, MAI, GLM /
  enz(pro)  "enzymes"      / LHC, PSET, CBM, LPB, RIB, GLM /
  mem(cmp)  "all membrane located components" / LHC, PSET, GLM, cpm, thy /
  cpmP(mem) "cytoplasmic membrane located proteins"   / GLM /
  thyP(mem) "thylakoid membrane located proteins"     / LHC, PSET /
  intP(pro) "intracellular proteins"    / LHC, PSET, CBM, LPB, RIB, MAI /
  conP(pro) "constrained protein mass fraction" / MAI, GLM /
  met(cmp)  "metabolites"               / hvi, atp, nadph, pre, lip /
  linreg    "linear regression parameters"   / interceptConst, slopeConst /


PARAMETERS

kcat(enz)  "kcat of enzymes = turnover number; PSET=55"
          /LHC     160
           PSET     31
           CBM       5
           LPB       7
           RIB      10
           GLM       5/

Km(enz)    "Km of enzymes = substrate affinity constant"
          /LHC      58
           PSET    108
           CBM     229
           LPB      13
           RIB     128
           GLM      50/

hc(enz)    "Hc of enzymes = Hill coefficient for cooperativity"
          /LHC     2.0043
           PSET    1.3989
           CBM     2.2477
           LPB     0.6744
           RIB     0.8659
           GLM     2.0000/


KmATP      "affinity constant of CBM for ATP"        /CBM      1/
KmNADPH    "affinity constant of CBM for NADPH"      /CBM      1/


sA(mem) "specific surface area of membrane located components"
          /LHC       1
           PSET      1
           GLM       1
           cpm       1
           thy       1/


hv      "extracellular light irradiation"
sub     "extracellular substrate concentration (CO2/HCO3-)"
glc     "extracellular glucose concentration"
Ki      "light inhibition constant for photosystems"


table linregTable(pro, linreg)  "constraints of linear regression for mu vs alpha under light limitation"
        interceptConst     slopeConst
LHC     0.3155446          -1.7828190
PSET    0.1488815          -0.2966529
CBM     0.0928786           0.7384860
LPB     0.02396478          0.01308973
RIB     0.1002482           0.9634276
MAI     0.2875577           0.3078048
GLM     0.03092466          0.05666372


table stoich(met, enz)  "reaction stoichiometry matrix"
      LHC  PSET    CBM    LPB   RIB   GLM
hvi     1    -1      0      0     0     0
atp     0     1     -1      0     0     1
nadph   0     1     -1      0     0     1
pre     0     0      1     -1    -1     0
lip     0     0      0      1     0     0;


POSITIVE VARIABLES

a(pro)           "fraction of ribosomes engaged in synthesis of protein x"
beta             "volume-to-surface ratio, increases with sphericity of a cell"
c(cmp)           "concentration of component"
v(enz)           "catalytic rate of enzyme"
intercept(pro)   "intercept of linear relation between mu and alpha"


FREE VARIABLE

* Since the objective variable has to be a "free" variable i.e. defined
* on the <-INF,+INF> interval, we transform mu to the logarithmic scale
* So, exp(logmu) = mu and maximizing logmu is equivalent to maximizing mu.
* The slope of linear regression constraints for alpha can also be 
* positive or negative
slope(pro)   "slope of linear relation between mu and alpha"
logmu        "natural logarithm of the specific growth rate";


EQUATIONS

volume       "cell volume is determined by beta and the cytoplasmic membrane surface. Constant."
alphaSum     "fractions of ribosomes engaged in synthesis of all proteins sum up to 1"
Pbal(pro)    "mass balance for proteins. Ribosomal synthesis - consumption by cell growth = 0"
Mbal(met)    "mass balance for metabolites. Production by enzyme reactions (rate * stoechiometry) - consumption by cell growth = 0"
cat_LHC      "catalytic rate of light harvesting and photosystems pool"
cat_PSET     "catalytic rate of NADPH and ATP generation by electron transport chain"
cat_CBM      "catalytic rate of carbon metabolism pool"
cat_LPB      "catalytic rate of lipid/membrane biosynthesis enzyme pool"
cat_RIB      "catalytic rate of ribosome pool"
cat_GLM      "catalytic rate of glucose import and catabolism"
maxP         "maximal intracellular protein concentration"
*maxM         "maximal total metabolite concentration"
LipBal       "lipid balance: thylakoid and cytoplasmic membrane equals total lipids"
int_cpm      "membrane integrity condition for cytoplasmic membrane"
int_thy      "membrane integrity condition for thylakoid membrane"
alph(pro)    "fraction of ribosomes engaged in protein synthesis"
;


volume..     beta*(sA('cpm')*c('cpm')+Sum(cpmP, sA(cpmP)*c(cpmP))) =E= 1;
alphaSum..   Sum(pro, a(pro)) =E= 1;
Pbal(pro)..  a(pro)*v('RIB') - exp(logmu)*c(pro) =E= 0;
Mbal(met)..  Sum(enz, stoich(met, enz)*v(enz)) - exp(logmu)*c(met) =E= 0;
*cat_PSET..   v('PSET') =E= kcat('PSET')*c('PSET')*c('hvi')/(c('hvi')*(1+c('hvi')/Ki) + Km('PSET'));
cat_LHC..    v('LHC') =E= kcat('LHC')*c('LHC')*rPower(hv, hc('LHC'))/(rPower(Km('LHC'), hc('LHC')) + rPower(hv, hc('LHC')));
cat_PSET..   v('PSET') =E= kcat('PSET')*c('PSET')*rPower(c('hvi'), hc('PSET')) /(rPower(c('hvi'), hc('PSET')) + rPower(Km('PSET'), hc('PSET')));
cat_CBM..    v('CBM') =E= kcat('CBM')*c('CBM')*c('nadph')*rPower(sub, hc('CBM'))*c('atp')/(c('nadph')*rPower(sub, hc('CBM'))*c('atp') + KmNADPH('CBM')*c('atp') + KmATP('CBM')*c('nadph') + KmATP('CBM')*rPower(sub, hc('CBM')) + rPower(Km('CBM'), hc('CBM'))*c('nadph'));
cat_LPB..    v('LPB') =E= kcat('LPB')*c('LPB')*rPower(c('pre'), hc('LPB'))/(rPower(Km('LPB'), hc('LPB')) + rPower(c('pre'), hc('LPB')));
cat_RIB..    v('RIB') =E= kcat('RIB')*c('RIB')*rPower(c('pre'), hc('RIB'))/(rPower(Km('RIB'), hc('RIB')) + rPower(c('pre'), hc('RIB')));
cat_GLM..    v('GLM') =E= kcat('GLM')*c('GLM')*rPower(glc, hc('GLM'))/(rPower(Km('GLM'), hc('GLM')) + rPower(glc, hc('GLM')));
maxP..       Sum(intP, c(intP)) =L= 1;
*maxM..       Sum(met, c(met)) =L= 1;
LipBal..     c('cpm') + c('thy') =E= c('lip');
int_cpm..    Sum(cpmP, c(cpmP))+0.1 =L= c('cpm');
int_thy..    Sum(thyP, c(thyP))+0.1 =L= c('thy');
alph(pro)..  a(pro) =E= intercept(pro)+(slope(pro)*exp(logmu));


MODEL CELL "the autotrophic cell model" /ALL/;


* ------------ STARTING CONCENTRATIONS ---------------------------------
*
hv  = 100;
sub = 100;
glc =   0;
*Ki  = 20;


* Initial values for cofactors (prevent division by zero, do not affect result)
c.l('atp') = 0.1;
c.l('nadph') = 0.1;


* ------------ OPTIONAL EXPERIMENTAL CONSTRAINTS -----------------------
*
* constraints for intercepts, representing maintenance expression at µ=0
intercept.lo(pro) = linregTable(pro, 'interceptConst');
intercept.up(pro) = linregTable(pro, 'interceptConst');

* constraints for slope, ONLY for selected proteins (e.g. to fix maintenance protein fraction)
slope.lo(pro)$conP(pro) = linregTable(pro, 'slopeConst')$conP(pro);
slope.up(pro)$conP(pro) = linregTable(pro, 'slopeConst')$conP(pro);


* ------------ ITERATIVE SOLVING OF MODEL ------------------------------
*
* Iterate over a FOR loop that tests different light or CO2 conditions
SET i                   "iteration driver" / 1*15 /;

* avoid zero values, that are not exported by gams
a.lo(pro) = 0.0001;

* report parameters of the model. '.l' is the current level, '.M' for 
* marginals is the shadow price = inrease in mu when constraint is relaxed by 1
PARAMETER report(*,*,*) "process level report" ;


LOOP (i,
    SOLVE CELL USING NLP MAXIMIZING logmu;
    report('model','sub',i) = sub;
    report('model','hv',i) = hv;
    report('model','glc',i) = glc;
    report('model','mu',i) = exp(logmu.l);
    report('model','beta',i) = beta.l;
    report('alpha',pro,i) = a.l(pro);
    report(' conc',cmp,i) = c.l(cmp);
    report(' rate',enz,i) = v.l(enz);
    report('inter',pro,i) = intercept.l(pro);
    report('slope',pro,i) = slope.l(pro);
    report(' shad',pro,i) = alph.M(pro);
    report(' kcat',enz,i) = kcat(enz);
    report('   Km',enz,i) = Km(enz);
    report('   hc',enz,i) = hc(enz);
    hv = hv/1.5;
);

DISPLAY report;


* ------- OPTIONAL SOLVING OF ALL COMBINATIONS OF TWO VARIABLES --------
*
* Iterate over two FOR loop that test different light and CO2 conditions
* For plotting 3D it's easier to have a linear instead of log distribution
*SET j                   "iteration driver" / 1*20 /;
*SET k                   "iteration driver" / 1*20 /;
*hv  = 100;
*sub = 100;
*
*
* this report needs to have 4 slots
*PARAMETER report2(*,*,*,*) "process level report" ;
*
*
*LOOP (j,
*    LOOP (k,
*        SOLVE CELL USING NLP MAXIMIZING logmu;
*        report2('model','sub',j,k) = sub;
*        report2('model','hv',j,k) = hv;
*        report2('model','mu',j,k) = exp(logmu.l);
*        report2('model','beta',j,k) = beta.l;
*        report2('alpha',pro,j,k) = a.l(pro);
*        report2(' conc',cmp,j,k) = c.l(cmp);
*        report2(' rate',enz,j,k) = v.l(enz);
*        hv = hv-5;
*    );
*    sub = sub-5;
*    hv  = 100;
*);
*
*DISPLAY report2;
