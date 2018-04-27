$ontext
PHOTOAUTOTROPHIC CELL MODEL
version: 1.0
subversion: 'random sampling of kinetic parameters'
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

kcat(enz)  "kcat of enzymes = turnover number"
Km(enz)    "Km of enzymes = substrate affinity constant"
hc(enz)    "Hc of enzymes = Hill coefficient for cooperativity"


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
glc = 0;

* Initial values for cofactors (prevent division by zero, do not affect result)
c.l('atp') = 0.1;
c.l('nadph') = 0.1;

PARAMETERS    
kcat_init(enz)    "Initial values for kcat"
          /LHC     160
           PSET     31
           CBM       5
           LPB       7
           RIB      10
           GLM       5/

Km_init(enz)      "Initial values for Km"
          /LHC      58
           PSET    108
           CBM     229
           LPB      13
           RIB     128
           GLM      50/
           
hc_init(enz)      "Initial values for Hill coefficient"
          /LHC     2.0043
           PSET    1.3989
           CBM     2.2477
           LPB     0.6744
           RIB     0.8659
           GLM     2.0000/
;

* ------------ OPTIONAL EXPERIMENTAL CONSTRAINTS -----------------------
*
* constraints for intercepts, representing maintenance expression at µ=0
intercept.lo(pro) = linregTable(pro, 'interceptConst');
intercept.up(pro) = linregTable(pro, 'interceptConst');

* constraints for slope, ONLY for selected proteins (e.g. to fix maintenance protein fraction)
slope.lo(pro)$conP(pro) = linregTable(pro, 'slopeConst')$conP(pro);
slope.up(pro)$conP(pro) = linregTable(pro, 'slopeConst')$conP(pro);


* Iterate over a FOR loop that sets random values for kcat and Km
* Iterate over a FOR loop that tests different conditions
SET j                   "iteration driver" / 1*1000 /;
SET i                   "iteration driver" / 1*10 /;
PARAMETER Res_mu        "Sum of Residuals between model-optimized and experimental growth rate";
PARAMETER Res_al        "Sum of Residuals between model-optimized and experimental protein fractions";
PARAMETER Res_mu_optim  "current optimal Sum of Residuals"; Res_mu_optim = 0.0209;
PARAMETER Res_al_optim  "current optimal Sum of Residuals"; Res_al_optim = 0.7081;


* report parameters of the model. '.l' is the current level, '.M' for 
* marginals is the shadow price = inrease in mu when constraint is relaxed by 1
PARAMETER report(*,*,*) "process level report";


* ------------ ITERATE OVER KINETIC PARAMETERS -------------------------
*
LOOP (j,
    hv = 50;
    Res_mu = 0;
    Res_al = 0;
    execseed = 1e8*(frac(jnow));
    LOOP(enz, kcat(enz) = uniformint(round(kcat_init(enz)*0.8), round(kcat_init(enz)*1.2)););
    LOOP(enz, Km(enz) = uniformint(round(Km_init(enz)*0.8), round(Km_init(enz)*1.2)););
    LOOP(enz, hc(enz) = uniform(hc_init(enz)*0.8, hc_init(enz)*1.2););
    
* ------------ ITERATE OVER CONDITION ----------------------------------
    
    LOOP (i,
        SOLVE CELL USING NLP MAXIMIZING logmu;
*       determine sum of residuals as a metric of divergence, either for alpha or mu
*       alpha residuals are calculated using experimentally determined proteome fractions
*       mu residuals are calculated using a fitted non-linear growth model
        Res_al = Res_al + Sum(enz, abs(a.l(enz)-(intercept.l(enz) + exp(logmu.l)*linregTable(enz, 'slopeConst'))));
        Res_mu = Res_mu + abs(exp(logmu.l)-(0.2154 * (1-hv/216.0161) * hv/(hv+(37.3387*(1-hv/216.0161)))));
        report('model','sub',i) = sub;
        report('model','hv',i) = hv;
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
        report('  Res','mu',i) = Res_mu;
        report('  Res','alpha',i) = Res_al;
        hv = hv/1.5;
    );
    if (Res_mu < Res_mu_optim AND Res_al < Res_al_optim, 
        DISPLAY report;
*        Res_mu_optim = Res_mu;
*        Res_al_optim = Res_al;
    );
);

