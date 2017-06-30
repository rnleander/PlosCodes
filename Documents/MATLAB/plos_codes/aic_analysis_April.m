
paramG1=[2 3 3 4 5];

%G1parts
g1_obs=266;
g1_1=-796.154;
g1_l=-762.659;
g1_2=-756.269;
g1_3=-756.016;
g1_e=-770.293;
g1_L=[g1_1 g1_e g1_l g1_2 g1_3];

[g1_aic_parts, g1_relp_parts]=akaikec(g1_L,paramG1,g1_obs)

paramG2=[2 3 3 5];

%G2parts
g2_obs=266;
g2_1=-596.015;
g2_l=-595.267;
g2_2=-595.267;
g2_3=-595.268;
g2_e=-597.567;
g2_L=[g2_1 g2_e g2_l g1_3];

[g2_aic_parts, g2_relp_parts]=akaikec(g2_L,paramG2,g2_obs)

paramIMT=[2 3 3 4 5];


%IMTparts
imt_obs=266;
imt_1=-892.381;
imt_l=-849.484;
imt_2=-848.715;
imt_3=-848.179;
imt_e=-852.060;
imt_L=[imt_1 imt_e imt_l imt_2 imt_3];

[imt_aic_parts, imt_relp_parts]=akaikec(imt_L,paramIMT,imt_obs)

param=[2 3 4 6];

% DMSO
dmso_obs = 379;
dmso1 = -847.544;
dmso2 = -757.852;
dmso3 = -757.859;
dmsoe = -775.211;
dmso_L=[dmso1 dmsoe dmso2 dmso3];

[dmso_aic,dmso_relp] = akaikec(dmso_L,param,dmso_obs);

% erlotinib
erlot_obs = 198;
erlot1 = -609.710;
erlot2 = -542.296;
erlot3 = -541.484;
erlote = -552.853;

erlot_L=[erlot1 erlote erlot2 erlot3];

[erlot_aic,erlot_relp] = akaikec(erlot_L,param,erlot_obs);

% CHX

chx_obs = 268;
chx1 = -761.073;
chx2 = -740.749;
chx3 = -740.750;
chxe = -742.745;

chx_L=[chx1 chxe chx2 chx3];

[chx_aic,chx_relp] = akaikec(chx_L,param,chx_obs);

%AT1

at1_obs = 188;
at11 = -423.389;
at12 = -390.634;
at13 = -390.618;
at1e = -396.004;

at1_L=[at11 at1e at12 at13];

[at1_aic,at1_relp] = akaikec(at1_L,param,at1_obs);

%MCF

mcf_obs = 106;
mcf1 = -323.614;
mcf2 = -299.723;
mcf3 = -299.723;
mcfe = -304.885;

mcf_L=[mcf1 mcfe mcf2 mcf3];

[mcf_aic,mcf_relp] = akaikec(mcf_L,param,mcf_obs);

%PC9

pc9_obs = 132;
pc91 = -454.058;
pc92 = -424.529;
pc93 = -424.529;
pc9e = -430.149;

pc9_L=[pc91 pc9e pc92 pc93];

[pc9_aic,pc9_relp] = akaikec(pc9_L,param,pc9_obs);
