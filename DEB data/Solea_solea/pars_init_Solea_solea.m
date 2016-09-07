function [par, metaPar, txtPar] = pars_init_Solea_solea(metaData)

metaPar.model = 'abj'; 

% reference parameter (not to be changed)
par.T_ref = C2K(20); free.T_ref = 0; units.T_ref = 'K';        label.T_ref = 'Reference temperature';

%% core primary parameters
par.z = 2.670;     free.z     = 1;   units.z     = '-';        label.z     = 'zoom factor';
par.F_m   = 6.5;   free.F_m   = 0;   units.F_m   = 'l/d.cm^2'; label.F_m   = '{F_m}, max spec searching rate';
par.kap_X = 0.8;   free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve';
par.kap_P = 0.1;   free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces';
par.v = 0.01124;   free.v     = 1;   units.v     = 'cm/d';     label.v     = 'energy conductance';
par.kap = 0.7353;  free.kap   = 1;   units.kap   = '-';        label.kap   = 'allocation fraction to soma';
par.kap_R = 0.95;  free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency';
par.p_M = 22.5;    free.p_M   = 1;   units.p_M   = 'J/d.cm^3'; label.p_M   = '[p_M], vol-spec somatic maint';
par.p_T   =  0;    free.p_T   = 0;   units.p_T   = 'J/d.cm^2'; label.p_T   = '{p_T}, surf-spec somatic maint';
par.k_J   = 0.002; free.k_J   = 0;   units.k_J   = '1/d';      label.k_J   = 'maturity maint rate coefficient';
par.E_G   = 5222;  free.E_G   = 1;   units.E_G   = 'J/cm^3';   label.E_G   = '[E_G], spec cost for structure';
par.E_Hb = 5.813e-2; free.E_Hb  = 1; units.E_Hb  = 'J';        label.E_Hb  = 'maturity at birth';
par.E_Hj = 3.916e0; free.E_Hj  = 1;  units.E_Hj  = 'J';        label.E_Hj  = 'maturity at metam';
par.E_Hp = 1.515e5; free.E_Hp  = 1;  units.E_Hp  = 'J';        label.E_Hp  = 'maturity at puberty';
par.h_a = 9.318e-9; free.h_a   = 1;  units.h_a   = '1/d^2';    label.h_a   = 'Weibull aging acceleration';
par.s_G   = 1e-4;  free.s_G   = 0;   units.s_G   = '-';        label.s_G   = 'Gompertz stress coefficient';

%% auxiliary parameters
par.T_A   = 6114;   free.T_A   = 0;    units.T_A = 'K';        label.T_A = 'Arrhenius temperature';
par.del_M = 0.1493; free.del_M = 1;    units.del_M = '-';      label.del_M = 'shape coefficient';
par.del_Me = 0.00843; free.del_Me = 1; units.del_Me = '-';     label.del_Me = 'shape coefficient for embryo';

%% environmental parameters (temperatures are in auxData)
par.f = 1.0;        free.f     = 0;    units.f = '-';          label.f    = 'scaled functional response for 0-var data';
par.f_tL = 0.6749;  free.f_tL  = 1;    units.f_tL = '-';       label.f_tL = 'scaled functional response for tL data';
par.f_Tab = 1.081;  free.f_Tab  = 1;   units.f_Tab = '-';      label.f_Tab = 'scaled functional response for Tab data';
par.f_Tabj = 1.117; free.f_Tabj  = 1; units.f_Tabj = '-';     label.f_Tabj = 'scaled functional response for Tabj data';
par.f_LN = 0.5286;  free.f_LN  = 1;    units.f_LN = '-';       label.f_LN = 'scaled functional response for LN data';

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);

%% Pack output:
txtPar.units = units; txtPar.label = label; par.free = free; 

