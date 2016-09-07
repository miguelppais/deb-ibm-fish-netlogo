function [prdData, info] = predict_Solea_solea(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
        
  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_tj = tempcorr(temp.tj, T_ref, T_A);
  TC_ap = tempcorr(temp.ap, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_R45 = tempcorr(temp.R45, T_ref, T_A);
  TC_tL = tempcorr(temp.tL, T_ref, T_A);
  TC_Tab = tempcorr(C2K(Tab(:,1)), T_ref, T_A);
  TC_Tabj = tempcorr(C2K(Tabj(:,1)), T_ref, T_A);  
  TC_LN = tempcorr(temp.LN, T_ref, T_A);

  % zero-variate data

  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_M;                % cm, total length at birth at f
  Ww_b = L_b^3 * (1 + f * w);       % g, wet weight at birth at f 
  a_b = t_b/ k_M; aT_b = a_b/ TC_ab; % d, age at birth at f and T

  % metam
  L_j = L_m * l_j;                  % cm, structural length at birth at f
  Lw_j = L_j/ del_M;                % cm, total length at birth at f
  tT_j = (t_j - t_b)/ k_M/ TC_tj;   % d, age at birth at metam

  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;                % cm, total length at puberty at f
  Ww_p = L_p^3 *(1 + f * w);        % g, wet weight at puberty 
  aT_p = t_p/ k_M/ TC_ap;           % d, age at birth at f and T

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate total length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 
 
  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
  RT_45 = TC_R45 * reprod_rate_j(45 * del_M, f, pars_R);  % #/d, reproduction rate for 45 cm

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
  % pack to output
  prdData.ab = aT_b;
  prdData.tj = tT_j;
  prdData.ap = aT_p;
  prdData.am = aT_m;
  prdData.Lb = Lw_b;
  prdData.Lj = Lw_j;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Wwb = Ww_b;
  prdData.Wwp = Ww_p;
  prdData.Wwi = Ww_i;
  prdData.R45 = RT_45;
  
  % uni-variate data
  
  % time-length for late juveniles and adults at f of 0-var data
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tL);
  Lw_i = l_i * L_m/ del_M; Lw_j = l_j * L_m/ del_M; % cm, total lengths
  rT_B = rho_B * k_M * TC_tL; tj = (t_j-t_b)/ k_M/ TC_tL; % 1/d, von Bert growth rate
  ELw = Lw_i - (Lw_i - Lw_j) * exp( - rT_B * (tL(:,1) - tj)); % cm, total length at time

  % Tab
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_Tab);
  Eab = t_b/ k_M ./ TC_Tab;              % d, age at birth     
  
  % Tabj
  [t_j t_p t_b l_j l_p l_b l_i r_j r_B info] = get_tj(pars_tj, f_Tabj);
  Eabj = (t_j - t_b)/ k_M./ TC_Tabj;     % d, time since birth at metam
  
  % LN
  EN = 365 * TC_LN * reprod_rate_j(LN(:,1) * del_M, f_LN, pars_R); % yearly reproductive output at length

  % pack to output
  prdData.tL = ELw;
  prdData.Tab = Eab;
  prdData.Tabj = Eabj;
  prdData.LN = EN;
  