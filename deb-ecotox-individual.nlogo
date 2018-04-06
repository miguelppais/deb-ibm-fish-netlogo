; script-file DEB-IBM with different DEBMoAs
; based on  "Limitations of extrapolating toxic effects on reproduction to the population level" Ecological Applications - by Ben Martin (btmarti25@gmail.com)
; ==========================================================================================================================================
; ========================== DEFINITION OF PARAMETERS AND STATE VARIABLES ==================================================================
; ==========================================================================================================================================
; global parameters: are accessible for patches and turtles
globals [
 repro-stress food-stress
 food-list   abundance-list  biomass-list  length-list ;book keeping variables for populations metrics
 mean-food-list mean-abundance-list  mean-biomass-list  mean-length-list
 ;volume ; volume of water being simulated (adjusted in startup to result in ~150 daphnia in control simulations)
 day
  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ; - - - - - - - - - - - - - - - EMBRYO (parameters used to calculate the costs for an egg / initial reserves) - - - - - - - - --
  embryo-timestep
  e_scaled_embryo
  e_ref
  U_E_embryo
  S_C_embryo
  U_H_embryo
  L_embryo
  dU_E_embryo
  dU_H_embryo
  dL_embryo
  L_0      ; cm, initial structural volume
  lower-bound ; lower boundary for shooting method
  upper-bound ; upper boundary for shooting method
  estimation  ; estimated value for the costs for an egg / initial reserve
  sim         ; this keeps track of how many times the calc-egg-size loop is r
  crit-mass
]
; ------------------------------------------------------------------------------------------------------------------------------------------
patches-own ; Resource state variables are patch variables
[
  R        ; # / cm^2, prey density
  d_X      ; change of prey density in time
]
; ------------------------------------------------------------------------------------------------------------------------------------------
; definition of parameters for the individuals:
; the notation follows the DEBtool-notation as far as possible
; deviation: rates are indicated with "_rate" rather than a dot
; each individual(turtle) in the model has the following parameters
turtles-own [
  ; - - - - - - - - - - - - - - - STATE VARIABLES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  L           ; cm, structural length
  dL          ; change of structural length in time
  U_H         ; t L^2, scaled maturity
  dU_H        ; change of scaled maturity in time
  U_E         ; t L^2, scaled reserves
  dU_E        ; change of scaled reserves in time
  e_scaled    ; - , scaled reserves per unit of structure
  U_R         ; t L^2, scaled energy in reproduction buffer (not standard DEB)
  dU_R        ; change of energy in reproduction buffer (reproduction rate)
  f           ; scaled functional response
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; - - - - - - - - - - - - - - - FLUXES (used by several submodels) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  S_A         ; assimilation flux
  S_C         ; mobilisation flux
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; - - - - - - - - - - - - - - - STANDARD DEB PARAMETERS (with dimension and name) - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  g           ; - , energy investment ratio
  v_rate      ; cm /t , energy conductance (velocity)
  kap         ; - , allocation fraction to soma
  kap_R       ; - , reproduction efficiency
  k_M_rate    ; 1/t, somatic maintenance rate coefficient
  k_J_rate    ; 1/t, maturity maintenance rate coefficient
  U_H^b       ; t L^2, scaled maturity at birth
  U_H^p       ; t L^2, scaled maturity at puberty
  scatter-multiplier    ; parameter that is used to randomize the input parameters
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; - - - - - - - - - - - - - - - PREY DYNAMICS (only relevant if prey-dynamics are not constant) - - - - - - - - - - - - - - - - - - - - -
  J_XAm_rate  ; # / (cm^2 t), surface-area-specific maximum ingestion rate
  H           ; # / cm^2, (half) saturation coefficient
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; - - - - - - - - - - - - - - - AGEING -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  q_acceleration  ; - , ageing acceleration
  dq_acceleration ; change of ageing acceleration in time
  h_rate          ; - , hazard rate
  dh_rate         ; change of hazard rate in time
;-------------------------------------------------- daphnia specific parameters -------------------------------------------------------------------
  molt-time ; state variable for time since last molt
  max-l  ; maximum length daphnia has reached
  juvenile ; 0 in embryo, else 1
  adult ; 0 if embryo or juvenile, else 1
  offspring-number ;number of offspring an individual has energy to produce
  develop-time ; state variable for number of timesteps til embryo hatching
]
; ==========================================================================================================================================
; ========================== SETUP PROCEDURE: SETTING INITIAL CONDITIONS ===================================================================
; =============================== ===========================================================================================================
to setup
  ca
  reset-ticks
  ifelse original-data? [fill-original-data] [choose-species]
  if add_my_pet? [convert-parameters]
  set abundance-list []
  set  biomass-list []
  set length-list []
  set food-list []
  set food-stress 1
  set repro-stress 1
  set crit-mass 0.4 ; if mass (L^3) falls below crit-mass of previous maximum max (max-L^3) then there is high mort probability in death? submodel
  if prey-dynamics = "chemostat" [
  if Rmax = 3170 [set volume 2100000 * .0145 / rho]
  if Rmax = 7925 [set volume 670000 * .0145 / rho]
  if Rmax = 15850 [set volume 320000 * .0145 / rho]
  if Rmax = 31700 [set volume 160000 * .0145 / rho]
  ]
  calc-embryo-reserve-investment
  crt initial_nr_embryos ;number of initial embryos to start simulations with
  ask  turtles  [
    newborn-initialization  ; first their individual variability in the parameter is set     ; then the initial energy is calculated for each
  ]
  ask patches [set R H_int]; set initial value of prey to their carrying capacity
  if prey-dynamics = "logistic" [
   ask patches [set R 1]
  ]
end
; ==========================================================================================================================================
; ========================== GO PROCEDURE: RUNNING THE MODEL ===============================================================================
; ==========================================================================================================================================
; the go statement below is the order in which all procedures are run each timestep

to go
  if ticks / timestep = start_stress_day [add-stress]  ;initate stress on this day
  ask turtles with [juvenile = 1] [
    calc-dU_E                       ; first all individuals calculate the change in their state variables based on the current conditions
    calc-dU_H
    calc-dU_R
    calc-dL
    calc-dq_acceleration
    calc-dh_rate
  ]
  if prey-dynamics = "chemostat" [
    ask patches  [calc-d_X-chemostat]
  ]
  if prey-dynamics = "logistic" [
    ask patches  [calc-d_X-logistic]
  ]
  ask turtles  [update-individuals]       ; the the state variables of the individuals updated based on the delta value
  update-environment  ; the the state variables of the environment are updated based on the delta value
  death?
  ask turtles with [U_H >= U_H^p] [
    set molt-time molt-time + (1 / timestep)
    if molt-time > t-molt [lay-eggs]
  ]

  tick
 ;----------------------------------
 ; updating and observation steps
 ;----------------------------------
  if count turtles < 1 [stop]
  set day ticks / timestep
  if ticks / timestep > 300  [list-calc] ;start recording population response at day 300
  if ticks / timestep = 600 ;calculate means for each variable over the 300 day observation period
  [
    set mean-abundance-list mean abundance-list
    set mean-biomass-list   mean biomass-list
    set mean-length-list mean length-list
    set mean-food-list mean food-list
  ]
  if ticks /  timestep > simulation_days  [stop] ;end simulation after day
end

; ==========================================================================================================================================
; ========================== SUBMODELS =====================================================================================================
; ==========================================================================================================================================

; ------------------------------------------------------------------------------------------------------------------------------------------
; ---------------- conversion of parameters: from add-my-pet to standard DEB parameters ----------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to convert-parameters
  set U_H^b_int precision (E_H^b / p_Am) 4
  set U_H^p_int precision (E_H^p / p_Am) 4
  set k_M_rate_int precision (p_M / E_G) 4
  set g_int precision ((E_G * v_rate_int / p_Am) / kap_int) 4
end



; ------------------------------------------------------------------------------------------------------------------------------------------
; ------------------------ newborn-initialization ------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to newborn-initialization
  ; individuals vary in their DEB paramters on a normal distribution with a mean on the input paramater and a coefficent of variation equal to the cv
  ; set cv to 0 for no variation
  set scatter-multiplier e ^ (random-normal 0 cv)
  set juvenile 1
  set L L_embryo * scatter-multiplier
  set U_E U_E_embryo * scatter-multiplier
  set max-l L
  set J_XAm_rate   J_XAm_rate_int * scatter-multiplier
  set g g_int / scatter-multiplier
  set U_H^b U_H^b_int / scatter-multiplier ;
  set U_H U_H^b
  set U_H^p U_H^p_int / scatter-multiplier ;
  set v_rate v_rate_int
  set kap kap_int
  set kap_R kap_R_int
  set k_M_rate k_M_rate_int
  set k_J_rate k_J_rate_int
   set H  H_int * scatter-multiplier
  if DEBMoA = "maintenance costs" and ticks / timestep > 150
  [ set k_M_rate k_M_rate_int  *   (1 + stress-level)
    set k_J_rate k_J_rate_int * (1 + stress-level)]
  if DEBMoA = "growth costs"  and ticks / timestep > 150
  [set k_m_rate k_m_rate_int / (1 + stress-level)
   set g (g_int * ( 1 + stress-level)) / scatter-multiplier ]
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- RESERVE DYNAMICS -------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; change in reserves: determined by the difference between assimilation (S_A) and mobilization (S_C) fluxes
; when food-dynamics are constant f = the value of f_scaled set in the user interface
; if food is set to  "logistic" f depends on prey density and the half-saturation coefficient (K)
; for embryos f = 0 because they do not feed exogenously

to calc-dU_E
  ifelse prey-dynamics = "constant" [set f f_scaled] [set f R / (H + R)]
  set e_scaled v_rate * (U_E / L ^ 3)
  set S_C L ^ 2 * (g * e_scaled / (g + e_scaled)) * (1 + (L / (g * (V_rate / ( g * K_M_rate)))))
  set S_A  food-stress * f * max-L ^ 2 ;
  set dU_E (S_A - S_C )
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- MATURITY AND REPRODUCTION  ---------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; change in maturity is calculated (for immature individuals only)
to calc-dU_H
  ifelse  adult = 0 ; they only invest into maturity until they reach puberty
  [set dU_H ((1 - kap) * S_C - k_J_rate * U_H) ]
  [ set dU_H 0 ]
end

; the following procedure calculates change in reproduction buffer if mature
to calc-dU_R
  if U_H >= U_H^p
  [ set dU_R  ((1 - kap) * S_C - k_J_rate * U_H^p)  ]
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- DYNAMICS OF STRUCTURAL LENGTH-------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; the following procedure calculates change in structural length, if growth in negative the individual does not have enough energy to pay somatic maintenance individuals shrink

to calc-dL
  set dL   ((1 / 3) * (((V_rate /( g * L ^ 2 )) * S_C) - k_M_rate * L))
end
;-------------------------------------------------------------------------------------------------------------------------------------------
;--------- LAY EGGS ------------------------------------------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------------------------------------------------------------
;the following procedure is run for mature individuals which reached their molt-time
to lay-eggs
  let offspring-number-1  floor ((U_R * kap_R )/ estimation )   ;calc how many offspring there is enough energy in the reproduction buffer for
  let bin_ct 0
  repeat offspring-number-1                                    ;for each embryo
  [if random-float 1 < repro-stress [set bin_ct bin_ct + 1]]   ;chance of dying due to direct reproduction stress
  set offspring-number bin_ct                                  ;number of embryos surviving stress
  set molt-time 0 ;reset
  if offspring-number-1 > 0
  [set U_R U_R - ((floor ((U_R ) / estimation)) * estimation )]
  hatch offspring-number ;reset state variables of offspring
  [
    set adult 0
    set juvenile 0
    set  max-l 0
    set offspring-number 0
    set L 0
    set U_E 0
    set U_H 0
    set U_R 0
    set dU_R  0
    set h_rate 0
    set dh_rate 0
    set q_acceleration 0
    set dq_acceleration 0
    set S_A 0
    set S_C 0
    set dL 0
    set develop-time floor(t-molt * timestep - 1) ;number of timesteps til embryo will hatch
  ]
end
; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- Resource dynamics ----------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; chemostat prey dynamics (rho is the dilution rate, Rmax is the maximum resource density)
to calc-d_X-chemostat
  set d_X  (Rmax - R) * (rho) - (sum [S_A * J_XAm_rate] of turtles-here / volume)
end

; logistic prey dynamics (rX is the growth rate, KX is the carrying capacity)
to calc-d_X-logistic
  set d_X (rX) * R * (1 - (R / KX)) - (sum [S_A * J_XAm_rate] of turtles-here / volume)
end
; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- AGEING -----------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; the following procedure calculates the change in damage enducing compounds of an individual
to calc-dq_acceleration
  set dq_acceleration (q_acceleration * (L ^ 3 / (v_rate / ( g * k_M_rate)) ^ 3) * sG + H_a)* f * (( v_rate / L) - ((3 / L)*  dL)) - ((3 / L ) * dL) * q_acceleration
end

; the following procedure calculates the change in damage in the individual
to calc-dh_rate
  set dh_rate q_acceleration - ((3 / L) * dL) * h_rate
end
; ------------------------------------------------------------------------------------------------------------------------------------------
; -----------------  update-individuals-----------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
to update-individuals
  set U_E U_E + dU_E / timestep
  set U_H U_H + dU_H / timestep
  set U_R U_R + dU_R / timestep
  if U_R < 0  [set U_R 0 ] ;repro buffer can't be negative
  set L L + dL / timestep
  if U_H > U_H^b
  [
    set q_acceleration q_acceleration + dq_acceleration  / timestep
    set h_rate h_rate + dh_rate  / timestep
  ]
  if l > max-l  [ set max-l l]
  if juvenile = 0 ; if an embryo
  [
    set develop-time develop-time - 1  ;start a embryo with development time = t-molt * timestep, each timestep develop-time is reduced an the embryo hatches when it = 0
    if develop-time = 0 [newborn-initialization]
  ]
  if adult = 0[ if U_H >= U_H^p [ set adult 1] ]
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; -----------------    mortality submodels    ----------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
to  death?
  ;Juvenile resource dependent mortality
  if any? turtles with [Juvenile = 1 and adult = 0]
    [
      ask turtles with [juvenile = 1 and adult = 0 ]
        [
          if random-float 1 < ( 1 - e_scaled) * ( 1 - (1 - juv-mort) ^ ( 1 / timestep))   [die]
        ]
    ]
  ;High mortality if individuals shink below 40% of their maximum size
  ask turtles with [   l ^ 3 < crit-mass  * max-l ^ 3 and random-float 1 < ( 1 - (1 - .35) ^ ( 1 / timestep))] [die]
  ;Ageing mortality
  ask turtles  [if random-float 1 <  ( 1 - (1 - h_rate) ^ ( 1 / timestep)) [die] ]
end
; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- record population metrics-----------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
to list-calc ;keep track of abundance, biomass, length, and resource density during obervation period (day 300-600)
  set abundance-list fput count turtles with [juvenile = 1] abundance-list
  set  biomass-list fput sum [L ^ 3] of turtles with [juvenile = 1] biomass-list
  if count turtles > 1[ set length-list fput mean [max-L] of turtles with [juvenile = 1] length-list]
  set food-list fput [R] of patch 0 0 food-list
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; -----------------  update-environment-----------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
to update-environment
  set day (ticks / timestep)
  ask patches
  [set R  R + (d_x / timestep) ]
end
; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- add stressor on a specific day (default 150) ---------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
to add-stress
  if DEBMoA = "maintenance costs" [ask turtles [ set k_m_rate k_m_rate * ( 1 + stress-level) set k_j_rate k_j_rate * (1 + stress-level)]]
  if DEBMoA = "growth costs"  [ask turtles [set k_m_rate k_m_rate / (1 + stress-level) set g g * ( 1 + stress-level)]]
  if DEBMoA = "feeding ability"[set food-stress 1 / (stress-level + 1)]
  if DEBMoA = "embryonic hazard"[set repro-stress exp(- stress-level)]
  if DEBMoA = "cost per egg" [ask turtles [set kap_R kap_R / (1 - stress-level)]]
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ------------------------ INITIAL ENERGY --------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; calculate the the amount of energy needed to create an offspring, and get initial state variables for newborns based on the DEB paramters entered in the interface.
;Because this is a time intensive process, we calculate these once at the begining of each simulation and use the values for all newborns
; estimation = the energy required to produce an offspring
; L_embryo is the length at birth
; U_E_embryo is the scaled reserve at birth
;these values are used as inputs to the newborn-initialization submodel
to calc-embryo-reserve-investment
  set L_0 .00001
  set embryo-timestep timestep * 100
  set lower-bound 0
  set upper-bound 1
  set sim 0
  loop
  [
    set sim sim + 1
    set estimation .5 * (lower-bound + upper-bound)
    set L_embryo  L_0
    set U_E_embryo estimation
    set U_H_embryo  0
    set e_scaled_embryo v_rate_int * (U_E_embryo / L_embryo  ^ 3)
    set e_ref 1
    while [U_H_embryo  < U_H^b_int and e_scaled_embryo > 1 ]
      [
        set e_scaled_embryo v_rate_int * (U_E_embryo / L_embryo  ^ 3)
        set S_C_embryo L_embryo  ^ 2 * (g_int * e_scaled_embryo / (g_int + e_scaled_embryo)) * (1 + (L_embryo  / (g_int * (v_rate_int / ( g_int * k_M_rate_int)))))
        set dU_E_embryo  ( -1 * S_C_embryo )
        set dU_H_embryo  ((1 - kap_int) * S_C_embryo - k_J_rate_int * U_H_embryo  )
        set dL_embryo   ((1 / 3) * (((V_rate_int /( g_int * L_embryo  ^ 2 )) * S_C_embryo) - k_M_rate_int * L_embryo ))
        set  U_E_embryo  U_E_embryo +  dU_E_embryo    / (embryo-timestep )
        set  U_H_embryo   U_H_embryo  +  dU_H_embryo   / (embryo-timestep )
        set  L_embryo   L_embryo  +  dL_embryo    / (embryo-timestep )
        set e_scaled_embryo v_rate_int * (U_E_embryo / L_embryo  ^ 3)
      ]
    if e_scaled_embryo <  (.01 +  e_ref) and e_scaled_embryo > (-.01 + e_ref) and U_H_embryo  >= U_H^b_int
      [  stop ]
    ifelse U_H_embryo  > U_H^b_int
      [ set upper-bound estimation]
      [  set lower-bound estimation ]
    if sim > 100
    [user-message ("Embryo submodel did not converge. Timestep may need to be smaller." )  stop ] ;if the timestep is too big relative to the speed of growth of species this will no converge
  ]
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- SPECIES PARAMETERS -----------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to fill-original-data
set add_my_pet? FALSE
set kap_int 0.678
set kap_R_int 0.95
set k_M_rate_int 0.3314
set k_J_rate_int 0.1921
set U_H^b_int 0.1108
set U_H^p_int 2.547
set v_rate_int 18.1
set g_int 10
set h_a 3.04E-6
set sG 0.019
set J_XAm_rate_int 380000
set H_int 1585
set juv-mort 0.09
set t-molt 2.8
set timestep 50
end



to choose-species
  if species = "Daphnia magna" [
    set p_Am 315.611
    set shape_factor 0.264
    set v_rate_int 0.1584
    set kap_int 0.61
    set kap_R_int 0.95
    set p_m 1453
    set E_G 4400
    set K_J_rate_int 0.002
    set E_H^b 0.01379
    set E_H^p 0.3211
    ;ageing
    set h_a 0.0003105
    set sG -0.3
    ;feeding
    set F_m  6.5
    set J_Xam_rate_int 380000
  ]

    if species = "Daphnia pulex" [
    set p_Am 520.825
    set shape_factor 0.37
    set v_rate_int 0.03627
    set kap_int 0.763
    set kap_R_int 0.95
    set p_m 1400
    set E_G 4400
    set K_J_rate_int 0.002
    set E_H^b 0.02251
    set E_H^p 0.6024
    ;ageing
    set h_a 0.0001116
    set sG 0.0001
    ;feeding
    set F_m  1
    set J_Xam_rate_int 380000
  ]

    if species = "Pomatoschistus minutus" [
    set p_Am 167.547
    set shape_factor 0.128
    set v_rate_int 0.01346
    set kap_int 0.9548
    set kap_R_int 0.95
    set p_m 128.7
    set E_G 5181
    set K_J_rate_int 0.002
    set E_H^b 0.004883
    set E_H^p 56.79
    ;ageing
    set h_a 8.868E-6
    set sG 0.0001
    ;feeding
    set F_m 6.5
    set J_Xam_rate_int 1
  ]

  if species = "Danio rerio" [
    set p_Am 246.32
    set shape_factor 0.1325
    set v_rate_int 0.0278
    set kap_int 0.4366
    set kap_R_int 0.95
    set p_m 500.9
    set E_G 4652
    set K_J_rate_int 0.01662
    set E_H^b 0.5402
    set E_H^p 2062
    ;ageing
    set h_a 1.96E-9
    set sG 0.0405
    ;feeding
    set F_m 1
    set J_Xam_rate_int 1
  ]


end
@#$#@#$#@
GRAPHICS-WINDOW
1410
10
1583
184
-1
-1
165.0
1
10
1
1
1
0
1
1
1
0
0
0
0
1
1
1
ticks
30.0

BUTTON
85
50
151
83
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
15
50
81
83
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
155
50
221
83
go-once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
10
520
91
580
v_rate_int
18.1
1
0
Number

INPUTBOX
10
160
91
220
kap_int
0.678
1
0
Number

INPUTBOX
10
220
91
280
kap_R_int
0.95
1
0
Number

INPUTBOX
10
280
91
340
k_M_rate_int
0.3314
1
0
Number

INPUTBOX
10
340
90
400
k_J_rate_int
0.1921
1
0
Number

INPUTBOX
10
580
90
640
g_int
10.0
1
0
Number

INPUTBOX
10
400
90
460
U_H^b_int
0.1108
1
0
Number

INPUTBOX
10
460
90
520
U_H^p_int
2.547
1
0
Number

INPUTBOX
255
220
348
282
H_int
1585.0
1
0
Number

INPUTBOX
255
160
348
220
J_XAm_rate_int
380000.0
1
0
Number

INPUTBOX
10
670
90
730
h_a
3.04E-6
1
0
Number

INPUTBOX
10
730
90
790
sG
0.019
1
0
Number

INPUTBOX
10
835
90
895
cv
0.05
1
0
Number

TEXTBOX
15
140
135
158
DEB-IBM parameters\n
11
0.0
1

TEXTBOX
255
142
405
160
Feeding parameters
11
0.0
1

TEXTBOX
15
650
165
668
Ageing parameters
11
0.0
1

TEXTBOX
99
853
249
888
Coefficient of variation of parameter values
11
0.0
1

INPUTBOX
255
375
350
435
t-molt
2.8
1
0
Number

INPUTBOX
255
315
350
375
juv-mort
0.09
1
0
Number

INPUTBOX
1160
380
1315
440
rho
0.05
1
0
Number

PLOT
535
160
1150
310
Population abundance
Day
Abundance
0.0
600.0
0.0
0.0
true
false
"" ""
PENS
"total" 1.0 0 -16777216 true "" "if ticks mod timestep = 0[\n  if any? turtles with [juvenile = 1]\n   [plot count turtles with [juvenile = 1]]]"

INPUTBOX
600
695
745
755
stress-level
0.0
1
0
Number

CHOOSER
600
649
745
694
DEBMoA
DEBMoA
"maintenance costs" "growth costs" "feeding ability" "cost per egg" "embryonic hazard"
2

MONITOR
1160
55
1370
100
Mean population abundance
mean-abundance-list
2
1
11

MONITOR
1160
103
1371
148
Mean population biomass (cubic mm)
mean-biomass-list
2
1
11

MONITOR
1160
150
1372
195
Mean length of individuals (mm)
mean-length-list
2
1
11

MONITOR
1160
194
1374
239
Mean resource density (fraction of H)
mean-food-list / H_int
2
1
11

CHOOSER
1160
440
1315
485
Rmax
Rmax
3170 7925 15850 31700
1

TEXTBOX
96
175
246
203
\"kappa\" Fraction of mobilized energy to soma
11
0.0
1

TEXTBOX
97
229
247
271
\"Kappa R\", Fraction of reproduction energy fixed in eggs
11
0.0
1

TEXTBOX
98
297
248
325
Somatic maintenance rate coefficient
11
0.0
1

TEXTBOX
95
354
245
382
Maturity maintenance rate coefficient
11
0.0
1

TEXTBOX
94
419
244
437
Scaled maturity at birth
11
0.0
1

TEXTBOX
96
480
246
498
Scaled maturity at puberty
11
0.0
1

TEXTBOX
97
539
247
557
Energy conductance
11
0.0
1

TEXTBOX
96
602
246
620
Energy investment ratio
11
0.0
1

TEXTBOX
96
697
246
715
Hazard rate
11
0.0
1

TEXTBOX
100
755
250
773
Gompertz stress coefficient
11
0.0
1

TEXTBOX
359
173
525
215
Maximum specific ingestion rate (nr/mm^2/day)
11
0.0
1

TEXTBOX
360
232
510
258
Half-saturation coefficent (nr)
11
0.0
1

TEXTBOX
1163
361
1313
379
Chemostat growth
11
0.0
1

TEXTBOX
1320
410
1470
428
Resource dilution rate
11
0.0
1

TEXTBOX
1319
460
1469
478
Maximum resource density
11
0.0
1

TEXTBOX
365
335
515
363
Juvenile resource dependent mortality rate
11
0.0
1

TEXTBOX
365
390
515
418
Days between reproductive molts
11
0.0
1

TEXTBOX
255
290
405
308
Daphnia specifc parameters
11
0.0
1

PLOT
534
10
1151
160
Population biomass
Day
Biomass
0.0
600.0
0.0
0.0
true
false
"" ""
PENS
"Biomass" 1.0 0 -16777216 true "" "if ticks mod timestep = 0[\nplot sum [L ^ 3] of turtles with [juvenile = 1]]"

PLOT
534
310
1151
471
Mean length of individuals in population
Day
Length (mm)
0.0
600.0
0.0
0.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "if ticks mod timestep = 0 [\nif count turtles with [juvenile = 1] > 0\n[plot mean [max-L] of turtles with [juvenile = 1]]]"

TEXTBOX
775
665
1095
790
Daphnia effects on reproduction (21 day tests)\n\nEffect level    Feeding    Maintenance    Growth     Embryonic\n      25%          0.09              0.22             0.36            0.29\n      50%          0.23              0.50             0.95            0.69\n      75%          0.48              0.92             2.30            1.39\n      90%          0.78              1.28             4.65            2.30\n      95%          0.97              1.44             6.80            3.00
11
0.0
1

PLOT
534
471
1152
621
Resource density (fraction of half-saturation coefficient)
Day
R / H
0.0
600.0
0.0
0.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "if ticks mod timestep = 0[\nplot mean [R] of patches / H_int]"

TEXTBOX
603
628
753
646
Stress parameters
11
0.0
1

TEXTBOX
610
790
760
808
Stress begins on this day
11
0.0
1

TEXTBOX
1162
7
1312
49
Population variables measured and averaged over days 300-600
11
0.0
1

SLIDER
15
85
220
118
initial_nr_embryos
initial_nr_embryos
1
300
150.0
1
1
NIL
HORIZONTAL

SLIDER
600
755
745
788
start_stress_day
start_stress_day
0
500
150.0
5
1
NIL
HORIZONTAL

MONITOR
250
25
307
70
Day
floor(day)
17
1
11

CHOOSER
315
25
480
70
species
species
"Daphnia magna" "Daphnia pulex" "Pomatoschistus minutus" "Danio rerio"
0

INPUTBOX
255
850
375
910
shape_factor
0.264
1
0
Number

INPUTBOX
255
610
375
670
p_M
1453.0
1
0
Number

INPUTBOX
255
670
375
730
E_G
4400.0
1
0
Number

INPUTBOX
255
730
375
790
E_H^b
0.01379
1
0
Number

INPUTBOX
255
790
375
850
E_H^p
0.3211
1
0
Number

INPUTBOX
255
490
375
550
p_Am
315.611
1
0
Number

CHOOSER
1160
245
1375
290
prey-dynamics
prey-dynamics
"constant" "chemostat" "logistic"
1

TEXTBOX
1160
495
1310
513
Logistic growth
11
0.0
1

INPUTBOX
255
550
375
610
F_m
6.5
1
0
Number

INPUTBOX
1160
515
1315
575
rX
10.0
1
0
Number

INPUTBOX
1160
575
1315
635
KX
250.0
1
0
Number

INPUTBOX
1160
290
1375
350
volume
194300.0
1
0
Number

SLIDER
250
75
480
108
simulation_days
simulation_days
100
5000
600.0
100
1
NIL
HORIZONTAL

SWITCH
255
455
375
488
add_my_pet?
add_my_pet?
1
1
-1000

TEXTBOX
385
510
535
536
Surface area-specific assimilation flux
11
0.0
1

TEXTBOX
385
565
535
591
Max surface area-specific search rate
11
0.0
1

TEXTBOX
385
625
535
651
Volume-specific somatic maintenance
11
0.0
1

TEXTBOX
385
685
535
711
Volume-specific cost for structure
11
0.0
1

TEXTBOX
385
755
535
773
Maturity at birth
11
0.0
1

TEXTBOX
385
810
535
828
Maturity at puberty
11
0.0
1

SLIDER
15
15
220
48
timestep
timestep
0
500
50.0
5
1
/day
HORIZONTAL

TEXTBOX
1320
605
1470
623
Carrying capacity
11
0.0
1

TEXTBOX
1320
545
1470
563
Growth rate
11
0.0
1

SLIDER
1160
665
1315
698
f_scaled
f_scaled
0
1
1.0
0.05
1
NIL
HORIZONTAL

TEXTBOX
1165
645
1315
663
Constant growth
11
0.0
1

TEXTBOX
1320
675
1470
693
Scaled functional response
11
0.0
1

PLOT
1530
365
1730
515
Individual growth
Day
L (cm)
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"Growth" 1.0 0 -16777216 true "" "if any? turtles with [who = 0] and ticks mod timestep = 0 [plot [L] of turtle 0]"

SWITCH
375
455
502
488
original-data?
original-data?
0
1
-1000

@#$#@#$#@
## WHAT IS IT?

This section could give a general understanding of what the model is trying to show or explain.

## HOW IT WORKS

This section could explain what rules the agents use to create the overall behavior of the model.

## HOW TO USE IT

This section could explain how to use the model, including a description of each of the items in the interface tab.

## THINGS TO NOTICE

This section could give some ideas of things for the user to notice while running the model.

## THINGS TO TRY

This section could give some ideas of things for the user to try to do (move sliders, switches, etc.) with the model.

## EXTENDING THE MODEL

This section could give some ideas of things to add or change in the procedures tab to make the model more complicated, detailed, accurate, etc.

## NETLOGO FEATURES

This section could point out any especially interesting or unusual features of NetLogo that the model makes use of, particularly in the Procedures tab.  It might also point out places where workarounds were needed because of missing features.

## RELATED MODELS

This section could give the names of models in the NetLogo Models Library or elsewhere which are of related interest.

## CREDITS AND REFERENCES

This section could contain a reference to the model's URL on the web if it has one, as well as any other necessary credits or references.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>ticks = 1100</exitCondition>
    <metric>defects</metric>
    <enumeratedValueSet variable="timestep">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim">
      <value value="10"/>
      <value value="25"/>
      <value value="50"/>
      <value value="75"/>
      <value value="100"/>
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bound-shift">
      <value value="0.05"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_P_H">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-t">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_U_b_H">
      <value value="9.993569821026685E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_B_H">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_g">
      <value value="0.15003750937734434"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="-0.5"/>
      <value value="0"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_U_p_H">
      <value value="4.000107169649555E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_K_J_rate">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_kappa">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0.5"/>
      <value value="0.75"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_kappa_r">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.00125"/>
      <value value="1.25E-5"/>
      <value value="1.25E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="arrhenius">
      <value value="6400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Int_J_X_Am_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_K_M_rate">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_v_rate">
      <value value="0.16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="survival" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="k_envir">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_s">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_g">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_t">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="6.413E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;flow-through&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="102"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.32"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016913"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="batch neos only" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults
get-embryos</setup>
    <go>go</go>
    <metric>count turtles with [U_H &gt;= U_H^b]</metric>
    <metric>count turtles with [U_H &gt;= U_H^b and l &lt; .13]</metric>
    <metric>count turtles with [l &gt;= .13 and l &lt; .23]</metric>
    <metric>count turtles with [l &gt;= .23]</metric>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.22"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_t">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016913"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_g">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="6.413E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_int">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_envir">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.066E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="individual" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos</setup>
    <go>go</go>
    <metric>mean [l / .1] of turtles</metric>
    <metric>mean [offspring-number] of turtles</metric>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_t">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="2564000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016913"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lambda">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_g">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_s_int">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="6.813E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_int">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_envir">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment new variabilty" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos</setup>
    <go>go</go>
    <metric>mean [l / .1] of turtles</metric>
    <metric>mean [offspring-number] of turtles</metric>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_t">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="2564000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016913"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lambda">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_g">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_s_int">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="6.813E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.47"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_int">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_envir">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment individual middle food cv.5" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos</setup>
    <go>go</go>
    <metric>mean [l / .1] of turtles</metric>
    <metric>mean [offspring-number] of turtles</metric>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_t">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="3696000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016913"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lambda">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_g">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_s_int">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="7.013E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_int">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_envir">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="survival" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="initial-neos">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_s">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_g">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_t">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="5000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_envir">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="7.013E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016913"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;flow-through&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="50"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_s">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.47597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="7.013E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_t">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_envir">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_g">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016913"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos</setup>
    <go>go</go>
    <metric>count turtles with [U_H &gt;= U_H^b]</metric>
    <metric>count turtles with [U_H &gt;= U_H^b and l &lt; .13]</metric>
    <metric>count turtles with [l &gt;= .13 and l &lt; .23]</metric>
    <metric>count turtles with [l &gt;= .23]</metric>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_s">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.47597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="7.013E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_t">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_envir">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_g">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016913"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment high food" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos</setup>
    <go>go</go>
    <metric>count turtles with [U_H &gt;= U_H^b]</metric>
    <metric>count turtles with [U_H &gt;= U_H^b and l &lt; .13]</metric>
    <metric>count turtles with [l &gt;= .13 and l &lt; .23]</metric>
    <metric>count turtles with [l &gt;= .23]</metric>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_s">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.47597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="7.013E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_t">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_envir">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_g">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016913"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="LowFoodNeos" repetitions="100" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults
get-embryos</setup>
    <go>go</go>
    <metric>count turtles with [U_H &gt;= U_H^b]</metric>
    <metric>count turtles with [U_H &gt;= U_H^b and l &lt; .13]</metric>
    <metric>count turtles with [l &gt;= .13 and l &lt; .23]</metric>
    <metric>count turtles with [l &gt;= .23]</metric>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="25000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_g">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_t">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016913"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c0_s">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_envir">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.066E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="7.013E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults
get-embryos</setup>
    <go>go</go>
    <metric>count turtles with [U_H &gt;= U_H^b]</metric>
    <metric>count turtles with [U_H &gt;= U_H^b and max-l &lt; .13]</metric>
    <metric>count turtles with [max-l &gt;= .13 and l &lt; .23]</metric>
    <metric>count turtles with [max-l &gt;= .23]</metric>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="7.013E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016999"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6729"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.471733"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="26000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="64640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.33866"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.471"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="no dd" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults
get-embryos</setup>
    <go>go</go>
    <metric>count turtles with [U_H &gt;= U_H^b]</metric>
    <metric>count turtles with [U_H &gt;= U_H^b and max-l &lt; .13]</metric>
    <metric>count turtles with [max-l &gt;= .13 and l &lt; .23]</metric>
    <metric>count turtles with [max-l &gt;= .23]</metric>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6729"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="26000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.471733"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016999"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="64640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="7.013E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.33866"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.471"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="highfood" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults
get-embryos</setup>
    <go>go</go>
    <metric>count turtles with [U_H &gt;= U_H^b]</metric>
    <metric>count turtles with [U_H &gt;= U_H^b and max-l &lt; .115]</metric>
    <metric>count turtles with [max-l &gt;= .115 and l &lt; .2]</metric>
    <metric>count turtles with [max-l &gt;= .2]</metric>
    <enumeratedValueSet variable="sm2">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="26000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.471"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016999"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="7.013E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_d">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.471733"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="46"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.33866"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6729"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="false">
    <setup>setup
get-neos</setup>
    <go>go</go>
    <metric>fitness</metric>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_d">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.33866"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.47133"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6729"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="7.013E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="26000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.566E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.016999"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.471"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <steppedValueSet variable="juv-mort" first="3.0E-4" step="1.0E-5" last="0.002"/>
    <enumeratedValueSet variable="f_scaled">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.354"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="fast" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos</setup>
    <go>go</go>
    <metric>count turtles with U_H &gt;= U_H^b</metric>
    <enumeratedValueSet variable="g_int">
      <value value="10.235"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.7367"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="c_d">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_e_rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.344"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25460000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.01888"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.957"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="26000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.566E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="9.04E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="5200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.344"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="7.5E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="timestepexp" repetitions="30" runMetricsEveryStep="false">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with[ experiment = 3]</metric>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.07597"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.4363"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
      <value value="100"/>
      <value value="200"/>
      <value value="500"/>
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.066E-9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.00137"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="1.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.4054"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.4054"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="26000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.544"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.006"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.01772"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="no dd mort final paper g 10" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="false">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>fitness</metric>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0.003" step="5.0E-4" last="0.015"/>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="27500000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.005"/>
      <value value="6.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="27500000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="realexp" repetitions="50" runMetricsEveryStep="false">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0" step="5.0E-4" last="0.02"/>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
      <value value="&quot;all&quot;"/>
      <value value="&quot;adult-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="27500000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="27000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="4.5E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="27000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neo
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="6.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="27000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="6.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="27000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="proportinal" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.0058"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.51"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="27000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="fixed food input" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.0058"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="27000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.0058"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="27000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="new" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="X_k">
      <value value="70000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="30000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="5.5E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="2.5E-6"/>
      <value value="5.0E-6"/>
      <value value="1.0E-5"/>
      <value value="2.5E-5"/>
      <value value="5.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="2.5E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="30000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;adult-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="6.5E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.0055"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="6.5E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="6.5E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="6.5E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="7.5E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="7.25E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="7.25E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="7.5E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="7.5E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="9.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="9.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="0.005"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0"/>
      <value value="0.005"/>
      <value value="0.01"/>
      <value value="0.015"/>
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;adult-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.1"/>
      <value value="0.02"/>
      <value value="0.03"/>
      <value value="0.04"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="juv-mort">
      <value value="8.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm2">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_p">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.005"/>
      <value value="0.01"/>
      <value value="0.015"/>
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ad-mortality-constant">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sm3">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crowding-mort">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;batch-advanced&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H_b">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="false">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-adults-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0" step="0.02" last="0.2"/>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6829"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-adults-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6829"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6829"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0.02" step="0.02" last="0.2"/>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0.04" step="0.01" last="0.15"/>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-adults-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0.02" step="0.02" last="0.3"/>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0.03" step="0.01" last="0.15"/>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-adult-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0.02" step="0.02" last="0.3"/>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-adult-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0.11" step="0.01" last="0.2"/>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="28000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="29000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-adult-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1405"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="false">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="29000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1405"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="false">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="29000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02719"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.328"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0" step="0.005" last="0.2"/>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001196"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
      <value value="&quot;e-all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1405"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02498"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6529"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="65640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles with [experiment = 1]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [experiment = 2]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [experiment = 3]</metric>
    <metric>count turtles with [max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro-int">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="35000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input4">
      <value value="66640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos4">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.19211"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.809"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="66640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001098"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6784"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults4">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02702"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [juvenile = 1 and experiment = 1]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [juvenile = 1 and [experiment = 2]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [juvenile = 1 and experiment = 3]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro-int">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="35000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input4">
      <value value="66640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos4">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.095"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.19211"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.809"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="66640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001098"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6784"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults4">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02702"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="false">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro-int">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="35000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-adult-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input4">
      <value value="66640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos4">
      <value value="5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality-constant" first="0.01" step="0.005" last="0.4"/>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.19211"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.809"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="66640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001098"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6784"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults4">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02702"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup
get-neos
get-adults</setup>
    <go>go</go>
    <metric>count turtles with [juvenile = 1 and experiment = 1]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; .12 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 1]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 1]</metric>
    <metric>count turtles with [juvenile = 1 and experiment = 2]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; .12 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 2]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 2]</metric>
    <metric>count turtles with [juvenile = 1 and experiment = 3]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; .12 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .12 and max-l &lt; .2 and experiment = 3]</metric>
    <metric>count turtles with [max-l &gt;= .2 and experiment = 3]</metric>
    <enumeratedValueSet variable="start-day">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos4">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant-ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp3?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos3">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults2">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp2?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exp1?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="water-change">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input4">
      <value value="66640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-day">
      <value value="43"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="1.809"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.19211"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.6784"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input2">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02702"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos2">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro-int">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.38"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="0.02547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="35000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.001098"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults3">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-adults4">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-neos">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input3">
      <value value="66640000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles with [juvenile = 1]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; 1.1 ]</metric>
    <metric>count turtles with [max-l &gt;= 1.1 and max-l &lt; 2 ]</metric>
    <metric>count turtles with [max-l &gt;= 2 ]</metric>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02702"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="days-between-repro-int">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles with [juvenile = 1]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; 1.1 ]</metric>
    <metric>count turtles with [max-l &gt;= 1.1 and max-l &lt; 2 ]</metric>
    <metric>count turtles with [max-l &gt;= 2 ]</metric>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-stress-int">
      <value value="1.37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.02702"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;0.5 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles with [juvenile = 1]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; 1.1 ]</metric>
    <metric>count turtles with [max-l &gt;= 1.1 and max-l &lt; 2 ]</metric>
    <metric>count turtles with [max-l &gt;= 2 ]</metric>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-stress-int">
      <value value="1.37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.4E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;0.5 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles with [juvenile = 1]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; 1.1 ]</metric>
    <metric>count turtles with [max-l &gt;= 1.1 and max-l &lt; 2 ]</metric>
    <metric>count turtles with [max-l &gt;= 2 ]</metric>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.595"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.4E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;0.5 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>consumed / ((ticks / timestep) * pred)</metric>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.2"/>
      <value value="0.6"/>
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.4E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;0.5 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pred">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles with [juvenile = 1]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; 1.1 ]</metric>
    <metric>count turtles with [max-l &gt;= 1.1 and max-l &lt; 2 ]</metric>
    <metric>count turtles with [max-l &gt;= 2 ]</metric>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="25640000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.4E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maint-stress-int">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;0.5 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pred">
      <value value="1.0E-9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="26670000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.4E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maint-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;1.3 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pred">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>sum [L ^ 3 ] of turtles with [juvenile = 1]</metric>
    <metric>sum [max-L ^ 3 ] of turtles with [juvenile = 1]</metric>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="26670000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.195"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.4E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maint-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;1.3 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pred">
      <value value="5.0E-12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles with [juvenile = 1]</metric>
    <metric>count turtles with [juvenile = 1 and max-l &lt; 1.1 ]</metric>
    <metric>count turtles with [max-l &gt;= 1.1 and max-l &lt; 2 ]</metric>
    <metric>count turtles with [max-l &gt;= 2 ]</metric>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;0.5 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="26670000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maint-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stress">
      <value value="0.159"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.4E-6"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean-abundance-list</metric>
    <metric>mean-biomass-list</metric>
    <metric>mean-adult-list</metric>
    <metric>mean-mat-list</metric>
    <metric>mean-length-list</metric>
    <metric>mean-food-list</metric>
    <metric>var-abundance-list</metric>
    <enumeratedValueSet variable="maint-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.04E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ss">
      <value value="0"/>
      <value value="0.22"/>
      <value value="0.5"/>
      <value value="0.92"/>
      <value value="1.28"/>
      <value value="1.44"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="66670000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="650000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_R">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stress">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_K">
      <value value="3170"/>
      <value value="7925"/>
      <value value="15850"/>
      <value value="31700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;1.3 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.019"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean-abundance-list</metric>
    <metric>mean-biomass-list</metric>
    <metric>mean-adult-list</metric>
    <metric>mean-mat-list</metric>
    <metric>mean-length-list</metric>
    <metric>mean-food-list</metric>
    <metric>var-abundance-list</metric>
    <enumeratedValueSet variable="maint-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.04E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ss">
      <value value="1.44"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="66670000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="155000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_R">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stress">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_K">
      <value value="31700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;1.3 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.019"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="food-stress" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean-abundance-list</metric>
    <metric>mean-biomass-list</metric>
    <metric>mean-adult-list</metric>
    <metric>mean-mat-list</metric>
    <metric>mean-length-list</metric>
    <metric>mean-food-list</metric>
    <metric>var-abundance-list</metric>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_R">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates, 3 adults&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_K">
      <value value="31700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;1.3 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stress">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.04E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.019"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maint-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="650000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ss">
      <value value="0.8125"/>
      <value value="0.677"/>
      <value value="0.562"/>
      <value value="0.508"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="66670000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="mort juv only e, L used for escaled calc" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean-abundance-list</metric>
    <metric>mean-biomass-list</metric>
    <metric>mean-adult-list</metric>
    <metric>mean-mat-list</metric>
    <metric>mean-length-list</metric>
    <metric>mean-food-list</metric>
    <metric>var-abundance-list</metric>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_R">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_K">
      <value value="31700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stress-type">
      <value value="&quot;maint-stress&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;1.3 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maint-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ss">
      <value value="0"/>
      <value value="1.28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stress">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.04E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.019"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="66670000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="155000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="10 rep xr xk" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean-abundance-list</metric>
    <metric>mean-biomass-list</metric>
    <metric>mean-adult-list</metric>
    <metric>mean-mat-list</metric>
    <metric>mean-length-list</metric>
    <metric>mean-food-list</metric>
    <metric>var-abundance-list</metric>
    <metric>var-adult-list</metric>
    <metric>maturation-flux</metric>
    <metric>reproduction-flux</metric>
    <metric>maturity-prob</metric>
    <enumeratedValueSet variable="K_int">
      <value value="1585"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-mass">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dd-mortality">
      <value value="&quot;e-juvenile-only&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="0.1108"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-conditions?">
      <value value="&quot;5 neonates&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="0.019"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_R">
      <value value="0.05"/>
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stress">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="input">
      <value value="66670000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="3.04E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.1921"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ss">
      <value value="0"/>
      <value value="0.22"/>
      <value value="0.5"/>
      <value value="0.92"/>
      <value value="1.28"/>
      <value value="1.44"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.678"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-constant">
      <value value="0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="380000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="9280"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.3314"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-between-molts">
      <value value="2.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stress-type">
      <value value="&quot;maint-stress&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_K">
      <value value="3170"/>
      <value value="7925"/>
      <value value="15850"/>
      <value value="31700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maint-stress-int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-level">
      <value value="&quot;1.3 mgC&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="2.547"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-l?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@
