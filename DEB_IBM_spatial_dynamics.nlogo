; script-file DEB-IBM (or however we are going to call it in the end)
; Authors: Ben Martin (btmarti25@gmail.com) and Elke Zimmer

; implementation of the standard DEB equations in an IBM
; check ... for the user-manual and the ODD
; published in "Dynamic Energy Budget theory meets individual-based modelling: a generic and accessible implementation", ... 2011

; ==========================================================================================================================================
; ========================== DEFINITION OF PARAMETERS AND STATE VARIABLES ==================================================================
; ==========================================================================================================================================

; global parameters: are accessible for patches and turtles
globals[
  U_E^0    ; t L^2, initial reserves of the embryos at the start of the simulation
  L_0      ; cm, initial structural volume
]  
; ------------------------------------------------------------------------------------------------------------------------------------------
; parameters for the environment: here only prey density

patches-own[
  X        ; # / cm^2, prey density
  d_X      ; change of prey density in time
]
; ------------------------------------------------------------------------------------------------------------------------------------------

; definition of parameters for the individuals: 
; the notation follows the DEBtool-notation as far as possible
; deviation: rates are indicated with "_rate" rather than a dot
; each individual(turtle) in the model has the following parameters
turtles-own[
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
  
  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                  
  ; - - - - - - - - - - - - - - - FLUXES (used by several submodels) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                     
 
  S_A         ; assimilation flux                                   
  S_C         ; mobilisation flux    
  
  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ; - - - - - - - - - - - - - - - EMBRYO (we use different state variable to not affect the state varibale of the mother) - - - - - - - - --
  e_scaled_embryo    
  e_ref                                       
  U_E_embryo                                                 
  S_C_embryo                                                 
  U_H_embryo                                                 
  L_embryo                                                   
  dU_E_embryo                                                 
  dU_H_embryo                                                 
  dL_embryo     
  ; parameters used to calculate the costs for an egg / initial reserves
  lower-bound ; lower boundary for shooting method                                                
  upper-bound ; upper boundary for shooting method                                                
  estimation  ; estimated value for the costs for an egg / initial reserve                                                 
  lay-egg?    ; parameter needed to hand over if an egg can be laid
  offspring-count ; with this parameter, the reproduction rate per turtle is shown on the interface
  sim         ; this keeps track of how many times the calc-egg-size loop is run
 
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
  ; parameter that is used to randomize the input parameters   
  scatter-multiplier    
 
  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ; - - - - - - - - - - - - - - - PREY DYNAMICS (only relevant if prey-dynamics are set to logistic) - - - - - - - - - - - - - - - - - - -                                                         
 
  J_XAm_rate  ; # / (cm^2 t), surface-area-specific maximum ingestion rate                                                    
  K           ; # / cm^2, (half) saturation coefficient                                                         
 
  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ; - - - - - - - - - - - - - - - AGEING -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                             
 
  q_acceleration  ; - , ageing acceleration
  dq_acceleration ; change of ageing acceleration in time                                                    
  h_rate          ; - , hazard rate
  dh_rate         ; change of hazard rate in time
  age-day         ; each turtle has a random whole number between 0 and timestep if the mod of ticks = the age day of a turtle is will check to see if it dies
                  ; based on the ageing submodel. This is because mortality probabilities are per day, and timesteps are smaller                                            

 f        ; - , scaled functional response
]

; ==========================================================================================================================================
; ========================== SETUP PROCEDURE: SETTING INITIAL CONDITIONS ===================================================================
; ==========================================================================================================================================

to setup
ca
reset-ticks 
 if add_my_pet? = "on" 
 [convert-parameters]
 
 set L_0 .00001
  
 crt 1000                   ; 10 turtles are created in the beginning    
 ask  turtles  [
  individual-variability  ; first their individual variability in the parameter is set
  calc-embryo-reserve-investment     ; then the initial energy is calculated for each
 setxy random-xcor random-ycor

 ]
  
 ask patches [ set X J_XAm_rate_int /   F_m ]; set initial value of prey to their carrying capacity
end

; ==========================================================================================================================================
; ========================== GO PROCEDURE: RUNNING THE MODEL ===============================================================================
; ==========================================================================================================================================
; the go statement below is the order in which all procedures are run each timestep

to go 
  
  ask turtles 
  [
    calc-dU_E                       ; first all individuals calculate the change in their state variables based on the current conditions
    calc-dU_H  
    calc-dU_R 
    calc-dL
  ] 
  if aging = "on"                   ; if the ageing submodel is turned on, the change in damage inducing compound and damage are calculated
  [ 
    ask turtles 
    [
      calc-dq_acceleration
      calc-dh_rate
    ]
  ]
  if food-dynamics = "logistic"     ; if prey dynamics are set to "logistic" the change in prey density is calculated
  [ask patches [calc-d_X]]
  
  update                           ; the the state variables of the individuals and prey are updated based on the delta value
  
  ask turtles 
  [ 
    if U_H >= U_H^p                 ; mature individual check if they have enough energy in their reproduction buffer to repdroduce
      [
        calc-lay-eggs
      ]
    if lay-egg? = 1  
      [
        calc-embryo-reserve-investment         ; if so, they calculate how much energy to invest in an embryo
        lay-eggs                    ; and they produce one offspring    
      ]  
  ]  
 movement-submodel
 
  ask patches [ set pcolor scale-color green X 2 0]
  tick
  do-plots                          ; then the plots are updated 
  if count turtles = 0 [stop]
end

; ==========================================================================================================================================
; ========================== SUBMODELS =====================================================================================================
; ==========================================================================================================================================

; ---------------- conversion of parameters: from add_my_pet to standard DEB parameters ----------------------------------------------------

to convert-parameters
  let p_am p_m * zoom / kap_int
  set U_H^b_int E_H_b / p_am
  set U_H^p_int E_H_p / p_am
  set k_M_rate_int p_m / E_G
  set g_int (E_G * v_rate_int / p_am) / kap_int
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ------------------------ INDIVIDUAL VARIABILITY ------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to individual-variability
  ; individuals vary in their DEB paramters on a normal distribution with a mean on the input paramater and a coefficent of variation equal to the cv
  ; set cv to 0 for no variation      
  set scatter-multiplier e ^ (random-normal 0 cv)
  set J_XAm_rate   J_XAm_rate_int * scatter-multiplier
  set g g_int / scatter-multiplier
  set U_H^b U_H^b_int / scatter-multiplier ; 
  set U_H^p U_H^p_int / scatter-multiplier ; 
 
  set v_rate v_rate_int   
  set kap kap_int
  set kap_R kap_R_int
  set k_M_rate k_M_rate_int                                           
  set k_J_rate k_J_rate_int
  set K  J_XAm_rate /   F_m 
  set age-day random timestep
end


; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- RESERVE DYNAMICS -------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; change in reserves: determined by the difference between assimilation (S_A) and mobilization (S_C) fluxes
; when food-dynamics are constant f = the value of f_scaled set in the user interface
; if food is set to  "logistic" f depends on prey density and the half-saturation coefficient (K)
; for embryos f = 0 because they do not feed exogenously

to calc-dU_E  
  
  if food-dynamics = "constant"
  [ ifelse U_H <= U_H^b
    [set f 0]
    [set f f_scaled] 
  ]
  if food-dynamics = "logistic"
  [ ifelse U_H <= U_H^b 
    [set f 0]
    [set f X / (K + X)]
  ]
  set e_scaled v_rate * (U_E / L ^ 3)
  set S_C L ^ 2 * (g * e_scaled / (g + e_scaled)) * (1 + (L / (g * (V_rate / ( g * K_M_rate)))))
  
  set S_A  f * L ^ 2 ;
  
  set dU_E (S_A - S_C )  
end 
; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- MATURITY AND REPRODUCTION  ---------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------ 
; change in maturity is calculated (for immature individuals only)

to calc-dU_H 
  
  ifelse U_H < U_H^p ; they only invest into maturity until they reach puberty 
    [set dU_H ((1 - kap) * S_C - k_J_rate * U_H) ]
    [set dU_H 0]
end

; the following procedure calculates change in reprobuffer if mature
to calc-dU_R  
  if U_H >= U_H^p
    [set dU_R  ((1 - kap) * S_C - k_J_rate * U_H^p) ]
end  

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- DYNAMICS OF STRUCTURAL LENGHT-------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; the following procedure calculates change in structural length, if growth in negative the individual does not have enough energy to pay somatic maintenance and the starvation submodel is run
; where growth is set to 0 and individuals divirt enough energy from development (for juveniles) or reprodution (for adults) to pay maintenance costs
to calc-dL 
  
  set dL   ((1 / 3) * (((V_rate /( g * L ^ 2 )) * S_C) - k_M_rate * L)) 
  
  if e_scaled < L / (V_rate / ( g * K_M_rate))  ; if growth is negative use starvation strategy 3 from the DEB book
    [
      set dl 0
      ifelse U_H < U_H^p
       [set dU_H (1 - kap) * e_scaled * L ^ 2 - K_J_rate * U_H^p - kap * L ^ 2 * ( L / (V_rate / ( g * K_M_rate)) - e_scaled)]
       [ set dU_R  (1 - kap) * e_scaled * L ^ 2 - K_J_rate * U_H^p - kap * L ^ 2 * ( L / (V_rate / ( g * K_M_rate)) - e_scaled)]
      set dU_E  S_A - e_scaled * L ^ 2
   ifelse U_H < U_H^p
 
 [  if dU_H < 0 [die]]
    
      [if U_R < 0 [die]]
    ]
 
end

;------------------------------------------------------------------------------------------------------------------------------------------
;---------- CHECK IF POSSIBLE TO LAY EGGS ------------------------------------------------------------------------------------------------- 
;------------------------------------------------------------------------------------------------------------------------------------------
; in the following, individuals determine if they have enough energy in their repro buffer to reproduce by creating an embryo with initial reserves set to the energy
; currently in their repro buffer * kap_R (conversion efficiancy of  reprobuffer to embryo) if the individual has enough energy to produce an offspring which will reach
; maturity and have a reserve density greater than the mothers when it hatches "lay-egg?" is set to one which will trigger the reproduction procedures "calc-egg-size" and "lay-eggs"
to calc-lay-eggs
  set L_embryo  L_0
  set U_E_embryo U_R * kap_R
  set U_H_embryo  0
  
  loop [   
    set e_scaled_embryo v_rate * (U_E_embryo / L_embryo  ^ 3)
    set S_C_embryo L_embryo  ^ 2 * (g * e_scaled_embryo / (g + e_scaled_embryo)) * (1 + (L_embryo  / (g * (V_rate / ( g * K_M_rate)))))
    
    set dU_E_embryo  ( -1 * S_C_embryo )  
    set dU_H_embryo  ((1 - kap) * S_C_embryo - k_J_rate * U_H_embryo ) 
    set dL_embryo  ((1 / 3) * (((V_rate /( g * L_embryo  ^ 2 )) * S_C_embryo) - k_M_rate * L_embryo )) 
    
    set  U_E_embryo  U_E_embryo +  dU_E_embryo  / timestep
    set  U_H_embryo  U_H_embryo  +  dU_H_embryo   / timestep
    set  L_embryo    L_embryo  +  dL_embryo   / timestep
    
    if U_H_embryo  > U_H^b * 1 [ set lay-egg? 1 stop]
    if e_scaled_embryo < e_scaled [stop]
    ] 
end 

; ------------------------------------------------------------------------------------------------------------------------------------------
; ------------------------ INITIAL ENERGY --------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; calculate the initial energy of the first individuals using a bisection method

to calc-embryo-reserve-investment
  set lower-bound 0
  ifelse ticks = 0
  [set upper-bound 100]
  [set upper-bound U_R * kap_R]
  set sim 0
 
  loop[
    set sim sim + 1
    
    set estimation .5 * (lower-bound + upper-bound)
    set L_embryo  L_0
    set U_E_embryo estimation
    set U_H_embryo  0
    set e_scaled_embryo v_rate * (U_E_embryo / L_embryo  ^ 3)
    
    ifelse ticks = 0[set e_ref 1][set e_ref e_scaled]  ; e_ref now determines which e_scaled_embryo to calculate: 1 for ticks = 0 (in the setup procedure), e_scaled otherwise
    
    while [U_H_embryo  < U_H^b and e_scaled_embryo > e_ref ] 
    ;     while [U_H_embryo  < U_H^b and e_scaled_embryo > 1 ] ; egg-size:  while [U_H_embryo  < U_H^b and e_scaled_embryo > e_scaled  ] ; I KEPT THIS LINE FOR NOW TO HAVE IT EASIER TO COMPARE
      [ set e_scaled_embryo v_rate * (U_E_embryo / L_embryo  ^ 3)
        set S_C_embryo L_embryo  ^ 2 * (g * e_scaled_embryo / (g + e_scaled_embryo)) * (1 + (L_embryo  / (g * (v_rate / ( g * k_M_rate)))))
        
        set dU_E_embryo  ( -1 * S_C_embryo )  
        set dU_H_embryo  ((1 - kap) * S_C_embryo - k_J_rate * U_H_embryo  ) 
        set dL_embryo   ((1 / 3) * (((V_rate /( g * L_embryo  ^ 2 )) * S_C_embryo) - k_M_rate * L_embryo )) 
        
        set  U_E_embryo  U_E_embryo +  dU_E_embryo    / (timestep )
        set  U_H_embryo   U_H_embryo  +  dU_H_embryo   / (timestep )
        set  L_embryo   L_embryo  +  dL_embryo    / (timestep )
      ] 
    
    if e_scaled_embryo <  .05 +  e_ref and e_scaled_embryo > -.05 + e_ref and U_H_embryo  >= U_H^b  [    
      ;      if e_scaled_embryo <  .05 +  1 and e_scaled_embryo > -.05 + 1 and U_H_embryo  >= U_H^b  [ ;egg-size: if e_scaled_embryo <  .05 +  e_scaled and e_scaled_embryo > -.05 + e_scaled and U_H_embryo  >= U_H^b  [stop]
      ;I KEPT THIS LINE FOR NOW TO HAVE IT EASIER TO COMPARE
      ifelse ticks = 0 ; 
      [set U_E^0 estimation
        set L L_0
        set U_E U_E^0
        set U_H 0
        set U_R 0
        set dU_R  0
        
        set age-day random timestep
        stop
      ][stop]]  
    
    ifelse U_H_embryo  > U_H^b  
      [ set upper-bound estimation ]
      [ set lower-bound estimation ]
    if sim > 100 [user-message ("Embryo submodel did not converge. Timestep may need to be smaller.") stop]
    ;if the timestep is too big relative to the speed of growth of species this will no converge
  ]
end

;-------------------------------------------------------------------------------------------------------------------------------------------
;--------- LAY EGGS ------------------------------------------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------------------------------------------------------------
;the following procedure is run for mature individuals which have enough energy to reproduce
; they create 1 offspring and give it the following state variables and DEB parameters
;the initial reserves is set to the value determined by the bisection method in "calc_egg_size"

to lay-eggs 
  hatch 1   
    [ 
      ;the following code give offspring varibility in their DEB paramters on a normal distribution with a mean on the input paramater and a coefficent of variation equal to the cv
      ; set cv to 0 for no variation  
  
      set scatter-multiplier e ^ (random-normal 0 cv)
      set J_XAm_rate   J_XAm_rate_int * scatter-multiplier
      set g g_int / scatter-multiplier
      set U_H^b    U_H^b_int / scatter-multiplier 
      set U_H^p    U_H^p_int / scatter-multiplier 
     
      set v_rate v_rate_int   
      
      set kap kap_int
      set kap_R kap_R_int
      set k_M_rate k_M_rate_int                                           
      set k_J_rate k_J_rate_int
      set  K J_XAm_rate /   F_m 
     
      set L L_0
      set U_E estimation
      set U_H 0
      set U_R 0
      set dU_R  0 
      set h_rate 0
      set dh_rate 0
      set q_acceleration 0
      set dq_acceleration 0
      set lay-egg? 0
      set age-day random timestep
    ]
  set lay-egg? 0  
  set U_R U_R - estimation / kap_R
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- LOGISTIC PREY ----------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
 ;the following procedure calculates change in prey density this procedure is only run when prey dynamics are set to "logistic" in the user interface

to calc-d_X 
   set d_X ((((X_r) * X * (1 - (X / X_k)) * volume))   - sum [ S_A * J_XAm_rate   ] of turtles-here) 
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- AGEING -----------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; the following procedure calculates the change in damage enducing compounds of an individual

to calc-dq_acceleration 
  set dq_acceleration (q_acceleration * (L ^ 3 / (v_rate / ( g * k_M_rate)) ^ 3) * sG + H_a) * e_scaled * (( v_rate / L) - ((3 / L)*  dL)) - ((3 / L ) * dL) * q_acceleration
end

; the following procedure calculates the change in damage in the individual
to calc-dh_rate 
  set dh_rate q_acceleration - ((3 / L) * dL) * h_rate
end


; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- UPDATE -----------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to update
; individuals update their state variables based on the calc_state variable proccesses
  ask turtles 
  [
    set U_E U_E + dU_E / timestep
    set U_H U_H + dU_H / timestep
    set U_R U_R + dU_R    / timestep
    set L L + dL    / timestep
    if U_H > U_H^b 
    [ set q_acceleration q_acceleration + dq_acceleration  / timestep
      set h_rate h_rate + dh_rate  / timestep
    ]
   if E_scaled <= 0.01 [die  ] ; starvation related mortality
    
   if aging = "on" [if ticks mod timestep = age-day [if random-float 1 < h_rate [die]] ] ;ageing related mortality
   if aging = "off" [if ticks mod timestep = age-day [if random-float 1 < background-mortality [die]] ]
 ]
  if food-dynamics = "logistic"[ ask patches [ set X X + ((d_X / timestep) / volume)]]
    ask patches [ if x < 0 [ set x .000001]]
     
end
; ------------------------------------------------------------------------------------------------------------------------------------------
; -------------- Movement submodel ---------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to movement-submodel
 ask turtles with [U_H > U_H^B] [if ticks mod 50 = 0[

let scale sum [x] of neighbors + [x] of patch-here
let p-a [x] of patch-at 0 1 / scale
let p-b [x] of patch-at 0 -1 / scale
let p-r [x] of patch-at 1 0 / scale
let p-l [x] of patch-at -1 0 / scale
let p-ar [x] of patch-at 1 1 / scale
let p-br [x] of patch-at 1 -1 / scale
let p-al [x] of patch-at -1 1 / scale
let p-bl [x] of patch-at -1 -1 / scale
let patch-h [x] of patch-at 0 0 / scale
let random-number random-float 1
if random-number < p-a [move-to patch-at 0 1]
if random-number >= p-a  and random-number < p-a + p-b [move-to patch-at 0 -1]
if random-number >= p-a + p-b and random-number < p-a + p-b + p-r [move-to patch-at 1 0]
if random-number >= p-a + p-b + p-r and  random-number < p-a + p-b + p-r + p-r [move-to patch-at -1 0]
if random-number >= p-a + p-b + p-r + p-l and  random-number < p-a + p-b + p-r + p-l + p-ar [move-to patch-at 1 1]
if random-number >= p-a + p-b + p-r + p-l + p-ar and  random-number < p-a + p-b + p-r + p-l + p-ar + p-br[move-to patch-at 1 -1]
if random-number >= p-a + p-b + p-r + p-l + p-ar + p-br and  random-number < p-a + p-b + p-r + p-l + p-ar + p-br + p-al [move-to patch-at -1 1]
if random-number >= p-a + p-b + p-r + p-l + p-ar + p-br + p-al and  random-number < p-a + p-b + p-r + p-l + p-ar + p-br + p-al + p-bl [move-to patch-at -1 -1]
if random-number >= p-a + p-b + p-r + p-l + p-ar + p-br + p-al + p-bl [move-to patch-at 0 0]    
  ]
  ] 
end
; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- PLOT -------------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to do-plots
  set-current-plot "stage class density"
 set-current-plot-pen "embryo"
    ifelse any? turtles with [U_H < U_H^b] [plot count turtles with [U_H < U_H^b]] 
    [plot 0]
   set-current-plot-pen "juvenile"
  ifelse any? turtles with [U_H > U_H^b and U_H < U_H^p] [plot count turtles with [U_H > U_H^b and U_H < U_H^p]]
  [plot 0]
  set-current-plot-pen "adult"
  ifelse any? turtles with [U_H >= U_H^p] [plot count turtles with [U_H >= U_H^p]]
  [plot 0]
  
 
  set-current-plot "food density"
  plot mean [X] of patches
  
  set-current-plot "population density"
 plot count turtles with [U_H > U_H^b]
  
  
  
  set-current-plot "size distribution"
  histogram [l / .054 ] of turtles with [U_H > U_H^b]

  
    set-current-plot "juv e distribution"
  histogram [e_scaled] of turtles with [U_H > U_H^b and U_H < U_H^p]
  
    set-current-plot "adult e distribution"
  histogram [e_scaled] of turtles with [U_H >= U_H^p]
end
  
@#$#@#$#@
GRAPHICS-WINDOW
1319
20
1809
531
-1
-1
5.93
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
80
0
80
1
1
1
ticks
30.0

SLIDER
279
575
391
608
f_scaled
f_scaled
0
1
0
.01
1
NIL
HORIZONTAL

BUTTON
42
70
108
103
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
42
35
108
68
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
42
104
108
137
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

PLOT
489
13
1311
294
stage class density
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"embryo" 1.0 0 -16777216 true "" ""
"juvenile" 1.0 0 -13345367 true "" ""
"adult" 1.0 0 -2674135 true "" ""

PLOT
489
532
773
727
size distribution
NIL
NIL
0.0
5.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.05 1 -16777216 true "" ""

SLIDER
110
35
243
68
timestep
timestep
0
1000
100
1
1
NIL
HORIZONTAL

MONITOR
358
34
415
79
days
ticks / timestep
1
1
11

MONITOR
275
34
359
79
NIL
count turtles\n
0
1
11

CHOOSER
279
530
458
575
food-dynamics
food-dynamics
"logistic" "constant"
0

INPUTBOX
39
467
120
527
v_rate_int
0.16
1
0
Number

INPUTBOX
39
527
120
587
kap_int
0.8
1
0
Number

INPUTBOX
39
586
120
646
kap_R_int
0.95
1
0
Number

INPUTBOX
41
203
122
263
k_M_rate_int
4
1
0
Number

INPUTBOX
40
646
120
706
k_J_rate_int
4
1
0
Number

INPUTBOX
41
263
121
323
g_int
0.15003750937734434
1
0
Number

INPUTBOX
42
323
121
383
U_H^b_int
9.993569821026685E-6
1
0
Number

INPUTBOX
42
383
121
443
U_H^p_int
4.000107169649555E-4
1
0
Number

INPUTBOX
279
680
391
755
F_m
1
1
0
Number

INPUTBOX
391
575
458
635
X_r
0.5
1
0
Number

INPUTBOX
391
635
458
695
X_k
2
1
0
Number

INPUTBOX
391
695
458
755
volume
0.01
1
0
Number

INPUTBOX
278
608
391
680
J_XAm_rate_int
1
1
0
Number

PLOT
489
292
1312
415
food density
NIL
NIL
0.0
10.0
0.0
0.0010
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

PLOT
489
414
1312
534
population density
NIL
NIL
0.0
10.0
0.0
0.0010
true
false
"" ""
PENS
"> 2.6" 1.0 0 -2674135 true "" ""

CHOOSER
300
252
417
297
aging
aging
"on" "off"
0

INPUTBOX
300
297
417
357
h_a
0.00125
1
0
Number

INPUTBOX
300
357
417
417
sG
-0.5
1
0
Number

PLOT
772
534
1053
727
juv e distribution
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.01 1 -16777216 true "" ""

PLOT
1051
534
1312
727
adult e distribution
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.01 1 -16777216 true "" ""

INPUTBOX
172
275
252
335
cv
0.1
1
0
Number

TEXTBOX
24
183
174
201
Standard DEB parameters
11
0.0
1

TEXTBOX
304
509
454
527
feeding related parameters
11
0.0
1

TEXTBOX
297
230
447
248
ageing related parameters
11
0.0
1

INPUTBOX
160
465
255
525
p_m
11200
1
0
Number

INPUTBOX
160
526
255
586
E_G
2800
1
0
Number

INPUTBOX
159
706
253
766
zoom
0.2666
1
0
Number

INPUTBOX
161
586
254
646
E_H_b
0.0373
1
0
Number

INPUTBOX
159
646
254
706
E_H_p
1.493
1
0
Number

CHOOSER
160
420
255
465
add_my_pet?
add_my_pet?
"on" "off"
0

TEXTBOX
164
256
314
274
intraspecific variation
11
0.0
1

TEXTBOX
159
183
309
201
temperature submodel
11
0.0
1

INPUTBOX
299
418
418
478
background-mortality
0.05
1
0
Number

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
NetLogo 5.0.5
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
0
@#$#@#$#@
