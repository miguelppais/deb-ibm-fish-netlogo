; Code pour simulation

extensions [profiler]

breed [ Juveniles Juvenile ]
breed [ Males Male ]
breed [ Females Female ]


patches-own  [  ]

turtles-own  [  sexe    W        puberty   Cohort                                        ; General characteristic
                L       Eh       En        lm.i
                M.t     M.c      M.n       M                                             ; Mortality
                L_mat.i  a_inh  K_dens.i   X                                              ; Maturity Length
                temp_L   PC

                 ; Individual DEB parameters
                f.i   size_food.i   R_Food_J.i R_Food_A.i F.ad_lib.i    factor.i         ; Food
                Ehp.i g.i Em.i Pam.i Phi.i Ehb.i   R
                Pam_t Phi_t Rad.i R_Ref.i

               ]

juveniles-own[ ]


females-own  [  R.max   clutch   A.clutch ]


males-own    [  Nest NbEggsMax  EggsQuantity    MyBreedingSeason  Territory
                Nb_built_nest   Eggs_duration   Time_building                            ; Nest building
                M.nid     M.nc                                                           ; Nest mortality
                Harvest Nb_F   Time_harvest  NumberFemales                               ; Egg harvest
                Nest.success   Time_Nest Time_Nest_fin   Days Time.Acc.i
                Nb_eggs]                                                                 ; End of the nest


 globals     [  ; Environmental variables
                Mesocosm    Superf.Meso       N.Cosme     Scale      file                ; Mesocosm
                Temperature Temperature_Input                                            ; Temperature
                Photoperiod Photoperiod_Input Photop.Thr                                 ; Photoperiod
                Food_J.t   Food_A.t    Food_Zoo.t                                                  ; Food
                Nauplii_Input    Copepode_Input    Cladocere_Input    Rotifere_Input     ; Zooplankton
                Gammare_inf5_Input   Gammare_supp5_Input   Aselle_inf5_Input             ; Macroinvertebrate
                Aselle_supp5_Input  Chironome_inf5_Input   Chironome_supp5_Input
                Time                                                                     ; Stepwise
                End_Experiment                                                           ; Duration of the experiment
                var.cosme   CV.c_M  CV.c_F                                               ; Inter-mesocosm variability

                ; Experimental constraint
                AdulteSeuil                                                              ; Size at sexing

                ; DEB parameters
                Phi                                                                      ; Feeding parameters
                L0 ShapeCoe alpha Kappa KappaM Pam Eg v PM Em g kM                       ; Growth parameters
                Eegg Ehb Ehp Kr Kj                                                       ; Reproduction parameters
                cv                                                                       ; Coefficient of variation
                PM_t Eg_t  Kj_t  v_t temp_Eh

                ; Food parameters
                a_food b_food  ab_food ratio_f                                           ; Kermoysan phD
                R_ref L_ref Qfood Bold.M
                Biomass K_dens                                                           ; self-inhibition parameters
                B.predator.CXX   a_Kdens Biomass.Fond.Init

                ; Temperature parameters
                Tc                                                                       ; DEB parameters
                CTM CTO CQ Zi Yi XXi Delta                                               ; Parameters from Hovel et al. 2015

                ; Mortality parameters
                Mr Mu bN  Mu_bN                                                          ; Basal mortality
                M.nn  P.Dest.nid                                                         ; Nest mortality
                a_t b_t  at_bt                                                           ; Mortality due to the temperature
                B.predator DP50_m m_dens P.Attaque m_meso                               ; Mortality density-dependant

                ; Reproduction parameters
                L_mat_founder                                                            ; Maturity Length
                a_R.max  b_R.max  A.clutch.max  Part_Rmax                                ; Clutch size parameters
                a_Nb.Egg  b_Nb.Egg                                                       ; Maximal number of eggs per male
                P.OL Breeding.Period                                                                    ; Egg mortality
                Time_harvest_max  Time_harvest_mean    Proba.Stop  Female_Select_Max
                R_min  R_diff Time_Dvpt_Eggs Time.Acc
                AreaPerMale  A.Territory.min Territoriality r.Territory Gamma.Compet     ; Male territory
                Nest.Success.Tot Nest.Built.Tot

                ;inter-year parameters
                rect_Kdens

                ;Initial physical length of founder sticklebacks
                TailleInitialeFem TailleInitialeMal                                      ; Standard length of females for each mesocosm
                TailleFem_C1 TailleFem_C2 TailleFem_C3 TailleFem_C4 TailleFem_C5
                TailleFem_C6 TailleFem_C7 TailleFem_C8 TailleFem_C9 TailleFem_C10
                TailleFem_C11 TailleFem_C12
                TailleMal_C1 TailleMal_C2 TailleMal_C3 TailleMal_C4 TailleMal_C5         ; Standard length of males for each mesocosm
                TailleMal_C6 TailleMal_C7 TailleMal_C8 TailleMal_C9 TailleMal_C10
                TailleMal_C11 TailleMal_C12

                ; Model output
                N.tot   N.tot.adult   SexRatio                                           ; Abundance and sex ratio
                F.M     F.F     F.J                                                      ; Frequency  of Males, Females and Juveniles
                N.F.C00 L.F.C00 CV.F.C00                                                 ; Number, Length and CV of female founders
                N.M.C00 L.M.C00 CV.M.C00                                                 ; Number, Length and CV of male founders
                N.F.CXX L.F.CXX CV.F.CXX F.F.CXX                                         ; Number, Length, CV and Frequency of females born in the mesocosm
                N.M.CXX L.M.CXX CV.M.CXX F.M.CXX                                         ; Number, Length, CV and Frequency of males born in the mesocosm
                N.M.m   N.M.im L.M.m   CV.M.m  F.M.m  F.M.im                             ; Number, Length, CV and Frequency of mature and immature males
                N.J     L.J     CV.J                                                     ; Number, Length and CV of  juveniles

                L.FreqF L.FreqM L.FreqJuv  L.FreqT

                 ; Toxic
                C_t C_t_input

                R_EC50 R_tox s_t.R
                gr.NF_EC50 gr.NF_tox  s_t.gr.NF   v.NF    v_t.NF
                gr.F_EC50 gr.F_tox  s_t.gr.F      v.F    v_t.F
                B.Period_EC50  B.Period_tox  s_t.B.Period   B.Period.t
                P.OL_EC50  P.OL_tox  s_t.P.OL   P.OL.t
]




;________________________________________________________________________________________________________________________________________________________________________
;                                                       Setup
;________________________________________________________________________________________________________________________________________________________________________

to Setup

  clear-all

  set-default-shape turtles "fish"

  Parameters
  Calc_Para

  Inputs

  Setup-Toxic-effect

  Setup-environment
  Setup-individual


  set NF.0 15
  set NM.0 10
  set pdt  1

  reset-ticks

end

;________________________________________________________________________________________________________________________________________________________________________
;                                                    Calculation of the parameters
;________________________________________________________________________________________________________________________________________________________________________
to Calc_Para
 ;Parameters initialization

  set CTM     Delta + CTO
  set Zi      ln( CQ ) * ( CTM - CTO )
  set Yi      ln( CQ ) * ( CTM - CTO + 2 )
  set XXi     ( Zi ^ 2 * (1 + ( 1 + 40 / Yi) ^ (0.5) ) ^ 2 ) / 400

  set Time_harvest_max (Time_harvest_mean / 100) * Time_Dvpt_Eggs                        ; in days à Topt
end

;________________________________________________________________________________________________________________________________________________________________________
;                                                    Simulation
;________________________________________________________________________________________________________________________________________________________________________
to Go

  Update-environment
  Update-Toxic-effect
 ; Move
  Update-DEB_parameters
  DEB_Juv
  DEB_F
  DEB_M
  BreedChange
  PubertyProcess
  ReproFemales
  ReproMales
  Mortality

  tick


  CalcOutputs


  if ticks = end_experiment [ stop ]

end


;________________________________________________________________________________________________________________________________________________________________________
;                                                    Procedures
;________________________________________________________________________________________________________________________________________________________________________


;_____________________________________________________________________________________________________________________________
;  I. Setup process
;_____________________________________________________________________________________________________________________________

to Setup-environment

   ; Mesocosm initialization
  ;ask patches [ set pcolor 89 + 5 + random 5 ]
  set Mesocosm     ( random N.Cosme ) + 1                                                                                                                ; Initialization of the mesocosm
  set var.cosme    random-normal 0 1                                                                                                                  ; Inter-mesocosm variability randomly initialized

 ; Environmental factors
  set Photoperiod   item Time Photoperiod_Input
  set Temperature   item Time Temperature_Input
  let Sum_Zoo       ( ( item Time Nauplii_Input ) + ( item Time Copepode_Input ) + ( item Time Cladocere_Input ) + ( item Time Rotifere_Input ) )
  let Sum_Food_J    ( ( item Time Gammare_inf5_Input ) + ( item Time Aselle_inf5_Input ) + ( item Time Chironome_inf5_Input ) )                          ; Food sum for juveniles
  let Sum_Food_A    ( ( item Time Gammare_supp5_Input ) + ( item Time Aselle_supp5_Input ) + ( item Time Chironome_supp5_Input ) )                       ; Food sum for adults
  set Food_Zoo.t    ( MiniS ( Sum_Zoo + var.cosme * (CV.c_F / 100 * Sum_Zoo) ) 0 )                                                                       ; Inter-mesocosm variability
  set Food_J.t      ( MiniS ( Sum_Food_J + var.cosme * (CV.c_F / 100 * Sum_Food_J) ) 0 )                                                                 ; Inter-mesocosm variability
  set Food_A.t      ( MiniS ( Sum_Food_A + var.cosme * (CV.c_F / 100 * Sum_Food_A) ) 0 )                                                                 ; Inter-mesocosm variability

  ; Tc temperature factor
  let Vi           ( CTM - Temperature ) / ( CTM - CTO )
  ifelse ( Temperature < CTM ) [ set Tc  Vi ^ XXi * exp( XXi * ( 1 - Vi )) ] [set Tc 0 ]


end


to Setup-Toxic-effect

  ; Toxic concentration in water
  set C_t           0

  set s_t.B.Period  0            ; stress function for nest mortality
  set s_t.P.OL      0            ; stress function for egg survival
  set s_t.R         0            ; stress function for energetic effects
  set s_t.gr.NF     0            ; stress function for juvenile growth
  set s_t.gr.F      0            ; stress function for adult growth

  set P.OL.t        P.OL                     ; Decrease of the probabiliy of eggs survival
  set B.Period.t    Breeding.Period          ; Increase of the nest mortality

  set v.NF          v                        ; Energy conductance for non founders
  set v.F           v                        ; Energy conductance for founders

end


to Setup-individual

  if(Mesocosm = 1) [ set TailleInitialeFem TailleFem_C1  set TailleInitialeMal TailleMal_C1 ]
  if(Mesocosm = 2) [ set TailleInitialeFem TailleFem_C2  set TailleInitialeMal TailleMal_C2 ]
  if(Mesocosm = 3) [ set TailleInitialeFem TailleFem_C3  set TailleInitialeMal TailleMal_C3 ]
  if(Mesocosm = 4) [ set TailleInitialeFem TailleFem_C4  set TailleInitialeMal TailleMal_C4 ]
  if(Mesocosm = 5) [ set TailleInitialeFem TailleFem_C5  set TailleInitialeMal TailleMal_C5 ]
  if(Mesocosm = 6) [ set TailleInitialeFem TailleFem_C6  set TailleInitialeMal TailleMal_C6 ]
  if(Mesocosm = 7) [ set TailleInitialeFem TailleFem_C7  set TailleInitialeMal TailleMal_C7 ]
  if(Mesocosm = 8) [ set TailleInitialeFem TailleFem_C8  set TailleInitialeMal TailleMal_C8 ]
  if(Mesocosm = 9) [ set TailleInitialeFem TailleFem_C9  set TailleInitialeMal TailleMal_C9 ]
  if(Mesocosm = 10)[ set TailleInitialeFem TailleFem_C10 set TailleInitialeMal TailleMal_C10]
  if(Mesocosm = 11)[ set TailleInitialeFem TailleFem_C11 set TailleInitialeMal TailleMal_C11]
  if(Mesocosm = 12)[ set TailleInitialeFem TailleFem_C12 set TailleInitialeMal TailleMal_C12]

  ; Individual length of the founder sticklebacks
  foreach ( n-values NF.0 [ ?1 -> ?1 ]) [ ?1 -> create-females 1 [
      set L item ?1 TailleInitialeFem
      Initialisation-Females-fondateur ]  ]

  foreach ( n-values NM.0 [ ?1 -> ?1 ] ) [ ?1 -> create-males 1 [
      set L item ?1 TailleInitialeMal
      Initialisation-Males-fondateur ] ]

   set Biomass.Fond.Init ( ( sum [W] of turtles with [(item 0 Cohort) = 0] ) / Superf.Meso )

end


;_____________________________________________________________________________________________________________________________
; II. Environment processes
;_____________________________________________________________________________________________________________________________


to Update-environment

  ; Environmental factors
  set Time         ticks                                                                                                                    ; Keep the last value of the scenario
  set Photoperiod  item Time Photoperiod_Input
  set Temperature  item Time Temperature_Input

  let Sum_Zoo       ( ( item Time Nauplii_Input ) + ( item Time Copepode_Input ) + ( item Time Cladocere_Input ) + ( item Time Rotifere_Input ) )
  let Sum_Food_J    ( ( item Time Gammare_inf5_Input ) + ( item Time Aselle_inf5_Input ) + ( item Time Chironome_inf5_Input ) )                          ; Food sum for juveniles
  let Sum_Food_A    ( ( item Time Gammare_supp5_Input ) + ( item Time Aselle_supp5_Input ) + ( item Time Chironome_supp5_Input ) )                       ; Food sum for adults
  set Food_Zoo.t    ( MiniS ( Sum_Zoo + var.cosme * (CV.c_F / 100 * Sum_Zoo) ) 0 )                                                                      ; Inter-mesocosm variability
  set Food_J.t      ( MiniS ( Sum_Food_J + var.cosme * (CV.c_F / 100 * Sum_Food_J) ) 0 )                                                                ; Inter-mesocosm variability
  set Food_A.t      ( MiniS ( Sum_Food_A + var.cosme * (CV.c_F / 100 * Sum_Food_A) ) 0 )                                                                ; Inter-mesocosm variability

  ; Tc temperature factor
  let Vi           ( CTM - Temperature ) / ( CTM - CTO )
  ifelse ( Temperature < CTM ) [ set Tc  Vi ^ XXi * exp( XXi * ( 1 - Vi )) ] [set Tc 0 ]

  ; Density-dependance for nest predation
  set B.predator      ( sum [ W ] of turtles with [ Eh >= Ehp.i ] / Superf.Meso  )                                                                      ; Biomass in mg / m^2
  set B.predator.CXX  (( sum [ W ] of turtles with [ L > 40 ] / Superf.Meso  ) - Biomass.Fond.Init)
  set P.Attaque       ( B.predator / ( DP50_m + B.predator ) )

  ; Density-dependance for the food
  set Biomass      (( ( sum [W] of turtles ) / Superf.Meso ) - Biomass.Fond.Init)                                                                      ; Biomass in mg / m^2

  ; Density-dependance for the mortality
  set N.tot count turtles
  set m_meso ( m_dens * N.tot )

  ; Reproduction : available territory
  ifelse ( any? Males with [R > R_min] ) [ set AreaPerMale (Superf.Meso / count Males with [R > R_min])][ set AreaPerMale Superf.Meso ]

  ; Acclimatation time
  ask males  with [(item 0 Cohort) = 0] [ set Time.Acc.i     ( Time.Acc.i + Tc ) ]

end



to Update-Toxic-effect

 ; Toxic concentration in water
 set C_t item Time C_t_input

ifelse ( C_t > 0 )[ set s_t.B.Period  ( 1 / (1 + exp(B.Period_tox * (ln(C_t) - (B.Period_EC50))))  )][set s_t.B.Period 0]          ; stress function for nest mortality
ifelse ( C_t > 0 )[ set s_t.P.OL  ( 1 / (1 + exp(P.OL_tox * (ln(C_t) - (P.OL_EC50))))  )][set s_t.P.OL 0]                          ; stress function for egg survival
ifelse ( C_t > 0 )[ set s_t.R     ( 1 / (1 + exp(R_tox * (ln(C_t) - (R_EC50))))  )][set s_t.R 0]                                   ; stress function for maturity effects
ifelse ( C_t > 0 )[ set s_t.gr.NF     ( 1 / (1 + exp(gr.NF_tox * (ln(C_t) - (gr.NF_EC50))))  )][set s_t.gr.NF 0]                   ; stress function for non founder growth
ifelse ( C_t > 0 )[ set s_t.gr.F     ( 1 / (1 + exp(gr.F_tox * (ln(C_t) - (gr.F_EC50))))  )][set s_t.gr.F 0]                       ; stress function for founder growth

set P.OL.t   ( P.OL / ( 1 +  s_t.P.OL ))                           ; Decrease of the probabiliy of eggs survival
set B.Period.t   Breeding.Period / ( 1 +  s_t.B.Period )           ; Increase of the nest mortality

set v.NF       v / ( 1 + s_t.gr.NF )                               ; Energy conductance for non founders
set v.F        v / ( 1 + s_t.gr.F )                                ; Energy conductance for founders

end


;_____________________________________________________________________________________________________________________________
; II. Fish processes
;_____________________________________________________________________________________________________________________________

;---------------------------------------------------INITIALISATION------------------------------------------------------------
 to Initialisation-Juveniles ;--------------------------------------------------------------------------

    ; Individual DEB parameters
    let scatter_multiplier e ^ (random-normal 0 cv)                                                     ; Log-normal distribution

    set Pam.i              ( Pam * scatter_multiplier )                                                 ; Maximum area specific assimilation rate  (J/ d / mm^2)
    set Phi.i              ( Phi * scatter_multiplier )                                                 ; Parameter to be fed ad libitum (J/ d / mm^2)
    set Ehb.i              ( Ehb * scatter_multiplier )                                                 ; Energy at hatching (J)
    set Ehp.i              ( Ehp * scatter_multiplier )                                                 ; Energy at maturity (J)
    set R_ref.i            ( R_ref * scatter_multiplier )                                               ;
    set Eh                 ( Ehb.i )                                                                    ; Initial energy for the maturity level (J)
    set Em.i               ( Pam.i / v )                                                                ; Maximum reserve density (J/mm^3)
    set g.i                Eg / ( Kappa * Em.i )                                                        ; Energy investment ratio (-)
    set En                 Em.i                                                                         ; Reserve density (J)
    set L                  L0                                                                           ; Length at hatchling (mm)
    set lm.i               Kappa * ( Pam.i / PM )                                                       ; Maximal length (mm)

    ; Other parameters
    set factor.i           ( Superf.Meso /(pi * ((L * R_ref.i / L_ref ) * 10 ^ (-3) ) ^ 2) )            ; Factor (foraging behavior)

    set size               L * 0.02                                                                     ; Proportion factor for graphing
    set W                  (ShapeCoe * L ) ^ 3                                                          ; Weight (mg)
    set M.n                Mu * W ^ bN                                                                  ; Basal mortality
    set R                  0                                                                            ; Reproduction
    set Puberty            0                                                                            ; Puberty level
    set color              grey                                                                         ; Color for graphic
    set sexe one-of        [ "male" "female" ]                                                          ; Sexe



 end


 to Initialisation-Males-fondateur ;--------------------------------------------------------------------

    ; Individual DEB parameters
    let scatter_multiplier e ^ (random-normal 0 cv)                                                     ; Log-normal distribution

    set Pam.i              ( Pam * scatter_multiplier )                                                 ; Maximum area specific assimilation rate  (J/ d / mm^2)
    set Phi.i              ( Phi * scatter_multiplier )                                                 ; Parameter to be fed ad libitum (J/ d / mm^2)
    set Eh                 ( Ehp  )                                                                     ; Initial energy for the maturity level (J)
    set Ehp.i              ( Ehp )                                                                      ; Energy at maturity (J)
    set Ehb.i              ( Ehb )                                                                      ; Energy at birth (J)
    set Em.i               ( Pam.i / v )                                                                ; Maximum reserve density (J/mm^3)
    set g.i                Eg / ( KappaM * Em.i )                                                       ; Energy investment ratio (-)
    set lm.i               KappaM * ( Pam.i / PM )                                                      ; Maximal Length (mm)
    set En                 Em.i                                                                         ; Reserve density (J)
    set f.i                1                                                                            ; Food functional response
    set R_ref.i            ( R_ref * scatter_multiplier )                                               ;

    ; Other parameters
    set Puberty            1                                                                            ; Puberty level
    set L_mat.i            L_mat_founder                                                                ; Maturity length
    set Time.Acc.i         ( Time.Acc * scatter_multiplier )
    set W                  (ShapeCoe * L) ^ 3                                                           ; Weight (mg)
    set M.n                Mu * W ^ bN                                                                  ; Basal mortality
    set factor.i           (pi * (R_ref.i * 10 ^ (-3) ) ^ 2) / Superf.Meso                              ; Factor (foraging behavior)
    set sexe               "male"                                                                       ; Sexe
    set color              blue                                                                         ; Color for gaphic
    set Cohort            [0 0 0]                                                                       ; Cohort notation for founders
    set MyBreedingSeason   1                                                                            ; Condition for male reproduction
    set Nest.success       0                                                                            ; Number of successful nest
    setxy random-xcor      random-ycor                                                                  ; Initial position in the mesocosm sAc1)

    M_Back_0                                                                                            ; Male re-initialization function for reproduction processes
    set R                  R_min                                                                        ; Reproduction parameter
    set Days               R_diff

 end


 to Initialisation-Females-fondateur ;------------------------------------------------------------------

    ; Individual DEB parameters
    let scatter_multiplier e ^ (random-normal 0 cv)                                                     ; Log-normal distribution

    set Pam.i              ( Pam * scatter_multiplier )                                                 ; Maximum area specific assimilation rate  (J/ d / mm^2)
    set Phi.i              ( Phi * scatter_multiplier )                                                 ; Parameter to be fed ad libitum (J/ d / mm^2)
    set Eh                 ( Ehp )                                                                      ; Initial energy for the maturity level (J)
    set Ehp.i              ( Ehp )                                                                      ; Energy at maturity (J)
    set Ehb.i              ( Ehb )                                                                      ; Energy at birth (J)
    set Em.i               ( Pam.i / v )                                                                ; Maximum reserve density (J/mm^3)
    set g.i                Eg / ( Kappa * Em.i )                                                        ; Energy investment ratio (-)
    set lm.i               Kappa * ( Pam.i / PM )                                                       ; Maximal Length (mm)
    set En                 Em.i                                                                         ; Reserve density (J)
    set f.i                1                                                                            ; Food functional response
    set R_ref.i            ( R_ref * scatter_multiplier )                                               ;

    ; Other parameters
    set Puberty            1                                                                            ; Puberty level
    set L_mat.i            L_mat_founder                                                                ; Maturity length
    set W                  (ShapeCoe * L) ^ 3                                                           ; Weight (mg)
    set M.n                Mu * W ^ bN                                                                  ; Basal mortality
    set factor.i           (pi * (R_ref.i * 10 ^ (-3) ) ^ 2) / Superf.Meso                              ; Factor (foraging behavior)
    set sexe               "female"                                                                     ; Sexe
    set color              orange                                                                       ; Color for gaphic
    set Cohort            [0 0 0]                                                                       ; Cohort notation for founders
    set R.max              a_R.max * ( L - L_mat.i )   + b_R.max                                        ; Maximum of eggs for females
    set R                  random-float R.max                                                           ; Initialisation of the number eggs at the beginning of the experiment
    setxy random-xcor      random-ycor                                                                  ; Initial position in the mesocosm

 end





;--------------------------------------------------------------DEB MODEL---------------------------------------------------------------------------------------
to Update-DEB_parameters ;------------------------------------------------------------------------------

  ; Parameters influencing by temperature
  set PM_t   PM * Tc                                                                                            ; Volume somatic maintenance costs (J/d/mm^3)

  set Kj_t   Kj * Tc                                                                                            ; Relative fraction of food for a maximal growth
  set kM     (PM_t / Eg)
 ; set v_t    v * Tc                                                                                            ; Energy conductance (mm /d)
  set v_t.NF   v.NF * Tc                                                                                        ; Energy conductance (mm /d)
  set v_t.F    v.F * Tc                                                                                         ; Energy conductance (mm /d)


  ask turtles  [
  set Pam_t  Pam.i * Tc                                                                                         ; Maximum area specific assimilation rate  (J/ d / mm^2)
  set Phi_t  Phi.i * Tc                                                                                         ; Parameter to be fed ad libitum (J/ d / mm^2)

  ; Available food per individual
  ifelse ( L < L_ref ) [ set Rad.i (L * R_ref.i / L_ref )  ]                                                    ; calculation of the predation distance
                       [ set Rad.i (R_ref.i) ]

  if (sexe = "male" and L >= L_mat.i)[ set Rad.i Rad.i * (1 + Bold.M) ]                                         ; aggressivity/boldness of males (King, 2013)
  set factor.i ((pi * (Rad.i * 10 ^ (-3) ) ^ 2) / Superf.Meso )                                                 ; Calculation of the factor for the water column

  ; Food density per fish
  set F.ad_lib.i    Phi_t * ( L * ShapeCoe) ^ 2                                                                 ; Food quantity to be fed ad libitum
  set size_food.i   ratio_f * ( ( b_food * ab_food ) * L + b_food )                                             ; maximum size to be eaten (Nettleships, 2012)
  ifelse ( size_food.i > 5 )[ set R_Food_J.i 1                                                                  ; food per size
                              set R_Food_A.i (( size_food.i - 5 ) / 5 )]
                            [ set R_Food_J.i ( size_food.i / 5 )
                              set R_Food_A.i 0 ]
  if ( R_Food_A.i > 1 )[ set R_Food_A.i 1 ]
  ]
                                                                                                                ; /!\ in general this case does not arrive /!\ : verification of the avalaible food in the environment
  let verif sum [( Food_Zoo.t + ( R_Food_J.i * Food_J.t ) + ( R_Food_A.i * Food_A.t )) * factor.i ] of turtles  ; Available food for all the fish
  ifelse (verif = 0) [set Qfood 0][set Qfood  (( Food_Zoo.t + Food_J.t + Food_A.t ) / ( verif ) )]              ; ratio betwwen the available food and the demand
  if  Qfood > 1 [ set Qfood 1 ]                                                                                 ; if there is enough food in the environment, there is no impact. If not, food amount per turtles is proportionnally reduced

  ask turtles  [
  let Food.t ( Food_Zoo.t + ( R_Food_J.i * Food_J.t ) + ( R_Food_A.i * Food_A.t ))                              ; All available food item
  set X (( Food.t * Qfood * factor.i ) / ( 0.5 * F.ad_lib.i ) )                                                 ; Calculation of the food density (J)
  set f.i ( X / ( 1 + X ))                                                                                      ; Calculation of the f functional response

  ; Density-dependance on growth - food competition
  set K_dens.i ( K_dens * ( W / Superf.Meso ) ^ a_Kdens )

  ifelse ((item 0  Cohort) != 0 )[set a_inh ( K_dens.i / (K_dens.i + Biomass ))][set a_inh (1)]
                                 ;[set a_inh ( K_dens.i / (K_dens.i + ( Biomass )))]

  ; set a_inh ( K_dens.i / (K_dens.i + Biomass ))
  set f.i (a_inh * f.i)
  ]


end





to DEB_Juv ;------ DEB model for juveniles -------------------------------------------------------------

 ask turtles with [ puberty = 0 ] [

   ; Initialisation
   let L_st  L * ShapeCoe                                                                              ; Structural length at the beginning of the time step

   let d_E  0
   let d_L  0
   let d_Eh 0

    ; Food
       ; Reserve density
       let nrj En / Em.i
       set d_E ( ( Pam_t / L_st ) * ( (f.i) - nrj ) ) * pdt

       ; Growth
       set temp_L  ((v_t.NF / (3 * (nrj + g.i))) * (nrj - (L_st / ( lm.i ) ))) * pdt
       ifelse( temp_L >= 0) [set d_L temp_L ] [ set d_L 0]

       ; Maturity
       let Nrj_Mat Kj_t * Eh
       set PC ((g.i * En) / (g.i + nrj)) * ((v_t.NF * L_st ^ 2) + (kM * L_st ^ 3))
       set temp_Eh  (  ( 1 - Kappa ) * PC  - Nrj_Mat  * pdt ) / ( 1 + s_t.R )
       ifelse ( Eh < Ehp.i) [ set d_Eh temp_Eh ] [ set d_Eh 0 ]

     ; Iterations
   set En    En  + d_E
   set L_st  L_st + d_L
   set Eh    Eh + d_Eh

   set L     L_st * ( 1 / ShapeCoe )
   set W     (ShapeCoe * L) ^ 3

   set size  L * 0.02                                                                                   ; Proportion factor for graphing
   ]

end



to DEB_F ;------ DEB model for females  ----------------------------------------------------------------

 ask females  with [ puberty = 1 ] [

   ; Initialisation
  let L_st  L * ShapeCoe                                                                               ; Structural length at the beginning of the time step

   let d_E  0
   let d_L  0
   let d_R  0

       ; Reserve density
       let nrj En / Em.i
       set d_E ( ( Pam_t / L_st ) * ( (f.i) - nrj ) ) * pdt

       ; Growth
       ifelse ( Cohort != [0 0 0] ) [ set temp_L  ((v_t.NF / (3 * (nrj + g.i))) * (nrj - (L_st / ( lm.i ) ))) * pdt ]
                                    [ set temp_L  ((v_t.F / (3 * (nrj + g.i))) * (nrj - (L_st / ( lm.i ) ))) * pdt ]
       ifelse( temp_L >= 0) [set d_L temp_L ] [ set d_L 0]

       ; Reproduction
       let Cost_Egg Kr / Eegg
       let Cost_Mat Kj_t * Ehp.i
       ifelse ( Cohort != [0 0 0] ) [ set PC ((g.i * En) / (g.i + nrj)) * ((v_t.NF * L_st ^ 2) + (kM * L_st ^ 3)) ]
                                    [ set PC ((g.i * En) / (g.i + nrj)) * ((v_t.F * L_st ^ 2) + (kM * L_st ^ 3)) ]
       if (photoperiod > Photop.Thr and clutch = 0)[
       ifelse ( Cohort != [0 0 0] ) [ set d_R  ( ( Cost_Egg *  ( ( 1 - Kappa ) * PC  - Cost_Mat ) )  * pdt ) / ( 1 + s_t.R ) ]
      [ set d_R  ( ( Cost_Egg *  ( ( 1 - Kappa ) * PC  - Cost_Mat ) )  * pdt ) ]
        ]

   ; Iterations
   set En    En  + d_E
   set L_st  L_st + d_L
   set R     R  + d_R

   set L     L_st * ( 1 / ShapeCoe )
   set W     ( ShapeCoe * L ) ^ 3

   set size  L * 0.02                                                                                   ; Proportion factor for graphing
   ]

end




to DEB_M ;------ DEB model for males  ------------------------------------------------------------------

 ask males with [ puberty = 1 ] [

   ; Initialisation
   let L_st  L * ShapeCoe                                                                              ; Structural length at the beginning of the time step

   let d_E  0
   let d_L  0
   let d_R  0

       ; Reserve density
       let nrj En / Em.i
       set d_E ( ( Pam_t / L_st ) * ( f.i  - nrj ) ) * pdt

       ; Growth
       ifelse ( Cohort != [0 0 0] ) [ set temp_L  ((v_t.NF / (3 * (nrj + g.i))) * (nrj - (L_st / ( lm.i ) ))) * pdt ]
                                    [ set temp_L  ((v_t.F / (3 * (nrj + g.i))) * (nrj - (L_st / ( lm.i ) ))) * pdt ]
       ifelse( temp_L >= 0) [set d_L temp_L ] [ set d_L 0]

       ; Reproduction
       let Cost_Mat Kj_t * Ehp.i
       ifelse ( Cohort != [0 0 0] ) [ set PC ((g.i * En) / (g.i + nrj)) * ((v_t.NF *  L_st ^ 2 ) + (kM * L_st ^ 3 )) ]
                                    [ set PC ((g.i * En) / (g.i + nrj)) * ((v_t.F *  L_st ^ 2 ) + (kM * L_st ^ 3 )) ]
       if (photoperiod > Photop.Thr) [
      ifelse ( Cohort != [0 0 0] ) [ set d_R  ( ( Kr *  ( ( 1 - KappaM ) * PC  - Cost_Mat ) )  * pdt ) / ( 1 + s_t.R ) ]
      [ set d_R  ( ( Kr *  ( ( 1 - KappaM ) * PC  - Cost_Mat ) )  * pdt ) ]
         ]

    ; Iterations
    set En    En  + d_E
    set L_st  L_st + d_L
    set R     R  + d_R

    set L      L_st * ( 1 / ShapeCoe )
    set W      (ShapeCoe * L) ^ 3

    set size  L * 0.02                                                                                   ; Proportion factor for graphing
 ]

end






;----------------------------------------- MOVEMENT AND PUBERTY PROCESSES --------------------------------------------------------------------
to Move;------------------------------------------------------------------------------------------------

  ask Females   [ setxy random-pxcor random-pycor ]
  ask Males     [ if (Territory = 0) [ setxy random-pxcor random-pycor  ] ]
  ask Juveniles [ setxy random-pxcor random-pycor ]

end


to BreedChange;-------Sexual classification due to the experimental constraints-------------------------

  ask Juveniles with [ L >= AdulteSeuil ][
   if(sexe = "male")[   set breed Males   set color blue   ]
   if(sexe = "female")[ set breed Females set color orange ]
   ]

end


to PubertyProcess;--------------------------------------------------------------------------------------

  ask turtles with [ L >= AdulteSeuil and Puberty = 0 and ( Eh >=  Ehp.i )  ][

      set Puberty 1
      set L_mat.i L

       if(sexe = "male") [
        set MyBreedingSeason 1
        set Nest.success     0
        M_Back_0
        set Time.Acc.i 0	
        ]

     if(sexe = "female") [
        F_Back_0	
        ]
  ]

    ask males with [lm.i = Kappa * ( Pam.i / PM ) and R >= R_min][
    set lm.i               KappaM * ( Pam.i / PM )

  ]

end



;----------------------------------------- SURVIVAL -----------------------------------------------------------------------------------------
to Mortality;------------------------------------------------------------------------------------------

   ; 1) Mortality of eggs and larvae : see the reproduction part


   ; 2) Individual mortality
   ask turtles [

     ; Basal mortality
     set M.n    ( bN * Mu_bN ) * W ^ bN ; Mu * W ^ bN
     if (breed = Males and Puberty = 1 and Nest.Success > 0) [ set M.n   M.n + Mr ]
    ]

     ; Mortality density-dependant (cannibalism and other conspecific attacks)
     ask turtles  with [(item 0 Cohort) != 0 ] [
     set M.n ( 1 - ( 1 - M.n ) / ( 1 + m_meso ))
    ]

     ask turtles [
     ; Mortality due to temperatures
     ifelse( Temperature >= 16 and Puberty = 0 )[   set M.t ( exp ( ( b_t * at_bt )* Temperature + b_t ) / ( 1 + exp ( ( b_t * at_bt )* Temperature + b_t ))) ] [ set M.t 0 ]
    ;[   set M.t ( exp ( a_t * Temperature + b_t ) / ( 1 + exp ( a_t * Temperature + b_t ))) ] [ set M.t 0 ]      ; The part due to the basal mortality is substracted

     ; Daily mortality and variability inter-mesocosm
     set M   M.n + M.t

     set M   bound ( M - var.cosme * (CV.c_M / 100 * M) ) 0 1

     ; Random draw to determine if the individual die
     if (random-float 1  < M  ) [ die ]
   ]

   ; 3) Predation on nests
   ask Males with [ Nest = 1 ] [
     set M.nc    P.Attaque * P.Dest.nid
     set M.nid   M.nn + M.nc
     if (random-float 1 < M.nid ) [ M_Back_0
    set Days R_diff ]
   ]

end




;-----------------------------------------REPRODUCTION----------------------------------------------------------------------------------
to F_Back_0;-------------------------------------------------------------------------------------------

    set A.clutch 0
    set clutch   0
    set R        0
    set R.max    ( MiniS ( a_R.max * ( L - L_mat.i )   + b_R.max ) 0 )

end



; 2.1. Reproduction des femelles-------------------------------------------------------------------------------------------------------
to ReproFemales ; -------------------------------------------------------------------------------------

  ; Reproduction des femelles
  ask Females with [ Puberty = 1] [
      if ( R >= R.max ) and (clutch = 0) [
        set clutch   1
        set A.clutch 0
        ]

      if ( clutch = 1 ) [                                                                                ; A.clutch in degree day. 1 day = 1 day * Tc
        set A.clutch A.clutch + Tc
        if ( A.clutch > A.clutch.max or R = 0 ) [ F_Back_0 ]
      ] ]

end



; 2.2. Reproduction des mâles-----------------------------------------------------------------------------------------------------------
to ReproMales ;----------------------------------------------------------------------------------------

 Get_Territory
 Build_Nest
 Harvest_Eggs
 Hatching

end

to M_Back_0 ;------------------------------------------------------------------------------------------

    set Time_building  0
    set Eggs_duration  0
    set Nest           0
    Set Time_Nest      0
    set Harvest        0
    set EggsQuantity   0
    set NumberFemales  0
    set Time_harvest   0
    set NbEggsMax      ( MiniS ( a_Nb.Egg * ( L - L_mat.i ) +  b_Nb.Egg ) 0 )
    set Days           0

end

to Get_Territory ;-------------------------------------------------------------------------------------

  if( photoperiod > Photop.Thr )  [

  ask Males with [ Puberty = 1 and Territory = 0 and MyBreedingSeason = 1 and R > R_min] [

    let List.Males.Compet Males with [ Territory = 1 and ( W <= [W] of myself - ( Gamma.Compet * [W] of myself )) ]

    ifelse (AreaPerMale >  A.Territory.min  )

        ; si le male a de la place pour faire son territoire
        [ set Territory 1  ]

        ; si le male n'a pas de la place pour faire son territoire
        [ if( any? List.Males.Compet )[
            ask one-of List.Males.Compet [ set Territory 0 ]
        set Territory 1 ] ]]

      ask Males with [ Puberty = 1 and Territory = 1 and MyBreedingSeason = 1 and R > R_min] [
         ask other males in-radius r.Territory [fd 1.2 * r.Territory ]
         set Days (Days + Tc )]                                ; Males which are too close from the territory have to step back
  ]

end

to Build_Nest;------------------------------------------------------------------------------------------

  ;; Construction du nid
  ask Males with [ Territory = 1  and Time.Acc.i >= 0 and Days > ( R_diff )] [
 ;; Si le mâle possède un territoire, contruit un nid
    if( Nest = 0 ) [
      set Nb_built_nest  ( Nb_built_nest + 1 )
            set Nest.Built.Tot ( Nest.Built.Tot + 1)
      set Nest           1
      set Time_building  ( ticks + 1 )

      set NbEggsMax      ( MiniS ( a_Nb.Egg * ( L - L_mat.i ) +  b_Nb.Egg ) 0 )   ]

  ]

;  ask Males [
;      if Nest = 1  [ ask patches in-radius  r.Territory [ set pcolor 68 ]]                              ; Color for the graph
;  ]

 end


to Harvest_Eggs;-----------------------------------------------------------------------------------------

  ask Males with [ Nest = 1 and ( ticks >= Time_building) and ( Harvest = 0)  ] [

    if ( EggsQuantity > 1 ) [ set Time_harvest Time_harvest + Tc ]                                       ; Relative age of the first harvested eggs ~ F(T°C)

     if any? (Females with [clutch = 1]) [

       set Nb_F  (( random Female_Select_Max  ) + 1)
       if( count Females with [clutch = 1] < Nb_F ) [set  Nb_F count Females with [clutch = 1] ]

       let Win max-n-of Nb_F Females with [clutch = 1] [R.max]
                                                                                                         ; Females with the maximum number of eggs are first selected
       foreach [who] of Win [ ?1 ->

         if ( EggsQuantity < NbEggsMax ) [

           set NumberFemales NumberFemales + 1
           let tmp.R Part_Rmax * [R.max] of Female ?1

           if (tmp.R > [R] of Female ?1)[ set tmp.R [R] of Female ?1 ]

           set EggsQuantity EggsQuantity + tmp.R

           ask female ?1 [ set R R - tmp.R ]
         ] ] ] ]

   ask Males with [ ( (Time_harvest >= Time_harvest_max ) or ( EggsQuantity >= NbEggsMax ) ) and ( Harvest = 0) ] [ set Harvest 1 ]  ; End of the harvest

   ask Males with [ Harvest = 1 ][ set Eggs_duration Eggs_duration + Tc ]                                ; Age of the eggs
   ask Males with [ Nest = 1 ][ set Time_Nest ( Time_Nest + Nest ) ]                                     ; Age of the nest in days

end



to Hatching;---------------------------------------------------------------------------------------------

 ask Males with [Nest = 1  and ( Eggs_duration >=  Time_Dvpt_Eggs  ) ] [

  if EggsQuantity > 0 [ set Nest.Success ( Nest.Success + 1 )
                        set Nest.Success.Tot ( Nest.Success.Tot + 1 )
     set Nb_eggs ( MiniS ( EggsQuantity * P.OL.t ) 0 )
     set EggsQuantity 0
    ]

 ifelse ( Nb_eggs > 1 )[
       let a Nest.success
       let Eggs_hatch ((random-float (0.5 * Nb_eggs )) + ( 0.5 * Nb_eggs  ))
      ;      print (who) print (Eggs_hatch) print(Nb_eggs)
      set Nb_eggs ( Nb_eggs - Eggs_hatch )
        hatch-juveniles Eggs_hatch [                                                                         ; Viability of the eggs in the nest
         Initialisation-Juveniles
         set xcor [xcor] of myself
         set ycor [ycor] of myself
         set Cohort  ( list ((item 0 Cohort) + 1) (item 2 Cohort) (a) )
    ]]
    [

      ;   Probabilité qu'un mâle arrête de fabriquer des nids
    if ( ticks > B.Period.t or (random-float 1 < ( bound ( Proba.Stop * Nest.Success ) 0 1 )) ) [
      set MyBreedingSeason 0
      set Territory        0
     ; ask patches in-radius r.Territory [ set pcolor 89 + 5 + random 5 ]
  ]


    set Time_Nest_fin Time_Nest

      M_Back_0

    ]
 ]
end


;_____________________________________________________________________________________________________________________________
; OUTPUTS : outputs of the model at the end of the simulation
;_____________________________________________________________________________________________________________________________

to CalcOutputs

  if ( ticks = End_experiment ) [

  ; Total number and frequence of fish
  set N.tot       count turtles
  set N.tot.adult count Males + count Females                                                   ; Total number of adults (Males + Females)

  ifelse(N.tot > 0)[set F.M ( count Males / N.tot ) ][set F.M "NA"]                             ; Males frequency
  ifelse(N.tot > 0)[set F.F ( count Females / N.tot ) ][set F.F "NA"]                           ; Females frequency
  ifelse(N.tot > 0)[set F.J ( count Juveniles / N.tot ) ][set F.J "NA"]                         ; Juveniles frequency

  ifelse(count Females  > 0)[set SexRatio count Males / count Females ][set SexRatio "NA"]      ; Sexratio

  ; Number and Length of Female C00
  set N.F.C00 count Females with [(item 0 Cohort) = 0]
  ifelse ( N.F.C00 = 0  )[set L.F.C00 "NA"][set L.F.C00 mean [L] of Females with [(item 0 Cohort) = 0]       ]
  ifelse ( N.F.C00 <= 1 )[set CV.F.C00 "NA"][set CV.F.C00 standard-deviation [L] of Females with [(item 0 Cohort) = 0]   / L.F.C00  ]

  ; Number and Length of Male C00
  set N.M.C00 count Males with [(item 0 Cohort) = 0]
  ifelse ( N.M.C00 = 0  )[set L.M.C00 "NA"][set L.M.C00 mean [L] of Males with [(item 0 Cohort) = 0]      ]
  ifelse ( N.M.C00 <= 1 )[set CV.M.C00 "NA"][set CV.M.C00 standard-deviation [L] of Males with [(item 0 Cohort) = 0]  / L.M.C00  ]

  ; Number, Frequency and Length of Female CXX
  set N.F.CXX count Females with [(item 0 Cohort) != 0]
  ifelse(count Females  > 0)[ set F.F.CXX N.F.CXX / count Females ][set F.F.CXX "NA"]
  ifelse ( N.F.CXX = 0  )[set L.F.CXX "NA"][set L.F.CXX mean [L] of Females with [(item 0 Cohort) != 0]       ]
  ifelse ( N.F.CXX <= 1 )[set CV.F.CXX "NA"][set CV.F.CXX standard-deviation [L] of Females with [(item 0 Cohort) != 0]   / L.F.CXX  ]

  ; Number, Frequency and Length of Male CXX
  set N.M.CXX count Males with [(item 0 Cohort) != 0]
  ifelse(count Males > 0)[set F.M.CXX N.M.CXX / count Males ][set F.M.CXX "NA"]
  ifelse ( N.M.CXX = 0  )[set L.M.CXX "NA"][set L.M.CXX mean [L] of Males with [(item 0 Cohort) != 0]       ]
  ifelse ( N.M.CXX <= 1 )[set CV.M.CXX "NA"][set CV.M.CXX standard-deviation [L] of Males with [(item 0 Cohort) != 0]    / L.M.CXX  ]

  ; Number and Length of Juveniles
  set N.J count Juveniles with [L >= 15]
  ifelse ( N.J = 0  )[set L.J "NA"  ][set L.J mean [L] of Juveniles with [L >= 15] ]  ; juveniles mean length
  ifelse ( N.J <= 1 )[set CV.J "NA" ][set CV.J standard-deviation [L] of Juveniles with [L >= 15] / mean [L] of Juveniles with [L >= 15] ]

  ; Number, Frequency and Length of mature Male
  set N.M.m  count Males with [Puberty = 1]
  set N.M.im count Males with [Puberty = 0]

  ifelse(count Males > 0)[set F.M.m N.M.m  / count Males ][set F.M.m  "NA"]
  ifelse(count Males > 0)[set F.M.im N.M.im  / count Males ][set F.M.im  "NA"]

  ifelse ( N.M.m = 0  )[set L.M.m "NA" ][set L.M.m mean [L] of Males with [Puberty = 1]      ]
  ifelse ( N.M.m <= 1 )[set CV.M.m "NA" ][set CV.M.m standard-deviation [L] of Males with [Puberty = 1]  / L.M.m ]


  set L.FreqF    Histo 80   Females 1    ; Call function Histo[Bsup Stage Step]
  set L.FreqM    Histo 80   Males 1      ; Call function Histo[Bsup Stage Step]
  set L.FreqJuv  Histo 26   Juveniles 1  ; Call function Histo[Bsup Stage Step]
  set L.FreqT    Histo 80   turtles 1    ; Call function Histo[Bsup Stage Step]


 ]

end



;_______________________________________________________________________________________________________________________________________________________________________________________________________________
; III. Model parameters
;_______________________________________________________________________________________________________________________________________________________________________________________________________________

to Inputs ; -----------------------------------------------------------------------------

    set N.Cosme 3
    set End_Experiment 211

; Condition 0 µg/L
  if ( Conc = 0 ) [

    set file "Input2012_C0.txt"

    set TailleMal_C1 [ 44.06 48.65 47.18 43.81 42.67 42.37 41.64 43.91 43.89 44.36 ] ; control mesocosm n°2
    set TailleMal_C2 [ 45.06 40.41 44.16 45.05 46.7  42.49 43.26 40.52 43.48 47.15 ] ; control mesocosm n°3
    set TailleMal_C3 [ 42.73 42.83 42.99 45.77 46.33 40.05 45.01 42.69 44.71 46.65 ] ; control mesocosm n°11

    set TailleFem_C1 [ 41.48 44.16 42.97 42.52 44.05 43.22 42.3  45    44.24 42.88 42.41 42.95 42.7  44.91 44.84 ]
    set TailleFem_C2 [ 43.41 42.54 41.07 41.67 42.11 43.08 41.42 42.46 44.76 41.49 44.52 42.29 43.59 43.89 45.61 ]
    set TailleFem_C3 [ 40    43.71 41.07 42.91 43.43 41.64 43.05 41.55 43.58 43.69 42    45.19 44.68 41.3  41.93 ]

 ]

; Condition 1 µg/L
  if ( Conc = 1 ) [
    set file "Input2012_C1.txt"

    set TailleMal_C1 [ 43.76 46.25	47.48	42.55	43.26	43.87	44.9864	44.72	43.73	45.02 ] ; mesocosm n°1
    set TailleMal_C2 [ 45.06 46.07	41.73	43.79	43.2	42.53	42.84	45.46	45.0121	48.75 ] ; mesocosm n°5
    set TailleMal_C3 [ 44.32 42.55	40.99	45.72	43.47	42.23	44.6884	51.26	48.55	49.24 ] ; mesocosm n°8

    set TailleFem_C1 [ 42.23 41.75	40.61	46.11	42.8	44.81	44.8	44.48	41.53	44.66	42.72	42.91	44.5	43.7	45.94 ]
    set TailleFem_C2 [ 41.33 40.25	41.3	41.96	43.64	44.3	42.99	44.03	41.15	42.06	42.17	43.71	40.61	41.03	42.84 ]
    set TailleFem_C3 [ 44.47 45.12	41.41	41.57	42.68	44.14	42.46	42.9	44.37	43.44	43.66	43.89	44.07	42.95	43.98 ]

 ]


; Condition 10 µg/L
    if ( Conc = 10 ) [

    set file "Input2012_C10.txt"

    set TailleMal_C1 [ 47.31 43.14	44.21	44.67	43.52	45.43	44.81	44.6686	50.22	49.61 ] ; mesocosm n°6
    set TailleMal_C2 [ 42.30 40.57	44.17	42.07	42.74	45.57	43.79	47.64	44.05	48.68 ]   ; mesocosm n°9
    set TailleMal_C3 [ 45.08	42.48	41.58	42.12	39.82	44.36	41.45	45.8232	48.52	49.24 ] ; mesocosm n°10

    set TailleFem_C1 [ 42.9	41.8	41.73	41.9	43.89	43.17	41.99	43.06	42.87	44.5	40.63	44.95	41.28	43.6	44.95   ]
    set TailleFem_C2 [ 42.43	42.75	43.31	43.46	41.24	43.23	43.08	45.31	44.64	45.37	44.71	44.43	42.89	44.45	42.91 ]
    set TailleFem_C3 [ 40.19	41.93	45.06	40.29	44.27	45.02	42.86	44.04	43.09	43.26	42.94	41.6	42.33	43.05	43.16 ]

 ]


; Condition 100 µg/L
    if ( Conc = 100 ) [

    set file "Input2012_C100.txt"

    set TailleMal_C1 [ 44.76	40.61	46.95	43.64	45.86	46.59	43.2	44.29	43.1891	44.77 ] ; mesocosm n°4
    set TailleMal_C2 [ 39.66	45.18	44.23	44.94	39.86	43.47	41.29	50.06	44.76	46.71 ]   ; mesocosm n°7
    set TailleMal_C3 [ 43.44	41.89	41.75	41.42	45.24	43.04	44.09	45.16	44.3219	47.52 ] ; mesocosm n°12

    set TailleFem_C1 [ 41.02	41.21	43.42	42.4	44.16	42.7	43.41	41.17	41.25	44.67	43.72	41.11	41.27	43.03	45.24 ]
    set TailleFem_C2 [ 42.98	41.09	42.29	41.11	44.55	41.11	42.41	40.95	43.39	44.92	42.96	41.24	42.99	41.92 42.42 ]
    set TailleFem_C3 [ 44.41	41.41	44.1	44.69	43.01	41.45	44.9	41.11	43.8	41.15	42.84	42.04	43.86	41.89	48.2  ]

]





; Temperature, photoperiod, food
     file-open file
              set Temperature_Input []  set Photoperiod_Input []
              set Nauplii_Input []  set Copepode_Input []  set Cladocere_Input []  set Rotifere_Input []                                                                        ; Zooplankton
              set Gammare_inf5_Input [] set Gammare_supp5_Input [] set Aselle_inf5_Input []  set Aselle_supp5_Input[] set Chironome_inf5_Input [] set Chironome_supp5_Input []  ; Macroinvertebrates
              set C_t_input []                                                                                                                                                  ; Concentration at time t

              while [ not file-at-end? ] [
              set Temperature_Input      sentence  Temperature_Input  (list file-read)
              set Photoperiod_Input      sentence  Photoperiod_Input (list file-read)
              set Nauplii_Input          sentence  Nauplii_Input (list file-read)
              set Copepode_Input         sentence  Copepode_Input (list file-read)
              set Cladocere_Input        sentence  Cladocere_Input (list file-read)
              set Rotifere_Input         sentence  Rotifere_Input (list file-read)
              set Gammare_inf5_Input     sentence  Gammare_inf5_Input (list file-read)
              set Gammare_supp5_Input    sentence  Gammare_supp5_Input (list file-read)
              set Aselle_inf5_Input      sentence  Aselle_inf5_Input (list file-read)
              set Aselle_supp5_Input     sentence  Aselle_supp5_Input (list file-read)
              set Chironome_inf5_Input   sentence  Chironome_inf5_Input (list file-read)
              set Chironome_supp5_Input  sentence  Chironome_supp5_Input (list file-read)
              set C_t_input              sentence  C_t_input (list file-read)

              ]
     file-close

end




to Parameters ; -------------------------------------------------------------------------

  ; Toxic parameters
;  set B.Period_EC50     14.14025                                                        ; Literature
;  set B.Period_tox     -0.12004
;  set P.OL_EC50         6.37678
;  set P.OL_tox          -2.16712
;  set R_EC50            6.66896                                                         ; Concentration without an effect
;  set R_tox             -0.7294                                                         ; Parameter for calculating the effect
;  set gr.NF_EC50         12.44867
;  set gr.NF_tox          -0.467161
;  set gr.F_EC50         10.09582
;  set gr.F_tox          -0.837434

  set B.Period_EC50     14.775                                                           ; Calibration
  set B.Period_tox      (- 0.109 )
  set P.OL_EC50         6.713
  set P.OL_tox          (- 1.996 )
  set R_EC50            0.501
  set R_tox             (-0.682)
  set gr.NF_EC50        7.024
  set gr.NF_tox         -0.431
  set gr.F_EC50         6.798
  set gr.F_tox         -0.894

  ; inter-year parameter
  set rect_Kdens        0.769

  ; Simulation parameters
  set AdulteSeuil       26                                                               ; Length at which individuals are sexed (mm) (normally it is 25 mm, but fixed to 26 mm to avoid error)
  set L0                5.72                                                             ; Standard length of juveniles at hatching (mm), INERIS

  ; Growth parameters
  set Kappa             0.756645                                                         ; DEB kappa parameter (-)
  set alpha             0.110948                                                         ; Kappa-male = kappa - alpha after puberty (-) calibrated with mesocosm data of 2010
  set Pam               2.42158                                                          ; Maximum area specific assimilation rate  (J/ d / mm^2)
  set v                 1.33408                                                          ; Energy conductance (mm /d)
  set PM                0.1186; 0.111384                                                      ; Volume somatic maintenance costs (J/d/mm^3)
  set Eg                1.09809                                                          ; Cost of synthesis of a unit of structure   (J/mm^3)
  set ShapeCoe          0.249668                                                         ; Weight (mg) W = ShapeCoe * L ^3
  set cv                0.061                                                            ; Variation coefficient of DEB parameters

  ; Reproduction parameters
  set Eegg              5.61405                                                          ; Energy in an eggs (J)
  set Ehb               1.32978                                                          ; Energy at hatching (J)
  set Ehp               374.79                                                           ; Energy at puberty (J)
  set Kj                0.00303578                                                       ; Maturity maintenance rate coefficient (d^-1)
  set Kr                0.978136                                                         ; Fraction of reproduction energy fixed in eggs
  set a_R.max           5.37                                                             ; Parameter to calculate the size of the clutch
  set b_R.max           16.8                                                             ; //
  set P.OL              0.8982                                                           ; Probability of egg survival
  set A.clutch.max      2                                                                ; Maximum period for which females can keep eggs before spawning (°C/J)
  set Proba.Stop        0.09                                                             ; Probability of stopping the reproduction processes in function of the number of successful eggs
  set Breeding.Period   117.23                                                           ; Duration time of the breeding period (d)
  set Part_Rmax         0.703                                                            ; Part of R.max for the female to spawn in a given nest
  set a_Nb.Egg          58.8                                                             ; Parameter to calculate the maximal number of harvested eggs by males
  set b_Nb.Egg          ( - 437 )                                                        ; //
  set Time_harvest_mean 12.7                                                             ; Maximal relative duration of the harvest of eggs by males (% cycle) : 2 - 3 days
  set Female_Select_Max 4                                                                ; Maximal number of chosen females by males to spawn in their nest
  set Time_Dvpt_Eggs    6.9                                                              ; Time to hatch at Topt (J)
  set Time.Acc          ( - 7 )                                                          ; Acclimitation time in mesocosms for males before starting the reproduction processes
  set R_min             1586.565                                                         ; Threshold for the male reproduction processes
  set R_diff            11.14                                                            ; Time between two nests (days)
  set L_mat_founder     33.01                                                            ; Length for a mature individual, value for the founders

  ; Temperature parameters
  set Delta             2
  set CTO               23
  set CQ                3

  ; Photoperiod
  set Photop.Thr        11.305                                                           ; Day duration EC50essary to start the reproduction cycle

  ; Food parameters
  set Phi               15.3058                                                          ; Relative fraction of food for a maximal growth
  set b_food            ( - 0.743 )                                                      ; Stickleback length-mouth (From De Kermoysan)
  set ab_food           ( - 0.400 )
  set a_food            b_food  * ab_food                                                ; Stickleback length-mouth (From De Kermoysan)
  set ratio_f           0.6                                                              ; Ratio between lengths of the prey and the mouth from Gill & Hart (1994)
  set R_ref             57.60                                                            ; Reference radius to find its prey (mm)
  set L_ref             23.54                                                            ; Reference length to find its prey (mm)
  set K_dens            ( 8702.16 * rect_Kdens )                                         ; Parameter of density-dependence to calculate a_inh
  let ab_Kdens          5.02E-05
  set a_Kdens           ( K_dens * ab_Kdens )                                            ; Parameter of density-dependence to calculate a_inh
  set Bold.M            0.169                                                            ; Aggressivity / boldness of males during foraging

  ; Environmntal parameters
  set Superf.Meso       20                                                               ; Size of a mesocosm (in m²)
  set scale             20 / 100                                                         ; Size of a patch
  set CV.c_M            23.92                                                            ; Coefficient of variation of the mortality ( inter-mesocosm variability in % )
  set CV.c_F            14.94                                                            ; Coefficient of variation of the food biomass inter-mesocosm ( inter-mesocosm variability in % )

  ; Territory parameters
  set A.Territory.min   0.218                                                            ; Minimal size of a territory for a male
  set Gamma.Compet      0.15                                                             ; Percentage of the length
  set r.Territory       (( A.Territory.min / pi ) ^ 0.5 ) * 1 / scale

  ; Mortality
  set Mr                0.000856                                                         ; Daily malus for males due to the reproduction processes
  set Mu                0.00323                                                          ; Natural mortality rate at unit weight
  set bN                ( - 0.0516 )                                                     ; Allometric scaling factor
  set Mu_bN             ( Mu / bN )
  set DP50_m            4000                                                             ; 50% of the density-dependance (mg/m²)
  set m_dens            8.60E-6                                                          ; Parameter for the density-dependant mortality (Median from Hazlerigg, 2012)
  set M.nn              0.0191                                                           ; Daily mortality of the nests (Data Louis)
  set P.Dest.nid        0.01                                                             ; Probability for a nest to be destroyed during an attack
  set a_t               3.64                                                             ; Parameter for the mortality due to the temperature
  set b_t               -93.7                                                            ; //
  set at_bt             ( a_t / b_t )

  set KappaM  Kappa - alpha
  set kM      (PM / Eg)

end







;_____________________________________________________________________________________________________________________________
; Reporter
;_____________________________________________________________________________________________________________________________

to-report L.Moy [Sx Co]
  ifelse count turtles with [breed = Co and sexe = Sx]  > 0
    [ report mean [L] of turtles with [breed = Co  and sexe = Sx]  ]
    [ report 0 ]
end

to-report W.Moy [Sx Co]
  ifelse count turtles with [breed = Co and sexe = Sx]  > 0
    [ report mean [W] of turtles with [breed = Co  and sexe = Sx]  ]
    [ report 0 ]
end

to-report W.Co [Sx Co]
  ifelse count turtles with [Cohort = Co and sexe = Sx]  > 0
    [ report mean [W] of turtles with [Cohort = Co  and sexe = Sx]  ]
    [ report 0 ]
end

to-report L.Co [Sx Co]
  ifelse count turtles with [Cohort = Co and sexe = Sx]  > 0
    [ report mean [L] of turtles with [Cohort = Co  and sexe = Sx]  ]
    [ report 0 ]
end

to-report bound [Nbr LowB UppB]                   ; Apply boundary to a number
  let tmp Nbr
  if( Nbr >= UppB ) [ set tmp UppB ]
  if( Nbr <= LowB ) [ set tmp LowB ]
  report tmp
end

to-report MiniS [Nbr LowBr]                       ; Apply boundary to a number
  let tmp Nbr
  if( Nbr < LowBr ) [ set tmp LowBr ]
  report tmp
end

to-report MaxiS [Nbr HigBr]                       ; Apply boundary to a number
  let tmp Nbr
  if( Nbr > HigBr ) [ set tmp HigBr ]
  report tmp
end


to-report Graph.Ind.R [Nb]
  ifelse any? turtles with [Who = Nb]
    [ report [r] of turtles with [Who = Nb]  ]
    [ report [0] ]
end

to-report Graph.Ind.F [Nb]
  ifelse any? turtles with [Who = Nb]
    [ report [f.i] of turtles with [Who = Nb]  ]
    [ report [0] ]
end

to-report Graph.Ind.M [Nb]
  ifelse any? turtles with [Who = Nb]
    [ report [1 - M] of turtles with [Who = Nb]  ]
    [ report [1.01] ]
end

to-report random-lognormal [moy sd]
  let Beta_ln    ln(1 + (sd ^ 2 / moy ^ 2))
  let Sigma_ln   sqrt Beta_ln
  let Mu_ln      ln(moy) - (Beta_ln / 2)

  report exp( random-normal Mu_ln Sigma_ln )
end

to-report Histo [UppB Stage Step]                 ; Histogram of individual length of values from 0 to Bsup for selected stade (Juvenile, Male, Female) with a step of "step" mm
  let temp []
  let Inter n-values UppB [ ?1 -> ?1 ]
  foreach Inter [ ?1 -> set temp lput (frequency (?1) ([L] of Stage) (Step) ) temp ]
  report temp
end

to-report frequency [val thelist Step]            ; Frequency of values between val and val+Step
  report length filter [i -> ( i >= val ) and (i < (val + Step))] thelist
end
@#$#@#$#@
GRAPHICS-WINDOW
37
325
1645
414
-1
-1
8.0
1
1
1
1
1
0
1
1
1
0
199
0
9
1
1
1
ticks
50.0

BUTTON
44
25
108
58
Setup
Setup
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
109
24
172
57
Go
Go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
52
183
224
216
NM.0
NM.0
0
100
10.0
1
1
NIL
HORIZONTAL

CHOOSER
99
74
237
119
Conc
Conc
0 1 10 100
3

SLIDER
50
241
222
274
pdt
pdt
0.1
1
1.0
0.05
1
NIL
HORIZONTAL

BUTTON
174
24
237
57
pdt
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

MONITOR
33
73
92
118
Mesocosm
Mesocosm
17
1
11

SLIDER
53
142
226
175
NF.0
NF.0
0
100
15.0
1
1
NIL
HORIZONTAL

BUTTON
277
29
349
63
Profiler
setup                  ;; set up the model\nprofiler:start         ;; start profiling\nrepeat End_experiment [ go ]       ;; run something you want to measure\nprofiler:stop          ;; stop profiling\nprint profiler:report  ;; view the results\nprofiler:reset         ;; clear the data
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
392
31
763
303
Length
NIL
NIL
0.0
70.0
0.0
100.0
true
false
"" ""
PENS
"default" 1.0 0 -8431303 true "" "histogram [L] of turtles"

PLOT
805
35
1179
301
Growth
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot [L] of turtle 1"
"pen-1" 1.0 0 -7500403 true "" "plot [L] of turtle 20"
"pen-2" 1.0 0 -2674135 true "" "plot [L] of turtle 45"
"pen-3" 1.0 0 -955883 true "" "plot [L] of turtle 50"
"pen-4" 1.0 0 -6459832 true "" "plot [L] of turtle 200"
"pen-5" 1.0 0 -14439633 true "" "plot [L] of turtle 250"

@#$#@#$#@
## WHAT IS IT?

blabla blablabla bla !!!!! Et bla !
## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
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
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

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

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Simulation" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>N.tot</metric>
    <metric>N.F.C00</metric>
    <metric>N.M.C00</metric>
    <metric>F.J</metric>
    <metric>SexRatio</metric>
    <metric>F.F</metric>
    <metric>F.M.m</metric>
    <metric>L.F.C00</metric>
    <metric>L.M.C00</metric>
    <metric>L.F.CXX</metric>
    <metric>L.M.CXX</metric>
    <metric>L.M.m</metric>
    <metric>L.J</metric>
    <metric>CV.F.C00</metric>
    <metric>CV.M.C00</metric>
    <metric>CV.F.CXX</metric>
    <metric>CV.M.CXX</metric>
    <metric>CV.M.m</metric>
    <metric>CV.J</metric>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Conc">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Devalaison" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
Calc_Para</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>count Juveniles with [ L &lt; 10 ]</metric>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2014"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Taille" repetitions="100" runMetricsEveryStep="false">
    <setup>Setup
Calc_Para</setup>
    <go>go</go>
    <metric>max [L] of turtles with [ Cohort != [0 0 0]]</metric>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Conc">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Nid" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>[ Nest] of males with [ Cohort = [ 0 0 0 ]]</metric>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Conc">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Nbr_temps" repetitions="50" runMetricsEveryStep="true">
    <setup>setup
Calc_Para</setup>
    <go>go</go>
    <metric>N.tot</metric>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2016"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Maturite" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>[L] of Males with [Puberty = 1]</metric>
    <metric>[L] of Males with [Puberty = 0]</metric>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2011"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Number_Time" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2011"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="SimulationBPA" repetitions="100" runMetricsEveryStep="false">
    <setup>setup
Calc_Para</setup>
    <go>go</go>
    <metric>N.tot</metric>
    <metric>count Males</metric>
    <metric>count Females</metric>
    <enumeratedValueSet variable="Year">
      <value value="2012"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Year">
      <value value="2012"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="11"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Frequence" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup
Calc_Para</setup>
    <go>go</go>
    <metric>L.FreqT</metric>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Conc">
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
