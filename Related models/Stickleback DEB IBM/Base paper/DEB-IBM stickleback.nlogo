extensions [profiler]

breed [ Juveniles Juvenile ]
breed [ Males Male ]
breed [ Females Female ]


patches-own  [  ]

turtles-own  [  sexe    W        puberty   Cohort                                        ; General characteristic
                L       Eh       En        lm.i
                M.t     M.c      M.n       M                                             ; Mortality
                L_mat.i  a_inh  K_dens.i   X                                              ; Maturity Length


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
                Nb_eggs Ratio_Eggs]                                                                 ; End of the nest


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
                PM_t Eg_t  Kj_t  v_t

                ; Food parameters
                a_food b_food  ab_food ratio_f                                           ; Kermoysan phD
                R_ref L_ref Qfood Bold.M
                Biomass K_dens                                                           ; self-inhibition parameters
                a_Kdens Biomass.Fond.Init

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
                P.OL   Breeding.Period                                                   ; Egg mortality
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
  Setup-environment
  Setup-individual

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
  Move
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
  ask patches [ set pcolor 89 + 5 + random 5 ]
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


;_____________________________________________________________________________________________________________________________
; II. Environment processes
;_____________________________________________________________________________________________________________________________


to Update-environment

  ; Environmental factors
  set Time         ticks
  let tmp length   Temperature_Input - 1
  if ( Time > tmp )[ set Time tmp ]                                                                                                                      ; Keep the last value of the scenario
                                                                                                                                                         ; Keep the last value of the scenario
  set Photoperiod  item Time Photoperiod_Input
  set Temperature  item Time Temperature_Input

  let Sum_Zoo      ( ( item Time Nauplii_Input ) + ( item Time Copepode_Input ) + ( item Time Cladocere_Input ) + ( item Time Rotifere_Input ) )
  let Sum_Food_J   ( ( item Time Gammare_inf5_Input ) + ( item Time Aselle_inf5_Input ) + ( item Time Chironome_inf5_Input ) )                          ; Food sum for juveniles
  let Sum_Food_A   ( ( item Time Gammare_supp5_Input ) + ( item Time Aselle_supp5_Input ) + ( item Time Chironome_supp5_Input ) )                       ; Food sum for adults
  set Food_Zoo.t    ( MiniS ( Sum_Zoo + var.cosme * (CV.c_F / 100 * Sum_Zoo) ) 0 )                                                                       ; Inter-mesocosm variability
  set Food_J.t      ( MiniS ( Sum_Food_J + var.cosme * (CV.c_F / 100 * Sum_Food_J) ) 0 )                                                                 ; Inter-mesocosm variability
  set Food_A.t      ( MiniS ( Sum_Food_A + var.cosme * (CV.c_F / 100 * Sum_Food_A) ) 0 )                                                                 ; Inter-mesocosm variability

  ; Tc temperature factor
  let Vi           ( CTM - Temperature ) / ( CTM - CTO )
  ifelse ( Temperature < CTM ) [ set Tc  Vi ^ XXi * exp( XXi * ( 1 - Vi )) ] [set Tc 0 ]

  ; Density-dependance for nest predation
  set B.predator  ( sum [ W ] of turtles with [ Eh >= Ehp.i ] / Superf.Meso  )                                                                          ; Biomass in mg / m^2
  set P.Attaque    ( B.predator / ( DP50_m + B.predator ) )

  ; Density-dependance for the food
  set Biomass      (( ( sum [W] of turtles ) / Superf.Meso )  - Biomass.Fond.Init)                                           ; Biomass in mg / m^2

  ; Density-dependance for the mortality
  set N.tot count turtles
  set m_meso ( m_dens * N.tot )

  ; Reproduction : available territory
    ifelse ( any? Males with [R > R_min] ) [ set AreaPerMale (Superf.Meso / count Males with [R > R_min])][ set AreaPerMale Superf.Meso ]

  ; Acclimatation time
  ask males  with [(item 0 Cohort) = 0] [set Time.Acc.i     ( Time.Acc.i + Tc )]
  ask males  with [(item 0 Cohort) = 0 and NbEggsMax > 0] [set Ratio_eggs (EggsQuantity / NbEggsMax)]

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
    set factor.i           ( Superf.Meso /(pi * ((L * R_ref.i / L_ref ) * 10 ^ (-3) ) ^ 2) )              ; Factor (foraging behavior)

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
    set factor.i           (pi * (R_ref.i * 10 ^ (-3) ) ^ 2) / Superf.Meso                                ; Factor (foraging behavior)
    set sexe               "male"                                                                       ; Sexe
    set color              blue                                                                         ; Color for gaphic
    set Cohort            [0 0 0]                                                                       ; Cohort notation for founders
    set MyBreedingSeason   1                                                                            ; Condition for male reproduction
    set Nest.success       0                                                                            ; Number of successful nest
    setxy random-xcor      random-ycor                                                                  ; Initial position in the mesocosm sAc1)

    M_Back_0                                                                                            ; Male re-initialization function for reproduction processes
    set R                  R_min                                                                        ; Reproduction parameter
    set Days               R_diff
    set Ratio_eggs            0
 end


 to Initialisation-Females-fondateur ;------------------------------------------------------------------

    ; Individual DEB parameters
    let scatter_multiplier               e ^ (random-normal 0 cv)                                       ; Log-normal distribution

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
    set factor.i           (pi * (R_ref * 10 ^ (-3) ) ^ 2) / Superf.Meso                                ; Factor (foraging behavior)
    set sexe               "female"                                                                     ; Sexe
    set color              orange                                                                       ; Color for gaphic
    set Cohort            [0 0 0]                                                                       ; Cohort notation for founders
    set R.max              a_R.max * ( L - L_mat.i )   + b_R.max                                          ; Maximum of eggs for females
    set R                  random-float R.max                                                           ; Initialisation of the number eggs at the beginning of the experiment
    setxy random-xcor      random-ycor                                                                  ; Initial position in the mesocosm

 end




;--------------------------------------------------------------DEB MODEL---------------------------------------------------------------------------------------
to Update-DEB_parameters ;------------------------------------------------------------------------------

  ; Parameters influencing by temperature
  set PM_t   PM * Tc                                                                                            ; Volume somatic maintenance costs (J/d/mm^3)
  set Kj_t   Kj * Tc                                                                                            ; Relative fraction of food for a maximal growth
  set kM     (PM_t / Eg)
  set v_t    v * Tc                                                                                             ; Energy conductance (mm /d)


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
  ifelse (verif = 0) [set Qfood 0][set Qfood (( Food_Zoo.t + Food_J.t + Food_A.t ) / ( verif ) )]               ; ratio betwwen the available food and the demand
  if  Qfood > 1 [ set Qfood 1 ]                                                                                 ; if there is enough food in the environment, there is no impact. If not, food amount per turtles is proportionnally reduced

  ask turtles  [
  let Food.t ( Food_Zoo.t + ( R_Food_J.i * Food_J.t ) + ( R_Food_A.i * Food_A.t ))                              ; All available food item
  set X (( Food.t * Qfood * factor.i ) / ( 0.5 * F.ad_lib.i ) )                                                 ; Calculation of the food density (J)
  set f.i ( X / ( 1 + X ))                                                                                      ; Calculation of the f functional response

  ; Density-dependance on growth - food competition
  set K_dens.i ( K_dens * ( W / Superf.Meso ) ^ a_Kdens )

  ifelse ((item 0  Cohort) != 0 )[set a_inh ( K_dens.i / (K_dens.i + Biomass ))][set a_inh (1)]

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
       let temp_L  ((v_t / (3 * (nrj + g.i))) * (nrj - (L_st / lm.i))) * pdt
       ifelse( temp_L >= 0) [set d_L temp_L ] [ set d_L 0]

       ; Maturity
       let Nrj_Mat Kj_t * Eh
       let PC ((g.i * En) / (g.i + nrj)) * ((v_t * L_st ^ 2) + (Km * L_st ^ 3))
       let temp_Eh  (  ( 1 - Kappa ) * PC  - Nrj_Mat ) * pdt
       ifelse ( Eh < Ehp.i) [ set d_Eh temp_Eh ] [ set d_Eh 0 ]

     ; Iterations
   set En  En  + d_E
   set L_st  L_st  + d_L
   set Eh Eh + d_Eh

   set L  L_st * ( 1 / ShapeCoe )
   set W       (ShapeCoe * L) ^ 3

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
       let temp_L  ((v_t / (3 * (nrj + g.i))) * (nrj - (L_st / lm.i))) * pdt
       ifelse( temp_L >= 0) [set d_L temp_L ] [ set d_L 0]

       ; Reproduction
       let Cost_Egg Kr / Eegg
       let Cost_Mat Kj_t * Ehp.i
       let PC ((g.i * En) / (g.i + nrj)) * ((v_t * L_st ^ 2) + (Km * L_st ^ 3))
       if (photoperiod > Photop.Thr and clutch = 0)[
           set d_R  ( ( Cost_Egg *  ( ( 1 - Kappa ) * PC  - Cost_Mat ) )  * pdt )
    ]

   ; Iterations
   set En  En  + d_E
   set L_st  L_st  + d_L
   set R  R  + d_R

   set L  L_st * ( 1 / ShapeCoe )
   set W       (ShapeCoe * L) ^ 3

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
       let temp_L  ((v_t / (3 * (nrj + g.i))) * (nrj - (L_st / lm.i ))) * pdt
       ifelse( temp_L >= 0) [set d_L temp_L ] [ set d_L 0]

       ; Reproduction
       ;let Cost_Egg Kr / Eegg
       let Cost_Mat Kj_t * Ehp.i
       let PC ((g.i * En) / (g.i + nrj)) * ((v_t *  L_st ^ 2 ) + (Km * L_st ^ 3 ))
       if (photoperiod > Photop.Thr) [
       set d_R  ( ( Kr *  ( ( 1 - KappaM ) * PC  - Cost_Mat ) )  * pdt )
         ]

    ; Iterations
    set En  En  + d_E
    set L_st  L_st  + d_L
    set R  R  + d_R

   set L  L_st * ( 1 / ShapeCoe )
   set W       (ShapeCoe * L) ^ 3

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
    			set lm.i             KappaM * ( Pam.i / PM )
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
     set M  bound ( M - var.cosme * (CV.c_M / 100 * M) ) 0 1

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

    ifelse ( AreaPerMale >  A.Territory.min )

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
  ask Males with [ Territory = 1  and Time.Acc.i >= 0 and Days > R_diff ] [
 ;; Si le mâle possède un territoire, contruit un nid
    if( Nest = 0 ) [
      set Nb_built_nest  ( Nb_built_nest + 1 )
            set Nest.Built.Tot ( Nest.Built.Tot + 1)
      set Nest           1
      set Time_building  ( ticks + 1 )

      set NbEggsMax      ( MiniS ( a_Nb.Egg * ( L - L_mat.i ) +  b_Nb.Egg ) 0 )   ]

  ]

  ask Males [
      if Nest = 1  [ ask patches in-radius  r.Territory [ set pcolor 68 ]]                              ; Color for the graph
  ]

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

   ask Males with [ ( ( Time_harvest >= Time_harvest_max ) or ( EggsQuantity >= NbEggsMax ) ) and ( Harvest = 0 ) ] [ set Harvest 1 ]  ; End of the harvest

   ask Males with [ Harvest = 1 ][ set Eggs_duration Eggs_duration + Tc ]                                ; Age of the eggs
   ask Males with [ Nest = 1 ][ set Time_Nest ( Time_Nest + Nest ) ]                                     ; Age of the nest in days

end



to Hatching;---------------------------------------------------------------------------------------------

 ask Males with [Nest = 1  and ( Eggs_duration >=  Time_Dvpt_Eggs  ) ] [

  if EggsQuantity > 0 [ set Nest.Success ( Nest.Success + 1 )
                        set Nest.Success.Tot ( Nest.Success.Tot + 1 )
     set Nb_eggs ( MiniS ( EggsQuantity * P.OL ) 0 )
     set EggsQuantity 0
    ]

 ifelse ( Nb_eggs > 1 )[
       let a Nest.success
       let Eggs_hatch ( (random-float (0.5 * Nb_eggs )) + ( 0.5 * Nb_eggs  ))
        ;      print (who) print (Eggs_hatch) print(Nb_eggs)
      set Nb_eggs ( Nb_eggs - Eggs_hatch )
        hatch-juveniles Eggs_hatch [                                                                         ; Viability of the eggs in the nest
         Initialisation-Juveniles
         set xcor [xcor] of myself
         set ycor [ycor] of myself
         set Cohort  ( list ((item 0 Cohort) + 1) (item 2 Cohort) (a) )
    ]]
    [

      ; Probabilité qu'un mâle arrête de fabriquer des nids
    if ( ticks > Breeding.Period or (random-float 1 < ( bound ( Proba.Stop * Nest.Success ) 0 1 )) ) [
      set MyBreedingSeason 0
      set Territory        0
      ask patches in-radius r.Territory [ set pcolor 89 + 5 + random 5 ]
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
  set N.J count Juveniles  with [L >= 15]
  ifelse ( N.J = 0  )[set L.J "NA"  ][set L.J mean [L] of Juveniles with [L >= 15] ]  ; juveniles mean length
  ifelse ( N.J <= 1 )[set CV.J "NA" ][set CV.J standard-deviation [L] of Juveniles with [L >= 15] / mean [L] of Juveniles with [L >= 15] ]

  ; Number, Frequency and Length of mature Male
  set N.M.m count Males with [Puberty = 1 and (item 0 Cohort) != 0]
  set N.M.im count Males with [Puberty = 0 and (item 0 Cohort) != 0]

  ifelse(count Males > 0)[set F.M.m N.M.m  / count Males ][set F.M.m  "NA"]
  ifelse(count Males > 0)[set F.M.im N.M.im  / count Males ][set F.M.im  "NA"]

  ifelse ( N.M.m = 0  )[set L.M.m "NA" ][set L.M.m mean [L] of Males with [Puberty = 1]      ]
  ifelse ( N.M.m <= 1 )[set CV.M.m "NA" ][set CV.M.m standard-deviation [L] of Males with [Puberty = 1]  / L.M.m ]


  set L.FreqF    Histo 80 Females 1    ; Call function Histo[Bsup Stage Step]
  set L.FreqM    Histo 80 Males 1      ; Call function Histo[Bsup Stage Step]
  set L.FreqJuv  Histo 26 Juveniles 1  ; Call function Histo[Bsup Stage Step]
  set L.FreqT    Histo 80 turtles 1     ; Call function Histo[Bsup Stage Step]



  ]
end



;_______________________________________________________________________________________________________________________________________________________________________________________________________________
; III. Model parameters
;_______________________________________________________________________________________________________________________________________________________________________________________________________________

to Inputs ; -----------------------------------------------------------------------------


  if ( Year = 2010) [                                                                     ; YEAR 2010

    set N.Cosme 11
    set file "Input2010.txt"

    set TailleMal_C1  [ 50.77 49.52 47.54 52.04 50.29 46.95 51.18 50.81 50.09 49.62 ]
    set TailleMal_C2  [ 49.11 52.62 52.13 45.67 51.71 53.07 49.92 52.47 51.16 47.04 ]
    set TailleMal_C3  [ 50.05 50.96 48.76 49.06 47.65 51.27 48.5  50.82 45.41 49.44 ]
    set TailleMal_C4  [ 46.38 50.05 54.56 48.37 48.67 47.65 49.24 47.95 48.64 47.47 ]
    set TailleMal_C5  [ 49.25 48.94 50.46 48.8  45.09 51.5  47.19 50.87 49    49.18 ]
    set TailleMal_C6  [ 46.92 49.65 48.69 49.72 50.91 45.83 52.6  50.56 48.51 45.72 ]
    set TailleMal_C7  [ 49    48.89 47.59 51.61 45.95 48.61 48.44 51.68 44.12 49.17 ]
    set TailleMal_C8  [ 48.68 48.66 46.92 49.49 46.33 45.81 51.27 49.2  50.67 44.89 ]
    set TailleMal_C9  [ 48.11 47.6  48.98 45.57 48.54 50.49 49    51.83 50.76 45.42 ]
    set TailleMal_C10 [ 48.68 51.15 47.4  51.71 49    49.2  48.85 49.96 56    47.36 ]
    set TailleMal_C11 [ 51.41 48.87 47.7  47.64 50.56 49.55 54.13 47.75 48.29 48.11 ]

    set TailleFem_C1  [ 46.76 48.92 50.38 46.48 49.4  51.09 46    57.29 57.33 55.89 54.11 57.88 57.99 53.25 56.26 ]
    set TailleFem_C2  [ 50.86 45.61 49.6  45.99 45.37 51.29 48.25 55.13 54.57 57.34 51.75 56.22 56.4  54.74 53.71 ]
    set TailleFem_C3  [ 50.32 50.69 46.66 48.54 46.85 48.17 51.98 50    58.16 53.28 52.87 51.66 53    55.28 57.34 ]
    set TailleFem_C4  [ 49.3  50.67 46.41 53.02 46.48 50.47 48.17 51.79 53.16 53.63 53.81 56.77 48.86 50.91 50.59 ]
    set TailleFem_C5  [ 47.28 51.81 53.16 46.83 47.06 51.63 47.13 51.45 49.87 55.12 51.7  58.79 52.28 53.53 51.2  ]
    set TailleFem_C6  [ 52.64 53.64 51.91 49.94 47.93 50.76 49.69 53.57 52.85 55.71 53.7  54.78 51.68 52.5  47.26 ]
    set TailleFem_C7  [ 47.33 49.67 50.21 48.42 48.95 48.93 46.87 52.69 51.63 57.16 51.14 55.52 55.24 52.61 55.64 ]
    set TailleFem_C8  [ 46.87 47.68 48.1  48.01 48.1  45.76 46.55 56.58 49.98 54.49 51.83 48.85 53.47 53.37 53.25 ]
    set TailleFem_C9  [ 49.37 48.94 45    48.84 48.25 49.31 48.35 55.9  56.44 54.29 58.1  57.94 58.04 51.3  52.02 ]
    set TailleFem_C10 [ 46.3  49.28 46.49 48.1  49.44 48.36 49.05 55.13 55.82 58.25 54.89 49.42 50    51.86 55.65 ]
    set TailleFem_C11 [ 46.01 51.98 54.84 46.79 47.79 50.12 46.74 54.71 51.43 50.27 51.17 54.97 55.28 55.03 56.44 ]

    set End_Experiment 211 ;; days

  ]


  if ( Year = 2011) [                                                                    ; YEAR 2011

    set N.Cosme 11
    set file "Input2011.txt"

    set TailleMal_C1  [ 41.54 44.34 42.46 41.63 41.46 43.74 42.75 44.51 39.18 36.35 ]
    set TailleMal_C2  [ 43.12 45.12 42.3  41.23 43.84 44.15 41.37 43.85 38.12 38.24 ]
    set TailleMal_C3  [ 43.72 39.79 41.26 41.28 42.54 43.51 41.29 42.12 38.59 33.76 ]
    set TailleMal_C4  [ 44.11 43.38 44.55 43.25 42.4  39.86 42.84 41.96 35.44 39.68 ]
    set TailleMal_C5  [ 42.38 41.94 40.32 42.63 43.16 44.82 42.55 44.38 34.38 33.65 ]
    set TailleMal_C6  [ 44.21 41.44 45.37 40.3  43.64 42.47 42.25 40.4  35.51 38.96 ]
    set TailleMal_C7  [ 43.82 43.33 44.78 45.12 43.91 44.99 45.47 41.67 38.54 36.58 ]
    set TailleMal_C8  [ 39.94 44.09 42.17 44.23 43.16 41.11 43.3  40.24 38.84 35.07 ]
    set TailleMal_C9  [ 43.75 44.07 41.09 41.07 40.2  39.46 44.77 44.38 37.28 37.48 ]
    set TailleMal_C10 [ 43.91 44.70 44.49 40.43 40.76 41.43 37.70 38.97 40.18 37.39 ]
    set TailleMal_C11 [ 44.3  42.42 43.71 42.75 44.53 44.76 40.87 35.46 44.06 42.77 ]


    set TailleFem_C1  [ 42.77 40.93 39.92 43.46 40.57 43.81 40.28 42.5  44.04 44.73 41.01 44.07 40.25 42.33 40.6 ]
    set TailleFem_C2  [ 45.05 44.69 41.12 39.56 43.54 44.65 39.66 43    42.9  44.24 40.36 43.75 41.93 41.14 41.84 ]
    set TailleFem_C3  [ 42.63 42.57 42.02 42.64 43.47 43.45 42.54 40.25 45.47 44.74 40.22 43.19 41.64 44.39 41.14 ]
    set TailleFem_C4  [ 40.56 39.01 39.01 43.75 43.76 41.12 40.04 43.58 42.03 42.58 41.86 42.61 39.62 41.18 42.79 ]
    set TailleFem_C5  [ 43.45 42.8  42.4  44.42 42.66 45.05 40.87 42.21 40.99 42.4  43.49 44.77 43.73 39.84 38.54 ]
    set TailleFem_C6  [ 39.82 40    41.31 44.49 41.33 44.14 41.55 42.42 41.84 44.96 43.81 42.99 41.59 42.26 39.34 ]
    set TailleFem_C7  [ 39.17 45.01 42.84 43.78 42.75 42.3  41.3  45.06 40.68 44.32 42.75 43.8  42.91 40.6  42.17 ]
    set TailleFem_C8  [ 43.77 44.41 41.07 40.75 43.73 45.13 44.96 45.05 44.19 44.85 42.42 42.87 42.94 44.08 42.63 ]
    set TailleFem_C9  [ 41.32 42.45 46.33 41.86 40.27 44.23 44.84 43.68 43.63 41.98 44.11 40.3  43.13 39.66 42.51 ]
    set TailleFem_C10 [ 41.89 43.07 42.35 44.27 43.40 41.34 42.67 43.19 42.30 40.86 44.68 44.79 41.83 40.12 39.78 ]
    set TailleFem_C11 [ 42.88 40.18 40.96 40.9  45.13 44.68 43.35 44.73 41.31 39.6  43.42 43.59 41.08 43.57 39.46 ]

    set End_Experiment 214 ;; days
 ]



  if ( Year = 2012) [                                                                    ; YEAR 2012

    set N.Cosme 3
    set file "Input2012.txt"

    set TailleMal_C1 [ 44.06 48.65 47.18 43.81 42.67 42.37 41.64 43.91 43.89 44.36 ] ; control mesocosm n°2
    set TailleMal_C2 [ 45.06 40.41 44.16 45.05 46.7  42.49 43.26 40.52 43.48 47.15 ] ; control mesocosm n°3
    set TailleMal_C3 [ 42.73 42.83 42.99 45.77 46.33 40.05 45.01 42.69 44.71 46.65 ] ; control mesocosm n°11

    set TailleFem_C1 [ 41.48 44.16 42.97 42.52 44.05 43.22 42.3  45    44.24 42.88 42.41 42.95 42.7  44.91 44.84 ]
    set TailleFem_C2 [ 43.41 42.54 41.07 41.67 42.11 43.08 41.42 42.46 44.76 41.49 44.52 42.29 43.59 43.89 45.61 ]
    set TailleFem_C3 [ 40    43.71 41.07 42.91 43.43 41.64 43.05 41.55 43.58 43.69 42    45.19 44.68 41.3  41.93 ]

    set End_Experiment 211 ;; days

     ]


  if ( Year = 2013) [                                                                    ; YEAR 2013

    set N.Cosme 3

    set file "Input2013.txt"

    set TailleMal_C1 [ 37.3   43.87  41.92  35.33  37.68  38    41  42  44  41.5] ; control mesocosm n°3
    set TailleMal_C2 [ 38.65  41.24  44.98  39.98  38.3   37    42  43  41  38  ] ; control mesocosm n°8
    set TailleMal_C3 [ 37.49  40.73  44.91  36.72  37     37    46  42  44  38  ] ; control mesocosm n°12

    set TailleFem_C1 [ 45.1   42.34  41.9   40.93  44.51  43.99  42.4   43.47  43.96  42.48  42.4   42.11  40.28  44.99  42.15 ]
    set TailleFem_C2 [ 42.25  42.65  44.97  41.64  42.15  42.52  44.24  43.54  41.73  42.84  44.71  44.92  44.62  41.35  41.7  ]
    set TailleFem_C3 [ 43.25  41.7   43.57  42.12  40.47  44.95  42.6   41.83  42.87  42.19  43.47  41.6   42.43  42.96  41.05 ]

    set End_Experiment 217 ;; days

    ]


  if ( Year = 2014) [                                                                    ; YEAR 2014

    set N.Cosme 3
    set file "Input2014.txt"

    set TailleMal_C1 [ 43.13  45.22  45.06  43.09  42.68  46.78  42.78  45.44  41.98  44.57 ] ; control mesocosm n°3
    set TailleMal_C2 [ 43.80  40.42  42.61  40.96  44.71  42.11  42.66  42.62  41.96  43.06 ] ; control mesocosm n°8
    set TailleMal_C3 [ 42.88  41.98  43.88  41.61  44.73  41.81  43.39  47.98  43.53  42.73 ] ; control mesocosm n°9

    set TailleFem_C1 [ 45.06  42.04  40.59  41.86  43.19  45.65  41.66  43.41  44.72  43.76  41.64  42.33  42.69  41.48  40.77 ] ; control mesocosm n°3
    set TailleFem_C2 [ 44.39  43.27  45.48  41.57  42.17  41.53  42.88  42.84  42.77  41.93  42.72  41.89  44.71  42.80  41.34 ] ; control mesocosm n°8
    set TailleFem_C3 [ 43.33  42.69  42.09  43.45  44.51  44.37  44.25  42.67  45.27  42.12  41.09  45.61  45.33  43.60  44.83 ] ; control mesocosm n°9

    set End_Experiment 215 ;; days
       ]


       if ( Year = 2016) [                                                                    ; YEAR 2016: the experiment stopped at 134 days, only one mesocosm

    set N.Cosme 1
    set file "Input2016.txt"

    set TailleMal_C1 [ 36.277  37.351  41.000  41.847  38.639  52.165  46.290  51.224  52.083  54.332 ] ; control mesocosm n°9
    set TailleFem_C1 [ 42.584  37.967  42.591  44.127  37.853  45.855  47.570  44.624  50.129  55.510  50.600  62.631  57.843  54.940  59.811 ] ; control mesocosm n°9

    set End_Experiment 139 ;; days
       ]





; Temperature, photoperiod, food
     file-open file
              set Temperature_Input []  set Photoperiod_Input []
              set Nauplii_Input []  set Copepode_Input []  set Cladocere_Input []  set Rotifere_Input []                                                                        ; Zooplankton
              set Gammare_inf5_Input [] set Gammare_supp5_Input [] set Aselle_inf5_Input []  set Aselle_supp5_Input[] set Chironome_inf5_Input [] set Chironome_supp5_Input []  ; Macroinvertebrates

              while [ not file-at-end? ] [
              set Temperature_Input      sentence Temperature_Input  (list file-read)
              set Photoperiod_Input      sentence Photoperiod_Input (list file-read)
              set Nauplii_Input     sentence Nauplii_Input (list file-read)
              set Copepode_Input     sentence Copepode_Input (list file-read)
              set Cladocere_Input     sentence Cladocere_Input (list file-read)
              set Rotifere_Input     sentence Rotifere_Input (list file-read)
              set Gammare_inf5_Input     sentence Gammare_inf5_Input (list file-read)
              set Gammare_supp5_Input     sentence Gammare_supp5_Input (list file-read)
              set Aselle_inf5_Input     sentence  Aselle_inf5_Input (list file-read)
              set Aselle_supp5_Input     sentence  Aselle_supp5_Input (list file-read)
              set Chironome_inf5_Input    sentence  Chironome_inf5_Input (list file-read)
              set Chironome_supp5_Input    sentence  Chironome_supp5_Input (list file-read)

              ]
     file-close

end




to Parameters ; -------------------------------------------------------------------------

; inter-year parameter
  if ( Year = 2010) [ set  rect_Kdens 0.974 ]
  if ( Year = 2011) [ set  rect_Kdens 1.026 ]

  if ( Year = 2012) [ set  rect_Kdens 0.759 ]
  if ( Year = 2013) [ set  rect_Kdens 0.628 ]
  if ( Year = 2014) [ set  rect_Kdens 0.724 ]

  ; Simulation parameters
  set AdulteSeuil       26                                                               ; Length at which individuals are sexed (mm) (normally it is 25 mm, but fixed to 26 mm to avoid error)
  set L0                5.72                                                             ; Standard length of juveniles at hatching (mm), INERIS

  ; Growth parameters
  set Kappa             0.756645                                                         ; DEB kappa parameter (-)
  set alpha             0.110948                                                         ; Kappa-male = kappa - alpha after puberty (-) calibrated with mesocosm data of 2010
  set Pam               2.42158                                                          ; Maximum area specific assimilation rate  (J/ d / mm^2)
  set v                 1.33408                                                          ; Energy conductance (mm /d)
  set PM                0.111384                                                         ; Volume somatic maintenance costs (J/d/mm^3)
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
  set Photop.Thr        11.305                                                           ; Day duration necessary to start the reproduction cycle

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
  set Bold.M            0.169                                                            ; Aggressivity/boldness of males during foraging

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
Year
Year
2010 2011 2012 2013 2014 2016
1

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
34
83
93
128
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
    <enumeratedValueSet variable="Year">
      <value value="2010"/>
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
    <metric>list ([ L ] of  turtles )</metric>
    <enumeratedValueSet variable="Year">
      <value value="2014"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Nid" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>Nest.Success.Tot</metric>
    <metric>Nest.Built.Tot</metric>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2011"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Frequence" repetitions="100" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>N.tot</metric>
    <metric>N.J</metric>
    <metric>F.J</metric>
    <metric>F.M</metric>
    <metric>F.F</metric>
    <metric>SexRatio</metric>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2011"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
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
  <experiment name="Maturite" repetitions="1000" runMetricsEveryStep="false">
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
  <experiment name="Frequence" repetitions="1000" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup
Calc_Para</setup>
    <go>go</go>
    <metric>L.FreqT</metric>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2014"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Repro" repetitions="1000" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>[Nest] of males with [Cohort = [0 0 0]]</metric>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2011"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Repro2" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean [Nest.Success] of males  with [Cohort = [0 0 0 ]]</metric>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2011"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Survie-Juv" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean [a_inh] of juveniles</metric>
    <metric>mean [a_inh] of juveniles with [L &lt; 15]</metric>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2011"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Year">
      <value value="2011"/>
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
  <experiment name="Maturite2" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count Males with [Puberty = 1]</metric>
    <metric>count Males with [Puberty = 0]</metric>
    <enumeratedValueSet variable="pdt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NF.0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Year">
      <value value="2010"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NM.0">
      <value value="10"/>
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
