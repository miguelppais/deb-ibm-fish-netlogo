;; Rémy Beaudouin 10/2013
;; Benoit Goussen 12/2013
;; Zebrafish agent-based model based on DEB model

; extensions [profiler]
; extensions [r]

;_____________________________________________________________________________________________________________________________
; Agents
;_____________________________________________________________________________________________________________________________
breed [ Juveniles Juvenile ]
breed [ Males Male ]
breed [ Females Female ]
breed [ EggMasses EggMass]

turtles-own   [Generation   Variability   Age  SR    M.R   M.M   M.Age S.rate Test]

Juveniles-own [Sexe       L       W     R     NRJ   dL    dR    dNRJ   food.i   Ec.i  SR.temperature t.puberty ]
Males-own     [Sexe       L       W     R     NRJ   dL    dR    dNRJ   food.i   Ec.i  territory    Engaged]
Females-own   [Sexe       L       W     R     NRJ   dL    dR    dNRJ   food.i   Ec.i  Rm.i         ]
EggMasses-own [Neggs      Sr.L   DD.i   H.rate ]
;_____________________________________________________________________________________________________________________________
; Cellule
;_____________________________________________________________________________________________________________________________
patches-own [Habitat Male.Tenant ]

;_____________________________________________________________________________________________________________________________
; Variables and parameter globals
;_____________________________________________________________________________________________________________________________
globals[
     ;------------------------------------  model variables   ------------------------------------
    Player       W.Tot      W.female  fish.D       Female.D   End-Experiment            pdt
    Food         f_N_P      f_I       f_T          f_h        FB             FBdead
    TIP          TIN        TNS       TPS          Neva       Nfix           Ec        HFc       Af      Hf
    Temperature  Jour.Temp  TempW     Photoperiod  JourPhoto  Photo          SSQ

  ;------------------------------------  model parametres ------------------------------------
    Tc      Ta    Tr    Wv      Wd    Nb.J   Nb.F   Nb.M   Np.b  Np.v FoodInput                               ;; Environment
    I_r     a_I   b_I   T_opta  k_T1  K_T2   Kn     Kpr    Kps   s                                            ;; Food (States   Parameters)
    r_max   K_r   K_ml  Kap     Kan   K_s    K_d    Khn   Khp                                                 ;; Food (States   Parameters)
    G2Kcal  h_n   h_p   Ksn     Knl                                                                           ;; Food (States   Parameters)
    Kappa   v     PAm   PM      Km    Eg     alpha  lf     sdi    PM.t   v.t   PAm.t                          ;; energy
    L.p     lp    Rm    P.t     H.max H.0.5  R.t    SR.mu  SR.sd  SR.a   SR.b   s.dmax T.t                    ;; Puberty and Reproduction
    M.a     M.b   M.c   M.d     M.p   M.e    M.g                                                              ;; Survival
    Aff     Zota  lb    L.b     Linf  Em     g      Fl.m   aW    bW    a.h   b.h                              ;; Growth

  ;------------------------------------ model outputs ------------------------------------
   N.tot    N.tot.adult  N.Fond                  ; Abundance
   F.M      F.F          F.J                     ; Frequency  (Males, Females, Juveniles)
   L.J      L.M          L.F                     ; Mean length
   CV.J     CV.M         CV.F                    ; CV length
   L.FreqM  L.FreqF      L.FreqJuv      L.FreqT  ; Length frequency distribution
   ]

;_____________________________________________________________________________________________________________________________
; Model initialization
;_____________________________________________________________________________________________________________________________

to SETUP
  clear-all
  set-default-shape Juveniles "fish"
  set-default-shape Males     "fish"
  set-default-shape Females   "fish"
  set-default-shape EggMasses "orbit 6"
  Inputs
  Model-parameters
  Setup-patches
  Setup-environment
  Setup-individual
  reset-ticks
end

;_____________________________________________________________________________________________________________________________
; Simulation
;_____________________________________________________________________________________________________________________________

to GO

  ;; Environment variables
  Update-environment
  Update-food

  ;; Agent behaviour
  Move
  Survive
  DEBmodel
  Hatching
  Puberty
  Spawn
  Aging

  ;; Time simulation
  tick
  CalcOutputs
  if ticks = End-Experiment [stop]

end

;_____________________________________________________________________________________________________________________________
; Simulation initialisation
;_____________________________________________________________________________________________________________________________

to Setup-environment

  set Photoperiod  select-list-T Jour.Temp Photo  0
  set Temperature  select-list-T Jour.Temp TempW  0
  set Tc           exp( Ta / Tr - Ta / (Temperature + 273 ) )

  set Af  2.5  * Wv               ; Kcal/pond
  set Hf  18.4 * Wv               ; Kcal/pond
  set TIN 0.97 * Wv               ; total dissolved inorganic nitrogen (g / pond)
  set TIP 0.02 * Wv               ; total dissolved inorganic phosphorus (g / pond)
  set TNS TIN                     ; total nitrogen in sediment (g / pond)
  set TPS TIP                     ; total phosporus in sediment (g / pond)

end


to Setup-patches

   ask patches[
     set Habitat "OpenWater"
     set pcolor  96
     set Male.Tenant 0 - 1 ]

   ask n-of Np.b patches with [ Habitat = "OpenWater" ][
     set Habitat "BreedingGrounds"
     set pcolor  94 ]

   ask n-of Np.v patches with [ Habitat = "OpenWater" ] [
     set Habitat "Vegetation"
     set pcolor  95 ]

end


to Setup-individual

  create-Juveniles Nb.J [
    Initialisation-Juveniles
    set SR normal-seuil SR.mu SR.sd 0 100
    set Generation 0
    set L ( random-float (L.p - L.b ) ) + L.b
    set Age 287 * L / Linf ]


  create-Males  Nb.M [
    set L random-float (Linf - L.p ) + L.p
    set W exp( aW * ln L  - bW )
    set Age  287 * L / Linf
    set Variability (random-normal 0 1 )
    set Generation 0
    Initialisation-Males
    setxy random-xcor random-ycor]

  create-females Nb.F [
    set L random-float (Linf - L.p ) + L.p
    set W exp( aW * ln L  - bW )
    set Age  287 * L / Linf
    set Variability (random-normal 0 1 )
    set Generation 0
    Initialisation-Females ]

  ;; Turtles juveniles, females, males
  set Player turtles with [ breed != eggmasses]

  ;; biomass
  set N.tot     count Player                            ; Total number of individuals (Juveniles + Males + Females)
  set W.Tot     (sum [W] of Player ) / 1000             ; Fish biomass (g)
  set W.female  (sum [W] of Females) / 1000             ; Females biomass (g)
  set fish.D    (N.tot  / Wv  )                          ; Fish density ( g / m3 )
  set Female.D  (W.female  / Wv  )                      ; Females density ( g / m3 )
  set FB        W.Tot * G2Kcal                          ; Fish biomass in Kcal/pond

  ;; Fish food consumption (HFCi / HFC )
  set  f_h ( 1 - exp( - s * ( Hf / (FB + 1E-5) )^ 2.2 ))                              ; food aviability, if FB= 0 --> 1E-5
  ask  Player [ set Ec.i   PAm * 1 * ( L * zota ) ^ 2 * (1 / 4186.8)  ]               ; nrj feeds ad-libitum in Kcal
  set  Ec  ( sum [Ec.i] of player)
  set  HFc Ec * f_h                                                                   ; global nrj feeds in Kcal

end

;_____________________________________________________________________________________________________________________________
; Non fish elements updates (Food, tempertaure, photoperiod)
;_____________________________________________________________________________________________________________________________

 to Update-environment

  set Photoperiod  select-list-T Jour.Temp Photo  ticks
  set Temperature  select-list-T Jour.Temp TempW  ticks

  ;; correct all parameters depending on time to temperature
  set Tc    exp( Ta / Tr - Ta / (Temperature + 273 ) )
  set v.t   v   * Tc                           ; Energy conductance (mm /d)
  set PAm.t PAm * Tc                           ; Maximum area specific assimilation rate  (J/ d / mm^2)
  set PM.t  PM  * Tc                           ; Volume somatic maintenance costs (J/d/mm^3)
  set Km    PM.t / Eg                          ; Somatic maintenance rate

  ;; Turtles juveniles, females, males
  set Player turtles with [ breed != eggmasses]

  ;; biommase
  set N.tot     count Player                            ; Total number of individuals (Juveniles + Males + Females)
  set W.Tot     (sum [W] of Player ) / 1000             ; Fish biomass (g)
  set W.female  (sum [W] of Females) / 1000             ; Females biomass (g)
  set fish.D    (N.tot / Wv     )                       ; Fish density ( g / m3 )
  set Female.D  (W.female  / Wv )                       ; Females density ( g / m3 )
  set FB        (W.Tot * G2Kcal )                       ; Fish biomass in Kcal/pond
  ask patches[   set Male.Tenant (- 1) ]

end

to Update-food

   ;; Temperature function
   ifelse( Temperature <=  T_opta ) [set  f_T exp( - k_T1 * ( Temperature - T_opta )^ 2 ) ][ set  f_T exp( - k_T2 * (T_opta - Temperature )^ 2 )  ]

   ;; solar radiation function
   set f_I ( I_r *  exp( -(a_I * AF + b_I * Hf) ) ) / I_r

   ;; heterotrophic food consumption and availability for fish (dimensionless)
   set  f_h ( 1 - exp( - s * ( Hf / (FB + 1E-5) )^ 2.2 ))                              ; Food aviability, if FB= 0 --> 1E-5
   set  Ec  sum [Ec.i] of Player                                                       ; Global fish energy demands (Kcal)
   set  HFc Ec * f_h                                                                   ; Global fish energy consumed (Kcal)

   ;; Variables computation
   let     INC    TIN /  Wv                                                     ;concentration in mg/L
   let     IPC    TIP /  Wv                                                     ;concentration in mg/L
   let     N.conc INC / ( INC + h_n )                                           ;total inorganic nitrogen mg/L
   let     P.conc IPC / ( IPC + h_p )                                           ;total inorganic phodporus mg/L
   ifelse (N.conc < P.conc )  [ set f_N_P N.conc ] [ set f_N_P P.conc ]
   let    AFR K_r * AF                                                          ;autotrpohic food loss due to phytoplankton loss
   let    AFG r_max * f_N_P * f_I * f_T * AF                                    ;autotrophic food growth due to phytoplankton growth (kcal day-1 pond-1)
   let    AFM K_ml * AF                                                         ;autotrophic food entering heterotrophic food pool (kcal day-1 pond-1) due to phytoplankton mortality and harvest by secondary producers
   let    HFD K_d * HF                                                          ;heterotrophic food loss rate (kcal day-1 pond-1) due to decomposition
   let    HFS K_s * HF                                                          ;heterotrophic food loss rate (kcal day-1 pond-1) due to sedimentation

   ;; Nitrogen dynamics
   set Neva  TIN * Knl
   set Nfix  Kn * AF * exp( - Kn * ( INC ^ 2 ) )
   let Naf   AFG * Kan - Nfix

   let dTIN AFR * Kan + HFD * Khn + TNS * Ksn  - Naf - Neva
   let dTNS HFS * Khn - TNS * Ksn

   set TIN  Positif ( TIN + dTIN ) 0
   set TNS  Positif ( TNS + dTNS  ) 0

   ;; Phosphorus dynamic
   let dTIP AFR * Kap + HFD * Khp + (TPS * Kpr / Wd )   - AFG * Kap  - (TIP * Kps / Wd)
   let dTPS HFS * Khp + (TIP * Kps / Wd) - (TPS * Kpr / Wd)
   set TIP  Positif ( dTIP + TIP ) 0
   set TPS  Positif ( dTPS + TPS ) 0

   ;; Autotrophic food dynamics
   let dAf (AFG  - AFR - AFM)
   set Af Positif ( Af + dAf ) 1E-6

   ;; Heterotrophic food dynamics
   let dHf FBdead + AFM - HFS  - HFD - HFC
   set Hf  Positif ( Hf + dHf ) 1E-6


   ;; nutriment input during monsoon period
     if( Photoperiod > P.t )[
       set TIN TIN + FoodInput * ( 0.97 * Wv )
       set TIP TIP + FoodInput * ( 0.01 * Wv )
       set Hf  Hf  + FoodInput * ( 18.4 * Wv )
       ]

   set FBdead 0  ; intialisation time step variable

 end


;_____________________________________________________________________________________________________________________________
; fish updates
;_____________________________________________________________________________________________________________________________

; Initialisation des agents
;-------------------------------------
 to Initialisation-Juveniles
    Set L L.b
    set Sexe"F"
    set NRJ 1
    setxy random-xcor random-ycor
    set color grey
    set size 0.5
    set W exp( aW * ln L  - bW )
    set Variability (random-normal 0 1 )
 end

 to Initialisation-Males
    set color 23
    set size 1
    set shape "fish"
    setxy random-xcor random-ycor
    set NRJ random-float 1
    set sexe "M"
 end

 to Initialisation-Females
    set color 43
    set size  1.1
    set shape "fish"
    setxy random-xcor random-ycor
    set NRJ random-float 1
    set R   random-float R.t
    set sexe "F"
 end

; Movement
;-------------------------------------
 to Move  ;; turtle procedure

  ask Juveniles [ let p.me one-of patches with [Habitat = "Vegetation" ]  move-to p.me ]

  ask Females [ rt random-float 360  fd random 2  ]

  let M.size sort-by [ [?1 ?2] -> [L] of ?1 > [L] of ?2 ] Males

  foreach M.size [ ?1 ->  ask ?1 [

    let p.me one-of patches with [Habitat = "BreedingGrounds" and Male.Tenant = -1 ]

    ifelse( p.me != nobody )[                          ; First case free breeding ground patches
       move-to p.me
       set territory 1
       ask p.me [ set Male.Tenant  [who] of myself]
    ][                                                 ; Second case No free breeding ground patches.
       lt random-float 360 fd random 2                 ; if no territory = move random
    ] ] ]

end

; Aging
;-------------------------------------
to Aging

  ask turtles [  set Age Age + 1  ]

end


; DEBmodel
;-------------------------------------
to DEBmodel

  ask Player [

    ;; individual fish heterotrophic food consumption (HFCi / HFC )
    set food.i bound ( f_h + Variability *  ( sdi * f_h) ) 0 1
    set Ec.i   PAm.t * 1 * ( L * zota ) ^ 2 * (1 / 4186.8)                   ; nrj feeds ad-libitum in Kcal

    ;;feeding motivation change after puberty
    if (breed = males)[ set food.i (Fl.m * food.i) ]

    ;; Eggs porduction rate related to female density
    if( breed = females )[  set Rm.i Rm * (1 - female.D /(H.0.5 + female.D) ) ]

    ;; intialisation
    let l.dyn ( L / Linf ) ; scaled length
    let r.dyn R
    let e.dyn NRJ
    let d.e 0
    let d.r 0
    let d.l 0
    let F.lim 0

   foreach (n-values (1 / pdt) [ ?1 -> ?1 ]) [

      ;; Stress alimentation
      set F.lim  ( 1 - ( ( 1 - ( 1 / ( 1 + ( lf ^ 3 /  l.dyn ^ 3 ))) ) * alpha )  )

      ;; Energie
      set d.e ( ( ( Km * g ) / l.dyn ) * ( F.lim *  Food.i - e.dyn ) * pdt )

      ;; Growth
      let rB ( Km  * g ) / ( 3 * ( e.dyn + g )  )
      set d.l ( rB * ( e.dyn - l.dyn ) * pdt )
      if ( d.l < 0 ) [set d.l 0]

      ;; Reproduction model
      if( breed = males   )[ set d.r ( ( 1    / ( 1 - lp ^ 3 ) *(   ( (g + l.dyn ) / (g + e.dyn) )  * e.dyn * l.dyn ^ 2  - lp ^ 3  ) ) * pdt )]
      if( breed = females )[ set d.r ( ( Rm.i / ( 1 - lp ^ 3 ) *(   ( (g + l.dyn ) / (g + e.dyn) )  * e.dyn * l.dyn ^ 2  - lp ^ 3  ) ) * pdt )]

      if( ( l.dyn < lp ) or ( dr < 0 ) )[ set d.r 0 ] ;; breed = juveniles

      ;; Environmental control on reproduction
      if( ( Photoperiod < P.t) or (Temperature < T.t) )[ set d.r 0 ]

      ;; iteration
      set e.dyn e.dyn + d.e
      set l.dyn l.dyn + d.l
      set r.dyn r.dyn + d.r
   ]

   set dNRJ e.dyn - NRJ             ; energy dynamic
   set NRJ  e.dyn                   ; Energy at t+1
   set dL   l.dyn - (L / Linf)
   set L    l.dyn * Linf            ; Length at t+1
   set dR   r.dyn - R
   set R    r.dyn                   ; Number of eggs at t + 1
   set W    exp(  aW * ln L  - bW ) ; mass (mg)
   ]

end


;; Death
;-------------------------------------
 to Survive

    ask EggMasses [
    set M.R    ( 1 -  (  M.c  / ( 1 + ( M.d * N.tot ) )  ) )                       ; Density-dependent background mortality probability
    set M.M    M.P                                                                 ; Daily predation mortality
    set S.rate 1 - ( M.M + M.R )                                                   ; Daily survival rate
    let tmp     ( binomial-drawn Neggs S.rate)
    set FBdead  FBdead +  (  0.16 / 1000 * G2Kcal * (Neggs - tmp) )                ; 0.16 mg larvae mass
    set Neggs   tmp ]                                                              ; Number of surviving eggs



   ask Player [
    set M.M    ( M.a * W ^ M.b )                                                   ; Daily predation mortality
    ifelse( Age < M.g)[ set M.Age  0 ][ set M.Age  M.e * (Age -  M.g)  ]           ; Ageing mortality

    set S.rate ( 1 -  ( M.M + M.Age))                                              ; Global daily survival rate

    if( S.rate * 100 ) <= (random-float 100) [
      set FBdead  FBdead +  ( W / 1000 * G2Kcal )
      die   ]
    ]

 end


;; Puberty
;-------------------------------------
 to Puberty

  ask Juveniles [
    set t.puberty t.puberty  + 1
    set SR.temperature SR.temperature + Temperature ]   ; cumualted temperature

  ask Juveniles with [L >= L.p][
   set SR.temperature SR.temperature / t.puberty        ; mean temperature during puberty
   set SR  SR + ( SR.a *( SR.temperature - SR.b ))
   set SR  bound SR 0 100                               ; Apply a boundary to SR
   if  SR >= (random-float 100) [ set sexe "M" ]        ; SR is male frequency

   if(sexe = "M")[  set breed Males       Initialisation-Males   ]
   if(sexe = "F")[  set breed females     Initialisation-females ]  ]

 end



;; spawning
;-------------------------------------
 to Spawn

  let F.size sort-by [ [?1 ?2] -> [L] of ?1 > [L] of ?2 ] Females with [R > R.t ]

  foreach F.size [ ?1 -> ask ?1[

     let Partner max-one-of males with [Territory = 1 and Engaged = 0 ] [L]

     ifelse ( Partner != nobody )[

       ; move-to Partner
       ask Partner[ set Engaged 1 ]

       hatch-EggMasses 1 [
          set H.rate H.max * ( 1 - female.D  / ( H.0.5 + female.D ) )
         set Generation ([Generation] of myself + 1)
         set Neggs [R] of myself
         set color yellow
         set SR normal-seuil SR.mu SR.sd 0 100
         set Age 0
         set DD.i 0
         setxy ([xcor] of Partner)   ([ycor] of Partner)
         set size 0.5  ]
     set R 0
     ][
     if (R > (  2 * R.t )  ) [ set R 0 ]  ; eggs are lost
     ]    ] ]

  ask males [ set Territory 0 set Engaged 0] ; initialisation reproduction

 end


;; Hatching
;-------------------------------------
 to Hatching
   ask EggMasses [

   set DD.i   DD.i + ( Temperature - b.h )

   if (  a.h <= DD.i ) [
      hatch-Juveniles NEggs * H.rate[
      Initialisation-Juveniles
      set SR ([SR] of myself)
      set Generation ([Generation] of myself )  ]
      die  ]
   ]
 end

;_____________________________________________________________________________________________________________________________
; OUTPUTS : Les sorties du modele en fin de simulation
;_____________________________________________________________________________________________________________________________

to CalcOutputs

 ;if( ticks >= 970 ) [ ; 3 * 365.25 / only the 3rd year saved

  set N.tot.adult count Males + count Females  ; Total number of adults (Males + Females)

  ifelse(N.tot > 0)[set F.m count Males / N.tot * 100 ][set F.M 0]      ; Males frequency
  ifelse(N.tot > 0)[set F.F count Females / N.tot * 100 ][set F.F 0]    ; Females frequency
  ifelse(N.tot > 0)[set F.J count Juveniles / N.tot * 100 ][set F.J 0]  ; Juveniles frequency

  set N.Fond count turtles with [Generation = 0]                            ; Number of fish introduced at the onset of experiment

  ifelse ( count Juveniles = 0)[set L.J 0 ][set L.J mean [L] of Juveniles]  ; juveniles mean length
  ifelse ( count Males = 0    )[set L.M 0 ][set L.M mean [L] of Males    ]  ; males mean length
  ifelse ( count Females = 0  )[set L.F 0 ][set L.F mean [L] of Females   ] ; females mean length

  ifelse ( count Juveniles <= 1 )[set CV.J 0 ][set CV.J 100 * (standard-deviation [L] of Juveniles) / mean [L] of Juveniles] ; juveniles length CV
  ifelse ( count Males     <= 1 )[set CV.M 0 ][set CV.M 100 * (standard-deviation [L] of Males    ) / mean [L] of Males    ] ; males length CV
  ifelse ( count Females   <= 1 )[set CV.F 0 ][set CV.F 100 * (standard-deviation [L] of Females  ) / mean [L] of Females  ] ; females length CV

  set L.FreqF    Histo 80 Females 1    ; Call function Histo[Bsup Stage Step]
  set L.FreqM    Histo 80 Males 1      ; Call function Histo[Bsup Stage Step]
  set L.FreqJuv  Histo 26 Juveniles 1  ; Call function Histo[Bsup Stage Step]
  set L.FreqT    Histo 80 Player 1     ; Call function Histo[Bsup Stage Step]


; ]

end

;_____________________________________________________________________________________________________________________________
; INPUTS : Les entrees du modele
;_____________________________________________________________________________________________________________________________

to Inputs ;; Donnees enregistrees pendant l'experience

  ;; la temperature et photopériode
  ;;----------------
             file-open "InputDataFromMay.txt"
              set Jour.Temp []    set TempW []  set Photo []
              while [ not file-at-end? ] [
              set Jour.Temp  sentence Jour.Temp  (list file-read)
              set TempW      sentence TempW  (list file-read)
              set Photo      sentence Photo (list file-read)]
            file-close

end

;_____________________________________________________________________________________________________________________________
; Modes parameters
;_____________________________________________________________________________________________________________________________

to Model-parameters

   ;; Model inputs
   ; --------------
   set pdt               0.1                         ; realtive time step of DEB model (pdt = DEB-pdt/IBM-pdt =  2h/24h = 0.083)
   set End-Experiment    experiment-duration           ; (days) experience duration (selected on interface (default 800))
   set Nb.F              initial-number-females        ; initial nb of males (selected on interface (default 10))
   set Nb.M              initial-number-males          ; initial nb of females (selected on interface (default 10))
   set Nb.J              initial-number-juveniles      ; initial nb of juveniles (selected on interface (default 0))

   ;; Environment
   ; --------------
   set Wv  18            ; water volume (m3) L 6m x l 6m x d 0.5m (selected on interface (default 1))
   set Wd  0.5           ; water depth (m)
   set Np.b 207          ; 23% breeding grounds
   set Np.v 207          ; 23% vegetation cover


   ;; Food
   ; --------------
   set FoodInput 0.0142
   set s         21.08                    ; proportionality coefficient of food nutrient quantity to fish biomass (dimensionless)

   ; nutriment dynamic
   set h_n 0.2                          ; half saturation Nitrogen (mgN/L)
   set h_p 0.02                         ; half saturation Phosporus (mgP/L)
   set Kn  0.01                         ; N-fixation coefficient of phytoplankton
   set Kan 0.0224                       ; N content of phytoplankton
   set Khn 0.0192                       ; N content of heterotrophic components
   set Ksn 0.003                        ; release coefficient of nitrogen in sediment
   set Knl 0.17                         ; coefficient of inorganic nitrogen loss to air
   set Kap 0.001                        ; P content of phytoplankton
   set Khp 0.001                        ; P content of heterotrophic components
   set Kpr 0.0006                       ; release coefficient of phosphorus in sediment (m.day-1)
   set Kps 0.28                         ; coefficient of inorganic phosphorus sedimentation to sediment (m.day-1)

   ; Heterotrophic food dynamics
   set K_s     0.14                     ; coefficient of heterotrophic food sedimentation  (day-1),
   set K_d     0.12                     ; oefficient of heterotrophic food  decomposition (day-1),

   ; Autotrohic food cpt
   set r_max   1.6                      ; maximum growth coefficient for phytoplankton growth (day-1)
   set K_r     0.1                      ; coefficient of phytoplankton respiration (day^-1),
   set K_ml    0.6                      ; coefficient of autotrophic food entering heterotrophic food pool (phytoplankton mortality and harvest; day^-1)

   ; limiting functions of solar radiation
   set I_r     6.547                    ; reference solar radiation 10^6 cal/ m^2/day
   set a_I     0.000017                 ; li et al. 2003 ;  pond/Kcal
   set b_I     0.000015                 ; li et al. 2003 ;  pond/Kcal

   ; Limiting functions of temperature
   set T_opta  30                       ; optimal temperature for phytoplankton growth (°C)
   set k_T1    0.004                    ; effects of temperature below Topta on growth (°C-2).
   set K_T2    0.008                    ; effects of temperature above Topta on growth (°C-2).

   ;; DEB parameter
   ;---------------------------
   set PAm     4.72                     ; Maximum area specific assimilation rate  (J/ d / mm^2)   # <- 2.463
   set v       0.60                     ; Energy conductance (mm /d)
   set Kappa   0.70                     ; fraction of energy to growth/somatic
   set sdi     0.235                    ; variability of the energy acquisition
   set lb      0.084                    ; Scaled length at birth
   set lf      0.156                    ; scaled length at half maximal assimilation
   set lp      0.58                     ; scaled length at puberty (-)
   set alpha   0.84                     ; fraction inaccessible
   set PM      0.44                     ; Volume somatic maintenance costs (J/d/mm^3)
   set Eg      2.35                     ; Cost of synthesis of a unit of structure  (J/mm^3)
   set Ta      3000                     ; Arrhenius temperature (k)
   set Tr      293                      ; Reference temperature (k)
   set Rm      406                      ; maximum number of egg/ day EatoFarl1974b T= 25.5 + 273 K
   set Zota    0.20                     ; Shape coefficient for adults V = zota * L ^3

   ;; Puberty and Reproduction
   ;---------------------------
   set P.t    12 / 24      ; Photoperiod threshold to the onset of the reproductive activity
   set T.t    22.5         ; Temperature threshold to the on/offset of the reproduction
   set H.0.5  24             ; Fish density inducing 50 % reduction fecundity
   set H.max  0.89           ; Hatching rate optimal
   set R.t    263            ; Number of eggs neccesary to spawn
   set SR.mu  50             ; sex-ratio genetic variability
   set SR.sd  23.1           ; sex-ratio genetic variability
   set SR.a   (- 0.0496)     ; slope sex-ratio genetic f(water temperature)
   set SR.b   27.9           ; sex-ratio genetic  °C
   set a.h    60.9           ; age at first feeding expressed in degree.day
   set b.h    10.3           ; Temperature threshold degree.day age at first feeding

   ;; Survival
   ;---------------------------
   set M.a  0.0292              ; natural mortality probability at unit weight
   set M.b  (- 0.382)           ; an allometric scaling factor
   set M.c  0.9576              ; density-independent mortality constant
   set M.d  0.0089              ; slope density-dependent mortality constant
   set M.p  0.025               ; the daily egg predation probability
   set M.g  550                 ; Age threshold of mortality due to aging
   set M.e  0.000002839         ; Ageing mortality

   ;; Growth
   ;----------------------------
   set Fl.m   0.91              ; Male appetite modified by male puberty
   set aW     3.205             ; allometric relation L/W
   set bW     5.1928            ; allometric relation L/W


   ;; Parameters functions of primary parameters
   ;---------------------------
   set Em   PAm / v                                 ; maximum reserve density (J/mm^3)
   set g    Eg / (Kappa * Em)                       ; energy investment ratio
   set Linf (v / ( ( PM / Eg ) * g ) ) * 1 / Zota   ; Maximum physical length (mm)
   set L.p  lp * Linf                               ; (mm) length at puberty
   set L.b  lb * Linf                               ; length at birth (mm)
   set G2Kcal  Em * 1000 / 4186                     ; zebrafish Kcal/g ; 1 Kcal = 4186 J ; d = 1 g/cm3 lindsey et al. 2010


end



;_____________________________________________________________________________________________________________________________
; Reporter
;_____________________________________________________________________________________________________________________________

to-report L.Moy [Sx Co]
  ifelse count  turtles with [breed = Co and sexe = Sx]  > 0
    [ report  mean [L] of turtles with [breed = Co  and sexe = Sx]  ]
    [ report 0 ]
end

to-report R-Ad-F
  ifelse count  females  > 0
    [ report  mean [R] of Females ]
    [ report 0 ]
end

to-report select-list-T [ListC ListD val2]
   let valtmp val2
   if( ( val2 / 365 )  > 1 ) [ set valtmp  ( val2 - ( 365 * floor (val2 / 365) ) )   ]
   let A (position valtmp ListC)
   report item A ListD
end

to-report normal-seuil [Moy Sd Tha Thb]
  let tirage random-normal(Moy) (Sd)
  while [(tirage < Tha) or (tirage > Thb )][ set tirage random-normal(Moy) (Sd) ]
  report tirage

end

to-report Histo [UppB Stage Step]                 ; Histogram of individual length of values from 0 to Bsup for selected stade (Juvenile, Male, Female) with a step of "step" mm
  let temp []
  let Inter n-values UppB [ ?1 -> ?1 ]
  foreach Inter [ ?1 -> set temp lput (frequency (?1) ([L] of Stage) (Step) ) temp ]
  report temp
end


to-report HistoLowerLT [LowerB UppB Stage Step]                 ; Histogram total lenght Ls * 1.25 = LT
  let temp []
  let Inter (n-values (  UppB - LowerB ) [ ?1 -> ?1 ] )
  foreach Inter [ ?1 -> set temp lput (frequency (?1 + LowerB ) (([L * 1.25] of Stage) ) (Step) ) temp ]
  report temp
end

to-report frequency [val thelist Step]            ; Frequency of values between val and val+Step
  report length filter [ ?1 -> (?1 >= val) and (?1 < (val + Step)) ] thelist
end

to-report bound [Nbr LowB UppB]                   ; Apply boundary to a number
  let tmp Nbr
  if( Nbr > UppB ) [ set tmp UppB ]
  if( Nbr < LowB ) [ set tmp LowB ]
  report tmp
end

to-report Positif [Nbr LowBr]                   ; Apply boundary to a number
  let tmp Nbr
  if( Nbr < LowBr ) [ set tmp LowBr ]
  report tmp
end

to-report binomial-drawn [n proba.binom ]
 let res.binom 0
 let l.n n-values n [ ?1 -> ?1 ]
 foreach l.n [  if random-float 1 <  proba.binom [ set res.binom res.binom + 1]   ]
 report res.binom
end
@#$#@#$#@
GRAPHICS-WINDOW
1069
10
1632
574
-1
-1
18.5
1
10
1
1
1
0
0
0
1
0
29
0
29
1
1
1
ticks
30.0

BUTTON
5
10
182
59
Setup
SETUP
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
5
63
182
108
Simulation
GO
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
7
476
355
638
Juvenile lenght
Time (d)
Lenght (mm)
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Juvenile M" 1.0 0 -13345367 true "" "plot L.Moy \"M\" juveniles"
"Juvenile F" 1.0 0 -2064490 true "" "plot L.Moy \"F\" juveniles"

BUTTON
5
110
182
149
Step by step
GO
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
4
643
355
804
Length distribution juveniles
Length (mm)
Count
0.0
15.0
0.0
10.0
true
false
"" ""
PENS
"Juveniles" 1.0 1 -16777216 false "" "histogram [L] of Juveniles"

MONITOR
1283
628
1347
673
Nb Juv
count Juveniles
17
1
11

PLOT
7
312
354
474
Abundances
Time (d)
Effective
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Juveniles" 1.0 0 -16777216 true "" "plot count Juveniles"
"Males" 1.0 0 -13345367 true "" "plot count Males"
"Females" 1.0 0 -2064490 true "" "plot count Females"
"total" 1.0 0 -10899396 true "" "plot count turtles with [breed != eggmasses]"

PLOT
356
476
705
639
Adult lenght
Time (d)
Length (mm)
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Males" 1.0 0 -13345367 true "" "plot L.Moy \"M\" Males"
"Femelles" 1.0 0 -2064490 true "" "plot L.Moy \"F\" Females"

MONITOR
1352
628
1405
673
Nb males
count Males
17
1
11

MONITOR
1409
628
1459
673
Nb Fem
count Females
17
1
11

PLOT
356
643
705
805
Length distribution males
Length (mm)
Count
18.0
45.0
0.0
10.0
true
false
"" ""
PENS
"Males" 1.0 1 -16777216 true "" "histogram [L] of player with [ L >= 18]"
"pen-1" 1.0 0 -13840069 true "" "plot-pen-reset\nplotxy 18  (0.056627185 * count player with [ L >= 18 ])\nplotxy 19  (0.049796353 * count player with [ L >= 18 ])\nplotxy 20  (0.06858616 * count player with [ L >= 18 ])\nplotxy 21  (0.102363165 * count player with [ L >= 18 ])\nplotxy 22  (0.10883644 * count player with [ L >= 18 ])\nplotxy 23  (0.067038776 * count player with [ L >= 18 ])\nplotxy 24  (0.049724555 * count player with [ L >= 18 ])\nplotxy 25  (0.078927574 * count player with [ L >= 18 ])\nplotxy 26  (0.103575812 * count player with [ L >= 18 ])\nplotxy 27  (0.052637194 * count player with [ L >= 18 ])\nplotxy 28  (0.043793243 * count player with [ L >= 18 ])\nplotxy 29  (0.080251814 * count player with [ L >= 18 ])\nplotxy 30  (0.053414717 * count player with [ L >= 18 ])\nplotxy 31  (0.056276341 * count player with [ L >= 18 ])\nplotxy 32  (0.016915638 * count player with [ L >= 18 ])\nplotxy 33  (0.003897849 * count player with [ L >= 18 ])\nplotxy 34  (0.006320148 * count player with [ L >= 18 ])\nplotxy 35  (0.001017035 * count player with [ L >= 18 ])\nplotxy 36  (0 * count player with [ L >= 18 ])\nplotxy 37  (0 * count player with [ L >= 18 ])\nplotxy 38  (0 * count player with [ L >= 18 ])\nplotxy 39  (0 * count player with [ L >= 18 ])"

PLOT
707
643
1061
805
Length distribution females
Length (mm)
Count
18.0
45.0
0.0
10.0
true
false
"" ""
PENS
"Females" 1.0 1 -2064490 true "" "histogram [L] of Females"
"Males" 1.0 0 -13791810 true "" "histogram [L] of Males"

PLOT
707
313
1062
475
Food
Time (d)
Food
0.0
10.0
0.0
1.1
true
true
"" ""
PENS
"Nourr" 1.0 0 -10899396 true "" "plot f_h"
"Sun" 1.0 0 -12345184 true "" "plot f_I"
"Tc" 1.0 0 -955883 true "" "plot f_T"
"Nutr" 1.0 0 -16777216 true "" "plot f_N_P"

PLOT
356
155
705
312
water temperature
Time (d)
Temperature (°C)
0.0
10.0
15.0
35.0
true
false
"" ""
PENS
"default" 1.0 0 -2674135 true "" "plot temperature"
"pen-1" 1.0 0 -7500403 false "" "plotxy 0 22.5\nplotxy End-Experiment 22.5"

PLOT
707
476
1062
640
Nutriments dynamics
Time (d)
Tot N/P (g/pond)
0.0
10.0
0.0
5.0
true
true
"" ""
PENS
"P water" 1.0 0 -8431303 true "" "plot TIP"
"N water" 1.0 0 -14439633 true "" "plot TIN"
"P sediment" 1.0 0 -3889007 true "" "plot TPS"
"N sediment" 1.0 0 -8330359 true "" "plot TNS"

PLOT
356
313
705
475
Biomass
Time (d)
Biomass (g/m3)
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"Total_1" 1.0 0 -10141563 true "" "plot ( ( sum [W] of Females + sum [W] of Males + sum [W] of Juveniles  ) / ( Wv * 1000 ) )"
"Adults" 1.0 0 -14730904 true "" "plot ( ( sum [W] of Females + sum [W] of Males ) / (Wv * 1000) )"
"Juveniles" 1.0 0 -14439633 true "" "plot ( ( sum [W] of Juveniles ) / (Wv * 1000) )"
"Females" 1.0 0 -10402772 true "" "plot female.D"
"Dead" 1.0 0 -16448764 true "" "plot FBdead / (Wv * 1000)"

PLOT
707
156
1061
313
Photoperiod
Time (d)
Photoperiod (%d)
0.0
10.0
0.4
0.6
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot Photoperiod"
"pen-1" 1.0 0 -7500403 false "" "plotxy 0 P.t\nplotxy End-Experiment P.t"

PLOT
7
154
354
311
Compartment energy density
Time (d)
Biomass (KCal/m3)
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Auto" 1.0 0 -8630108 true "" "plot AF / Wv"
"Hetero" 1.0 0 -8275240 true "" "plot HF / Wv"
"Fish" 1.0 0 -955883 true "" "plot FB / Wv"

SLIDER
189
47
464
80
initial-number-males
initial-number-males
0
100
35.0
1
1
NIL
HORIZONTAL

SLIDER
189
83
465
116
initial-number-females
initial-number-females
0
100
35.0
1
1
NIL
HORIZONTAL

TEXTBOX
248
23
415
43
Initial number of individuals
11
0.0
1

SLIDER
189
118
466
151
initial-number-juveniles
initial-number-juveniles
0
500
300.0
1
1
NIL
HORIZONTAL

TEXTBOX
548
25
715
45
Environment setup
11
0.0
1

SLIDER
490
43
697
76
experiment-duration
experiment-duration
0
5 * 364.5
1822.0
1
1
days
HORIZONTAL

MONITOR
821
57
934
102
Fish biomass dead
FBdead
3
1
11

MONITOR
708
10
816
55
Fish consumption
HFC
2
1
11

MONITOR
708
57
816
102
food aviability
f_h
3
1
11

MONITOR
708
103
816
148
Fish energy needs
sum ( [Ec.i] of player )
2
1
11

BUTTON
842
107
913
140
profiler
;setup                  ;; set up the model\n;profiler:start         ;; start profiling\n;repeat 10 [ go ]       ;; run something you want to measure\n;profiler:stop          ;; stop profiling\n;print profiler:report  ;; view the results\n;profiler:reset         ;; clear the data
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
175
112
220
Kcal System
HF + AF + FB
2
1
11

MONITOR
1491
628
1548
673
SSQ
SSQ
3
1
11

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

orbit 6
true
0
Circle -7500403 true true 116 11 67
Circle -7500403 true true 26 176 67
Circle -7500403 true true 206 176 67
Circle -7500403 false true 45 45 210
Circle -7500403 true true 26 58 67
Circle -7500403 true true 206 58 67
Circle -7500403 true true 116 221 67

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
  <experiment name="Spence2007" repetitions="1000" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>L.FreqT</metric>
    <enumeratedValueSet variable="initial-number-juveniles">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experiment-duration">
      <value value="1348"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number-females">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number-males">
      <value value="35"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Hazlerigg2014" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>L.FreqT</metric>
    <enumeratedValueSet variable="initial-number-juveniles">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experiment-duration">
      <value value="1110"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number-females">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number-males">
      <value value="35"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="SimFigArt" repetitions="1000" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>N.tot</metric>
    <metric>N.tot.adult</metric>
    <metric>F.m</metric>
    <metric>F.F</metric>
    <metric>F.J</metric>
    <metric>L.J</metric>
    <metric>L.M</metric>
    <metric>L.F</metric>
    <metric>CV.J</metric>
    <metric>CV.M</metric>
    <metric>CV.F</metric>
    <metric>( ( sum [W] of Females + sum [W] of Males + sum [W] of Juveniles  ) / ( Wv * 1000 ) )</metric>
    <metric>( ( sum [W] of Females + sum [W] of Males ) / (Wv * 1000) )</metric>
    <metric>( ( sum [W] of Juveniles ) / (Wv * 1000) )</metric>
    <metric>AF / Wv</metric>
    <metric>HF / Wv</metric>
    <metric>FB / Wv</metric>
    <enumeratedValueSet variable="initial-number-males">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number-juveniles">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experiment-duration">
      <value value="1110"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number-females">
      <value value="35"/>
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
