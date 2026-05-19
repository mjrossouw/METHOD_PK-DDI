$SIZES      PD=-1000 LVR=-150 LTH=-200 MAXFCN=10000000  LNP4=-150000

$PROBLEM    INH_METHOD_test
;; 1. Based on: run077
;; 2. Description: 2 CMPT - Base model (exl DATA-INH DATA V2) - Hepatic extraction+INCL PART+MET_BIO
;; x1. Author: Jacques

; Settings for the memory of NONMEM

;;2020-02-01

$INPUT      
  NROW
  ID
  USUBJID=DROP
  WHAT=DROP
  NTIME
  AMT
  RATE=DROP
  DV_ORIG=DROP
  DV
  DAT2=DROP
  TIME
  EVID
  MDV
  OCC
  BLQ
  CENS
  VPCTIME
  INTENSIVE
  SEX
  AGE
  HT_M
  WEIGHT_KG
  FFM_KG
  RACE 
  STRATUM
  CIG_DAY
  SMOK_CAT
  ALC_G_DAY
  CREAT_SER_UMMOL_L
  AST
  ALT
  BILI
  GFR_ML_MIN
  CRCL_ML_MIN
  FED_STATE
  METFORMIN
  EFV_ETV_TFV
  DOLUTEGRAVIR
  TLD
  FLAG
  UNOBSERVED
  NAT


; IGNORE=@ will skip any line starting with any non-numerical character			
$DATA INH_METHOD_FINAL_V2.csv IGNORE=@ IGNORE=(FLAG.EQ.1)       	
	  IGNORE=(ID.EQ.4)
      IGNORE=(ID.EQ.5)
      IGNORE=(ID.EQ.7)

$ABB DERIV2=NO ; Prevents  the computation of second derivatives, which are needed only for the Laplacian method.

$SUBROUTINE ADVAN4 TRANS1 ; 2 compartments
; $SUBROUTINE ADVAN12 TRANS1 ; 3 compartments

$MIX
NSPOP=2               ; 2 POPULATIONS (SLOW AND FAST ACETYLATORS)
P(1)=THETA(13)         ; PROBABILITY OF BEING A FAST ACETYLATOR
P(2)=1-P(1)           ; PROBABILITY OF NOT BEING A FAST ACETYLATOR (THEREFORE SLOW)	

;----------------------------------------------------------------------------------------------------------------------------------------------
$PK  

; ------- BSV
BSVCL   = ETA(1)
BSVV   = ETA(2)
BSVKA   = ETA(3)
BSVBIO = ETA(4)
BSVV3 = ETA(5)
BSVQ = ETA(6)
BSVV4 = ETA(7)
BSVQ2 = ETA(8)
BSVLAG = ETA(9)

; ---------- BOV   
BOVCL =  0 
IF (OCC.EQ.1)BOVCL = ETA(10)
IF (OCC.EQ.2)BOVCL = ETA(11)
BOVBIO = 0 
IF (OCC.EQ.1)BOVBIO = ETA(12)
IF (OCC.EQ.2)BOVBIO = ETA(13)
BOVKA =  0 
IF (OCC.EQ.1)BOVKA = ETA(14)
IF (OCC.EQ.2)BOVKA = ETA(15)
BOVLAG = 0 
IF (OCC.EQ.1)BOVLAG = ETA(16)
IF (OCC.EQ.2)BOVLAG = ETA(17)

;-------------------------------- Allometric scaling ---------------------------------
 ;------------------------------------------------------------------
;   Body-size covariates (values already calculated in R)
;   HT_M        : height in metres
;   WEIGHT_KG   : total body weight in kg
;   FFM_KG      : fat-free mass in kg
;------------------------------------------------------------------
TVWT = 60 ; 61.3kg in METHOD TRIAL -> Use rounded numbers
TVFAT = 15; 15.7kg in METHOD TRIAL
TVFFM = 45; 45.6kg in METHOD TRIAL

FAT = WEIGHT_KG - FFM_KG              ; adipose mass (kg)
IF (FAT < 0) FAT = 0                  ; safety clamp

;--------- Allometric scaling and covariates
ALLMCL_WT  = (WEIGHT_KG / TVWT )**0.75   ; CL ~ WT^0.75
ALLMV_WT   = (WEIGHT_KG / TVWT )         ; V  ~ WT^1.0

ALLMCL_FFM = (FFM_KG    / TVFFM)**0.75   ; CL ~ FFM^0.75
ALLMV_FFM  = (FFM_KG    / TVFFM)         ; V  ~ FFM^1.0

ALLMCL_FAT = (FAT       / TVFAT)**0.75   ; optional fat-based scaling
ALLMV_FAT  = (FAT       / TVFAT)

; --- proportion fat, handy for exploratory plots -----------------
PER_FAT = FAT / WEIGHT_KG

;--------- Allometric scaling for liver
ALLMCL_WT_HEP = (WEIGHT_KG/70)**0.75  ; CL ~ WT^0.75
ALLMV_WT_HEP = (WEIGHT_KG/70)			; V  ~ WT^1.0
 
ALLMCL_FFM_HEP = (FFM_KG/56.1)**0.75	; CL ~ FFM^0.75
ALLMV_FFM_HEP = (FFM_KG/56.1)			; V  ~ FFM^1.0

 
;---------Typical values--------------------------------------------------------------------------------------------------------------------------------------------------------

;------------------------------Mixture Modeling------------------------
EST = MIXEST

IF(MIXNUM.EQ.1) IMP_NAT = 1 ;; CL FOR FAST ACETYLATORS

IF(MIXNUM.EQ.2) IMP_NAT = 2 ;; CL FOR SLOW ACETYLATORS

IF(NAT.LT.0)NAT2_MIX= IMP_NAT ;-1 = Missing genotype info therefore use the imputed


IF(NAT2_MIX.EQ.1) THEN
      TVCL = THETA(1)*ALLMCL_FFM ; CL FOR FAST ACETYLATORS
      ;TVBIO = THETA(4)
   
ENDIF
      
IF (NAT2_MIX.EQ.2) THEN
      TVCL = THETA(2)*ALLMCL_FFM ; CL FOR SLOW ACETYLATORS
      ;TVBIO = THETA(17)
     
ENDIF

;IF (NAT2_MIX.EQ.2) THEN
       ;TVCL = THETA(3)*ALLMCL_FFM ; fast
       ;TVBIO = THETA(18)
;ENDIF   

;-----------------FED STATE------------------------------
FED = 0
IF (FED_STATE.EQ.1) FED = 1      ; 0 = fasted, 1 = fed

; ---------- FED effect on KA (multiplicative) ----------
 
 FED_KA  = 1                      ; fasted reference
 IF(FED.EQ.0) FED_KA = 1 ; reference
IF (FED.EQ.1) THEN
  FED_KA = THETA(14)             ; fold-change when fed
ENDIF                            ; (θ14 lower bound ≥ 0.01)

; ---------- FED effect on LAG (multiplicative) ----------

 FED_LAG = 1                      ; fasted reference
 IF(FED.EQ.0) FED_LAG = 1 ; reference
IF (FED.EQ.1) THEN
  FED_LAG = THETA(15)            ; fold-change when fed
ENDIF                            ; (θ15 lower bound ≥ 0.01)
; ---------- METFORMIN effect on BIO (multiplicative) ----------
IF(METFORMIN.EQ.0) MET_BIO = 1 ; reference
IF (METFORMIN.EQ.1) THEN
  MET_BIO = THETA(18)             ; fold-change when on METFORMIN
ENDIF                            ; (θ16 lower bound ≥ 0.01)

;----------------------------------------------------------

TVV = THETA(3)*ALLMV_FFM
TVKA = THETA(4)*FED_KA
TVBIO = THETA(5)*MET_BIO
TVLAG = THETA(8)*FED_LAG
TVV3 = THETA(9)*ALLMV_FFM
TVQ = THETA(10)*ALLMCL_FFM
TVV4 = THETA(11)*ALLMV_FFM
TVQ2 = THETA(12)*ALLMCL_FFM

;----------------------------------------------------------------------------------------------------------------------------------------------------

;--------------------------HEPATIC CL----------------------------------------;
TVQH=THETA(16)*ALLMCL_FFM_HEP   ; PLASMA FLOW RATE
TVFU=THETA(17)            ; UNBOUND PLASMA FRACTION OF INH



;-----------Define parameters------------------------------------------------------------------------------------------------------------------------------------------
CLINT  = TVCL*EXP(BSVCL+BOVCL) ; CLEARANCE 
V   = TVV*EXP(BSVV) ; CENTRAL VOL. 
KA  = TVKA*EXP(BSVKA+BOVKA) ; ABS. RATE CONSTANT
BIO = TVBIO*EXP(BSVBIO+BOVBIO) ; BIOAVAILABILITY
LAG =TVLAG*EXP(BOVLAG) ; LAG TIME
V3 = TVV3*EXP(BSVV3) ; PERIPH VOL
Q = TVQ*EXP(BSVQ) ; INTER COMPT CL
V4 = TVV4*EXP(BSVV4) ; PERIPH VOL
Q2 = TVQ2*EXP(BSVQ2) ; INTER COMPT CL
QH = TVQH
FU = TVFU

;-----------------------------------------------------------------------------------------------------------------------------------------------------
CL = CLINT ;+ CLr
; re-parameterization

;--------  Transfer constants for liver model ---------

; define hepatic extraction
EH  = (CLINT*FU)/((CLINT*FU)+QH) ; fraction undergoing first pass extraction
FH  = 1 - EH ; fraction available after 1st pass to go to systemic circulation
CLH = QH*EH

; K = CL/V ;(rate constant of elimination)
K = CLH /V
K23 = Q/V   ; (rate constant from central to peripheral 1)
K32 = Q/V3 ;(rate constant from peripheral 1 to central)
K24 = Q2/V ;(rate constant from central to peripheral 2)
K42 = Q2/V4 ; (rate constant from peripheral 2 to central)
;      KA (rate constant of absorption)
ALAG1 = LAG
F1 = FH*BIO
KA = KA

; necessary for these ADVANs (1-4,12,12)
S2 = V

A_0(1) = 1E-6
A_0(2) = 1E-6
A_0(3) = 1E-6

;---------------------------------------------------------------------------------------------------------------------------------------------------------
;================================================================;

$ERROR
IPRED=A(2)/V


; DEFINE LLOQ VALUE 

LLOQ = 0.105 ; 

CENS_THR = LLOQ 

PROP = IPRED*THETA(6)

ADD = THETA(7)+(CENS_THR*0.2) 


IF (ICALL/=4.AND.BLQ==1) THEN
	ADD = ADD + LLOQ
ENDIF

W = SQRT(ADD**2+PROP**2)

; Protective code
IF (W.LE.0.000001) W=0.000001

IRES=DV-IPRED
IWRES=IRES/W

Y = IPRED + W*ERR(1)

; To prevent simulation (ICALL==4) of negative values. It set a positive lower bound for Y, so that VPCs in the log-scale can be plotted
IF (ICALL==4.AND.Y<=CENS_THR) Y = 0


; To calculate time after dose.
IF(AMT>0) THEN
 TIMEDOSE = TIME
 AMOUNTDOSE = AMT
ENDIF

TAD = TIME-TIMEDOSE

VARCL = BSVCL + BOVCL
VARBIO = BSVBIO + BOVBIO
VARAUC = BSVBIO + BOVBIO - BSVCL - BOVCL
VARABS = BOVKA + BSVKA - BSVLAG - BOVLAG ;-BSVMTT - BOVMTT

;------------------------------------------RETRIEVE AMOUNT IN EACH COMPARTMENT---------------------------------------------------------------------------------------
AA1 = A(1)
AA2 = A(2)
AA3 = A(3)
; AA4 = A(4)

;--------------------------------------------------------------------------------------------------------------------------------------------------------------------
$THETA
(0, 31.4,90) ; 1 CL [L/h] fast
(0, 11.7,50) ; 2 CL [L/h] slow
(0, 45.1,250) ; 3 V [L]
(0, 3.04,8) ; 4 KA [1/h] - Fasted(baseline)
(1) FIX ; 5 BIO
(0, 0.178,0.5) ; 6 PROP []
(0) FIX ; 7 ADD [mg/L]
(0) FIX ; 8 ALAG1_fasting(Baseline) (0.01 to 2 h)
(0, 34.4,800) ; 9 V3 [L]
(0, 3.14,90) ; 10 Q [L/h]
(0, 0,800) FIX ; 11 V4 [L]
(0, 0,90) FIX ; 12 Q2 [L/h]
(0.05, 0.691,0.95) ; 13 CL_PROP_FAST (PROPORTION WITH fast(0.26)/intermed(0.55) CL-FAST/INTERMED ACETYLATORS)
(1) FIX ; 14 FED_KA  fold-change (1 = no effect)
(1) FIX ; 15 FED_LAG fold-change (1 = no effect)
(90) FIX ; 16 QH (L/h) -> NB Mathematically: Avg Blood flow to liver in a 70kg man (90 L/H) * (FFM avg of METHOD cohort/Avg FFM of a 70kg man)**0.75 = 76.3L/H
(0.95) FIX ; 17 FU(%)
(0.01, 1,10) ; 18 MET_BIO fold-change (1 = no effect)

$OMEGA BLOCK(1) 0.0411 ; 1 BSV CL
$OMEGA BLOCK(1) 0 FIX  ; 2 BSV V
$OMEGA BLOCK(1) 0 FIX  ; 3 BSV KA
$OMEGA BLOCK(1) 0 FIX  ; 4 BSV BIO
$OMEGA BLOCK(1) 0 FIX  ; 5 BSVV3
$OMEGA BLOCK(1) 0 FIX  ; 6 BSVQ 
$OMEGA BLOCK(1) 0 FIX  ; 7 BSVV4
$OMEGA BLOCK(1) 0 FIX  ; 8 BSVQ2 
$OMEGA BLOCK(1) 0 FIX  ; 9 BSVLAG
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0 FIX  ; 10 BOVCL (OCC1)
$OMEGA BLOCK(1) 0 FIX  ; 11 BOVCL (OCC2)
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0.0515 ; 12 BOVBIO (OCC1)
$OMEGA  BLOCK(1) SAME  ; 13 BOVBIO (OCC2)
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 2.26 ; 14 BOVKA (OCC1)
$OMEGA  BLOCK(1) SAME ; 15 BOVKA (OCC2)
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0 FIX  ; 16 BOVLAG (OCC1)
$OMEGA  BLOCK(1) SAME ; 17 BOVLAG (OCC2)
;--------------------------------------------------------------------------------------------------------------------------------------------------------

$SIGMA 1 FIX 
;-------------------------------------------------------------------------------------------------------------------------------------------------------

$ESTIMATION MSFO=run079.msf MAXEVAL=9999 PRINT=1 METHOD=1 INTER NOABORT
NSIG=3 NONINFETA=1 ETASTYPE=1 ; REPEAT

;As the model becomes more complex, you can use MATRIX=S and then remove the $COVARIANCE step completely when the model is too complex to obtain precisions
$COVARIANCE PRINT=E UNCONDITIONAL; MATRIX=S 

;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE WRESCHOL NROW ID OCC
TIME TAD VPCTIME AA1 AA2 INTENSIVE AA3; AA4
Y DV PRED RES WRES IPRED IRES IWRES CWRESI OBJI
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=sdtab079.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE ID OCC
CL V KA BIO ALAG1 V3 Q V4 Q2
BSVCL BSVV BSVKA BSVBIO BSVV3 BSVQ BSVV4 BSVQ2 BSVLAG
BOVLAG BOVCL BOVKA BOVBIO
VARCL VARBIO VARAUC
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=patab079.csv ;parameter
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE ID OCC
WEIGHT_KG HT_M AGE FFM_KG FAT PER_FAT
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=cotab079.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE ID OCC
SEX NAT2_MIX
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=catab079.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE WRESCHOL NROW ID OCC AMT ;
TIME TAD VPCTIME AA1 AA2 AA3; AA4
Y DV MDV PRED RES WRES IPRED IRES IWRES CWRESI CWRES OBJI INTENSIVE
CL V KA BIO ALAG1 V3 Q V4 Q2
BSVCL BSVV BSVKA BSVBIO BSVV3 BSVQ BSVV4 BSVQ2
BSVLAG BOVCL BOVKA BOVBIO BOVLAG NAT2_MIX
VARCL VARBIO VARAUC
WEIGHT_KG HT_M AGE FFM_KG FAT PER_FAT
SEX UNOBSERVED BLQ EVID METFORMIN
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=mytab079.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE
NROW ID NTIME AMT DV TIME EVID MDV OCC BLQ CENS VPCTIME INTENSIVE
SEX AGE HT_M WEIGHT_KG FFM_KG RACE STRATUM CIG_DAY SMOK_CAT ALC_G_DAY
CREAT_SER_UMMOL_L AST ALT BILI GFR_ML_MIN CRCL_ML_MIN FED_STATE
METFORMIN EFV_ETV_TFV DOLUTEGRAVIR TLD FLAG UNOBSERVED NAT NAT2_MIX
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=vpctab079.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------------------------------------------------------------------------


