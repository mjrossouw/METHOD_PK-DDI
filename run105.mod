;$SIZES      PD=-1000 LVR=-150 LTH=-200 MAXFCN=10000000  LNP4=-150000

$PROBLEM    RIF_METHOD
;; 1. Based on: run102
;; 2. Description: 1C(EXL)-A2=MET_BIO-FIXADD
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
  NAT=DROP
  rec_dose=DROP
  INC_DOSE=DROP
  
  
  
 ; IGNORE=@ will skip any line starting with any non-numerical character			
$DATA RIF_METHOD_FINAL.csv IGNORE=@ IGNORE=(FLAG.EQ.1)
	  IGNORE=(ID.EQ.5)
      IGNORE=(ID.EQ.8)
	  IGNORE=(ID.EQ.15)
	  IGNORE=(ID.EQ.38)

$ABB DERIV2=NO ; Prevents  the computation of second derivatives, which are needed only for the Laplacian method.


$SUBROUTINE ADVAN2 TRANS1 ; 1 compartment with Lag

;----------------------------------------------------------------------------------------------------------------------------------------------
$PK  

; ------- BSV
BSVCL   = ETA(1)
BSVV    = ETA(2)
BSVKA   = ETA(3)
BSVBIO  = ETA(4)
BSVV3   = ETA(5)
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
 
;----------------------------------------------------Adding covariates-------------------------------------------------------------------------------------------------------------

IF(METFORMIN.EQ.0) MET_BIO = 1
IF(METFORMIN.EQ.1) MET_BIO = (1+THETA(12))


;---------Typical values--------------------------------------------------------------------------------------------------------------------------------------------------------

;----------------------------------------------------------
TVCL = THETA(1)*ALLMCL_FFM;*CL_RIF ; 
TVV = THETA(2)*ALLMV_FFM
TVKA = THETA(3);
TVBIO = THETA(4)*MET_BIO
TVV3 = THETA(7);*ALLMV_WT
TVQ = THETA(8);*ALLMCL_WT
TVV4 = THETA(9);*ALLMV_WT
TVQ2 = THETA(10);*ALLMCL_WT
TVLAG = THETA(11);

;----------------------------------------------------------------------------------------------------------------------------------------------------

;-----------Define parameters------------------------------------------------------------------------------------------------------------------------------------------
CL  = TVCL*EXP(BSVCL+BOVCL) ; CLEARANCE ; add CLINT when using liver compartment
V   = TVV*EXP(BSVV) ; CENTRAL VOL. 
KA  = TVKA*EXP(BSVKA+BOVKA) ; ABS. RATE CONSTANT
BIO = TVBIO*EXP(BSVBIO+BOVBIO) ; BIOAVAILABILITY
V3 = TVV3*EXP(BSVV3) ; PERIPH VOL
Q = TVQ*EXP(BSVQ) ; INTER COMPT CL
V4 = TVV4*EXP(BSVV4) ; PERIPH VOL2
Q2 = TVQ2*EXP(BSVQ2) ; INTER COMPT CL2
LAG = TVLAG*EXP(BSVLAG)
;-----------------------------------------------------------------------------------------------------------------------------------------------------

; re-parameterization

K = CL/V ;(rate constant of elimination)
ALAG1 = LAG
F1 = BIO


; necessary for these ADVANs (1-4,12,12)
S2 = V

A_0(1) = 1E-9
A_0(2) = 1E-9
;A_0(3) = 1E-6


$ERROR

IPRED=A(2)/V


LLOQ = 0.117 ; DEFINE YOUR OWN LLOQ HERE


CENS_THR = LLOQ 

PROP = IPRED*THETA(5)


ADD = THETA(6)+(CENS_THR*0.2) 


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
;AA3 = A(3)
; AA4 = A(4)

;--------------------------------------------------------------------------------------------------------------------------------------------------------------------
$THETA
(0, 15.8,50) ; 1 CL [L/h] 
(0, 51.1,100) ; 2 V [L]
(0.01, 1.05,10) ; 3 KA [1/h] 
(1) FIX ; 4 BIO
(0.01, 0.321,1) ; 5 PROP []
(0) FIX ; 6 ADD [mg/L]
(0, 0,800) FIX ; 7 V3 [L]
(0, 0,90) FIX ; 8 Q [L/h]
(0, 0,800) FIX ; 9 V4 [L]
(0, 0,90) FIX ; 10 Q2 [L/h]
(0) FIX ; 11 LAG (0.01 to 2 h)
(-0.99, -0.241,10) ; 12 MET_BIO SCALER

$OMEGA BLOCK(1) 0.162 ; 1 BSV CL
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
$OMEGA BLOCK(1) SAME  ; 11 BOVCL (OCC2)
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0.0922 ; 12 BOVBIO (OCC1)
$OMEGA  BLOCK(1) SAME  ; 13 BOVBIO (OCC2)
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0.0925 ; 14 BOVKA (OCC1)
$OMEGA  BLOCK(1) SAME ; 15 BOVKA (OCC2)
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0 FIX  ; 16 BOVLAG (OCC1)
$OMEGA  BLOCK(1) SAME ;  17 BOVLAG (OCC2)
;--------------------------------------------------------------------------------------------------------------------------------------------------------

$SIGMA 1 FIX 
;-------------------------------------------------------------------------------------------------------------------------------------------------------

$ESTIMATION MSFO=run105.msf MAXEVAL=9999 PRINT=1 METHOD=1 INTER NOABORT
NSIG=3 NONINFETA=1 ETASTYPE=1 ; REPEAT

;As the model becomes more complex, you can use MATRIX=S and then remove the $COVARIANCE step completely when the model is too complex to obtain precisions
$COVARIANCE PRINT=E UNCONDITIONAL; MATRIX=S

;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE WRESCHOL NROW ID OCC
TIME TAD VPCTIME AA1 AA2 INTENSIVE ;AA3; AA4
Y DV PRED RES WRES IPRED IRES IWRES CWRESI OBJI
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=sdtab105.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE ID OCC
CL V KA BIO ALAG1 ;V3 Q V4 Q2
BSVCL BSVV BSVKA BSVBIO BSVLAG ;BSVV3 BSVQ BSVV4 BSVQ2
BOVCL BOVKA BOVBIO BOVLAG
VARCL VARBIO VARAUC
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=patab105.csv ;parameter
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE ID OCC
WEIGHT_KG HT_M AGE FFM_KG FAT PER_FAT
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=cotab105.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE ID OCC
SEX
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=catab105.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE WRESCHOL NROW ID OCC EVID AMT ;
TIME TAD VPCTIME AA1 AA2 ;AA3; AA4
Y DV MDV PRED RES WRES IPRED IRES IWRES CWRESI CWRES OBJI INTENSIVE
CL V KA BIO ALAG1 ;V3 Q V4 Q2
BSVCL BSVV BSVKA BSVBIO ;BSVV3 BSVQ BSVV4 BSVQ2 BSVLAG
BOVCL BOVKA BOVBIO BOVLAG
VARCL VARBIO VARAUC
WEIGHT_KG HT_M AGE FFM_KG FAT PER_FAT
SEX METFORMIN
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=mytab105.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE
NROW ID NTIME AMT DV TIME EVID MDV OCC BLQ CENS VPCTIME
INTENSIVE SEX AGE HT_M WEIGHT_KG FFM_KG RACE STRATUM CIG_DAY
SMOK_CAT ALC_G_DAY CREAT_SER_UMMOL_L GFR_ML_MIN CRCL_ML_MIN
FED_STATE METFORMIN EFV_ETV_TFV DOLUTEGRAVIR TLD FLAG
UNOBSERVED
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=vpctab105.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------

