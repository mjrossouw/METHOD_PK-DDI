$PROBLEM    PZA_METHOD
;; 1. Based on: run033B
;; 2. Description: 1 CMPT(ALL DATA) - PZA Base model+ADVAN5-TRANSIT+FIXED NN+BSVMTT
;; x1. Author: Jacques
;$SIZES      PD=-1000 LVR=-150 LTH=-200 MAXFCN=10000000  LNP4=-150000


; Settings for the memory of NONMEM

;;2020-02-01

$INPUT      
  NROW ID USUBJID=DROP WHAT=DROP NTIME AMT	RATE=DROP DV_ORIG=DROP	DV	DAT2=DROP TIME EVID	
  MDV OCC BLQ CENS VPCTIME INTENSIVE SEX AGE HT_M WEIGHT_KG	FFM_KG RACE	STRATUM	CIG_DAY	
  SMOK_CAT ALC_G_DAY CREAT_SER_UMMOL_L GFR_ML_MIN CRCL_ML_MIN FED_STATE	METFORMIN EFV_ETV_TFV	
  DOLUTEGRAVIR TLD FLAG	UNOBSERVED NAT=DROP	rec_dose=DROP INC_DOSE=DROP

  
  
 ; IGNORE=@ will skip any line starting with any non-numerical character			
$DATA PZA_METHOD.csv IGNORE=@ IGNORE=(FLAG.EQ.1)
			;IGNORE=(INTENSIVE.EQ.0)
			;IGNORE=(OCC.GT.2) 

$ABB DERIV2=NO ; Prevents  the computation of second derivatives, which are needed only for the Laplacian method.


$SUBROUTINE ADVAN5 TRANS1 ;TOL=9 ATOL=9 SSTOL=6 SSATOL=6

;------------------------------------------------------------------------------
$MODEL        NCOMPS=6 ; NUMBER OF COMPARTMENTS (ABSORPTION COMPATMENT (DEFINED AS FIRST ONE) AND CENTRAL COMPARTMENT DEFIEND AS 2ND COMPARTMENT
              COMP=(TRANSIT1,DEFDOSE)				;1 GUT TRANIST 1 (F1 is associated with first compartment)
			  COMP=(TRANSIT2)						    ;2 GUT TRANIST 2
			  COMP=(TRANSIT3)						    ;3 GUT TRANIST 3
			  COMP=(TRANSIT4)						    ;4 GUT TRANIST 4
              COMP=(ABS) 							      ;5 GUT ABS
              COMP=("CENTRAL",DEFOBS) 			;6 CENTRAL CMT




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
BSVMTT = ETA(18)

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
BOVMTT = 0 
IF (OCC.EQ.1)BOVMTT = ETA(19)
IF (OCC.EQ.2)BOVMTT = ETA(20)

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
 
;---------Typical values--------------------------------------------------------------------------------------------------------------------------------------------------------
IF(METFORMIN.EQ.0) MET_KA = 1
IF(METFORMIN.EQ.1) MET_KA = (1+THETA(14))

;----------------------------------------------------------
TVCL = THETA(1)*ALLMCL_FFM 
TVV = THETA(2)*ALLMV_FFM
TVKA = THETA(3)*MET_KA;*FED_KA* SPARSE_EFFECT
TVBIO = THETA(4)
TVLAG = THETA(7);*FED_LAG
TVV3 = THETA(8)*ALLMV_FFM
TVQ = THETA(9)*ALLMCL_FFM
TVV4 = THETA(10)*ALLMV_FFM
TVQ2 = THETA(11)*ALLMCL_FFM
TVMTT = THETA(12)
TVNN  = THETA(13)


;----------------------------------------------------------------------------------------------------------------------------------------------------

;-----------Define parameters------------------------------------------------------------------------------------------------------------------------------------------
CL  = TVCL*EXP(BSVCL+BOVCL) ; CLEARANCE 
V   = TVV*EXP(BSVV) ; CENTRAL VOL. 
KA  = TVKA*EXP(BSVKA+BOVKA) ; ABS. RATE CONSTANT
BIO = TVBIO*EXP(BSVBIO+BOVBIO) ; BIOAVAILABILITY
LAG =TVLAG*EXP(BOVLAG) ; LAG TIME
V3 = TVV3*EXP(BSVV3) ; PERIPH VOL
Q = TVQ*EXP(BSVQ) ; INTER COMPT CL
V4 = TVV4*EXP(BSVV4) ; PERIPH VOL
Q2 = TVQ2*EXP(BSVQ2) ; INTER COMPT CL
MTT = TVMTT*EXP(BSVMTT+BOVMTT)  ; MTT TIME
NN  = TVNN                      ; Number of transit compartments 
;-----------------------------------------------------------------------------------------------------------------------------------------------------

; re-parameterization
F1 		= BIO
KTR = (NN+1)/MTT
K12	=	KTR  ;Rate between transit CMT
K23	=	KTR   ;Rate between transit CMT
K34	=	KTR   ;Rate between transit CMT
K45	=	KTR   ;Rate between transit CMT
K56 =	KA
K60	=  CL/V
S6  = 	V 							;CENETRAL COMPARTMENT SCALAR (based on numbering in $MODEL)

ALAG1 = LAG


IF (NEWIND/=2.OR.EVID>=3) THEN 
   	TNXD=TIME 
	PNXD=AMT 
	TIMEDOSE = TIME
	AMOUNTDOSE = AMT
ENDIF

TDOS=TNXD 
PD=PNXD 

IF(AMT>0) THEN 
	TNXD=TIME
	PNXD=AMT
ENDIF


;================================================================;

$ERROR

IPRED=A(6)/V

LLOQ = 0.2 ; DEFINE YOUR OWN LLOQ HERE


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

AUC_INF=AMOUNTDOSE*BIO/CL

VARCL = BSVCL + BOVCL
VARBIO = BSVBIO + BOVBIO
VARAUC = BSVBIO + BOVBIO - BSVCL - BOVCL
VARABS = BOVKA + BSVKA -BSVMTT - BOVMTT ;BSVLAG - BOVLAG

;------------------------------------------RETRIEVE AMOUNT IN EACH COMPARTMENT---------------------------------------------------------------------------------------
AMOUNT_1 = A(5) ; ABS CMT
AMOUNT_2 = A(6) ; CENTRAL CMT
CONC_MOD = A(6)/V


;--------------------------------------------------------------------------------------------------------------------------------------------------------------------
$THETA
(0, 3.85,10) ; 1 CL [L/h] 
(0, 39.6,50) ; 2 V [L]
(0, 3.25,30) ; 3 KA [1/h] 
(1) FIX ; 4 BIO
(0, 0.0729,1) ; 5 PROP []
(0, 0.867,10) ; 6 ADD [mg/L]
(0) FIX ; 7 ALAG1_fasting(Baseline) (0.01 to 2 h)
(0, 0,800) FIX ; 8 V3 [L]
(0, 0,90) FIX ; 9 Q [L/h]
(0, 0,800) FIX ; 10 V4 [L]
(0, 0,90) FIX ; 11 Q2 [L/h]
(0.01, 0.242) ; 12 MTT
(4) FIX ; 13 NN
(-0.99, 0,10)FIX ; 14 MET_Ka SCALER

$OMEGA BLOCK(1) 0.1 ; 1 BSV CL
$OMEGA BLOCK(1) 0 FIX  ; 2 BSV V
$OMEGA BLOCK(1) 0 FIX  ; 3 BSV KA
$OMEGA BLOCK(1) 0 FIX  ; 4 BSV BIO
$OMEGA BLOCK(1) 0 FIX  ; 5 BSVV3
$OMEGA BLOCK(1) 0 FIX  ; 6 BSVQ 
$OMEGA BLOCK(1) 0 FIX  ; 7 BSVV4
$OMEGA BLOCK(1) 0 FIX  ; 8 BSVQ2 
$OMEGA BLOCK(1) 0 FIX  ; 9 BSVLAG
$OMEGA BLOCK(1) 0 FIX  ; 18 BSVMTT
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0 FIX  ; 10 BOVCL (OCC1)
$OMEGA BLOCK(1) SAME  ; 11 BOVCL (OCC2)
;----------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0.3 ; 12 BOVBIO (OCC1)
$OMEGA  BLOCK(1) SAME  ; 13 BOVBIO (OCC2)
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0.3 ; 14 BOVKA (OCC1)
$OMEGA  BLOCK(1) SAME ; 15 BOVKA (OCC2)
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0 FIX  ; 16 BOVLAG (OCC1)
$OMEGA  BLOCK(1) SAME ;  17 BOVLAG (OCC2)
;--------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA BLOCK(1) 0.3 ; 19 BOVMTT (OCC1)
$OMEGA  BLOCK(1) SAME ;  20 BOVMTT (OCC2)
;--------------------------------------------------------------------------------------------------------------------------------------------------------

$SIGMA 1 FIX 
;-------------------------------------------------------------------------------------------------------------------------------------------------------

$ESTIMATION MSFO=run033C.msf MAXEVAL=9999 PRINT=1 METHOD=1 INTER NOABORT
NSIG=3 SIGL=9 NONINFETA=1 ETASTYPE=1 ; REPEAT

;As the model becomes more complex, you can use MATRIX=S and then remove the $COVARIANCE step completely when the model is too complex to obtain precisions
$COVARIANCE PRINT=E UNCONDITIONAL; MATRIX=S

;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE WRESCHOL NROW ID OCC
TIME TAD VPCTIME INTENSIVE ;AA1 AA2 ;AA3; AA4
Y DV PRED RES WRES IPRED IRES IWRES CWRESI OBJI
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=sdtab033C.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE ID OCC
CL V KA BIO ALAG1 V3 Q V4 Q2
BSVCL BSVV BSVKA BSVBIO BSVV3 BSVQ BSVV4 BSVQ2 BSVLAG
BOVCL BOVKA BOVBIO BOVLAG
VARCL VARBIO VARAUC
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=patab033C.csv ;parameter
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE ID OCC
WEIGHT_KG HT_M AGE FFM_KG FAT PER_FAT
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=cotab033C.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE ID OCC
SEX
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=catab033C.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE WRESCHOL NROW ID OCC EVID MTT ;
TIME TAD VPCTIME ;AA1 AA2 ;AA3; AA4
Y DV MDV PRED RES WRES IPRED IRES IWRES CWRESI CWRES OBJI INTENSIVE
CL V KA BIO ALAG1 V3 Q V4 Q2
TVCL TVKA TVBIO TVV
BSVCL BSVV BSVKA BSVBIO BSVV3 BSVQ BSVV4 BSVQ2 BSVLAG BSVMTT
BOVCL BOVKA BOVBIO BOVLAG BOVMTT
VARCL VARBIO VARAUC AUC_INF AMOUNT_1 AMOUNT_2 CONC_MOD
WEIGHT_KG HT_M AGE FFM_KG FAT PER_FAT UNOBSERVED BLQ
METFORMIN AMT
NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=mytab033C.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$TABLE
NROW ID NTIME AMT DV TIME EVID MDV OCC BLQ CENS VPCTIME INTENSIVE
SEX AGE HT_M WEIGHT_KG FFM_KG RACE STRATUM CIG_DAY SMOK_CAT ALC_G_DAY
CREAT_SER_UMMOL_L GFR_ML_MIN CRCL_ML_MIN FED_STATE METFORMIN
EFV_ETV_TFV DOLUTEGRAVIR TLD FLAG UNOBSERVED

NOPRINT NOAPPEND ONEHEADER FORMAT=, FILE=vpctab033C.csv
;-------------------------------------------------------------------------------------------------------------------------------------------------------

