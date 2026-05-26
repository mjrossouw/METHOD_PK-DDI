$PROBLEM    PZA_METHOD
$INPUT      NROW ID AMT DV TIME EVID MDV OCC BLQ CENS 
  FFM_KG METFORMIN FLAG
$DATA       PZA_METHOD.csv IGNORE=@
$ABBR DERIV2=NO
$SUBROUTINE ADVAN5 TRANS1 
;-------------------------------------------------------------
$MODEL        
NCOMPS=6 ; 
COMP=(TRANSIT1,DEFDOSE)				
COMP=(TRANSIT2)						    
COMP=(TRANSIT3)						    
COMP=(TRANSIT4)						    
COMP=(ABS) 							   
COMP=("CENTRAL",DEFOBS) 			
;-------------------------------------------------------------
$PK  
; ------- BSV
BSVCL   = ETA(1)
; ---------- BOV   
BOVBIO = 0 
IF (OCC.EQ.1)BOVBIO = ETA(2)
IF (OCC.EQ.2)BOVBIO = ETA(3)
BOVKA =  0 
IF (OCC.EQ.1)BOVKA = ETA(4)
IF (OCC.EQ.2)BOVKA = ETA(5)
BOVMTT = 0 
IF (OCC.EQ.1)BOVMTT = ETA(6)
IF (OCC.EQ.2)BOVMTT = ETA(7)
;------------------- Allometric scaling ----------------------
;   Body-size covariates 
;   FFM_KG      : fat-free mass in kg
;-------------------------------------------------------------
TVFFM = 45; 45.6kg in METHOD TRIAL
;--------- Allometric scaling and covariates------------------
ALLMCL_FFM = (FFM_KG    / TVFFM)**0.75   ; CL ~ FFM^0.75
ALLMV_FFM  = (FFM_KG    / TVFFM)         ; V  ~ FFM^1.0
;---------Typical values--------------------------------------
IF(METFORMIN.EQ.0) MET_KA = 1 IF(METFORMIN.EQ.1) MET_KA = (1+THETA(9))
;----------------------------------------------------------
TVCL = THETA(1)*ALLMCL_FFM 
TVV = THETA(2)*ALLMV_FFM
TVKA = THETA(3)
TVBIO = THETA(4)
TVMTT = THETA(7)
TVNN  = THETA(8)
;-----------Define parameters---------------------------------
CL  = TVCL*EXP(BSVCL)   ; CLEARANCE 
V   = TVV)              ; CENTRAL VOL. 
KA  = TVKA*EXP(BOVKA)   ; ABS. RATE CONSTANT
BIO = TVBIO*EXP(BOVBIO) ; BIOAVAILABILITY
MTT = TVMTT*EXP(BOVMTT) ; MEAN TRANSIT TIME (MTT)
NN  = TVNN              ; Number of transit compartments 
;-------------------------------------------------------------
; re-parameterization
F1 		= BIO
KTR = (NN+1)/MTT
K12	=	KTR  ;Rate between transit CMT
K23	=	KTR   ;Rate between transit CMT
K34	=	KTR   ;Rate between transit CMT
K45	=	KTR   ;Rate between transit CMT
K56 =	KA
K60	=  CL/V
S6  = 	V ;CENETRAL COMPARTMENT SCALAR 
;============================================================;
$ERROR
IPRED=A(6)/V
LLOQ = 0.2 
CENS_THR = LLOQ 
PROP = IPRED*THETA(5)
ADD = THETA(6)+(CENS_THR*0.2) 

; BLQ handling method M7+, Wijk et al. CPT:PSP 2025 14(6)1042:1049
IF (ICALL/=4.AND.BLQ==1) ADD = ADD + LLOQ

W = SQRT(ADD**2+PROP**2)
IF (W.LE.0.000001) W=0.000001

IRES=DV-IPRED
IWRES=IRES/W

Y = IPRED + W*ERR(1)
IF (ICALL==4.AND.Y<=CENS_THR) Y = 0

;----------- RETRIEVE AMOUNT IN EACH COMPARTMENT--------------
AMOUNT_1 = A(5) ; ABS CMT
AMOUNT_2 = A(6) ; CENTRAL CMT
CONC_MOD = A(6)/V
;-------------------------------------------------------------
$THETA
(3.88) ; 1 CL [L/h] 
(40.1) ; 2 V [L]
(3.45) ; 3 KA [1/h] 
(1) FIX ; 4 BIO
(0.074) ; 5 PROP []
(0.895) ; 6 ADD [mg/L]
(0.424) ; 7 MTT
(4) FIX ; 8 NN

$OMEGA BLOCK(1) 0.0940 ; 1 BSV CL
;-------------------------------------------------------------
$OMEGA BLOCK(1) 0.0282 ; 2 BOVBIO (OCC1)
$OMEGA  BLOCK(1) SAME  ; 3 BOVBIO (OCC2)
;-------------------------------------------------------------
$OMEGA BLOCK(1) 0.653 ; 4 BOVKA (OCC1)
$OMEGA  BLOCK(1) SAME ; 5 BOVKA (OCC2)
;-------------------------------------------------------------
$OMEGA BLOCK(1) 2.10 ; 6 BOVMTT (OCC1)
$OMEGA  BLOCK(1) SAME ;  7 BOVMTT (OCC2)
;-------------------------------------------------------------
$SIGMA 1 FIX 
;-------------------------------------------------------------
$ESTIMATION MAXEVAL=9999 PRINT=1 METHOD=1 INTER
