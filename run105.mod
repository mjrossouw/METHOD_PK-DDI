$PROBLEM    RIF_METHOD
$INPUT      NROW ID AMT DV TIME EVID MDV OCC BLQ CENS FFM_KG    
  METFORMIN FLAG UNOBSERVED
$DATA       RIF_METHOD_FINAL.csv IGNORE=@
$ABBR DERIV2=NO
$SUBROUTINE ADVAN2 TRANS1 ; 
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
;------------------- Allometric scaling ----------------------
TVFFM = 45; 45.6kg in METHOD TRIAL
;--------- Allometric scaling and covariates ----------------
ALLMCL_FFM = (FFM_KG    / TVFFM)**0.75   ; CL ~ FFM^0.75
ALLMV_FFM  = (FFM_KG    / TVFFM)         ; V  ~ FFM^1.0
;---------------------Adding covariates-----------------------
IF(METFORMIN.EQ.0) MET_BIO = 1
IF(METFORMIN.EQ.1) MET_BIO = (1+THETA(7))
;----------------------Typical values-------------------------
TVCL = THETA(1)*ALLMCL_FFM;
TVV = THETA(2)*ALLMV_FFM
TVKA = THETA(3);
TVBIO = THETA(4)*MET_BIO
;---------------------Define parameters-----------------------
CL  = TVCL*EXP(BSVCL) ; CLEARANCE ; 
V   = TVV; CENTRAL VOL. 
KA  = TVKA*EXP(BOVKA) ; ABS. RATE CONSTANT
BIO = TVBIO*EXP(BOVBIO) ; BIOAVAILABILITY
;-------------------------------------------------------------
; re-parameterization
K = CL/V ;(rate constant of elimination)
F1 = BIO
S2 = V

A_0(1) = 1E-9
A_0(2) = 1E-9
;-------------------------------------------------------------
$ERROR
IPRED=A(2)/V
LLOQ = 0.117 ; 
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
;------------RETRIEVE AMOUNT IN EACH COMPARTMENT-------------
AA1 = A(1)
AA2 = A(2)
;-------------------------------------------------------------
$THETA
(15.8) ; 1 CL [L/h] 
(50.8) ; 2 V [L]
(1.03) ; 3 KA [1/h] 
(1) FIX ; 4 BIO
(0.321) ; 5 PROP []
(0) FIX ; 6 ADD [mg/L]
(-0.240) ; 7 MET_BIO SCALER
;-------------------------------------------------------------
$OMEGA BLOCK(1) 0.164 ; 1 BSV CL
;-------------------------------------------------------------
$OMEGA BLOCK(1) 0.0930 ; 2 BOVBIO (OCC1)
$OMEGA  BLOCK(1) SAME  ; 3 BOVBIO (OCC2)
;-------------------------------------------------------------
$OMEGA BLOCK(1) 0.0936 ; 4 BOVKA (OCC1)
$OMEGA  BLOCK(1) SAME ; 5 BOVKA (OCC2)
;-------------------------------------------------------------
$SIGMA 1 FIX 
;-------------------------------------------------------------
$ESTIMATION MAXEVAL=9999 PRINT=1 METHOD=1 INTER
;-------------------------------------------------------------
