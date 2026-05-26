$PROBLEM    INH_METHOD
$INPUT      NROW ID AMT DV TIME EVID MDV OCC BLQ CENS FFM_KG 
  METFORMIN FLAG NAT
$DATA       INH_METHOD_FINAL_V2.csv IGNORE=@
$ABBR DERIV2=NO
$SUBROUTINE ADVAN4 TRANS1
;-------------------------------------------------------------
$MIX
NSPOP=2               ; 2 POPULATIONS (SLOW AND FAST/INTERMEDIATE ACETYLATORS)
P(1)=THETA(10)         ; PROBABILITY OF BEING A FAST/INTERMEDIATE ACETYLATOR
P(2)=1-P(1)           ; PROBABILITY OF BEING A SLOW ACETYLATOR
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
; -----------------------Allometric scaling ------------------
;   Body-size covariates 
TVFFM = 45; 45.6kg in METHOD TRIAL
;--------- Allometric scaling and covariates
ALLMCL_FFM = (FFM_KG    / TVFFM)**0.75   ; CL ~ FFM^0.75
ALLMV_FFM  = (FFM_KG    / TVFFM)         ; V  ~ FFM^1.0

;--------- Allometric scaling for liver
ALLMCL_FFM_HEP = (FFM_KG/56.1)**0.75	; CL ~ FFM^0.75
ALLMV_FFM_HEP = (FFM_KG/56.1)			; V  ~ FFM^1.0
;---------Typical values--------------------------------------
;---------------------Mixture Model---------------------------
EST = MIXEST
IF(MIXNUM.EQ.1) IMP_NAT = 1 ; CL FOR FAST/INTERMEDIATE ACETYLATORS

IF(MIXNUM.EQ.2) IMP_NAT = 2 ;; CL FOR SLOW ACETYLATORS

IF(NAT.LT.0)NAT2_MIX= IMP_NAT ;-1 = Missing genotype info therefore use the imputed

IF(NAT2_MIX.EQ.1) THEN
      TVCL = THETA(1)*ALLMCL_FFM ; CL FOR FAST/INTERMEDIATE ACETYLATORS   
ENDIF
      
IF (NAT2_MIX.EQ.2) THEN
      TVCL = THETA(2)*ALLMCL_FFM ; CL FOR SLOW ACETYLATORS
ENDIF

; -------- METFORMIN effect on BIO (multiplicative) ----------
IF(METFORMIN.EQ.0) MET_BIO = 1 IF(METFORMIN.EQ.1)   MET_BIO = THETA(13)                             
;------------------------------------------------------------
TVV = THETA(3)*ALLMV_FFM
TVKA = THETA(4)
TVBIO = THETA(5)*MET_BIO
TVV3 = THETA(8)*ALLMV_FFM
TVQ = THETA(9)*ALLMCL_FFM
;--------------------------HEPATIC CL-------------------------
TVQH=THETA(11)*ALLMCL_FFM_HEP; PLASMA FLOW RATE
TVFU=THETA(12)               ; UNBOUND PLASMA FRACTION OF INH
;-----------Define parameters---------------------------------
CLINT  = TVCL*EXP(BSVCL) ; INTRINSIC CLEARANCE 
V   = TVV ; CENTRAL VOL. 
KA  = TVKA*EXP(BOVKA) ; ABS. RATE CONSTANT
BIO = TVBIO*EXP(BOVBIO) ; BIOAVAILABILITY
V3 = TVV3*; PERIPHERAL VOLUME
Q = TVQ; INTERCOMPARTMNTAL CL
QH = TVQH ; BLOOD FLOW TO THE LIVER
FU = TVFU ; FRACTION UNBOUND
;-------------------------------------------------------------
CL = CLINT ;

; re-parameterization
;--------  Transfer constants for liver model ---------
EH  = (CLINT*FU)/((CLINT*FU)+QH) ; fraction undergoing first pass extraction
FH  = 1 - EH ; fraction available after 1st pass to go to systemic circulation
CLH = QH*EH

; K = CL/V ;(rate constant of elimination)
K = CLH /V
K23 = Q/V   ; (rate constant from central to peripheral 1)
K32 = Q/V3 ;(rate constant from peripheral 1 to central)

F1 = FH*BIO
S2 = V

A_0(1) = 1E-6
A_0(2) = 1E-6
A_0(3) = 1E-6
;============================================================;

$ERROR
IPRED=A(2)/V
LLOQ = 0.105 ; 
CENS_THR = LLOQ 
PROP = IPRED*THETA(6)
ADD = THETA(7)+(CENS_THR*0.2) 

; BLQ handling method M7+, Wijk et al. CPT:PSP 2025 14(6)1042:1049
IF (ICALL/=4.AND.BLQ==1) ADD = ADD + LLOQ

W = SQRT(ADD**2+PROP**2)
IF (W.LE.0.000001) W=0.000001

IRES=DV-IPRED
IWRES=IRES/W

Y = IPRED + W*ERR(1)

;-------------------------------------------------------------
$THETA
(28.8) ; 1 CL [L/h] FAST/INTERMEDIATE ACETYLATORS
(10.5) ; 2 CL [L/h] SLOW ACETYLATORS
(41.3) ; 3 V [L]
(3.14) ; 4 KA [1/h] 
(1) FIX ; 5 BIO
(0.179) ; 6 PROP []
(0) FIX ; 7 ADD [mg/L]
(30.2) ; 8 V3 [L]
(2.97) ; 9 Q [L/h]
(0.695) ; 10 PROPORTION WITH FAST/INTERMED CL
(90) FIX ; 11 QH (L/h) -> AVERAGE BLOODFLOW TO THE LIVER FOR A 70kg man (90 L/H) * (FFM AVERAGE OF METHOD COHORT/AVERAGE FFM OF A 70kg MAN)**0.75 = 76.3L/H
(0.95) FIX ; 12 FRACTION UNBOUND
(0.857) ; 13 MET_BIO fold-change (1 = no effect)

$OMEGA BLOCK(1) 0.0282 ; 1 BSV CL
;-------------------------------------------------------------
$OMEGA BLOCK(1) 0.0464 ; 2 BOVBIO (OCC1)
$OMEGA  BLOCK(1) SAME  ; 3BOVBIO (OCC2)
;-------------------------------------------------------------
$OMEGA BLOCK(1) 1.24 ; 4 BOVKA (OCC1)
$OMEGA  BLOCK(1) SAME ; 5 BOVKA (OCC2)
;-------------------------------------------------------------
$SIGMA 1 FIX 
;-------------------------------------------------------------
$ESTIMATION MAXEVAL=9999 PRINT=1 METHOD=1 INTER 
 
