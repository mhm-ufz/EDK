C********************************************************************************
C
C    GRGSUB:  Nonlinear Optimization Subroutine
C
C    METHOD and REMARCKS:
C    GRG2 uses an implementation of the generalized reduced gradient (GRG)
C    algorithm. It seeks a feasible solution first (if one is not provided)
C    and then retains feasibility as the objective is improved. It uses a
C    robust implementation of the BFGS quasi-Newton algorithm as its default
C    choice for determining a search direction. A limited-memory conjugate
C    gradient method is also available, permitting solutions of problems with
C    hundreds or thousands of variables. The problem Jacobian is stored and
C    manipulated as a dense matrix, so the effective size limit is one to two
C    hundred active constraints (excluding simple bounds on the variables,
C    which are handled implicitly). 
C
C    GRG2 is a program that solves nonlinear optimization problems of the
C    following form:
C
C      Minimize or maximize    gp(X),
C
C      subject to:
C                  glbi <= gi(X) <= gubi     for i=1,...,m, i # p
C                  xlbj <= xj    <= xubj	   for j=1,...,n.
C
C      X is a vector on n variables, x1 ,...,xn, and the functions
C      g1 ,...,gm all depend on X.
C
C    Any of these functions may be nonlinear. Any of the bounds may be infinite
C    and any of the constraints may be absent. If there are no constraints,
C    the problem is solved as an unconstrained optimization problem. Upper and
C    lower bounds on the variables are optional and, if present, are not treated
C    as additional constraints but are handled separately. The program solves
C    problems of this form by the Generalized Reduced Gradient Methods
C                                                                 
C*******************************************************************************
C
      SUBROUTINE GRGSUB(INPRNT,OTPRNT,NCORE,NNVARS,NFUN,MAXBAS,         GRG25260
     1 MAXHES,NNOBJ,TTITLE,BLVAR,BUVAR,BLCON,BUCON,DEFAUL,FPNEWT,FPINIT,GRG25270
     2 FPSTOP,FPSPIV,PPH1EP,NNSTOP,IITLIM,LLMSER,IIPR,IIPN4,IIPN5,      GRG25280
     3 IIPN6,IIPER,IIDUMP,IIQUAD,LDERIV,MMODCG,                         GRG25290
     4 RAMCON,RAMVAR,XX,FCNS,INBIND,RMULTS,NONBAS,REDGR,                GRG25300
     5 NBIND,NNONB,INFORM,Z)                                            GRG25310
C                                                                       GRG25320
C     THIS IS THE SUBROUTINE INTERFACE FOR GRG2. DESPITE THE LONG ARGUMEGRG25330
C     LIST, IT IS ACTUALLY EASY TO USE. MANY OF THE INPUT ARGUMENTS HAVEGRG25340
C     DEFAULTS AND NEED NOT BE SPECIFIED.                               GRG25350
C                                                                       GRG25360
C     THE ARGUMENTS ARE DIVIDED INTO 5 SECTIONS:                        GRG25370
C                                                                       GRG25380
C     1. ARGUMENTS {INPRNT} THROUGH {DEFAUL} -                          GRG25390
C        VALUES FOR THESE MUST BE PROVIDED                              GRG25400
C                                                                       GRG25410
C     2. ARGUMENTS {FPNEWT} THROUGH {RAMVAR}-                           GRG25420
C        THESE ARE ALSO INPUT ARGUMENTS. VALUES FOR THESE ARE           GRG25430
C        OPTIONAL. THE ARRAY {DEFAUL} TELLS GRGSUB WHICH ONES           GRG25440
C        YOU WANT TO USE THE GRG2 DEFAULT VALUES FOR, AS IS EXPLAINED   GRG25450
C        BELOW                                                          GRG25460
C                                                                       GRG25470
C     3. ARGUMENT XX - THIS IS BOTH AN INPUT AND OUTPUT ARGUMENT.       GRG25480
C        VALUES MUST BE PROVIDED PRIOR TO CALLING GRGSUB. ANY SUCH      GRG25490
C        VALUES WHICH DO NOT SATISFY THE BOUNDS ON THE VARIABLE WILL    GRG25500
C        BE CHANGED TO THE BOUND NEAREST THE VALUE. ON TERMINATION      GRG25510
C        XX GIVES GRGSUB-S ESTIMATE OF THE OPTIMUM                      GRG25520
C                                                                       GRG25530
C     4. ARGUMENTS {FCNS} THROUGH {INFORM} - THESE ARE OUTPUT ARGUMENTS GRG25540
C        SEE BELOW FOR DESCRIPTIONS.                                    GRG25550
C                                                                       GRG25560
C     5. Z - THIS IS AN INPUT ARGUMENT, THE WORK ARRAY USED BY GRG2.    GRG25570
C        IT MUST BE DIMENSIONED IN THE USER-S CALLING PROGRAM, AND {NCORGRG25580
C        MUST BE SET TO IT-S DIMENSION THERE. IF THE DIMENSION IS NOT   GRG25590
C        LARGE ENOUGH, AN ERROR MESSAGE WILL BE PRINTED, TELLING YOU    GRG25600
C        HOW MUCH CORE IS NEEDED.                                       GRG25610
C                                                                       GRG25620
C        ++++++++++++++++++++++++++++++++++                             GRG25630
C        + INPUT VARIABLES AND PARAMETERS +                             GRG25640
C        ++++++++++++++++++++++++++++++++++                             GRG25650
C                                                                       GRG25660
C                                                                       GRG25670
C     -DEFAUL(1<=I<=19)                                                 GRG25680
C              .NE. 1.0 - USE PARAMETER VALUES ASSIGNED BY USER         GRG25690
C              .EQ. 1.0 - USE DEFAULT PARAMETERS VALUES                 GRG25700
C                                                                       GRG25710
C                                                                       GRG25720
C     INPUT PROGRAM PARAMETERS TO BE PROVIDED IF DEFAUL(I) .NE. 1.0     GRG25730
C     -------------------------------------------------------------     GRG25740
C                                                                       GRG25750
C     (I)-(INPUT VARIABLE)--(INTERNAL GRG2 VARIABLE)--(DEFAULT VALUE)--(GRG25760
C                                                                       GRG25770
C     1   FPNEWT--EPNEWT---1.0E-04--- A CONSTRAINT IS ASSUMED TO BE     GRG25780
C                      BINDING IF IT IS WITHIN THIS EPSILON             GRG25790
C                      OF ONE OF ITS BOUNDS.                            GRG25800
C     2   FPINIT--EPINIT---1.0E-04--- IF IT IS DESIRED TO RUN THE       GRG25810
C                      PROBLEM WITH {EPNEWT} INITIALLY SET FAIRLY       GRG25820
C                      LARGE AND THEN TIGHTENED AT THE END OF THE       GRG25830
C                      OPTIMIZATION THEN THIS IS ACCOMPLISHED BY        GRG25840
C                      ASSIGNING {EPINIT} THE INITIAL TOLERANCE         GRG25850
C                      AND {EPNEWT} THE FINAL ONE.                      GRG25860
C     3   FPSTOP--EPSTOP---1.0E-04--- IF THE FRACTIONAL CHANGE IN THE   GRG25870
C                      OBJECTIVE IS LESS THAN {EPSTOP} FOR {NSTOP}      GRG25880
C                      CONSECUTIVE ITERATIONS, THE PROGRAM WILL         GRG25890
C                      STOP. THE PROGRAM WILL ALSO STOP IF KUHN-TUCKER  GRG25900
C                      OPTIMALITY CONDITIONS ARE SATISFIED TO WITHIN    GRG25910
C                      {EPSTOP}.                                        GRG25920
C     4   FPSPIV--EPSPIV---10.0E-3--- IF, IN CONSTRUCTING THE BASIS     GRG25930
C                      INVERSE, THE ABSOLUTE VALUE OF A PROSPECTIVE     GRG25940
C                      PIVOT ELEMENT IS LESS THAN {EPSPIV}, THE         GRG25950
C                      PIVOT WILL BE REJECTED AND ANOTHER PIVOT         GRG25960
C                      ELEMENT WILL BE SOUGHT.                          GRG25970
C     5   PPH1EP--PH1EPS--- 0.0   --- IF NONZERO, THE PHASE 1 OBJECTIVE GRG25980
C                      IS AUGMENTED BY A MULTIPLE OF THE TRUE           GRG25990
C                      OBJECTIVE.  THE MULTIPLE IS SELECTED SO THAT,    GRG26000
C                      AT THE INITIAL POINT, THE RATIO OF THE TRUE      GRG26010
C                      OBJECTIVE AND SUM OF THE INFEASIBILITIES IS      GRG26020
C                      {PH1EPS}.                                        GRG26030
C     6   NNSTOP--NSTOP --- 3     --- IF THE FRACTIONAL CHANGE IN THE   GRG26040
C                      OBJECTIVE IS LESS THAN {EPSTOP} FOR {NSTOP}      GRG26050
C                      CONSECUTIVE ITERATIONS, THE PROGRAM WILL         GRG26060
C                      STOP.                                            GRG26070
C     7   IITLIM--ITLIM --- 10    --- IF SUBROUTINE NEWTON TAKES        GRG26080
C                      {ITLIM} ITERATIONS WITHOUT CONVERGING            GRG26090
C                      SATISFACTORILY, THE ITERATIONS ARE STOPPED       GRG26100
C                      AND CORRECTIVE ACTION IS TAKEN.                  GRG26110
C     8   LLMSER--LIMSER---10,000 --- IF THE NUMBER OF COMPLETED ONE    GRG26120
C                      DIMENSIONAL SEARCHES EQUALS {LIMSER},            GRG26130
C                      OPTIMIZATION WILL TERMINATE.                     GRG26140
C     9   IIPR  --IPR   --- 0 - SUPPRESS ALL OUTPUT PRINTING EXCEPT     GRG26150
C                           INITIAL AND FINAL REPORTS.                  GRG26160
C                   --- 1 - PRINT ONE LINE OF OUTPUT FOR EACH ONE       GRG26170
C                           DIMENSIONAL SEARCH.                         GRG26180
C                   --- 2 - PROVIDE MORE DETAILED INFORMATION ON        GRG26190
C                           THE PROGRESS OF EACH ONE DIMENSIONAL        GRG26200
C                           SEARCH.                                     GRG26210
C                   --- 3 - EXPAND THE OUTPUT TO INCLUDE THE PROBLEM    GRG26220
C                           FUNCTION VALUES AND VARIABLE VALUES AT      GRG26230
C                           EACH ITERATION AS WELL AS THE SEPARATION    GRG26240
C                           OF CONSTRAINTS INTO NONBINDING AND          GRG26250
C                           BINDING AND VARIABLES INTO BASIC,           GRG26260
C                           SUPERBASIC AND NONBASIC.                    GRG26270
C                   --- 4 - AT EACH ITERATION THE REDUCED GRADIENT,     GRG26280
C                           THE SEARCH DIRECTION AND THE TANGENT        GRG26290
C                           VECTOR ARE PRINTED.                         GRG26300
C                   --- 5 - PROVIDES DETAILS OF THE BASIS INVERSION     GRG26310
C                           PROCESS INCLUDING THE INITIAL BASIS AND     GRG26320
C                           ITS INVERSE.  ALSO DISPLAYS THE VARIABLE    GRG26330
C                           VALUES AND CONSTRAINT ERRORS FOR EACH       GRG26340
C                           NEWTON ITERATION.                           GRG26350
C                   --- 6 - THIS IS THE MAXIMUM LEVEL OF PRINT          GRG26360
C                           AVAILABLE AND INCLUDES ALL OF THE ABOVE     GRG26370
C                           ALONG WITH DETAILED PROGRESS OF THE         GRG26380
C                           BASIS CONSTRUCTION PHASE, INCLUDING         GRG26390
C                           THE BASIS INVERSE AT EACH PIVOT.            GRG26400
C     10  IIPN# --IPN#  ---   0   --- IF IIPN# IS GREATER THAN ZERO THENGRG26410
C     11                             WILL BE SET TO # AFTER IIPN# ITERATGRG26420
C     12                                                                GRG26430
C     13  IIPER --IPER  ---   0   --- IF IIPER IS GREATER THAN ZERO THENGRG26440
C                                FOR EVERY IIPER-TH ITERATION, PRINT    GRG26450
C                                USING THE CURRENT VALUE OF {IPR}       GRG26460
C                                OTHERWISE USE IPR=1 .                  GRG26470
C     14  IIDUMP--IDUMP ---   0   --- IF IIDUMP IS GREATER THAN ZERO THEGRG26480
C                                AFTER IIDUMP ITERATIONS, PROVIDE AN    GRG26490
C                                OUTPUT DUMP ON OUTPUT UNIT 7 AND       GRG26500
C                                TERMINATE COMPUTATION.                 GRG26510
C                                                                       GRG26520
C     15  IIQUAD--IQUAD---- 0 - METHOD FOR INITIAL ESTIMATES OF BASIC   GRG26530
C                           VARIABLES FOR EACH ONE DIMENSIONAL          GRG26540
C                           SEARCH                                      GRG26550
C                   --- 0 - TANGENT VECTOR AND LINEAR EXTRAPOLATION     GRG26560
C                           WILL BE USED.                               GRG26570
C                   --- 1 - QUADRATIC EXTRAPOLATION WILL BE USED.       GRG26580
C                                                                       GRG26590
C     16  LDERIV--KDERIV--- 0 - METHOD FOR OBTAINING PARTIAL DERIVATIVE GRG26600
C                   --- 0 - FORWARD DIFFERENCE APPROXIMATION            GRG26610
C                   --- 1 - CENTRAL DIFFERENCE APPROXIMATION            GRG26620
C                   --- 2 - USER SUPPLIED SUBROUTINE {PARSH} IS USED    GRG26630
C                                                                       GRG26640
C     17  MMODCG--MODCG --- 0 - {MODCG} AND {MAXHES} (SEE BELOW) CONTROLGRG26650
C                           USE OF A CONJUGATE GRADIENT ( CG )          GRG26660
C                           METHOD.  IF THE NUMBER OF SUPERBASIC        GRG26670
C                           VARIABLES EXCEEDS {MAXHES}, THE CG          GRG26680
C                           METHOD INDICATED BY {MODCG} IS USED.        GRG26690
C                           DEFAULT VALUE OF MODCG=1 .  TO USE A        GRG26700
C                           CG METHOD AT EACH ITERATION, SET            GRG26710
C                           MAXHES=0 .                                  GRG26720
C                   --- 1 - USES FLETCHER-REEVES FORMULA.               GRG26730
C                   --- 2 - USES POLAK-RIBIERE FORMULA.                 GRG26740
C                   --- 3 - USES PERRY'S FORMULA.                       GRG26750
C                   --- 4 - USES 1 STEP VERSION OF DFP.                 GRG26760
C                   --- 5 - USES 1 STEP VERSION OF BFS.                 GRG26770
C                                                                       GRG26780
C     18  RAMCON AND RAMVAR--NAMCON AND NAMVAR--BLANKS                  GRG26790
C                   REAL ARRAYS OF DIMENSION NFUN AND NNVARS RESPECTIVELGRG26800
C                   FOR ENTERING NAMES OF FUNCTIONS AND VARIABLES.  OBJEGRG26810
C                   AND CONSTRAINT FUNCTION NAMES GO IN {RAMCON}, VARIABGRG26820
C                   NAMES GO IN {RAMVAR}.                               GRG26830
C                                                                       GRG26840
C     19 {NO GRGSUB VARIABLE}--MAXIM--FALSE                             GRG26850
C                   IF DEFAUL(19)=1.0, OBJECTIVE WILL BE MINIMIZED      GRG26860
C                   IF NOT EQUAL  1.0, OBJECTIVE WILL BE MAXIMIZED      GRG26870
C                                                                       GRG26880
C                                                                       GRG26890
C     OTHER INPUT VARIABLES --- THESE HAVE NO DEFAULT VALUES, SO        GRG26900
C     +++++++++++++++++++++     VALUES MUST BE PROVIDED BY THE USER     GRG26910
C                                                                       GRG26920
C     LOGICALS                                                          GRG26930
C     ========                                                          GRG26940
C                                                                       GRG26950
C     -INPRNT   FALSE:  DO NOT PRINT ANY ECHO BACK OF INPUT DATA        GRG26960
C              TRUE:   PRINT INPUT DATA                                 GRG26970
C                                                                       GRG26980
C     -OTPRNT   FALSE:  DO NOT PRINT ANY FINAL RESULTS                  GRG26990
C              TRUE:   PRINT FINAL RESULTS                              GRG27000
C                                                                       GRG27010
C                                                                       GRG27020
C     SCALARS                                                           GRG27030
C     =======                                                           GRG27040
C                                                                       GRG27050
C     (INPUT VARIABLE)--(GRG2 INTERNAL VARIABLE)--(MEANING)             GRG27060
C                                                                       GRG27070
C     NCORE --NCORE -- DIMENSION OF THE {Z} ARRAY                       GRG27080
C     NNVARS--NVARS -- NUMBER OF VARIABLES                              GRG27090
C     NFUN  --NROWS -- NUMBER OF FUNCTIONS INCLUDING OBJECTIVE          GRG27100
C     MAXBAS--MAXB  -- UPPER LIMIT ON NUMBER OF BINDING                 GRG27110
C                      CONSTRAINTS - USE {NFUN} IF UNSURE OF SMALLER LIMGRG27120
C     MAXHES--MAXR  -- MAXIMUM ALLOWABLE SIZE OF APPROXIMATE            GRG27130
C                      HESSIAN - USE {NVARS} IF YOU WANT A QUASI-NEWTON GRG27140
C                      TO BE USED AT EVERY ITERATION (FASTEST METHOD IF GRG27150
C                      NOT TOO MANY VARIABLES)                          GRG27160
C     NNOBJ --NOBJ  -- INDEX OF COMPONENT OF VECTOR {G} IN              GRG27170
C                      SUBROUTINE GCOMP CORRESPONDING TO                GRG27180
C                      OBJECTIVE FUNCTION                               GRG27190
C                                                                       GRG27200
C     ARRAYS                                                            GRG27210
C     ======                                                            GRG27220
C                                                                       GRG27230
C     BLVAR --ALB   -- REAL ARRAY CONTAINING LOWER BOUNDS OF VARIABLES. GRG27240
C                      MUST HAVE DIMENSION EQUAL TO {NNVARS}. IF        GRG27250
C                      A VARIABLE XX(I) HAS NO LOWER BOUND, SET BLVAR(I)GRG27260
C                      TO -1.0E31                                       GRG27270
C                                                                       GRG27280
C     BUVAR --UB    -- REAL ARRAY CONTAINING UPPER BOUNDS OF VARIABLES. GRG27290
C                      MUST HAVE DIMENSION EQUAL TO VALUE OF {NNVARS}.  GRG27300
C                      IF XX(I) HAS NO UPPER BOUND, SET BUVAR(I) TO     GRG27310
C                      1.0E31                                           GRG27320
C                                                                       GRG27330
C     BLCON --ALB   -- REAL ARRAY CONTAINING LOWER BOUNDS OF ALL FUNCTIOGRG27340
C                      DEFINED IN {GCOMP}. MUST HAVE DIMENSION EQUAL TO GRG27350
C                      VALUE OF {NFUN}. IF A FUNCTION G(I) HAS NO LOWER GRG27360
C                      BOUND, SET BLCON(I) TO -1.0E31                   GRG27370
C                                                                       GRG27380
C     BUCON --UB    -- REAL ARRAY CONTAINING UPPER BOUNDS OF ALL FUNCTIOGRG27390
C                      DEFINED IN {GCOMP}. MUST HAVE DIMENSION EQUAL TO GRG27400
C                      VALUE OF {NFUN}. IF A FUNCTION G(I) HAS NO UPPER GRG27410
C                      BOUND, SET BUCON(I) TO 1.0E31                    GRG27420
C                                                                       GRG27430
C                                                                       GRG27440
C     NOTE 1: IT DOES NOT MATTER WHAT YOU USE FOR THE BOUNDS OF         GRG27450
C             THE OBJECTIVE FUNCTION IN {BUCON} AND {BLCON}, BUT        GRG27460
C             SOMETHING MUST BE THERE                                   GRG27470
C                                                                       GRG27480
C     NOTE 2: IF YOU WISH TO FIX A VARIABLE AT A CERTAIN VALUE AND      GRG27490
C             HAVE GRGSUB LEAVE IT UNCHANGED, SET ITS ENTRIES IN BLVAR  GRG27500
C             AND BUVAR TO THAT VALUE                                   GRG27510
C                                                                       GRG27520
C     NOTE 3: IF G(I) IS AN EQUALITY CONSTRAINT, EQUAL TO, SAY, B, SET  GRG27530
C             BLCON(I)=BUCON(I)=B                                       GRG27540
C                                                                       GRG27550
C     NOTE 4: IF A CONSTRAINT G(I) IS TO BE IGNORED IN THE CURRENT CALL GRG27560
C             TO GRGSUB, SET BLCON(I) TO -1.0E31 AND BUCON(I) TO 1.0E31 GRG27570
C                                                                       GRG27580
C     XX    --X     -- ARRAY CONTAINING INITIAL VALUES OF               GRG27590
C                      VARIABLES. MUST HAVE DIMENSION EQUAL TO VALUE    GRG27600
C                      OF {NNVARS}.                                     GRG27610
C                                                                       GRG27620
C                                                                       GRG27630
C        +++++++++++++++++++++++++++++++++++                            GRG27640
C        + OUTPUT VARIABLES AND PARAMETERS +                            GRG27650
C        +++++++++++++++++++++++++++++++++++                            GRG27660
C                                                                       GRG27670
C     XX     -- REAL ARRAY OF DIMENSION EQUAL TO VALUE OF {NNVARS}.     GRG27680
C               XX(I) CONTAINS FINAL VALUE OF X(I), WHERE X IS ARRAY    GRG27690
C               OF VARIABLES IN GCOMP.                                  GRG27700
C                                                                       GRG27710
C     FCNS   -- REAL ARRAY OF DIMENSION EQUAL TO VALUE OF {NFUN}.       GRG27720
C               FCNS(I) CONTAINS FINAL VALUE OF G(I) WHERE G IS ARRAY   GRG27730
C               OF FUNCTIONS IN GCOMP.                                  GRG27740
C                                                                       GRG27750
C     INBIND -- INTEGER ARRAY OF DIMENSION EQUAL TO VALUE OF NFUN. THE FGRG27760
C               NBIND POSITIONS OF INBIND CONTAIN THE INDICES OF THOSE  GRG27770
C               CONSTRAINTS (COMPONENTS OF THE G VECTOR) WHICH EQUAL EITGRG27780
C               THEIR LOWER OR UPPER BOUNDS AT TERMINATION.             GRG27790
C                                                                       GRG27800
C     RMULTS -- REAL ARRAY OF DIMENSION NFUN. THE FIRST NBIND POSITIONS GRG27810
C               RMULTS CONTAIN THE LAGRANGE MULTIPLIERS OF THE BINDING  GRG27820
C               CONSTRAINTS, IN THE SAME ORDER AS THE INDICES IN INBIND.GRG27830
C                                                                       GRG27840
C     NONBAS -- INTEGER ARRAY OF DIMENSION EQUAL TO VALUE OF {NNVARS}.  GRG27850
C               ITS FIRST {NNONB} POSITIONS CONTAIN THE INDICES OF THOSEGRG27860
C               COMPONENTS OF XX WHICH ARE NOT BASIC (I.E. EITHER       GRG27870
C               SUPERBASIC OR NONBASIC) AT TERMINATION. THE REMAINING   GRG27880
C               POSITIONS CONTAIN NO USEFUL INFORMATION.                GRG27890
C                                                                       GRG27900
C     REDGR  -- REAL ARRAY OF SAME DIMENSION AS NONBAS. ITS FIRST {NNONBGRG27910
C               POSITIONS CONTAIN REDUCED GRADIENTS OF THE VARIABLES WHOGRG27920
C               INDICES ARE IN NONBAS, ORDERED IN THE SAME WAY AS NONBASGRG27930
C                                                                       GRG27940
C     NBIND  -- INTEGER SCALAR. NUMBER OF BINDING CONSTRAINTS. SEE DESCRGRG27950
C               OF INBIND AND RMULTS ABOVE.                             GRG27960
C                                                                       GRG27970
C     NNONB  -- INTEGER SCALAR. SEE NONBAS AND REDGR EXPLNANATION ABOVE GRG27980
C                                                                       GRG27990
C                                                                       GRG28000
C     TERMINATION MESSAGES                                              GRG28010
C     ++++++++++++++++++++                                              GRG28020
C                                                                       GRG28030
C                                                                       GRG28040
C     INFORM - 0 - KUHN-TUCKER CONDITIONS SATISFIED                     GRG28050
C            - 1 - FRACTIONAL CHANGE IN OBJECTIVE LESS THAN             GRG28060
C                  {EPSTOP} FOR {NSTOP} CONSECUTIVE ITERATION           GRG28070
C            - 2 - ALL REMEDIES HAVE FAILED TO FIND A BETTER            GRG28080
C                  POINT                                                GRG28090
C            - 3 - NUMBER OF COMPLETED ONE DIMENSIONAL SEARCHES         GRG28100
C                  EQUAL TO {LIMSER}                                    GRG28110
C            - 4 - SOLUTION UNBOUNDED                                   GRG28120
C            - 5 - FEASIBLE POINT NOT FOUND                             GRG28130
C                                                                       GRG28140
      IMPLICIT  REAL*8(A-H,O-Z)                                         GRG28150
      INTEGER*4 NCORE,NNVARS,NFUN,MAXBAS,MAXHES                         GRG28160
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT          GRG28170
      LOGICAL MAXIM,INPRNT,OTPRNT                                       GRG28180
      DIMENSION BLVAR(NNVARS),BUVAR(NNVARS),BLCON(NFUN),BUCON(NFUN)     GRG28190
      DIMENSION XX(NNVARS),FCNS(NFUN),RMULTS(MAXBAS)                    GRG28200
      DIMENSION RAMCON(NFUN),RAMVAR(NNVARS),INBIND(MAXBAS),Z(NCORE)     GRG28210
      DIMENSION NONBAS(NNVARS),REDGR(NNVARS)                            GRG28220
      COMMON /INOUT/ TITLE(19)                                          GRG28230
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT       GRG28240
      COMMON/GOBK/NVARS,NROWS,MAXR,MAXB                                 GRG28250
      COMMON/DYNAM1/KX,KG,KALB,KUB,KICAND,KISTAT,KIFIX                  GRG28260
      COMMON/DYNAM2/KU,KGRADF,KINBC,KIBC,KIBV,KINBV,KIUB                GRG28270
      COMMON/DYNAM3/KR,KV,KD,KGBEST,KXBEST,KXB1,KXB2,KXB3               GRG28280
      COMMON/DYNAM4/KDBND,KCNORM,KXSTAT,KGG,KRR,KY,KGRADP               GRG28290
      COMMON/DYNAM5/KROWB,KCOLB,KBINV,KGRAD,KICOLS,KINORM               GRG28300
      COMMON /DYNAM6/ KX0,KG0,KCON,KVAR,KINBVP                          GRG28310
      COMMON /INITBK/  INIT,LASTCL                                      GRG28320
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR             GRG28330
      COMMON/REV/IREST,IPER,IDUMP                                       GRG28340
      COMMON/MAXBK/MAXIM                                                GRG28350
      COMMON/INFBK/INFO                                                 GRG28360
      DIMENSION DEFAUL(19)                                              GRG28370
      DIMENSION TTITLE(19)                                              GRG28380
      CALL INITLZ                                                       GRG28390
C                                                                       GRG28400
C     VALIDITY CHECK FOR VARIABLES                                      GRG28410
C                                                                       GRG28420
      IF (NNVARS.GE.1) GO TO 20                                         GRG28430
      WRITE (6,10) NNVARS                                               GRG28440
   10 FORMAT ( 28H0FATAL ERROR - NNVARS LT 1 (,I5,  1H))                GRG28450
      GO TO 570                                                         GRG28460
   20 IF (NFUN.GE.1) GO TO 40                                           GRG28470
      WRITE (6,30) NFUN                                                 GRG28480
   30 FORMAT ( 26H0FATAL ERROR - NFUN LT 1 (,I5,  1H))                  GRG28490
      GO TO 570                                                         GRG28500
   40 IF (MAXBAS.GE.1) GO TO 60                                         GRG28510
      WRITE (6,50) MAXBAS                                               GRG28520
   50 FORMAT ( 28H0FATAL ERROR - MAXBAS LT 1 (,I5,  1H))                GRG28530
      GO TO 570                                                         GRG28540
   60 IF (MAXHES.GE.0) GO TO 80                                         GRG28550
      WRITE (6,70) MAXHES                                               GRG28560
   70 FORMAT ( 28H0FATAL ERROR - MAXHES LT 0 (,I5,  1H))                GRG28570
      GO TO 570                                                         GRG28580
   80 IF (1.LE.NNOBJ.AND.NNOBJ.LE.NFUN) GO TO 100                       GRG28590
      WRITE (6,90) NNOBJ                                                GRG28600
   90 FORMAT ( 38H0FATAL ERROR - NNOBJ LT 1 OR GT NFUN (,I5,  1H))      GRG28610
      GO TO 570                                                         GRG28620
  100 DO 120 I=1,NNVARS                                                 GRG28630
         IF (BUVAR(I).GE.BLVAR(I)) GO TO 120                            GRG28640
         WRITE (6,110) I,BUVAR(I),BLVAR(I)                              GRG28650
  110 FORMAT ( 35H0FATAL ERROR - BUVAR(I) LT BLVAR(I),/,  5H I = ,I5, 12GRG28660
     1H BUVAR(I) = ,E15.8, 12H BLVAR(I) = ,E15.8)                       GRG28670
         GO TO 570                                                      GRG28680
  120 CONTINUE                                                          GRG28690
      DO 140 I=1,NFUN                                                   GRG28700
         IF (BUCON(I).GE.BLCON(I).OR.I.EQ.NOBJ) GO TO 140               GRG28710
         WRITE (6,130) I,BUCON(I),BLCON(I)                              GRG28720
  130 FORMAT ( 35H0FATAL ERROR - BUCON(I) LT BLCON(I),/,  5H I = ,I5, 12GRG28730
     1H BUCON(I) = ,E15.8, 12H BLVAR(I) = ,E15.8)                       GRG28740
         GO TO 570                                                      GRG28750
  140 CONTINUE                                                          GRG28760
      IF (DEFAUL(1).EQ.1.0) GO TO 160                                   GRG28770
      IF (FPNEWT.GT.0.0) GO TO 160                                      GRG28780
      WRITE (6,150) FPNEWT                                              GRG28790
  150 FORMAT ( 28H0FATAL ERROR - FPNEWT LE 0 (,E15.8,  1H))             GRG28800
      GO TO 570                                                         GRG28810
  160 IF (DEFAUL(2).EQ.1.0) GO TO 180                                   GRG28820
      IF (FPINIT.GE.0.0) GO TO 180                                      GRG28830
      WRITE (6,170) FPINIT                                              GRG28840
  170 FORMAT ( 28H0FATAL ERROR - FPINIT LT 0 (,E15.8,  1H))             GRG28850
      GO TO 570                                                         GRG28860
  180 IF (DEFAUL(3).EQ.1.0) GO TO 200                                   GRG28870
      IF (FPSTOP.GT.0.0) GO TO 200                                      GRG28880
      WRITE (6,190) FPSTOP                                              GRG28890
  190 FORMAT ( 28H0FATAL ERROR - FPSTOP LE 0 (,E15.8,  1H))             GRG28900
      GO TO 570                                                         GRG28910
  200 IF (DEFAUL(4).EQ.1.0) GO TO 220                                   GRG28920
      IF (FPSPIV.GT.0.0) GO TO 220                                      GRG28930
      WRITE (6,210) FPSPIV                                              GRG28940
  210 FORMAT ( 28H0FATAL ERROR - FPSPIV LE 0 (,E15.8,  1H))             GRG28950
      GO TO 570                                                         GRG28960
  220 IF (DEFAUL(5).EQ.1.0) GO TO 240                                   GRG28970
      IF (PPH1EP.GE.0.0) GO TO 240                                      GRG28980
      WRITE (6,230) PPH1EP                                              GRG28990
  230 FORMAT ( 28H0FATAL ERROR - PPH1EP LT 0 (,E15.8,  1H))             GRG29000
      GO TO 570                                                         GRG29010
  240 IF (DEFAUL(6).EQ.1.0) GO TO 260                                   GRG29020
      IF (NNSTOP.GE.1) GO TO 260                                        GRG29030
      WRITE (6,250) NNSTOP                                              GRG29040
  250 FORMAT ( 28H0FATAL ERROR - NNSTOP LT 1 (,I5,  1H))                GRG29050
      GO TO 570                                                         GRG29060
  260 IF (DEFAUL(7).EQ.1.0) GO TO 280                                   GRG29070
      IF (IITLIM.GE.1) GO TO 280                                        GRG29080
      WRITE (6,270) IITLIM                                              GRG29090
  270 FORMAT ( 28H0FATAL ERROR - IITLIM LT 1 (,I5,  1H))                GRG29100
      GO TO 570                                                         GRG29110
  280 IF (DEFAUL(8).EQ.1.0) GO TO 300                                   GRG29120
      IF (LLMSER.GE.1) GO TO 300                                        GRG29130
      WRITE (6,290) LLMSER                                              GRG29140
  290 FORMAT ( 28H0FATAL ERROR - LLMSER LT 1 (,I5,  1H))                GRG29150
      GO TO 570                                                         GRG29160
  300 IF (DEFAUL(9).EQ.1.0) GO TO 320                                   GRG29170
      IF (0.LE.IIPR.AND.IIPR.LE.6) GO TO 320                            GRG29180
      WRITE (6,310) IIPR                                                GRG29190
  310 FORMAT ( 34H0FATAL ERROR - IIPR LT 0 OR GT 6 (,I5,  1H))          GRG29200
      GO TO 570                                                         GRG29210
  320 IF (DEFAUL(10).EQ.1.0) GO TO 340                                  GRG29220
  340 IF (DEFAUL(11).EQ.1.0) GO TO 360                                  GRG29230
      IF (IIPN5.GE.0) GO TO 360                                         GRG29240
      WRITE (6,350) IIPN5                                               GRG29250
  350 FORMAT ( 27H0FATAL ERROR - IIPN5 LT 0 (,I5,  1H))                 GRG29260
      GO TO 570                                                         GRG29270
  360 IF (DEFAUL(12).EQ.1.0) GO TO 380                                  GRG29280
      IF (IIPN6.GE.0) GO TO 380                                         GRG29290
      WRITE (6,370) IIPN6                                               GRG29300
  370 FORMAT ( 27H0FATAL ERROR - IIPN6 LT 0 (,I5,  1H))                 GRG29310
      GO TO 570                                                         GRG29320
  380 IF (DEFAUL(13).EQ.1.0) GO TO 400                                  GRG29330
      IF (IIPER.GE.0) GO TO 400                                         GRG29340
      WRITE (6,390) IIPER                                               GRG29350
  390 FORMAT ( 27H0FATAL ERROR - IIPER LT 0 (,I5,  1H))                 GRG29360
      GO TO 570                                                         GRG29370
  400 IF (DEFAUL(14).EQ.1.0) GO TO 420                                  GRG29380
      IF (IIDUMP.GE.0) GO TO 420                                        GRG29390
      WRITE (6,410) IIDUMP                                              GRG29400
  410 FORMAT ( 28H0FATAL ERROR - IIDUMP LT 0 (,I5,  1H))                GRG29410
      GO TO 570                                                         GRG29420
  420 IF (DEFAUL(15).EQ.1.0) GO TO 440                                  GRG29430
      IF (IIQUAD.EQ.0.OR.IIQUAD.EQ.1) GO TO 440                         GRG29440
      WRITE (6,430) IIQUAD                                              GRG29450
  430 FORMAT ( 33H0FATAL ERROR - IIQUAD NE 0 OR 1 (,I5,  1H))           GRG29460
      GO TO 570                                                         GRG29470
  440 IF (DEFAUL(16).EQ.1.0) GO TO 460                                  GRG29480
      IF (LDERIV.EQ.0.OR.LDERIV.EQ.1.OR.LDERIV.EQ.2) GO TO 460          GRG29490
      WRITE (6,450) LDERIV                                              GRG29500
  450 FORMAT ( 38H0FATAL ERROR - LDERIV NE 0 OR 1 OR 2 (,I5,  1H))      GRG29510
      GO TO 570                                                         GRG29520
  460 IF (DEFAUL(17).EQ.1.0) GO TO 480                                  GRG29530
      IF (MMODCG.GE.1.AND.MMODCG.LE.5) GO TO 480                        GRG29540
      WRITE (6,470) MMODCG                                              GRG29550
  470 FORMAT ( 36H0FATAL ERROR - MMODCG GT 5 OR LT 1 (,I5,  1H))        GRG29560
      GO TO 570                                                         GRG29570
  480 CONTINUE                                                          GRG29580
      DO 490 I=1,19                                                     GRG29590
  490 TITLE(I)=TTITLE(I)                                                GRG29600
      NVARS=NNVARS                                                      GRG29610
      NROWS=NFUN                                                        GRG29620
      MAXR=MAXHES                                                       GRG29630
      MAXB=MAXBAS                                                       GRG29640
      IPER=IIPER                                                        GRG29650
      IDUMP=IIDUMP                                                      GRG29660
      MAXIM=.TRUE.                                                      GRG29670
      IF (DEFAUL(19).EQ.1.0) MAXIM=.FALSE.                              GRG29680
      CALL SETUP (NCORE)                                                GRG29690
      INIT=1                                                            GRG29700
      NOBJ=NNOBJ                                                        GRG29710
      CALL SETDAT (NNVARS,NFUN,BLVAR,BUVAR,BLCON,BUCON,XX,DEFAUL,FPNEWT,GRG29720
     1FPINIT,FPSTOP,FPSPIV,PPH1EP,NNSTOP,IITLIM,LLMSER,IIPR,IIPN4,IIPN5,GRG29730
     2IIPN6,IIPER,IIDUMP,IIQUAD,LDERIV,MMODCG,RAMCON,RAMVAR,Z(KX),Z(KALBGRG29740
     3),Z(KUB),Z(KISTAT),Z(KIFIX),Z(KCON),Z(KVAR),Z(KX0),Z,NCORE,INPRNT,GRG29750
     4MAXHES)                                                           GRG29760
      IF (INPRNT) GO TO 500                                             GRG29770
      ITIPR3=IPR3                                                       GRG29780
      IPR3=-2                                                           GRG29790
  500 CONTINUE                                                          GRG29800
      CALL TABLIN (Z(KX),Z(KALB),Z(KUB),Z(KISTAT),Z(KIFIX),Z(KCON),Z(KVAGRG29810
     1R),Z(KGBEST),Z(KG),Z(KX0),Z(KG0)  )                               GRG29820
      IF (.NOT.INPRNT) IPR3=ITIPR3                                      GRG29830
      CALL REPORT (Z(KG),Z(KX),MP1,N,Z(KCON),Z(KVAR),Z(KX0))            GRG29840
      INIT=0                                                            GRG29850
      CALL GRGITN (Z(KBINV),Z(KGRAD),Z(KR),Z(KALB),Z(KUB),Z(KX),Z(KGRADFGRG29860
     1),Z(KG),Z(KV),Z(KD),Z(KU),Z(KGBEST),Z(KXBEST),Z(KINBV),Z(KIUB),Z(KGRG29870
     2INBC),Z(KROWB),Z(KCOLB),Z(KIBC),Z(KIBV),Z(KISTAT),Z(KXSTAT),Z(KIFIGRG29880
     3X),Z(KXB1),Z(KXB2),Z(KXB3),Z(KGRADP),Z(KDBND),Z(KCNORM),Z(KGG),Z(KGRG29890
     4RR),Z(KICOLS),Z(KINORM),Z(KY),Z(KINBVP),Z(KX0),Z,NCORE,Z(KICAND) )GRG29900
      IF (.NOT.OTPRNT) GO TO 510                                        GRG29910
      CALL OUTRES (Z(KG),Z(KX),Z(KINBV),Z(KIBV),Z(KALB),Z(KIUB),Z(KUB),ZGRG29920
     1(KIBC),Z(KU),Z(KINBC),Z(KGRADF),Z(KISTAT),Z(KGG),Z(KCON),Z(KVAR),ZGRG29930
     2(KIFIX),Z(KGRADP),Z(KX0),Z(KG0),Z(KGBEST),Z(KD),                  GRG29940
     3 Z(KDBND),Z(KICOLS)  )                                            GRG29950
      CALL REPORT(Z(KG),Z(KX),MP1,N,Z(KCON),Z(KVAR),Z(KX0) )            GRG29960
  510 CONTINUE                                                          GRG29970
C                                                                       GRG29980
C     SET OUTPUT ARGUMENTS TO FINAL VALUES FROM Z ARRAY                 GRG29990
C                                                                       GRG30000
      DO 520 I=1,NNVARS                                                 GRG30010
         XX(I)=Z(I)                                                     GRG30020
  520 CONTINUE                                                          GRG30030
      J=0                                                               GRG30040
      IEND=KALB-1                                                       GRG30050
      DO 530 I=KG,IEND                                                  GRG30060
         J=J+1                                                          GRG30070
         FCNS(J)=Z(I)                                                   GRG30080
  530 CONTINUE                                                          GRG30090
      DO 540 I=1,NB                                                     GRG30100
         J=KU+I-1                                                       GRG30110
         RMULTS(I)=Z(J)                                                 GRG30120
         IF (MAXIM) RMULTS(I)=-RMULTS(I)                                GRG30130
         J=KIBC+I-1                                                     GRG30140
         INBIND(I)=IGTSUB(Z(KIBC),I)                                    GRG30150
  540 CONTINUE                                                          GRG30160
      DO 550 I=1,N                                                      GRG30170
         NONBAS(I)=0                                                    GRG30180
         REDGR(I)=0.0                                                   GRG30190
  550 CONTINUE                                                          GRG30200
      J=0                                                               GRG30210
      DO 560 I=1,N                                                      GRG30220
         NINDEX=IGTSUB(Z(KINBV),I)                                      GRG30230
         IF (NINDEX.GT.NVARS) GO TO 560                                 GRG30240
         J=J+1                                                          GRG30250
         NONBAS(J)=NINDEX                                               GRG30260
         L=KGRADF+I-1                                                   GRG30270
         REDGR(J)=Z(L)                                                  GRG30280
         IF (MAXIM) REDGR(J)=-REDGR(J)                                  GRG30290
  560 CONTINUE                                                          GRG30300
      NNONB=J                                                           GRG30310
      NBIND=NB                                                          GRG30320
      INFORM=INFO                                                       GRG30330
      RETURN                                                            GRG30340
  570 CONTINUE                                                          GRG30350
      STOP                                                              GRG30360
C                                                                       GRG30370
C++++++++++++++++++++     END OF GRGSUB   +++++++++++++++++             GRG30380
C                                                                       GRG30390
      END                                                               GRG30400
      FUNCTION IGTSUB(IA,INDX)                                          GRG30410
      DIMENSION IA(1)                                                   GRG30420
      IGTSUB = IA(INDX)                                                 GRG30430
      RETURN                                                            GRG30440
C                                                                       GRG30450
      END                                                               GRG30460
      SUBROUTINE SETDAT(NNVARS,NFUN,BLVAR,BUVAR,BLCON,BUCON,XX,         GRG30470
     1 DEFAUL,FPNEWT,FPINIT,FPSTOP,FPSPIV,PPH1EP,NNSTOP,IITLIM,         GRG30480
     2 LLMSER,IIPR,IIPN4,IIPN5,IIPN6,IIPER,IIDUMP,IIQUAD,LDERIV,        GRG30490
     3 MMODCG,RAMCON,RAMVAR,X,ALB,UB,ISTAT,IFIX,CON,VAR,X0,Z,           GRG30500
     4 NCORE,INPRNT,MAXHES)                                             GRG30510
      IMPLICIT  REAL*8(A-H,O-Z)                                         GRG30520
      LOGICAL MAXIM,INPRNT                                              GRG30530
      INTEGER*4 NCORE,NNVARS,NFUN,MAXBAS,MAXHES                         GRG30540
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT,LASTZ    GRG30550
      DIMENSION BLVAR(NNVARS),BUVAR(NNVARS),BLCON(NFUN),BUCON(NFUN)     GRG30560
      DIMENSION RAMCON(NFUN),RAMVAR(NNVARS),XX(NNVARS)                  GRG30570
      DIMENSION Z(NCORE)                                                GRG30580
      DIMENSION DEFAUL(19)                                              GRG30590
      COMMON /DIMEN/ M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT     GRG30600
      DIMENSION CON(1),VAR(1),X0(1)                                     GRG30610
      DIMENSION X(1),ALB(1),UB(1)                                       GRG30620
      DIMENSION ISTAT(1),IFIX(1)                                        GRG30630
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR             GRG30640
      COMMON /TOLS/ EPS,PLINFY,PLZERO,TOLX,TOLZ                         GRG30650
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM                        GRG30660
      COMMON/MNGRG/EPSTOP,LIMSER,NSTOP,IERR,IPN4,IPN5,IPN6              GRG30670
      COMMON /QUADBK/ A1,A2,A3,ICON,IQUAD                               GRG30680
      COMMON /INOUT/ TITLE(19)                                          GRG30690
      COMMON /PARDAT/ KDERIV                                            GRG30700
      COMMON /MISC/ MAXR,NSEAR,JP,LV,JQQ                                GRG30710
      COMMON /GOBK/ NVARS,NROWS,MAXRR,MAXB                              GRG30720
      COMMON /CGBK/ MODCG                                               GRG30730
      COMMON /SETIN/ LASTZ                                              GRG30740
      COMMON/REV/IREST,IPER,IDUMP                                       GRG30750
      COMMON/PH1BK/PHMULT,PH1EPS,INITPH                                 GRG30760
      COMMON/INGRG/EPINIT                                               GRG30770
      COMMON/MAXBK/MAXIM                                                GRG30780
      DATA BLANK/8H        /                                            GRG30790
      IREST=0                                                           GRG30800
      MAXB=NBMAX                                                        GRG30810
      IF (.NOT.INPRNT) GO TO 10                                         GRG30820
      WRITE (6,200) NVARS,NROWS,MAXRR,MAXB                              GRG30830
      WRITE (6,280) LASTZ                                               GRG30840
   10 CONTINUE                                                          GRG30850
C                                                                       GRG30860
C     INITIALIZE BOUNDS ET AL. TO DEFAULT VALUES                        GRG30870
C                                                                       GRG30880
      DO 20 I=1,N                                                       GRG30890
         IFIX(I)=0                                                      GRG30900
         X(I)=0.                                                        GRG30910
         ALB(I)=-PLINFY                                                 GRG30920
         UB(I)=PLINFY                                                   GRG30930
   20 CONTINUE                                                          GRG30940
C                                                                       GRG30950
C     DEFAULT CONSTRAINT STATUS IS INEQUALITY GREATER THAN 0             GRG30960
C                                                                       GRG30970
      DO 30 I=1,MP1                                                     GRG30980
         ISTAT(I)=2                                                     GRG30990
         ALB(N+I)=0.                                                    GRG31000
         UB(N+I)=PLINFY                                                 GRG31010
   30 CONTINUE                                                          GRG31020
C                                                                       GRG31030
C     CHECK DEFAULT PROGRAM PARAMETERS                                  GRG31040
C                                                                       GRG31050
      EPBOUN=1.0E-4                                                     GRG31060
      MODCG=1                                                           GRG31070
      EPNEWT=1.0E-4                                                     GRG31080
      EPSTOP=1.0E-4                                                     GRG31090
      EPSPIV=0.001                                                      GRG31100
      EPINIT=0.0                                                        GRG31110
      PH1EPS=0.                                                         GRG31120
      NSTOP=3                                                           GRG31130
      ITLIM=10                                                          GRG31140
      LIMSER=10000                                                      GRG31150
      IPR=1                                                             GRG31160
      IPN5=0                                                            GRG31170
      IPN4=0                                                            GRG31180
      IPN6=0                                                            GRG31190
      IDUMP=0                                                           GRG31200
      IPER=0                                                            GRG31210
      IQUAD=0                                                           GRG31220
      KDERIV=0                                                          GRG31230
      MAXR=MAXRR                                                        GRG31240
      IF (DEFAUL(17).NE.1.0) MODCG=MMODCG                               GRG31250
      IF (DEFAUL(1).NE.1.0) EPNEWT=FPNEWT                               GRG31260
      IF (DEFAUL(3).NE.1.0) EPSTOP=FPSTOP                               GRG31270
      IF (DEFAUL(4).NE.1.0) EPSPIV=FPSPIV                               GRG31280
      IF (DEFAUL(2).NE.1.0) EPINIT=FPINIT                               GRG31290
      IF (DEFAUL(5).NE.1.0) PH1EPS=PPH1EP                               GRG31300
      IF (DEFAUL(6).NE.1.0) NSTOP=NNSTOP                                GRG31310
      IF (DEFAUL(7).NE.1.0) ITLIM=IITLIM                                GRG31320
      IF (DEFAUL(8).NE.1.0) LIMSER=LLMSER                               GRG31330
      IF (DEFAUL(9).NE.1.0) IPR=IIPR                                    GRG31340
      IF (DEFAUL(10).NE.1.0) IPN4=IIPN4                                 GRG31350
      IF (DEFAUL(11).NE.1.0) IPN5=IIPN5                                 GRG31360
      IF (DEFAUL(12).NE.1.0) IPN6=IIPN6                                 GRG31370
      IF (DEFAUL(14).NE.1.0) IDUMP=IIDUMP                               GRG31380
      IF (DEFAUL(13).NE.1.0) IPER=IIPER                                 GRG31390
      IF (DEFAUL(15).NE.1.0) IQUAD=IIQUAD                               GRG31400
      IF (DEFAUL(16).NE.1.0) KDERIV=LDERIV                              GRG31410
      IPR3=IPR-1                                                        GRG31420
C                                                                       GRG31430
C     ASSIGN BOUNDS OF VARIABLES                                        GRG31440
C                                                                       GRG31450
      DO 40 I=1,N                                                       GRG31460
         ALB(I)=BLVAR(I)                                                GRG31470
         UB(I)=BUVAR(I)                                                 GRG31480
         IFIX(I)=0                                                      GRG31490
         IF (ALB(I).EQ.UB(I)) IFIX(I)=1                                 GRG31500
   40 CONTINUE                                                          GRG31510
C                                                                       GRG31520
C     ASSIGN BOUNDS OF CONSTRAINTS                                      GRG31530
C                                                                       GRG31540
      DO 50 I=1,MP1                                                     GRG31550
         IF (I.EQ.NOBJ) GO TO 50                                        GRG31560
         ALB(N+I)=BLCON(I)                                              GRG31570
         UB(N+I)=BUCON(I)                                               GRG31580
         ISTAT(I)=2                                                     GRG31590
         IF (ALB(N+I).EQ.-PLINFY.AND.UB(N+I).EQ.PLINFY) ISTAT(I)=0      GRG31600
         IF (ALB(N+I).EQ.UB(N+I)) ISTAT(I)=1                            GRG31610
   50 CONTINUE                                                          GRG31620
C                                                                       GRG31630
C     SET MAXR=0, IF CONJUGATE GRADIENT METHOD IS USED                  GRG31640
C                                                                       GRG31650
      IF (MAXHES.EQ.0) MAXR=0                                           GRG31660
C                                                                       GRG31670
C     ASSIGN INITIAL VARIABLE VALUES                                    GRG31680
C                                                                       GRG31690
      DO 60 I=1,N                                                       GRG31700
   60 X(I)=XX(I)                                                        GRG31710
      IF (INPRNT) WRITE (6,190) (X(I),I=1,N)                            GRG31720
C                                                                       GRG31730
C     CHECK NAMES OF VARIABLES AND FUNTIONS                             GRG31740
C                                                                       GRG31750
      IF (DEFAUL(18).EQ.1.0) GO TO 90                                   GRG31760
      DO 70 I=1,N                                                       GRG31770
   70 VAR(I)=BLANK                                                      GRG31780
      DO 80 I=1,MP1                                                     GRG31790
   80 CON(I)=BLANK                                                      GRG31800
      GO TO 120                                                         GRG31810
   90 CONTINUE                                                          GRG31820
      DO 100 I=1,N                                                      GRG31830
  100 VAR(I)=RAMVAR(I)                                                  GRG31840
      DO 110 I=1,MP1                                                    GRG31850
  110 CON(I)=RAMCON(I)                                                  GRG31860
  120 CONTINUE                                                          GRG31870
C                                                                       GRG31880
C     CHECK TO SEE THAT VARIABLES ARE WITHIN BOUNDS                     GRG31890
C                                                                       GRG31900
      DO 140 I=1,N                                                      GRG31910
         X0(I)=X(I)                                                     GRG31920
         IF ((X(I).GE.ALB(I)).AND.(X(I).LE.UB(I))) GO TO 140            GRG31930
         IF (X(I).GT.UB(I)) GO TO 130                                   GRG31940
         WRITE (6,230) I,X(I),ALB(I)                                    GRG31950
         X(I)=ALB(I)                                                    GRG31960
         GO TO 140                                                      GRG31970
  130    WRITE (6,240) I,X(I),UB(I)                                     GRG31980
         X(I)=UB(I)                                                     GRG31990
  140 CONTINUE                                                          GRG32000
      ISTAT(NOBJ)=0                                                     GRG32010
      ALB(N+NOBJ)=-PLINFY                                               GRG32020
      UB(N+NOBJ)=PLINFY                                                 GRG32030
      IF (EPINIT.EQ.0.0) EPINIT=EPNEWT                                  GRG32040
      IF (.NOT.INPRNT) GO TO 150                                        GRG32050
      WRITE (6,170) EPNEWT,EPINIT,EPSTOP,EPSPIV,PH1EPS                  GRG32060
      WRITE (6,210) NSTOP,ITLIM,LIMSER                                  GRG32070
      WRITE (6,220) IPR,IPN4,IPN5,IPN6,IPER,IDUMP                       GRG32080
      IF (IQUAD.EQ.1) WRITE (6,160)                                     GRG32090
      IF (IQUAD.EQ.0) WRITE (6,180)                                     GRG32100
      IF (KDERIV.EQ.0) WRITE (6,270)                                    GRG32110
      IF (KDERIV.EQ.1) WRITE (6,260)                                    GRG32120
      IF (KDERIV.EQ.2) WRITE (6,250)                                    GRG32130
      IF (MAXIM) WRITE (6,290)                                          GRG32140
      IF (.NOT.MAXIM) WRITE (6,300)                                     GRG32150
      WRITE (6,310) MAXR                                                GRG32160
  150 CONTINUE                                                          GRG32170
  160 FORMAT (80H0QUADRATIC EXTRAPOLATION FOR INITIAL ESTIMATES OF BA   GRG32180
     1 VARIABLES WILL BE USED. )                                        GRG32190
  170 FORMAT (9H0EPNEWT =,E12.4,2X,8HEPINIT =,E12.4,2X,8HEPSTOP =,E12.4,GRG32200
     12X,7HEPPIV =,E12.4,2X,8HPH1EPS =,E12.5)                           GRG32210
  180 FORMAT (71H0TANGENT VECTORS WILL BE USED FOR INITIAL ESTIMATES    GRG32220
     1BASIC VARIABLES )                                                 GRG32230
  190 FORMAT (1X,8(E15.7))                                              GRG32240
  200 FORMAT (24H0NUMBER OF VARIABLES IS ,I5/24H0NUMBER OF FUNCTIONS IS GRG32250
     1,I5/42H0SPACE RESERVED FOR HESSIAN HAS DIMENSION ,I5/33H0LIMIT ON GRG32260
     2BINDING CONSTRAINTS IS ,I5)                                       GRG32270
  210 FORMAT (8H0NSTOP =,I5,2X,7HITLIM =,I5,2X,8HLIMSER =,I10)          GRG32280
  220 FORMAT (6H0IPR =,I5,2X,5HPN4 =,I5,2X,5HPN5 =,I5,2X,5HPN6 =,I5,2X,5GRG32290
     1HPER =,I5,2X,6HDUMP =,I5)                                         GRG32300
  230 FORMAT (17H FOR SUBSCRIPT = ,I5,27H INITIAL VARIABLE VALUE OF ,E15GRG32310
     1.8,30H WAS CHANGED TO LOWER BOUND = ,E15.8)                       GRG32320
  240 FORMAT (17H FOR SUBSCRIPT = ,I5,27H INITIAL VARIABLE VALUE OF ,E15GRG32330
     1.8,30H WAS CHANGED TO UPPER BOUND = ,E15.8)                       GRG32340
  250 FORMAT (47H0THE USER'S OWN PARSH SUBROUTINE WILL BE USED. )       GRG32350
  260 FORMAT (68H0THE FINITE DIFFERENCE PARSH USING CENTRAL DIFFERENC   GRG32360
     1WILL BE USED )                                                    GRG32370
  270 FORMAT (68H0THE FINITE DIFFERENCE PARSH USING FORWARD DIFFERENC   GRG32380
     1WILL BE USED )                                                    GRG32390
  280 FORMAT (29H0ACTUAL LENGTH OF Z ARRAY IS ,I6)                      GRG32400
  290 FORMAT (39H0OBJECTIVE FUNCTION WILL BE MAXIMIZED. )               GRG32410
  300 FORMAT (39H0OBJECTIVE FUNCTION WILL BE MINIMIZED. )               GRG32420
  310 FORMAT (21H0LIMIT ON HESSIAN IS ,I5)                              GRG32430
      RETURN                                                            GRG32440
C                                                                       GRG32450
C     END OF SETDAT                                                     GRG32460
C                                                                       GRG32470
      END                                                               GRG32480
C
C
C     GRG2-1
C
C
      SUBROUTINE CHUZQ(GRAD,BINV,V,IBC,D,INBV,ALB,UB,COLB,X,IBV,IUB,
     1 ICAND)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 GRAD
      INTEGER *4M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON /DIMEN/ M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      DIMENSION V(NBMAX), IBC(NBMAX), D(N), INBV(N), ALB(NPMP1), UB(NPMP
     11), COLB(NBMAX), X(NPMP1)
      DIMENSION IBV(M), IUB(N)
      DIMENSION ICAND(NPNBMX)
      DIMENSION BINV(NBMAX,NNBMAX),GRAD(MP1,N)
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON/EPSCOM/EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      COMMON/BESTBK/STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON/MISC/MAXR,NSEAR,JP,LV,JQQ
      COMMON/TOLS/EPS,PLINFY,PLZERO,TOLX,TOLZ
C
C     CHOOSE MAXIMUM PIVOT FROM COLUMNS OF B2.
C     D(1),...,D(NSUPER) WILL HOLD THE PIVOT VECTOR.
C
      IF (IPR.GE.5) WRITE (6,210)
      JQ=1
      IF (NSUPER.EQ.1) GO TO 70
C
C     COMPUTE PIVOT ROW IN B2
C
      DO 10 I=1,NB
10        V(I)=BINV(LV,I)
      PIVOT=0.0D0
      DO 50 I=1,NSUPER
          SUM=0.0D0
          II=INBV(I)
          DO 40 J=1,NB
              JJ=IBC(J)
              IF (II.GT.N) GO TO 20
              TS=GRAD(JJ,II)
              GO TO 30
C
C     SLACK COLUMN
C
20            CONTINUE
              TS=-1.0D0
              IF (II-N.NE.JJ) GO TO 40
30            SUM=SUM+V(J)*TS
40        CONTINUE
          D(I)=SUM
          IF (PIVOT.GE.DABS(SUM)) GO TO 50
          PIVOT=DABS(SUM)
          JQ=I
50    CONTINUE
      IF (IPR.GE.5) WRITE (6,190) JQ,PIVOT
C
C     CHOOSE ONE AWAY FROM ITS BOUNDS IF POSSIBLE.
C
      TOL=0.1D0*PIVOT
      DMAX=0.0D0
      JQ2=0
      DO 60 J=1,NSUPER
          IF (DABS(D(J)).LT.TOL) GO TO 60
          K=INBV(J)
          XJ=X(K)
          D1=DABS(XJ-ALB(K))
          D2=DABS(UB(K)-XJ)
          IF (D1.GT.D2) D1=D2
          IF (DMAX.GT.D1) GO TO 60
          DMAX=D1
          JQ2=J
 60    CONTINUE
      IF (JQ2.GT.0) JQ=JQ2
C
C     NOW PIVOT
C
70    ICOL=INBV(JQ)
      IF (IPR.GE.3) WRITE (6,200) ICOL,JQ
      IF (ICOL.GT.N) GO TO 90
C
C     SELECT COLUMN FROM GRAD ARRAY
C
      DO 80 I=1,NB
          II=IBC(I)
80        V(I)=GRAD(II,ICOL)
      GO TO 110
90    CONTINUE
C
C     SLACK COLUMN
C
      K=ICOL-N
      DO 100 I=1,NB
          II=IBC(I)
          V(I)=0.0D0
          IF (II.EQ.K) V(I)=-1.0D0
100   CONTINUE
110   CONTINUE
      DO 130 I=1,NB
          SUM=0.0D0
          DO 120 J=1,NB
120           SUM=SUM+V(J)*BINV(I,J)
130       COLB(I)=SUM
      PIVOT=1.0D0/COLB(LV)
      IF (IPR.GE.5) WRITE (6,230) PIVOT
      DO 140 I=1,NB
140       BINV(LV,I)=BINV(LV,I)*PIVOT
      DO 160 I=1,NB
          IF (I.EQ.LV) GO TO 160
          CONST=COLB(I)
          DO 150 J=1,NB
              BINV(I,J)=BINV(I,J)-CONST*BINV(LV,J)
150       CONTINUE
160   CONTINUE
C
C     UPDATE INDEX SETS OF BASIC AND NONBASIC VARIABLES AND IUB
C
      JQ2=IBV(LV)
      JQ3=ICAND(LV)
      IBV(LV)=ICOL
      IF (ICOL.LE.N) GO TO 163
      ICOL = ICOL-N
      DO 162 I = 1,NB
         IF ( IBC(I) .NE. ICOL ) GO TO 162
         ICOL = N + I
         GO TO 163
  162 CONTINUE
  163 CONTINUE
      ICAND(LV)=ICOL
      IF (JQ.EQ.NSUPER) GO TO 168
      K=NSUPER-1
      DO 165 I = JQ,K
      II=I+NB
      ICAND(II)=ICAND(II+1)
165       INBV(I)=INBV(I+1)
168   CONTINUE
      IUB(NSUPER)=0
      ICAND(NB+NSUPER)=JQ3
      IF (DABS(X(JQ2)-UB(JQ2)).LE.EPNEWT) IUB(NSUPER)=1
      IF (DABS(X(JQ2)-ALB(JQ2)).LE.EPNEWT) IUB(NSUPER)=-1
      INBV(NSUPER)=JQ2
      IF (IUB(NSUPER).NE.0) NSUPER=NSUPER-1
      IF (IPR.LT.4) RETURN
      WRITE (6,240) (INBV(I),I=1,N)
      WRITE (6,250) (IUB(I),I=1,N)
	WRITE (6,260)(IBV(I),I=1,NB)
      WRITE (6,270) (ICAND(I),I=1,II)
      IF (IPR.GE.5) WRITE (6,220)
      RETURN
C
180   FORMAT (29H0PIVOT IN CHUZQ TOO SMALL , =,E13.6)
190   FORMAT(12H CHUZQ JQ =   ,I5,8H PIVOT = ,E13.6)
200   FORMAT (9H VARIABLE,I4,33H ENTERING BASIS - SUPERBASIC NO.   ,I4)
210   FORMAT (15H CHUZQ ENTERED   )
220   FORMAT (17H CHUZQ COMPLETED   )
230   FORMAT (13H BINV PIVOT =,E13.6)
240   FORMAT (8H INBV IS,25I4/(1X,25I4))
250   FORMAT (8H IUB IS ,25I4/(1X,25I4))
260   FORMAT (8H IBV IS ,25I4/(1X,25I4))
270   FORMAT (9H ICAND IS,25I4/(9X,25I4))
C
C     END OF CHUZQ
C
      END
       SUBROUTINE CHUZR (ALB,UB,X,DX,G,IBV,JP,MOVE)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL MOVE
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON/TOLS/EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON/EPSCOM/EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      DIMENSION ALB(NPMP1), UB(NPMP1), X(NPMP1), IBV(M), DX(NBMAX)
      DIMENSION  G(MP1)
      IF (IPR.GE.5) WRITE (6,200)
      MOVE=.TRUE.
      THETA=PLINFY
      PSI=PLINFY
      PERTBN = EPNEWT
      JP=0
      DO 100 J = 1,NB
          T=-DX(J)
          IF (DABS(T).LE.EPS) GO TO 100
          K=IBV(J)
          IF (T.LT.0.0D0) GO TO 50
          D=X(K)-ALB(K) + PERTBN
          GO TO 60
50        D=X(K)-UB(K) - PERTBN
60        CONTINUE
          IF(DABS(D) .GT. 1.0D20) D=DSIGN(1.0D20,D)
          T = D/T
          IF (PSI.GT.T) PSI=T
          JP=J
100   CONTINUE
      IF (JP.EQ.0) GO TO 160
      IF (IPR.GE.5) WRITE (6,220) PSI,(DX(I),I=1,NB)
C
C  SECOND PASS OF HARRIS
C
      TMAX=0.0D0
      DO 150 J = 1,NB
          T=-DX(J)
          IF (DABS(T).LT.EPS) GO TO 150
          K=IBV(J)
          IF (T.LT.0.0D0) GO TO 120
          D = X(K)-ALB(K)
          GO TO 130
120       D=X(K)-UB(K)
130       CONTINUE
          IF(DABS(D) .GT. 1.0D20) D=DSIGN(1.0D20,D)
          TR=D/T
          IF (IPR.GE.5) WRITE (6,230) K,TR
          IF (TR.GT.PSI) GO TO 150
          IF (TMAX.GT.DABS(T)) GO TO 150
          TMAX=DABS(T)
          JP=J
      THETA = TR
          IF (IPR.GE.5) WRITE (6,240) JP,TMAX
150   CONTINUE
      JR=IBV(JP)
      IF(JR.LE.N) GO TO 2000
      BNDJR=ALB(JR)
      PIVOT = DX(JP)
      IF (PIVOT.GT.0.0D0) BNDJR=UB(JR)
      XJR=X(JR)
      IF (JR.GT.N) XJR=G(JR-N)
      THETA=(BNDJR-XJR)/PIVOT
2000  CONTINUE
      MOVE = (THETA .GT. TOLX )
160   IF (IPR.GE.5) WRITE (6,210)
      RETURN
200   FORMAT (15H CHUZR ENTERED   )
210   FORMAT (17H CHUZR COMPLETED   )
220   FORMAT (6H PSI =,E13.6,6H DX IS/(1X,10E13.6))
230   FORMAT (4H K =,I5,5H TR =,E13.6)
240   FORMAT (5H JP =,I5,7H TMAX =,E13.6)
C
C     END OF CHUZR
C
      END
      SUBROUTINE ADDCOL
     1 (R,Y,NY)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 NY
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      DIMENSION R(NRTOT)
      DIMENSION Y(NY)
      NSPREV=NSUPER-1
      DIAG=1.0D0
      DO 100 I = 1,NSUPER
100       Y(I)=0.0D0
C
C  INSERT NEW COLUMN OF R
C
      Y(NSUPER)=DIAG
      K=NSPREV*NSUPER/2
      DO 200 I = 1,NSUPER
200       R(K+I)=Y(I)
      RETURN
C
C     END OF ADDCOL
C
      END
      SUBROUTINE DELCOL(R,JQ)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,Mp1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      DIMENSION R(NRTOT)
      COMMON/TOLS/EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON/EPSCOM/EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
C
C  DELETE THE JQ'TH COLUMN OF UPPER TRIANGULAR R
C
      NS=NSUPER
      IF (JQ.GE.NS) GO TO 60
      K=JQ*(JQ+1)/2
      I1=JQ+1
      DO 50 I = I1,NS
          K=K+I
          T1=R(K-1)
          T2=R(K)
          D=DSQRT(T1*T1 + T2*T2)
          R(K-1) = D
          IF (I.EQ.NS) GO TO 20
          CS=T1/D
          SN=T2/D
          J1=I+1
          K1=K+I
          DO 10 J = J1,NS
              T1=R(K1-1)
              T2=R(K1)
              R(K1-1) = CS*T1 + SN*T2
              R(K1) = SN*T1 - CS*T2
              K1=K1+J
10        CONTINUE
20        K1=I-1
          J2=K-I
          J1=J2-I+2
          DO 30 J = J1,J2
30            R(J)=R(J+K1)
50    CONTINUE
60    NSUPER=NSUPER-1
      RETURN
C
C
C     END OF DELCOL
C
      END
      SUBROUTINE CONSBS(BINV,GRAD,IBV,X,ALB,UB,INBV,IUB,IBC,INBC,
     1 G,IFIK,ISTAT,DBND,CNORM,GG,RR,ICOLS,INORM,IWORK,ICAND)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL *4 GRAD
      INTEGER *4M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      LOGICAL BSCHNG,EXBAS,MAXIM
      LOGICAL UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES
      COMMON /DIMEN/ M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      DIMENSION BINV(NBMAX,NNBMAX), GRAD(MP1,N)
      DIMENSION IBV(M),X(NPMP1),ALB(NPMP1),UB(NPMP1),INBV(N),IUB(N),
     1 IBC(NBMAX)
      DIMENSION INBC(M),G(MP1),IFIK(N),ISTAT(MP1),DBND(NPNBMX)
      DIMENSION CNORM(NPNBMX),IWORK(MP1),ICAND(NPNBMX)
      DIMENSION GG(MP1),RR(NBMAX),ICOLS(NPNBMX),INORM(NPNBMX)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
C***********************************************************************
C
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON /BESTBK/ STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON /MISC/ MAXR,NSEAR,JP,LV,JQQ
      COMMON /PARDAT/ KDERIV
      COMMON /EPSCOM/ EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      COMMON /TOLS/ EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON/MNGRG/EPSTOP,LIMSER,NSTOP,IERR,IPN4,IPN5,IPN6
      COMMON /COUNTS/ NFTN,NGRAD,NMINV,NNFAIL,NCALLS,NIT,NBS,NSTEPC,NDUB
      COMMON/MAXBK/MAXIM
      COMMON/BASBK/BSCHNG
      COMMON/CONGRG/NSEAR0,LVLAST
      COMMON /SRCHLG/ UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES,UNCONP
C     ******************************************************************
C
C     SUBROUTINE CONSBS CONSTRUCTS A BASIS FROM THE JACOBIAN OF
C     BINDING CONSTRAINTS USING A MODIFIED COMPLETE PIVOTING PROCEDURE.
C
C     INPUT VARIABLES ARE -
C      ISTAT  ... STATUS ARRAY FOR CONSTRAINTS
C
C     OUTPUT VARIABLES ARE -
C      NB     ... NO. OF BASIC VARIABLES
C      IBV    ... INDEX SET OF BASIC VARIABLES
C      INBV   ... INDEX SET OF NON-BASIC VARIABLES
C      IBC    ... INDEX SET OF BINDING CONSTRAINTS
C      INBC   ... INDEX SET OF NONBINDING CONSTRAINTS
C      BINV   ... INVERSE OF BASIS
C      IUB    ... BOUND INDICATOR ARRAY FOR NONBASICS
C
C     ******************************************************************
C
C     INITIALIZATIONS
C
      NPIV=0
      NBP=NB
      IF (NBP.EQ.0) GO TO 7
      DO 4 I =1,NBP
4         IWORK(I)=IBC(I)
7     CONTINUE
      IF (BSCHNG.OR.NB.EQ.0) IBASCT=0
      IF (.NOT.BSCHNG) IBASCT=IBASCT+1
      MAXCT=2*(NSUPER+1)
      BSCHNG=.FALSE.
      EXBAS=.FALSE.
      NB=0
      NNBC=0
      IF (IPR.GE.5) WRITE (6,740)
C         **************************************************************
C         THIS DETERMINES INDEX SETS OF BINDING AND NONBINDING CONSTRAIN
C     IBC AND INBC AND SETS SLACKS OF BINDING CONSTRAINTS
C         **************************************************************
      DO 50 I=1,MP1
          II=ISTAT(I)
C
C     IF IGNORED CONSTRAINT OR OBJECTIVE, SKIP
C
          IF (II.EQ.0) GO TO 50
          IF (II.EQ.2) GO TO 20
C
C     EQUALITY CONSTRAINT
C
          BL=ALB(N+I)
          IF (DABS(G(I)-BL).GE.EPNEWT) GO TO 10
          NB=NB+1
          IBC(NB)=I
          X(N+I)=BL
          GO TO 50
10        NNBC=NNBC+1
          INBC(NNBC)=I
          GO TO 50
C
C     INEQUALITY CONSTRAINT
C
20        BU=UB(N+I)
          IF (DABS(G(I)-BU).GE.EPNEWT) GO TO 30
          NB=NB+1
          IBC(NB)=I
          X(N+I)=BU
          GO TO 50
30        BL=ALB(N+I)
          IF (DABS(G(I)-BL).GE.EPNEWT) GO TO 40
          NB=NB+1
          IBC(NB)=I
          X(N+I)=BL
          GO TO 50
40        NNBC=NNBC+1
          INBC(NNBC)=I
50    CONTINUE
      IF (NB.GT.NBMAX) GO TO 500
      IF (IPR.LT.3) GO TO 70
      IF (NB.GT.0) WRITE (6,580) (IBC(I),I=1,NB)
      IF (NNBC.GT.0) WRITE (6,590) (INBC(I),I=1,NNBC)
70    CONTINUE
      IF (NB.NE.NBP) BSCHNG=.TRUE.
      IF (IBASCT.GE.MAXCT) BSCHNG=.TRUE.
      IF (LV.NE.0) BSCHNG=.TRUE.
      IF (JSTFES) BSCHNG=.TRUE.
      IF (NSEAR.EQ.NSEAR0) BSCHNG=.TRUE.
      IF (BSCHNG.OR.NB.EQ.0) GO TO 77
      DO 74 I =1,NB
          IF (IBC(I).NE.IWORK(I)) BSCHNG=.TRUE.
74    CONTINUE
77    CONTINUE
C     ******************************************************************
C     110      COMPUTE GRADIENTS OF CONSTRAINTS
      NGRAD=NGRAD+1
      OBJBST=G(NOBJ)
C
C***********************************************************************
      IF (PENLTY) SAVE = PENTRM
C***********************************************************************
C
C     ******************************************************************
      IF (KDERIV.EQ.0) CALL PARSHF(X,G,GRAD,IFIK,GG,UB)
      IF (KDERIV.EQ.1) CALL PARSHC(X,G,GRAD,IFIK,GG,IWORK)
      IF (KDERIV.EQ.2) CALL PARSH(G,X,MP1,N,GRAD)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF (.NOT. PENLTY) GO TO 1000
      PENTRM = SAVE
C
C  PERMITS USER TO EXCLUDE HARD CONSTRAINTS FROM G VECTOR
C
      IF (NHARD .LE. 0) GO TO 1000
C
C   COMPUTE GRADIENT OF OBJECTIVE INCLUDING PENALTY TERMS
C
      CONST = 1.0D0
      IF (MAXIM) CONST = -CONST
      DO 1010 I = 1,N
         SUM = 0.0D0
         DO 1020 J = 1,NHARD
             SUM = SUM + GRAD(J,I)/G(J)**2
1020     CONTINUE
         GRAD(NOBJ,I) = GRAD(NOBJ,I) - RPEN*CONST*SUM
         GRDPEN(I) = -SUM
1010  CONTINUE
1000  CONTINUE
C***********************************************************************
C
      IF (NINF.NE.0) G(NOBJ)=OBJBST
C     ******************************************************************
C
C     PRINT INITIAL GRAD ARRAY IF IPN4 < 0
C
      IF (IPN4.GE.0) GO TO 90
      WRITE (6,540)
      DO 80 I=1,MP1
80        WRITE(6,550) I,(GRAD(I,J),J=1,N)
      IPN4=0
90    CONTINUE
      IF (.NOT.MAXIM) GO TO 100
      DO 95 I = 1,N
95        GRAD(NOBJ,I)=-GRAD(NOBJ,I)
C
C     MORE INITIALIZATIONS
C
100   DO 110 I = 1,N
          ICOLS(I)=I
C
C     SET DISTANCES OF VARIABLES FROM NEAREST BOUND
C
          DD=DMIN1(X(I)-ALB(I),UB(I)-X(I))
          IF (DD.LT.EPNEWT) DD=0.0D0
          IF (DD.GT.0.1D0*PLINFY) DD=PLINFY
          DBND(I)=DD
110   CONTINUE
      IF (IPR.GE.5) WRITE (6,730) (DBND(I),I=1,N)
      NPNB=N+NB
      IF (NB.EQ.0) GO TO 430
      DO 115 I = 1,NB
          IWORK(I)=IBV(I)
          IBV(I)=0
          J=N+I
          DBND(J)=0.0D0
          CNORM(J)=1.0D0
          INORM(J)=I
115       ICOLS(J)=J
      IF (.NOT.BSCHNG) GO TO 120
      NCAND=N+NB
      DO 117 I =1,NCAND
117       ICAND(I)=I
C
C     FILL IN BINV WITH GRADIENTS OF BINDING CONSTRAINTS
C
120   CONTINUE
      IF (IPR.GE.5) WRITE (6,690) (ICOLS(J),J=1,NPNB)
      IF (IPR.GE.5) WRITE (6,530)
      DO 140 I=1,NB
          II=IBC(I)
          DO 130 J=1,N
130           BINV(I,J)=GRAD(II,J)
          IF (IPR.GE.5) WRITE (6,520) (BINV(I,J),J=1,N)
140   CONTINUE
C     UU IS RELATIVE PIVOT TOLERANCE
      UU=0.01D0
C
C     COUNT NUMBER OF SLACKS IN OLD BASIS IF THERE IS NO BASIS CHANGE
C
      NS=0
      IF (NCAND.GT.NB) GO TO 150
      DO 145 I=1,NB
          IF (ICAND(I).GT.N) NS=NS+1
145   CONTINUE
C
C     START MAIN PIVOT LOOP.  NO RETURNS TO ABOVE HERE
C
150   CONTINUE
C
C     IF NO SLACKS LEFT TO PIVOT ON IN OLD BASIS THEN GO TO REGULAR PIVO
C
      IF (NS.EQ.0) GO TO 160
C
C     OTHERWISE PIVOT ON SLACKS FIRST SINCE NO CHANGE TO BINV.
C
      DO 155 I =1,NB
          JPIV=ICAND(I)
          IF (JPIV.LE.N) GO TO 155
          NS=NS-1
          ICAND(I)=0
      GO TO 270
155   CONTINUE
160   CONTINUE
C
C     CASE 1 WHEN NCAND=NB, REGULAR PIVOTS
C
      IF (NCAND.NE.NB) GO TO 185
C
C     FOR UNPIVOTED ROWS AND UNPIVOTED COLUMNS OF BINV, FIND LARGEST
C     ABSOLUTE ELEMENT.
C
      ELTMAX=0.0D0
      DO 180 JJJ=1,NCAND
          J=ICAND(JJJ)
          IF (J.EQ.0) GO TO 180
          DO 170 I=1,NB
              IF (IBV(I).NE.0) GO TO 170
              C=DABS(BINV(I,J))
              IF (C.LT.ELTMAX) GO TO 170
              ELTMAX=C
              JPIV=J
              INORM(J)=I
              JJJT=JJJ
170       CONTINUE
180   CONTINUE
C
C     PIVOT ON ELTMAX IF LARGE ENOUGH
C
      IF (ELTMAX.LT.EPSPIV) GO TO 183
      ICAND(JJJT)=0
      GO TO 270
183   CONTINUE
C
C     OTHERWISE SET NCAND TO N+NB AND GO TO 185
C
      NCAND=N+NB
      BSCHNG=.TRUE.
      EXBAS=.TRUE.
      IF (IPR.GE.4) WRITE (6,760) ELTMAX
185   CONTINUE
      ELTMAX=0.0D0
C
C     COMPUTE MAX ABS ELEMENT FOR UNPIVOTED COLUMNS
C
      DO 200 J = 1,N
C
C     DETERMINE MAX ABS ELEMENT IN UNPIVOTED ROWS OF JCOL
C
          IMAX=0
          CMAX=0.0D0
          DO 190 I=1,NB
C
C     PICK UNPIVOTED ROWS
C
              IF (IBV(I).NE.0) GO TO 190
              C=DABS(BINV(I,J))
              IF (C.LT.CMAX) GO TO 190
              CMAX=C
              IMAX=I
190       CONTINUE
          CNORM(J)=CMAX
          INORM(J)=IMAX
          IF (CMAX.LT.ELTMAX) GO TO 200
          ELTMAX=CMAX
200   CONTINUE
      IF (IPR.LT.5) GO TO 210
      WRITE (6,700) (CNORM(J),J=1,NPNB)
      WRITE (6,710) (INORM(J),J=1,NPNB)
C
C     CHOOSE PIVOT COLUMN, JPIV
C     PICK COLUMN FARTHEST FROM BOUND PROVIDED COLUMN NORM IS NOT TOO
C     SMALL.
C
210   CONTINUE
      DMAX=0.0D0
      CMAX=0.0D0
      JMAX=0
      IF (ELTMAX .LT. EPSPIV) GO TO 235
      DO 230 JCOL=1,N
          C=CNORM(JCOL)
C
C     SKIP IF COLUMN NORM TOO SMALL
C
          IF (C.LT.UU*ELTMAX) GO TO 230
          DD=DBND(JCOL)
          IF (DD.LT.DMAX) GO TO 230
          IF (DD.GT.DMAX) GO TO 220
C
C     DISTANCE FROM BOUNDS EQUAL, PICK LARGEST COLUMN NORM
C
          IF (C.LE.CMAX) GO TO 230
          CMAX=C
          JMAX=JCOL
          GO TO 230
C
C     NEW COLUMN FARTHER FROM BOUND
C
220       DMAX=DD
          CMAX=C
          JMAX=JCOL
230   CONTINUE
235   CONTINUE
      IF (IPR.GE.5) WRITE (6,720) JMAX,CMAX,ELTMAX
      JPIV=JMAX
C
C     IF ALL COLUMNS TOO SMALL PICK SLACK IN UNPIVOTED ROW.
      IF (JPIV.EQ.0) JPIV=INORM(1)+N
240   CONTINUE
C
C     IF NOT AT BOUND HAVE SATISFACTORY PIVOT COLUMN
C
      IF (DMAX.NE.0.0D0) GO TO 270
C
C
C     PIVOT COLUMN AT BOUND.  COMPUTE LARGEST FREE COLUMN.
C
      CMAX=0.0D0
      JMAX=0
      DO 260 J=1,N
          IF (DBND(J).EQ.0.0D0) GO TO 260
          C=CNORM(J)
          IF (C.LT.CMAX) GO TO 260
          CMAX=C
          JMAX=J
260   CONTINUE
C
C     PIVOT ON ABOVE FREE COLUMN IF NORM GREATER THAN PIVOT TOLERANCE
C
      IF (CMAX.LT.EPSPIV) GO TO 265
C
C     FOUND FREE COLUMN SO USE IT
      JPIV=JMAX
      GO TO 270
C
C     NO FREE COLUMN.  PICK BEST OF UNPIVOTED SLACK OR CURRENT COLUMN
C     AT BOUND
C
265   IF (CNORM(JPIV).GE.1.0D0) GO TO 270
C
C     SLACK COLUMN BIGGER SO CHOOSE IT
C
      IF (JPIV.GT.N) GO TO 266
      IF (N+INORM(JPIV).EQ.LVLAST) GO TO 270
      JPIV=N+INORM(JPIV)
      GO TO 270
266   CONTINUE
      DO 268 I = 1,NB
          IF (IBV(I).NE.0) GO TO 268
          JPIV=I+N
          IF (JPIV.NE.LVLAST) GO TO 270
268   CONTINUE
270   CONTINUE
C
C     PIVOT TO UPDATE BINV
C
      IROW=INORM(JPIV)
      JCOL=JPIV
      IF (JPIV.GT.N) GO TO 340
C
C     PIVOT COLUMN IS NON-SLACK OR UPDATED SLACK
C
290   CONTINUE
      PIV=BINV(IROW,JPIV)
      IF (IPR.GE.5) WRITE (6,620) IROW
      IF (IPR.GE.5) WRITE (6,630) JPIV
      IF (IPR.GE.5) WRITE (6,560) PIV
      DO 300 I=1,NB
          RR(I)=BINV(I,JPIV)
          BINV(I,JPIV)=0.0D0
300   CONTINUE
      IF (IPR.GE.5) WRITE (6,570) (RR(I),I=1,NB)
      BINV(IROW,JPIV)=1.0D0
      DO 310 J=1,N
310       BINV(IROW,J)=BINV(IROW,J)/PIV
      DO 330 I=1,NB
          IF (I.EQ.IROW) GO TO 330
          R=RR(I)
          DO 320 J=1,N
              BINV(I,J)=BINV(I,J)-R*BINV(IROW,J)
320       CONTINUE
330   CONTINUE
C
C     UPDATE COLUMN STATUS BY SWAPPING PIVOT COLUMN WITH UPDATED SLACK
C     PLUS OTHER BOOKKEEPING FOR CASE OF PIVOT IN NON-SLACK OR UPDATED
C     SLACK COLUMN
C
      DBND(JPIV)=0.0D0
      J=ICOLS(JPIV)
      IBV(IROW)=J
      I=N+IROW
      ICOLS(JPIV)=ICOLS(I)
      ICOLS(I)=J
      GO TO 360
340   CONTINUE
C
C     NON-UPDATED SLACK COLUMN PIVOT LOGIC
C
      DO 350 J=1,N
          BINV(IROW,J)=-BINV(IROW,J)
350   CONTINUE
      IBV(IROW)=JPIV
360   CONTINUE
      IF (IPR.LT.6) GO TO 368
      WRITE (6,690) (ICOLS(J),J=1,NPNB)
      WRITE (6,680) (IBV(J),J=1,NB)
      WRITE (6,530)
      DO 365 I = 1,NB
365       WRITE (6,520) (BINV(I,J),J=1,N)
368   CONTINUE
      NPIV=NPIV+1
C
C     IF HAVE NOT MADE NB PIVOTS GET ANOTHER PIVOT COLUMN
C
      IF (NPIV.LT.NB) GO TO 150
C
C     ALL PIVOTS DONE.  NOW REARRANGE COLUMNS TO GET INVERSE
C
C
C     NOW LOOK FOR SLACK COLUMNS WHICH WERE UPDATED BY PIVOT IN ROW J,
C     J=1,...,NB SINCE THESE WILL FORM BINV.
C
      DO 381 J = 1,NB
          JJ=N+J
C
C     PROCESS NON SLACK COLUMNS FIRST
C
          DO 380 I = 1,N
              II=ICOLS(I)
              IF (II.NE.JJ) GO TO 380
C
C     NOW SWAP COLUMN I AND J
C
              DO 375 K = 1,NB
                  TMP=BINV(K,J)
                  BINV(K,J) = BINV(K,I)
                  BINV(K,I)=TMP
375           CONTINUE
              ICOLS(I)=ICOLS(J)
              ICOLS(J)=0
              GO TO 381
380       CONTINUE
381   CONTINUE
      IF (IPR3.GE.4) WRITE (6,690) (ICOLS(K),K=1,NPNB)
C
C     NOW TO PROCESS SLACKS
C
      DO 400 J = 1,NB
          JJ=N+J
          DO 390 I = 1,NB
              II=N+I
              II=ICOLS(II)
              IF (II.NE.JJ) GO TO 390
C
C     NOW HAVE FOUND SLACK. GET + OR - UNIT VECTOR.
C
              DO 382 K = 1,NB
                  BINV(K,J)=0.0D0
382           CONTINUE
C
C     IF UPDATED SLACK, GET NEGATIVE UNIT VECTOR.
C
              IF (I.NE.J) GO TO 385
              BINV(J,J)=-1.0D0
              GO TO 400
385           CONTINUE
C
C     UPDATED SLACK -- GET A UNIT VECTOR
C
              BINV(I,J)=1.0D0
390       CONTINUE
400   CONTINUE
      IF (IPR.LT.5) GO TO 430
      WRITE (6,530)
      DO 425 I = 1,NB
425   WRITE (6,520) (BINV(I,J),J=1,NB)
C
C     SET UP INDEX SET OF NONBASIC VARIABLES, SUPERBASICS FIRST AND
C     NONBASIC SLACKS LAST
C
430   CONTINUE
      IBCTR=N
      DO 440 I=1,NPNB
440       ICOLS(I)=0
      IF (NB.EQ.0) GO TO 448
      DO 442 I = 1,NB
          J=IBV(I)
442       ICOLS(J)=1
      NCAND=NB
      DO 444 I =1,NCAND
444       ICAND(I)=IBV(I)
C
C     REINDEX SLACKS IN IBV
C
      DO 446 I = 1,NB
          J=IBV(I)
          IF (J.GT.N) J=N+IBC(J-N)
          IBV(I)=J
446   CONTINUE
448   CONTINUE
C
C     DO NONBASIC SLACKS
C
      ICTRB=N+NB
      IF (NB.EQ.0) GO TO 470
      DO 460 I=1,NB
          K=N+I
          IF (ICOLS(K).EQ.1) GO TO 460
          ICAND(ICTRB)=K
          ICTRB=ICTRB-1
          K=N+IBC(I)
          INBV(IBCTR)=K
          IBCTR=IBCTR-1
460   CONTINUE
C
C     HANDLE REST OF NONBASICS
C
470   CONTINUE
      ICSUPE=NB
      NSUPER=0
      DO 490 I=1,N
          IF (ICOLS(I).EQ.1) GO TO 490
          IF (DBND(I).GT.0.0D0) GO TO 480
          INBV(IBCTR)=I
          IBCTR=IBCTR-1
          ICAND(ICTRB)=I
          ICTRB=ICTRB-1
          GO TO 490
480       CONTINUE
C
C     SUPERBASIC VARIABLES
C
          ICSUPE=ICSUPE+1
          ICAND(ICSUPE)=I
          NSUPER=NSUPER+1
          INBV(NSUPER)=I
490   CONTINUE
492   CONTINUE
      DO 495 I = 1,N
          IUB(I)=0
          K=INBV(I)
          IF (DABS(X(K)-UB(K)).LT.EPNEWT) IUB(I)=1
          IF (DABS(X(K)-ALB(K)).LT.EPNEWT) IUB(I)=-1
495   CONTINUE
      LVLAST=0
      IF (.NOT.EXBAS) GO TO 498
      DO 497 I=1,NB
          II=IBV(I)
          DO 496 J=1,NB
              IF (IWORK(J).EQ.II) GO TO 497
496       CONTINUE
          GO TO 498
497   CONTINUE
      BSCHNG=.FALSE.
      IF (IPR.GE.4) WRITE (6,770)
498   CONTINUE
      IF (IPR3.LT.2) RETURN
      WRITE (6,660) (INBV(I),I=1,N)
      WRITE (6,670)(IUB(I),I=1,N)
      IF (NB.GT.0) WRITE (6,680) (IBV(I),I=1,NB)
      IF (IPR.GE.5) WRITE (6,750)
      RETURN
C
C     NUMBER OF BINDING CONSTRAINTS TOO LARGE
C
500   WRITE (6,640) NB
      STOP
C
520   FORMAT (1X,8D16.7)
530   FORMAT (8H BINV IS)
540   FORMAT (15H0GRAD ARRAY IS )
550   FORMAT (7H0FOR G(,I3,1H)/(1X,1P10E13.5))
560   FORMAT (6H PIV =,E13.6)
570   FORMAT (5H RR =,10E12.5)
580   FORMAT (25H BINDING CONSTRAINTS ARE   ,25I4/(1X,30I4))
590   FORMAT(28H NONBINDING CONSTRAINTS ARE   ,20I4/(1X,30I4))
600   FORMAT (I5,1X,20HBINDING CONSTRAINTS )
610   FORMAT (I5,24H NONBINDING CONSTRAINTS   )
620   FORMAT (14H PIVOT ROW IS ,I5)
630   FORMAT (17H PIVOT COLUMN IS ,I5)
640   FORMAT (64H0DIMENSIONS WILL BE EXCEEDED BY NUMBER OF BINDING CONST
     1RAINTS = ,I5)
660   FORMAT(8H INBV IS,30I4/(1X,30I4))
670   FORMAT(8H IUB IS ,30I4/(1X,30I4))
680   FORMAT(7H IBV IS ,30I4/(1X,30I4))
690   FORMAT(9H ICOLS IS ,25I4/(1X,30I4))
700   FORMAT (10H CNORM IS ,6E13.6/(1X,10E13.6))
710   FORMAT(9H INORM IS  ,25I4/(1X,30I4))
720   FORMAT (7H JMAX =,I5,2X,6HCMAX =,E13.6,2X,8HELTMAX =,E13.6)
730   FORMAT(7H DBND =,8E15.7/(7X,8E15.7))
740   FORMAT (16H CONSBS ENTERED   )
750   FORMAT (18H CONSBS COMPLETED   )
760   FORMAT (48H BASIS CANDIDATE LIST EXPANDED BECAUSE ELTMAX =  ,
     1 E12.4)
770   FORMAT (22H BASIS SAME AFTER ALL   )
C
C     END OF CONSBS
C
      END
      SUBROUTINE DIREC (GRADF,R,GRADFP,IUB,INBV,D,ISTAT,Y,U,
     1 X,ALB,UB,INBVP)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES
      LOGICAL MOVE,RESTRT,DROP,VARMET,CONJGR
      LOGICAL MODR,BSCHNG,UPDATE
      LOGICAL DFAIL
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      DIMENSION GRADF(N), R(NRTOT), GRADFP(N), IUB(N), INBV(N), D(N),
     1 ISTAT(MP1),Y(N)
      DIMENSION X(NPMP1),ALB(NPMP1),UB(NPMP1)
      DIMENSION U(NBMAX)
      DIMENSION INBVP(N)
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON /EPSCOM/ EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      COMMON /TOLS/ EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON /SRCHLG/ UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES,UNCONP
      COMMON /LOGBLK/ MOVE,RESTRT,DROP,VARMET,CONJGR
      COMMON/MISC/MAXR,NSEAR,JP,LV,JQQ
      COMMON/MNGRG/EPSTOP,LIMSER,NSTOP,IERR,IPN4,IPN5,IPN6
      COMMON/BESTBK/STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON/DIRGRG/COND,UPDATE,NSUPP
      COMMON/BASBK/BSCHNG
      COMMON/DFBLK/DFAIL
      IF (IPR.GE.5) WRITE (6,340)
      UPDATE=.FALSE.
      DFAIL=.FALSE.
      VARMET=(NSUPER.LE.MAXR)
      CONJGR=(.NOT.VARMET)
      IF (JSTFES) RESTRT=.TRUE.
      IF (NSUPER.EQ.0) GO TO 145
      IF (CONJGR) GO TO 110
C
C     COMPLEMENTARY DFP VARIABLE METRIC ALGORITHM
C
      IF (NSEAR.EQ.0) GO TO 60
      IF (BSCHNG) GO TO 50
      IF (RESTRT) GO TO 60
      MOVE=(STEP.GT.EPS)
C
C     MSGCG CAUSES CG SUBROUTINE TO PRINT ON NEXT CALL
C
      IF (IPR3.GE.0) MSGCG=1
      MODR=.FALSE.
      IF (.NOT.MOVE) GO TO 70
      YTP=0.0D0
      GTP=0.0D0
      GTP1=0.0D0
      DO 10 I=1,NSUPER
          P=D(I)
          Y(I)=GRADF(I)-GRADFP(I)
          YTP=Y(I)*P+YTP
          GTP=GRADFP(I)*P+GTP
          GTP1=GRADF(I)*P+GTP1
10    CONTINUE
C
C     USE COMDFP UPDATE ONLY IF GTP1/GTP < 0.9(SAY).  (NOTE THAT GTP<0)
C
      MODR=(GTP1.GT.0.99D0*GTP)
      IF (.NOT.MODR.AND.IPR.GE.4) WRITE (6,290) GTP1,GTP
      GO TO 70
C
C     HESSIAN UPDATE WHEN BASIC VARIABLE HITS BOUND
C
50    CONTINUE
      IF (NSUPP.GT.MAXR.AND.IPR3.GE.0) WRITE (6,310)
      CALL RESETR (R,COND)
      GO TO 80
C
C     RESET HESSIAN
C
60    CONTINUE
      CALL RESETR (R,COND)
      GO TO 80
C
C     HESSIAN UPDATE WHEN NO VARIABLE HITS BOUND
C
70    CONTINUE
      IF (.NOT.MODR) GO TO 80
      KADD=4
      KSUB=4
      CALL COMDFP (R,Y,GRADFP,NSUPER,STEP,YTP,GTP,KADD,KSUB)
      UPDATE=.TRUE.
C
C     COMPUTE SEARCH DIRECTION,D
C
80    CONTINUE
      DO 90 J=1,NSUPER
90        D(J)=-GRADF(J)
      CALL RTRSOL (R,D,NSUPER)
C
C     COMPUTE CONDITION NUMBER OF DIAGONAL OF R
C
      DMIN=PLINFY
      DMAX=0.0D0
      K=0
      DO 100 I=1,NSUPER
          K=K+I
          T=DABS(R(K))
          IF (DMIN.GT.T) DMIN=T
          IF (DMAX.LT.T) DMAX=T
100   CONTINUE
      COND=PLINFY
      IF (DMIN.LT.EPS) GO TO 120
      COND=(DMAX/DMIN)**2
      GO TO 140
C
C     CONJUGATE GRADIENT METHOD
C
110   CONTINUE
      IF (.NOT.UNCON.OR.BSCHNG) RESTRT=.TRUE.
      CALL CG(GRADFP,GRADF,D,MSGCG,INBV,IUB)
C
C     CHECK IF DIRECTION IS DOWNHILL
C
120   CONTINUE
      SUM=0.0D0
      DO 130 I=1,NSUPER
130       SUM=SUM+D(I)*GRADF(I)
      IF (SUM.LT.-EPS) GO TO 145
C
C     BAD DIRECTION.  RESET
C
      IF (RESTRT) GO TO 235
      IF (IPR3.GE.2) WRITE (6,300)
      RESTRT=.TRUE.
      IF (VARMET) GO TO 60
      GO TO 110
140   CONTINUE
      IF (IPR3.LT.5) GO TO 145
      K=NSUPER*(NSUPER+1)/2
      WRITE (6,320) (R(I),I=1,K)
C
C     THIS CODE DECIDES IF ANY VARIABLES AT BOUNDS ARE TO BE RELEASED
C     FROM THEM.
C
145   CONTINUE
      DNORM=0.0D0
      IF (NSUPER.EQ.0) GO TO 155
      DO 150 I=1,NSUPER
150       DNORM=DNORM+D(I)**2
155   UNORM=1.0D0
      IF (NB.EQ.0) GO TO 170
      UNORM=0.0D0
      DO 160 I=1,NB
160       UNORM=UNORM+U(I)**2
      UNORM=DSQRT(UNORM/NB)
      IF (UNORM.EQ.0.0D0) UNORM=1.0D0
170   CONTINUE
      TOLD=EPSTOP*UNORM
      K=NSUPER+1
      IF (K.GT.N) GO TO 220
      AXMULT=0.0D0
      IDROP=0
      DO 200 I=K,N
          II=INBV(I)
          IF (II.LE.N) GO TO 190
C
C     SLACK VARIABLE.  SEE IF SLACK IS FOR EQUALITY ROW
C
          J=II-N
          IF (ISTAT(J).EQ.1) GO TO 200
190       CONTINUE
          TST=-GRADF(I)
          IF (IUB(I).EQ.1) TST=GRADF(I)
          IF (TST.LE.TOLD) GO TO 200
C
C     PROFITABLE COLUMN.  SEE IF BEST SO FAR
C
          TST=DABS(TST)
          IF (TST.LT.AXMULT) GO TO 200
          AXMULT=TST
          IDROP=I
200   CONTINUE
      IF (IPR.GE.5) WRITE (6,240) AXMULT,IDROP
C
C     DROP SET IN GRG TO FORCE DROPPING A CONSTRAINT
C
      IF ((IDROP.NE.0).AND.(DROP)) GO TO 210
      IF ((IDROP.EQ.0).OR.(DNORM.GT.AXMULT**2/4.0D0)) GO TO 215
C
C     ADD NONBASIC VARIABLE TO SUPERBASIC SET.
C     INTERCHANGE NONBASICS IN POSITION NSUPER+1 AND IDROP
C
210   SV=GRADF(IDROP)
      NSUPER=NSUPER+1
      ISV1=INBV(NSUPER)
      ISV2=IUB(NSUPER)
      II=INBV(IDROP)
      IF (IPR3.GE.2) WRITE (6,260) II,NSUPER
      SV3=GRADF(NSUPER)
      INBV(NSUPER)=INBV(IDROP)
      GRADF(NSUPER)=SV
      INBV(IDROP)=ISV1
      IUB(IDROP)=ISV2
      IUB(NSUPER)=0
      GRADF(IDROP)=SV3
      IF (IPR3.GE.2) WRITE (6,270) (INBV(I),I=1,N)
      IF (IPR3.GE.2) WRITE (6,280)(IUB(I),I=1,N)
C
C     UPDATE SEARCH DIRECTION
C
      D(NSUPER)=-SV
C
C     UPDATE HESSIAN OR RESET CG METHOD
C
      IF (VARMET) CALL ADDCOL (R,Y,N)
      IF (CONJGR) RESTRT=.TRUE.
      DFAIL=.FALSE.
      GO TO 220
215   DROP=.FALSE.
220   DO 230 I=1,N
230       GRADFP(I)=GRADF(I)
      NSUPP=NSUPER
      IF (DFAIL) GO TO 237
      IF (IPR.GE.5) WRITE (6,350)
      RETURN
C
C     NEGATIVE GRADIENT DIRECTION NOT DOWNHILL
C
235   DO 236 I=1,NSUPER
          IF (D(I).NE.0.0D0) GO TO 239
236   CONTINUE
C     DIRECTION VECTOR ZERO.  TRY DROPPING A CONSTRAINT.
      DROP=.TRUE.
      DFAIL=.TRUE.
      IF (IPR.GE.1) WRITE (6,380)
      GO TO 145
237   WRITE (6,370)
      RETURN
239   WRITE (6,330)
      DFAIL=.TRUE.
      RETURN
C
C
C
240   FORMAT (9H AXMULT =,E13.6,2X,7HIDROP =,I5)
250   FORMAT(71H IN DIREC, INDICATORS SHOW CONSTRAINED, BUT VARIABLE DID
     1 NOT HIT BOUND   )
260   FORMAT(23H BOUND ON VARIABLE NO.   ,I5,20H RELAXED.  NSUPER = ,I5)
270   FORMAT(8H INBV IS,30I4)
280   FORMAT(8H IUB IS ,30I4)
290   FORMAT(36H MODR FALSE, SKIP UPDATE OF HESSIAN   /
     1 7H GTP1 =,E13.6,6H GTP =   ,E13.6)
300   FORMAT (31H DIRECTION NOT DOWNHILL, RESET   )
310   FORMAT (27H SWITCH TO VARIABLE METRIC   )
320   FORMAT (5H R IS ,8E15.7/(1X,8E15.7))
330   FORMAT(81H0NEGATIVE GRADIENT DIRECTION NOT DOWNHILL.  CHECK DERIVA
     1TIVES AND/OR TOLERANCES.   )
340   FORMAT (15H DIREC ENTERED )
350   FORMAT (17H DIREC COMPLETED   )
360   FORMAT(45H MODR FAILS 10% TEST BUT RESET TO TRUE.GTP1= ,E10.2,
     1  7H GTP = ,E10.2)
370   FORMAT (90H0DIRECTION VECTOR ZERO AND NO CONSTRAINT COULD BE DROPP
     1ED.  KUHN-TUCKER CONDITION IMPLIED.   )
380   FORMAT (52H DIRECTION VECTOR ZERO.  TRY DROPPING A CONSTRAINT.   )
C
C     END OF DIREC
      END
      SUBROUTINE GRGITN (BINV,GRAD,R,ALB,UB,X,GRADF,G,V,D,U,GBEST,XBEST,
     1INBV,IUB,INBC,ROWB,COLB,IBC,IBV,ISTAT,XSTAT,IFIK,XB1,XB2,XB3,GRADF
     2P,DBND,CNORM,GG,RR,ICOLS,INORM,Y,INBVP,X0,Z,NCORE,ICAND)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 GRAD
      INTEGER *4IREST,NCORE,LASTZ
      LOGICAL UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES
      LOGICAL MOVE,RESTRT,DROP,VARMET,CONJGR
      LOGICAL DEGEN,DRFLAG,MAXIM,BSCHNG,UPDATE
      LOGICAL DFAIL,UNCONP
      INTEGER *4M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON /DIMEN/ M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON /NINTBK/ NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      DIMENSION ALB(NPMP1), UB(NPMP1), X(NPMP1), GRADF(N), G(MP1), V(NBM
     1AX), D(N), U(NBMAX), GBEST(MP1), XBEST(MPNBMX), INBV(N), IUB(N), I
     2NBC(M), IBC(NBMAX), IBV(M), ISTAT(MP1)
      DIMENSION XSTAT(N), IFIK(N), XB1(NBMAX), XB2(NBMAX), XB3(NBMAX), G
     1RADFP(N), Y(N)
      DIMENSION ROWB(NBMAX), COLB(NBMAX)
      DIMENSION GRAD(MP1,N), BINV(NBMAX,NNBMAX), R(NRTOT)
      DIMENSION DBND(NPNBMX), CNORM(NPNBMX), GG(MP1), RR(NBMAX), ICOLS(N
     1PNBMX), INORM(NPNBMX)
      DIMENSION INBVP(N), X0(N)
      DIMENSION Z(NCORE), ICAND(NPNBMX)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
      COMMON /USOBJ/ USEOBJ,USBEST,PENBST
C***********************************************************************
C
      COMMON /SETIN/ LASTZ
      COMMON /INOUT/ TITLE(19)
      COMMON /CGBK/ MODCG
      COMMON /PARDAT/ KDERIV
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON /EPSCOM/ EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      COMMON /TOLS/ EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON /MNGRG/ EPSTOP,LIMSER,NSTOP,IERR,IPN4,IPN5,IPN6
      COMMON /MISC/ MAXR,NSEAR,JP,LV,JQQ
      COMMON /COUNTS/ NFTN,NGRAD,NMINV,NNFAIL,NCALLS,NIT,NBS,NSTEPC,NDUB
      COMMON /SRCHLG/ UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES,UNCONP
      COMMON /LOGBLK/ MOVE,RESTRT,DROP,VARMET,CONJGR
      COMMON /BESTBK/ STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON /REV/ IREST,IPER,IDUMP
      COMMON /INGRG/ EPINIT
      COMMON /PH1BK/ PHMULT,PH1EPS,INITPH
      COMMON /DIRGRG/ COND,UPDATE,NSUPP
      COMMON /MAXBK/ MAXIM
      COMMON /BASBK/ BSCHNG
      COMMON /CONGRG/ NSEAR0,LVLAST
      COMMON/DFBLK/DFAIL
      COMMON /INFBK/ INFO
C         **************************************************************
C
C
C     INITIALIZE PERFORMANCE COUNTERS
      NCALLS=0
C     100  NCALLS =  NUMBER OF NEWTON CALLS
      NIT=0
C     110  NIT    =  CUMULATIVE NO. OF NEWTON ITERATIONS
      NFTN=1
C     120  NFTN   =  TOTAL NO. OF GCOMP CALLS
      NGRAD=0
C     130  NGRAD  =  NO OF GRADIENT CALLS
      NSEAR=0
C     140  NSEAR  =  NO OF ONE DIMENSIONAL SEARCHES
      NSEAR0=0
      ISTOP=0
C     150  ISTOP  =  NO OF CONSECUTIVE TIMES RELATIVE CHANGE IN FUNCTION
C                LESS THAN EPSTOP
      NBS=0
C     160  NBS    =  NO. OF TIMES BASIC VARIABLE VIOLATES BOUND
      NNFAIL=0
C     170  NNFAIL =  NO. TIMES NEWTON FAILED TO CONVERGE
      NSTEPC=0
C     180  NSTEPC =  NO. TIMES STEP SIZE CUT BACK WHEN NEWTON FAILS
C
C     ADJUSTMENTS FOR USING TWO CONSTRAINT TOLERANCES
C
      EPLAST=EPNEWT
      IF (EPINIT.NE.EPNEWT) EPNEWT=EPINIT
      IF (EPLAST.NE.EPNEWT) EPSTOP=10.0D0*EPSTOP
C
      STEP=0.0D0
      IPR=IPR3+1
      LINECT=60
      IPRHLD=IPR
      IPRHD3=IPR3
C
C     OTHER INITIALIZATIONS
C
10    CONTINUE
      DROP=.FALSE.
      MOVE=.TRUE.
      RESTRT=.FALSE.
      UNBD=.FALSE.
      JSTFES=.FALSE.
      DEGEN=.FALSE.
      UNCON=.FALSE.
      UNCONP = .FALSE.
      DRFLAG=.FALSE.
      BSCHNG=.TRUE.
      IDEGCT=0
      NSUPER=0
      TRUBST=0.0D0
      NSUPP=0
      IERR=0
      NFAIL=0
      MSGCG=1
      STPBST=1.0D0
      LV=0
      ISTOP=0
      OBJTST=PLINFY
      STEP=0.0D0
      COND=1.0D0
      NINF=0
      NB=0
      IF (IPER.NE.0) IPR=1
      IF (IPER.NE.0) IPR3=0
      DO 20 I=1,MP1
20        X(N+I)=G(I)
      DO 30 I=1,NBMAX
          DO 30 J=1,NNBMAX
30            BINV(I,J)=0.0D0
C
C
C         **************************************************************
C     THIS IS RETURN POINT FOR MAIN LOOP.
C     NO TRANSFERS ARE MADE TO STATEMENTS BEFORE
C
C     ------------------------------------------------------------------
40    CONTINUE
C
C     COMPUTE BASIS INVERSE, EXCEPT WHEN DEGENERATE
C
C     ------------------------------------------------------------------
      CALL CONSBS (BINV,GRAD,IBV,X,ALB,UB,INBV,IUB,IBC,INBC,G,IFIK,ISTAT
     1,DBND,CNORM,GG,RR,ICOLS,INORM,GBEST,ICAND)
C     ------------------------------------------------------------------
      IF (NINF.EQ.0.OR.PH1EPS.EQ.0.0D0) GO TO 50
      INITPH=1
      CALL PH1OBJ (INBC,G,UB,ALB)
      INITPH=0
50    CONTINUE
      IF (NSEAR.NE.NSEAR0) GO TO 100
C
C     INITIALIZATIONS THAT MUST BE DONE AFTER FIRST CONSBS CALL
C
      INITPH=2
      CALL PH1OBJ (INBC,G,UB,ALB)
      INITPH=0
      IF (NB.EQ.0) GO TO 70
      DO 60 I=1,NB
          K=IBV(I)
60        XBEST(I)=X(K)
70    IF (NNBC.EQ.0) GO TO 100
      DO 80 I=1,NNBC
          K=INBC(I)
80        XBEST(NB+I)=G(K)
      DO 90 I=1,MP1
90        GBEST(I)=G(I)
100   CONTINUE
C
C     COMPUTE REDUCED GRADIENT
C


C
C     ------------------------------------------------------------------
      CALL REDGRA (GRADF,BINV,GRAD,INBC,IBV,U,INBV,IBC,G,ALB,UB)
C     ------------------------------------------------------------------
      IF (IPR.LT.4) GO TO 140
      DO 110 I=1,N
110       XSTAT(I)=GRADF(I)
      IF (.NOT.MAXIM.OR.NINF.NE.0) GO TO 130
      DO 120 I=1,N
120       XSTAT(I)=-XSTAT(I)
130   WRITE (6,600) (XSTAT(I),I=1,N)
140   IF (NSEAR .EQ. 0 .AND..NOT. DEGEN) GO TO 320
C
C     STOP CRITERIA.  CHECK IF CURRENT POINT OPTIMAL
C     TEST IF KUHN-TUCKER CONDITIONS SATISFIED
C
150   CONTINUE
      DO 190 I=1,N
          II=INBV(I)
          TST=GRADF(I)
          IF (II.LE.N) GO TO 160
          IF (ISTAT(II-N).EQ.1) GO TO 190
160       CONTINUE
          IF (IUB(I).EQ.0) GO TO 180
          IF (IUB(I).EQ.1) GO TO 170
          IF (TST.LT.-EPSTOP) GO TO 200
          GO TO 190
170       CONTINUE
          IF (TST.GT.EPSTOP) GO TO 200
          GO TO 190
180       IF (DABS(TST).GT.EPSTOP) GO TO 200
190   CONTINUE
C     250  OPTIMALITY TEST MET.GO TO OUTPUT PHASE.
      GO TO 450
C     ------------------------------------------------------------------
C     CHECKS IF RELATIVE CHANGE IN OBJECTIVE IS LESS THAN EPSTOP FOR
C     NSTOP CONSECUTIVE ITERATIONS.
200   CONTINUE
      IF (DEGEN) GO TO 250
      IF (DABS(G(NOBJ)-OBJTST).GT.DABS(OBJTST*EPSTOP)) GO TO 210
      ISTOP=ISTOP+1
      IF (ISTOP.GE.NSTOP) GO TO 460
      IBS=0
      IF (BSCHNG) IBS=1
      GO TO 220
C     ------------------------------------------------------------------
210   ISTOP=0
      OBJTST=G(NOBJ)
220   CONTINUE
C
C     STOP IF TOO MANY LINEAR SEARCHES
C
      IF (NSEAR.GE.LIMSER) GO TO 490
C
C     COMPUTE SEARCH DIRECTION FOR SUPERBASICS
C
250   CONTINUE
      CALL DIREC (GRADF,R,GRADFP,IUB,INBV,D,ISTAT,Y,U,X,ALB,UB,INBVP)
      IF (DFAIL) GO TO 520
      IF (IPR.GE.4) WRITE (6,660) (D(I),I=1,NSUPER)
      IF ((.NOT.DRFLAG).OR.DROP) GO TO 260
      IF (NFAIL.EQ.2) GO TO 410
      GO TO 480
260   CONTINUE
      IF (NB.EQ.0) GO TO 300
C
C     COMPUTE TANGENT VECTOR V
C
      CALL TANG (BINV,GRAD,IBC,INBV,V,D,RR)
      IF (IPR.GE.4) WRITE (6,640) (V(I),I=1,NB)
C
C     FIND JP, INDEX OF FIRST BASIC VARIABLE TO HIT A BOUND
C
       CALL  CHUZR(ALB,UB,X,V,G,IBV,JP,MOVE)
      IF (MOVE) GO TO 300
C
C     DEGENERATE AND NO MOVE IN BASICS IS POSSIBLE
C
      JR=IBV(JP)
      IF (IPR.GE.3) WRITE (6,670) JR
      LV=JP
      LVLAST=LV
      DEGEN=.TRUE.
      IDEGCT=IDEGCT+1
      IF (IDEGCT.LT.15) GO TO 320
      WRITE (6,800) IDEGCT
      GO TO 480
C
C     GO TO PRINT SECTION AND RETURN
C
290   CONTINUE
C
C     EXCHANGE BASIC WITH SOME SUPERBASIC AND UPDATE BINV
C
      CALL CHUZQ (GRAD,BINV,V,IBC,D,INBV,ALB,UB,COLB,X,IBV,IUB,ICAND)
C
C     SET LOGICALS FOR USE BY DIREC
C
      RESTRT=.FALSE.
      UNCON=.FALSE.
      BSCHNG=.TRUE.
      MXSTEP=.TRUE.
C
C     NOW GO TO BEGIN NEW ITERATION FOR DEGENERATE CASE
C
      GO TO 100
300   CONTINUE
      DEGEN=.FALSE.
      DRFLAG=.FALSE.
      IDEGCT=0
C     ------------------------------------------------------------------
      CALL SEARCH (BINV,GRAD,D,X,G,IBV,V,XB1,XB2,XB3,INBV,ALB,UB,XSTAT,I
     1STAT,INBC,GBEST,XBEST,IBC,ROWB,COLB,RR,DBND)
C     ------------------------------------------------------------------
C     IF ABSOLUTE VALUE OF X'S IS VERY SMALL, CHANGE TO 0 TO AVOID
C     UNDERFLOW.
      DO 310 I=1,N
          IF (DABS(X(I)).LT.EPS) X(I)=0.0D0
310   CONTINUE
      NSEAR=NSEAR+1
      IF (NSEAR.EQ.IPN4) IPR=4
      IF (NSEAR.EQ.IPN5) IPR=5
      IF (NSEAR.EQ.IPN6) IPR=6
      IPR3=IPR-1
C
C     PRINT SECTION
C
320   CONTINUE
      IF (IPER.NE.0) GO TO 330
      IF (IPR.LT.1) GO TO 380
      LINECT=LINECT+1
      IF (LINECT.LT.48.AND.IPR3.EQ.0) GO TO 340
      IF (IPR3.EQ.0.AND.NSEAR.GT.1) WRITE (6,760)
      WRITE (6,770)
      LINECT=0
      GO TO 340
330   CONTINUE
      K=NSEAR/IPER*IPER
      IF (K.NE.NSEAR.AND.K.NE.NSEAR-1) GO TO 340
      IF (NSEAR.EQ.0) WRITE (6,770)
      IF (NSEAR.LT.2) GO TO 340
      IF (K.EQ.NSEAR-1) WRITE (6,770)
      IPR=IPRHLD
      IPR3=IPRHD3
340   GRNORM=0.0D0
      IF (NSUPER.EQ.0) GO TO 360
      DO 350 I=1,NSUPER
          IF (DABS(GRADF(I)).GT.GRNORM) GRNORM=DABS(GRADF(I))
350   CONTINUE
360   CONTINUE
      IF (NSUPER.GT.MAXR) COND=0.0D0
      IF (DEGEN) STEP=0.0D0
      IF (MAXIM.AND.NINF.EQ.0) G(NOBJ)=-G(NOBJ)
      IF (DEGEN) WRITE (6,780) NSEAR,G(NOBJ),NB,NSUPER,NINF,GRNORM,COND,
     1UPDATE,STEP,DEGEN
      IF (.NOT.DEGEN) WRITE (6,780) NSEAR,G(NOBJ),NB,NSUPER,NINF,GRNORM,
     1COND,UPDATE,STEP
C
C***********************************************************************
      IF ( PENLTY ) WRITE (6,900) USEOBJ,PENTRM
900   FORMAT(5X,15HTRUE OBJECTIVE= ,E20.6,5X,13HPENALTY TERM= ,E20.6)
C***********************************************************************
C
      IF (IPR.LT.3) GO TO 370
      WRITE (6,580)
      WRITE (6,740) (I,G(I),I=1,MP1)
      WRITE (6,590)
      WRITE (6,740) (I,X(I),I=1,N)
370   CONTINUE
      IF (MAXIM.AND.NINF.EQ.0) G(NOBJ)=-G(NOBJ)
      IF (IPER.EQ.0) GO TO 380
      IF (K.NE.NSEAR-1) GO TO 380
      IF (NSEAR.NE.1) WRITE (6,770)
      IPR=1
      IPR3=0
380   IF (DEGEN) GO TO 290
      IF (NSEAR .EQ. 0) GO TO 150
C
C     IF SUPERBASIC HAS HIT BOUND, DELETE APPROPRIATE COLUMNS OF HESSIAN
C
      IF (.NOT.MXSTEP) GO TO 400
      IF (NSUPER .EQ. 0) GO TO 400
      IF (NSUPER .EQ. 0) GO TO 400
      II=NSUPER
      III=NSUPER
      DO 390 KK=1,II
          I=II+1-KK
          J=INBV(I)
          IF (DABS(X(J)-ALB(J)).GT.EPNEWT.AND.DABS(X(J)-UB(J)).GT.EPNEWT
     1    ) GO TO 390
          III=III-1
          IF (VARMET) CALL DELCOL(R,I)
          DO 385 K=I,III
385           GRADFP(K)=GRADFP(K+1)
390   CONTINUE
400   CONTINUE
      IF (SUCCES) GO TO 440
C
C     TROUBLE--NO FUNCTION DECREASE
C     FORCE BASIS CHANGE
C
      IF (NB.EQ.0) GO TO 410
      IF (BSCHNG) GO TO 410
      NB=0
      IF (IPR.GE.1) WRITE (6,790)
      LINECT=LINECT+2
      ISTOP=ISTOP-1
      GO TO 40
C
C     TRY STEEPEST DESCENT STEP, DROPPING A CONSTRAINT, THEN GIVE UP
C
410   NFAIL=NFAIL+1
      GO TO (420,430,500), NFAIL
420   IF (RESTRT) GO TO 410
      IF (NSUPER.LT.2) GO TO 410
      RESTRT=.TRUE.
      IF (IPR.GE.1) WRITE (6,700)
      LINECT=LINECT+1
      GO TO 250
430   IF (DROP) GO TO 410
      IF (NSUPER.EQ.N) GO TO 410
      IF (IPR.GE.1) WRITE (6,690)
      LINECT=LINECT+1
      DROP=.TRUE.
      DRFLAG=.TRUE.
      GO TO 250
440   CONTINUE
      IF (UNBD) GO TO 510
      NFAIL=0
      RESTRT=.FALSE.
      DROP=.FALSE.
      GO TO 40
C
C     KUHN TUCKER CONDITIONS SATISFIED
C
450   WRITE (6,610) EPSTOP
      INFO = 0
      GO TO 480
C
C     FRACTIONAL CHANGE CRITERION SATISFIED
C
460   IF (IBS.EQ.1.OR.NB.EQ.0) GO TO 470
      IF (BSCHNG) IBS=1
      IF (BSCHNG) GO TO 220
      NB=0
      ISTOP=ISTOP-1
      GO TO 40
470   WRITE (6,630) EPSTOP,ISTOP
      INFO = 1
      LINECT=LINECT+2
      IERR=1
      IF (NSUPER.EQ.N.OR.ISTOP.GE.NSTOP+2) GO TO 480
      DRFLAG=.TRUE.
      DROP=.TRUE.
      IF (IPR.GE.1) WRITE (6,690)
      LINECT=LINECT+1
      GO TO 250
480   IF (EPNEWT.EQ.EPLAST) GO TO 520
      EPNEWT=EPLAST
      WRITE (6,750) EPNEWT
      LINECT=LINECT+2
      EPSTOP=0.1D0*EPSTOP
      NSEAR0=NSEAR
      PHHOLD=PH1EPS
      PH1EPS=0.2D0
      DO 485 I =1,N
          IF (IFIK(I).NE.0) GO TO 485
          TS=UB(I)+EPNEWT
          IF (X(I).GT.TS) X(I)=TS
          TS=ALB(I)-EPNEWT
          IF (X(I).LT.TS) X(I)=TS
485   CONTINUE
      CALL GCOMP(G,X)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
C
      IF ( .NOT. PENLTY ) GO TO 487
      CONST = 1.0D0
      IF (MAXIM) CONST = -CONST
      USEOBJ = G(NOBJ)
      G(NOBJ) = G(NOBJ) + RPEN*CONST*PENTRM
487   CONTINUE
C***********************************************************************
C
      IF (MAXIM) G(NOBJ)=-G(NOBJ)
      NFTN=NFTN+1
      GO TO 10
C
C     TOO MANY LINEAR SEARCHES
C
490   WRITE (6,650) NSEAR
      INFO = 3
      IERR=11
      GO TO 520
C
C     NO IMPROVEMENT IN LINESEARCH.  ALL REMEDIES FAILED.
C
500   IERR=2
      INFO = 2
      WRITE (6,570)
      LINECT=LINECT+2
      GO TO 480
C
C     UNBOUNDED SOLUTION
C
510   IERR=20
      INFO = 4
      WRITE (6,620) NDUB
C
C     NORMAL TERMINATION STEPS
C
520   CONTINUE
      IF (NINF.EQ.0) GO TO 540
C
C     SOLUTION INFEASIBLE
C
      INFO = 5
      WRITE (6,550) TRUOBJ,NINF
      IERR=9
      DO 530 I=1,NNBC
          J=INBC(I)
          IF ((G(J).GT.ALB(N+J)).AND.(G(J).LT.UB(N+J))) GO TO 530
          WRITE (6,560) J,G(J)
530   CONTINUE
      G(NOBJ)=TRUOBJ
540   CONTINUE
      IF (EPNEWT.NE.EPLAST) EPSTOP=0.1D0*EPSTOP
      EPNEWT=EPLAST
      PH1EPS=PHHOLD
      RETURN
C
C
C
550   FORMAT (69H0FEASIBLE POINT NOT FOUND.  VALUE OF TRUE OBJECTIVE AT
     1TERMINATION = ,1PE13.6/14H THE FOLLOWING,I4,32H CONSTRAINTS WERE I
     2N VIOLATION: )
560   FORMAT (I5,1PE13.6)
570   FORMAT (71H0ALL REMEDIES HAVE FAILED TO FIND A BETTER POINT.  PROG
     1RAM TERMINATED. )
580   FORMAT (5H G IS)
590   FORMAT (5H X IS)
600   FORMAT (21H REDUCED GRADIENT IS /,(1X,1P8D16.8))
610   FORMAT (72H0TERMINATION CRITERION MET.  KUHN-TUCKER CONDITIONS SAT
     1ISFIED TO WITHIN ,1PE12.5,1X,17HAT CURRENT POINT )
620   FORMAT (60H0SOLUTION UNBOUNDED--FUNCTION IMPROVING AFTER DOUBLING
     1STEP ,I4,6H TIMES)
630   FORMAT (48H0TOTAL FRACTIONAL CHANGE IN OBJECTIVE LESS THAN ,1PE12.
     15,1X,3HFOR,I4,1X,23HCONSECUTIVE ITERATIONS )
640   FORMAT (19H TANGENT VECTOR IS /(1X,1P10E13.5))
650   FORMAT (56H0NUMBER OF COMPLETED ONE-DIMENSIONAL SEARCHES = LIMSER
     1=,I5,28H.  OPTIMIZATION TERMINATED. )
670   FORMAT (29H0BASIS DEGENERATE--VARIABLE  ,I5,15H LEAVING BASIS )
660   FORMAT (21H DIRECTION VECTOR IS /(1X,1P10E13.5))
680   FORMAT (6H JP = ,I4,2X,7HVNORM =,1PE13.5)
690   FORMAT (27H TRY DROPPING A CONSTRAINT )
700   FORMAT (27H TRY STEEPEST DESCENT STEP )
720   FORMAT (66H0A RESTART FROM A PREVIOUS DUMP HAS BEEN MADE WITH X'S
     1 AS FOLLOWS )
740   FORMAT (5(I4,1PE14.7,3X),I4,1PE14.7)
750   FORMAT( 63H0CONSTRAINT TOLERANCE HAS BEEN TIGHTENED TO ITS FINAL V
     1ALUE OF ,1PE12.5)
760   FORMAT (1H1)
770   FORMAT (10H0ITERATION,2X,9HOBJECTIVE,3X,10HNO.BINDING,2X,9HNO.SUPE
     1R-,4X,6HNUMBER,5X,9HNORM RED.,4X,7HHESSIAN,3X,7HHESSIAN,13X,10HDEG
     2ENERATE/3X,6HNUMBER,4X,8HFUNCTION,3X,11HCONSTRAINTS,3X,6HBASICS,4X
     3,10HINFEASIBLE,3X,8HGRADIENT,3X,9HCONDITION,3X,6HUPDATE,4X,8HSTEPS
     4IZE,4X,4HSTEP)
780   FORMAT (I6,1PE18.6,I7,I12,I11,1PE15.3,1PE12.3,L5,1PE15.3,L7)
790   FORMAT (36H0TRY EXPANDING BASIS CANDIDATE LIST )
800   FORMAT (15H0DEGENERATE FOR ,I5,27H STEPS.  PROBABLY CYCLING.   )
C
C     END OF GRGITN
C
      END
      SUBROUTINE OUTRES(G,X,INBV,IBV,ALB,IUB,UB,IBC,U,INBC,GRADF,ISTAT,
     1 GG,CON,VAR,IFIK,GRADFP,X0,XBEST,GBEST,D,DBND,ICOLS)
C
C     THIS SUBROUTINE PRINTS RESULTS OF OPTIMIZATION.
C
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL MAXIM
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
C
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
      COMMON /USOBJ/ USEOBJ,USBEST,PENBST
C***********************************************************************
C
      INTEGER *4M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON /DIMEN/ M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/MAXBK/MAXIM
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON /COUNTS/ NFTN,NGRAD,NMINV,NNFAIL,NCALLS,NIT,NBS,NSTEPC,NDUB
      COMMON /BESTBK/ STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON/MNGRG/EPSTOP,LIMSER,NSTOP,IERR,IPN4,IPN5,IPN6
      COMMON /MISC/ MAXR,NSEAR,JP,LV,JQQ
      COMMON /INOUT/ TITLE(19)
      COMMON /TOLS/ EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON /PARDAT/ KDERIV
      DIMENSION G(MP1), X(NPMP1), IUB(N), INBV(N), IBV(M), ALB(NPMP1), U
     1B(NPMP1), IBC(NBMAX), U(NBMAX), INBC(MP1), GRADF(N), ISTAT(MP1),
     2 GG(MP1), ICOLS(NPNBMX), X0(N)
      DIMENSION CON(MP1), VAR(N), GBEST(MP1)
      DIMENSION IFIK(N),GRADFP(N),XBEST(MP1),D(N),DBND(NPNBMX)
      DIMENSION SYMBOL(12)
      DATA SYMBOL/8H  OBJ   ,8HIGNORED ,8HLOWERBND,8HUPPERBND,
     1 8HEQUALITY,8H  FREE  ,8H FIXED  ,8HNO BOUND,8HSUPBASIC,8HNONBASIC
     2 ,8H BASIC  ,8HVIOLATED/
      DATA IBL/2H  /,IU/2H:U/,IL/2H:L/
      IF (MAXIM) G(NOBJ)=-G(NOBJ)
      IF (IPR3.GE.-1) GO TO 10
      WRITE (6,270) G(NOBJ)
      WRITE (6,280) (X(I),I=1,N)
      GO TO 210
10    CONTINUE
      IF (.NOT.MAXIM.OR.NINF.NE.0) GO TO 18
      DO 12 I = 1,N
12        GRADF(I)=-GRADF(I)
      IF (NB.EQ.0) GO TO 18
      DO 14 I = 1,NB
14        U(I)=-U(I)
18    CONTINUE
      WRITE (6,350)
      WRITE (6,360) TITLE
      WRITE (6,370)
      DO 20 I=1,MP1
          GBEST(I)=SYMBOL(6)
20        GG(I)=PLINFY
25    CONTINUE
      IF (NB.EQ.0) GO TO 40
      DO 30 I=1,NB
          K=IBC(I)
30        GG(K)=U(I)
40    CONTINUE
      DO 70 I = 1,MP1
          NI=N+I
          TS1=UB(NI)-G(I)
          TS2=G(I)-ALB(NI)
          TS=TS1
          INBC(I)=IU
          IF (TS2.GT.TS1) GO TO 45
          TS=TS2
          INBC(I)=IL
45        CONTINUE
          X(NI)=TS
          IF (TS.LT.-EPNEWT) GO TO  65
          IF (TS.GT.EPNEWT) GO TO 60
          IF (UB(NI).EQ.ALB(NI)) GO TO 50
          IF (DABS(G(I)-UB(NI)).LT.EPNEWT) GBEST(I)=SYMBOL(4)
          IF (DABS(G(I)-ALB(NI)).LT.EPNEWT) GBEST(I)=SYMBOL(3)
          GO TO 70
50        GBEST(I)=SYMBOL(5)
          INBC(I)=IBL
          GO TO 70
60        IF (ISTAT(I).NE.0) GO TO 70
          INBC(I)=IBL
          GBEST(I)=SYMBOL(2)
          IF (I.EQ.NOBJ) GBEST(I)=SYMBOL(1)
          GO TO 70
65        GBEST(I)=SYMBOL(12)
70    CONTINUE
      DO 90 I = 1,MP1
          IF ((I/40)*40.NE.I) GO TO 75
          WRITE (6,350)
          WRITE (6,360) TITLE
          WRITE (6,370)
75        CONTINUE
          IF (I.EQ.NOBJ) GO TO 85
          IF (GG(I).NE.PLINFY) GO TO 80
          WRITE (6,380) I,CON(I),XBEST(I),G(I),GBEST(I),X(N+I),INBC(I)
          GO TO 90
80        WRITE (6,380) I,CON(I),XBEST(I),G(I),GBEST(I),X(N+I),INBC(I),
     1 GG(I)
          GO TO 90
85        WRITE (6,380) I,CON(I),XBEST(I),G(I),GBEST(I)
90    CONTINUE
C
C***********************************************************************
      IF ( PENLTY ) WRITE (6,500) USEOBJ
500   FORMAT (18H TRUE OBJECTIVE IS ,E12.5)
C***********************************************************************
C
      DO 100 I=1,N
          GRADFP(I)=PLINFY
          ICOLS(I)=IBL
          D(I)=SYMBOL(10)
          IF (IFIK(I).NE.0) GO TO 94
          TS1=UB(I)-X(I)
          TS2=X(I)-ALB(I)
          TS=TS1
          K=IU
          IF (TS2.GT.TS1) GO TO 92
          TS=TS2
          K=IL
92        CONTINUE
          IF (TS.GT.1.0D20) GO TO 96
          IF (TS.LT.EPNEWT) GO TO 98
          DBND(I)=TS
          ICOLS(I)=K
          GO TO 100
94        DBND(I)=SYMBOL(7)
          GO TO 100
96        DBND(I)=SYMBOL(8)
          GO TO 100
98        IF (DABS(UB(I)-X(I)).LT.EPNEWT) DBND(I)=SYMBOL(4)
          IF (DABS(ALB(I)-X(I)).LT.EPNEWT) DBND(I)=SYMBOL(3)
100   CONTINUE
      IF (NSUPER.EQ.0) GO TO 120
      DO 110 I=1,NSUPER
          K=INBV(I)
110       IF (K.LE.N) D(K)=SYMBOL(9)
120   CONTINUE
      IF (NB.EQ.0) GO TO 140
      DO 130 I = 1,NB
          K=IBV(I)
130       IF (K.LE.N) D(K)=SYMBOL(11)
140   CONTINUE
      DO 150 I=1,N
          K=INBV(I)
150       IF (K.LE.N) GRADFP(K)=GRADF(I)
      IF (MP1+N.LE.34) GO TO 155
      WRITE(6,350)
      WRITE (6,360) TITLE
155   WRITE (6,390)
      DO 190 I=1,N
          IF ((I/40)*40.NE.I) GO TO 158
          WRITE (6,350)
          WRITE (6,360) TITLE
          WRITE (6,390)
158       CONTINUE
          IF (GRADFP(I).EQ.PLINFY) GO TO 170
          IF (ICOLS(I).NE.IBL) GO TO 160
          WRITE (6,400) I,VAR(I),X0(I),X(I),D(I),DBND(I),GRADFP(I)
          GO TO 190
160       WRITE (6,410) I,VAR(I),X0(I),X(I),D(I),DBND(I),ICOLS(I),GRADFP
     1    (I)
          GO TO 190
170       IF (ICOLS(I).NE.IBL) GO TO 180
          WRITE (6,400) I,VAR(I),X0(I),X(I),D(I),DBND(I)
          GO TO 190
180       WRITE (6,410) I,VAR(I),X0(I),X(I),D(I),DBND(I),ICOLS(I)
190   CONTINUE
      WRITE (6,330) TITLE
210   WRITE (6,240) NSEAR
      TS=0.0D0
      IF (NCALLS.EQ.0) GO TO 220
      XNIT=NIT
      XNCALL=NCALLS
      TS=XNIT/XNCALL
220   CONTINUE
      WRITE (6,250) NCALLS,NIT,TS
      I=N*(KDERIV+1)
      IF (KDERIV.EQ.2) I=0
      I=I*NGRAD+NFTN
      WRITE (6,260) NFTN,NGRAD,I
      WRITE (6,290) NBS
      WRITE (6,300) NNFAIL,NSTEPC
      RETURN
C
C
240   FORMAT (39H0NUMBER OF ONE-DIMENSIONAL SEARCHES =  ,I5)
250   FORMAT (15H0NEWTON CALLS =,I5,2X,19HNEWTON ITERATIONS =,I5,2X,9HAV
     1ERAGE =,F8.2)
260   FORMAT (17H0FUNCTION CALLS =,I5,2X,17HGRADIENT CALLS = ,I5/46H  AC
     1TUAL FUNCTION CALLS (INC. FOR GRADIENT) = ,I6)
270   FORMAT (26H0FINAL OBJECTIVE VALUE IS ,E15.8)
280   FORMAT (31H0FINAL VALUES OF VARIABLES ARE /(1X,10E13.6))
290   FORMAT (51H0NUMBER OF TIMES BASIC VARIABLE VIOLATED A BOUND = ,I5)
300   FORMAT (45H0NUMBER OF TIMES NEWTON FAILED TO CONVERGE = ,I5/49H0TI
     1MES STEPSIZE CUT BACK DUE TO NEWTON FAILURE = ,I5)
330   FORMAT (16H1RUN STATISTICS    /1H0,19A4)
350   FORMAT (1H1,15X,13HFINAL RESULTS)
360   FORMAT (1H0,19A4)
370   FORMAT (26H0SECTION 1 -- FUNCTIONS   //49X,8HDISTANCE/
     1 16X,7HINITIAL,6X,5HFINAL,17X,4HFROM,6X,8HLAGRANGE/
     2 4H NO.,3X,4HNAME,6X,5HVALUE,7X,5HVALUE,5X,6HSTATUS,4X,
     3 7HNEAREST,4X,10HMULTIPLIER/50X,5HBOUND)
380   FORMAT (I4,1X,A8,1P2E12.5,1X,A8,1PE10.3,A2,1PE12.5)
390   FORMAT(/23H0SECTION 2 -- VARIABLES //49X,8HDISTANCE/
     1 16X,7HINITIAL,6X,5HFINAL,17X,4HFROM,6X,7HREDUCED/
     2 4H NO.,3X,4HNAME,6X,5HVALUE,7X,5HVALUE,5X,6HSTATUS,4X,
     3 7HNEAREST,5X,8HGRADIENT/50X,5HBOUND)
400   FORMAT (I4,1X,A8,1P2E12.5,1X,A8,2X,A8,2X,1PE12.5)
410   FORMAT (I4,1X,A8,1P2E12.5,1X,A8,1PE10.3,A2,1PE12.5)
C
C     END OF OUTRES
C
      END
      SUBROUTINE REPORT(G,X,MP1,N,CON,VAR,X0)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 N,MP1
      DIMENSION X(N),G(MP1)
      DIMENSION VAR(N),CON(MP1),X0(N)
      RETURN
C
C     END OF REPORT
C
      END
      SUBROUTINE TABLIN(X,ALB,UB,ISTAT,IFIK,CON,VAR,GBEST,G,X0,G0)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER *4 GBEST,NONE,LITST,ITYPE
      INTEGER*4 II
      LOGICAL MAXIM
      INTEGER *4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON /DIMEN/ M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      DIMENSION GBEST(MP1),LITST(7),ITYPE(6)
      DIMENSION CON(MP1), VAR(N)
      DIMENSION G(1),X0(1),G0(1)
      DIMENSION X(NPMP1), ALB(NPMP1), UB(NPMP1)
      DIMENSION ISTAT(MP1), IFIK(N)
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
      COMMON /USOBJ/ USEOBJ,USBEST,PENBST
C***********************************************************************
C
      COMMON /TOLS/ EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON /INOUT/ TITLE(19)
      COMMON/MAXBK/MAXIM
      DATA LITST/3HUL ,3HLL ,3HEQ ,3H***,3H   ,4HFREE,3HFX /
      DATA ITYPE/4HEQ  ,4HLE  ,4HGE  ,4HRNGE,4HOBJ ,4HNA  /
      DATA NONE/4HNONE/
      CALL GCOMP (G,X)
      DO 10 I = 1,MP1
         G0(I) = G(I)
   10 CONTINUE
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( .NOT. PENLTY ) GO TO 100
      IF ( VIOL ) RETURN
      CONST = 1.0D0
      IF (MAXIM) CONST = -CONST
      USEOBJ = G(NOBJ)
      G(NOBJ) = G(NOBJ) + RPEN*CONST*PENTRM
100   CONTINUE
C***********************************************************************
C
      IF (IPR3.LT.-1) GO TO 935
C
C     PRINT INITIAL TABLES
C
      WRITE (6,1480) TITLE
      WRITE (6,1490)
      DO 825 I = 1,MP1
          NI=N+I
          GBEST(I)=ITYPE(4)
          IF (ALB(NI).EQ.UB(NI)) GBEST(I)=ITYPE(1)
          IF (ALB(NI).NE.-PLINFY.AND.UB(NI).NE.PLINFY) GO TO 825
          IF (ALB(NI).EQ.-PLINFY.AND.UB(NI).EQ.PLINFY) GO TO 820
          IF (ALB(NI).NE.-PLINFY) GBEST(I)=ITYPE(3)
          IF (UB(NI).NE.PLINFY) GBEST(I)=ITYPE(2)
          GO TO 825
820       GBEST(I)=ITYPE(6)
825   CONTINUE
      GBEST(NOBJ)=ITYPE(5)
      DO 890 I=1,MP1
          NI=N+I
          GI=G(I)
          TEMP=CON(I)
          II=LITST(5)
          IF (ISTAT(I).EQ.1) GO TO 840
          IF (ISTAT(I).EQ.0) GO TO 830
          IF (DABS(GI-UB(NI)).LT.EPNEWT) II=LITST(1)
          IF (DABS(ALB(NI)-GI).LT.EPNEWT) II=LITST(2)
          IF (GI.GT.UB(NI)+EPNEWT.OR.GI.LT.ALB(NI)-EPNEWT) II=LITST(4)
          GO TO 860
830       WRITE (6,1500) I,TEMP,II,GBEST(I),G(I)
          GO TO 890
840       IF (DABS(GI-ALB(NI)).LT.EPNEWT) GO TO 850
          II=LITST(4)
          GO TO 860
850       II=LITST(3)
860       CONTINUE
          IF (UB(NI).EQ.PLINFY) GO TO 870
          IF (ALB(NI).EQ.-PLINFY) GO TO 880
          WRITE (6,1500) I,TEMP,II,GBEST(I),G(I),ALB(NI),UB(NI)
          GO TO 890
870       WRITE (6,1520) I,TEMP,II,GBEST(I),G(I),ALB(NI),NONE
          GO TO 890
880       WRITE (6,1510) I,TEMP,II,GBEST(I),G(I),NONE,UB(NI)
890   CONTINUE
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( PENLTY ) WRITE(6,1600) USEOBJ
1600  FORMAT(18H TRUE OBJECTIVE IS  ,E12.5 )
C***********************************************************************
C
      WRITE (6,1530)
      DO 930 I=1,N
          TEMP=VAR(I)
          II=LITST(5)
          IF (ALB(I).EQ.X(I)) II=LITST(2)
          IF (UB(I).EQ.X(I)) II=LITST(1)
          IF (ALB(I).EQ.UB(I)) II=LITST(7)
          IF (ALB(I).EQ.-PLINFY.AND.UB(I).EQ.PLINFY) GO TO 900
          IF (ALB(I).EQ.-PLINFY) GO TO 910
          IF (UB(I).EQ.PLINFY) GO TO 920
          WRITE (6,1540) I,TEMP,II,X(I),ALB(I),UB(I)
          GO TO 930
900       II=LITST(6)
          WRITE (6,1550) I,TEMP,II,X(I),NONE,NONE
          GO TO 930
910       WRITE (6,1560) I,TEMP,II,X(I),NONE,UB(I)
          GO TO 930
920       WRITE (6,1570) I,TEMP,II,X(I),ALB(I),NONE
930   CONTINUE
      WRITE (6,1580)
935   IF (MAXIM) G(NOBJ)=-G(NOBJ)
      RETURN
1480  FORMAT (1H1,30X,25HOUTPUT OF INITIAL VALUES   ///1X,19A4//
     1 26H0SECTION 1 -- FUNCTIONS             )
1490  FORMAT (5X,10H FUNCTION ,19X,7HINITIAL,8X,5HLOWER,9X,5HUPPER/ 4H N
     1O.,3X,4HNAME,5X,6HSTATUS,2X,4HTYPE,7X,5HVALUE,9X,5HLIMIT,9X, 5HLIM
     2IT)
1500  FORMAT (I4,2X,A8,3X,A3,4X,A4,1P3E15.7)
1510  FORMAT (I4,2X,A8,3X,A3,4X,A4,1PE15.7,8X,A4,1PE18.7)
1520  FORMAT (I4,2X,A8,3X,A3,4X,A4,1P2E15.7,6X,A4)
1530  FORMAT (////24H0SECTION 2 -- VARIABLES /5X,8HVARIABLE,13X,7HINITIA
     1L,9X,5HLOWER,11X,5HUPPER/4H NO.,3X,4HNAME,5X,6HSTATUS,5X,5HVALUE,
     2 10X,5HLIMIT,11X,5HLIMIT)
1540  FORMAT (I4,2X,A8,3X,A3,2X,1P3E15.7)
1550  FORMAT (I4,2X,A8,3X,A4,1PE16.7,6X,A4,11X,A4)
1560  FORMAT (I4,2X,A8,3X,A3,1PE17.7,6X,A4,1PE18.7)
1570  FORMAT (I4,2X,A8,3X,A3,2X,1P2E15.7,6X,A4)
1580  FORMAT (1H1)
C
C     END OF TABLIN
C
      END

C
C
C     GRG2-2
C
C
      SUBROUTINE NEWTON (BINV,X,INBV,G,XSTAT,D,IBC,IBV,ROWB,COLB,INBC,
     1 ERROR,DBND)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES,MAXIM
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
C
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
      COMMON /USOBJ/ USEOBJ,USBEST,PENBST
C***********************************************************************
C
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON/TOLS/EPS,PLINFY,PLZERO,TOLX,TOLZ
      DIMENSION BINV(NBMAX,NNBMAX), X(NPMP1), INBV(N), G(MP1), XSTAT(N),
     1 D(N), IBC(NBMAX), IBV(M), ROWB(NBMAX), COLB(NBMAX)
      DIMENSION INBC(M),ERROR(NBMAX)
      DIMENSION DBND(NPNBMX)
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON /COUNTS/ NFTN,NGRAD,NMINV,NNFAIL,NCALLS,NIT,NBS,NSTEPC,NDUB
      COMMON /SRCHLG/ UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES,UNCONP
      COMMON /REDNEW/ CORNER,XB,XSTEP
      COMMON/MISC/MAXR,NSEAR,JP,LV,JQQ
      COMMON/NEWSRC/ITER
      COMMON/MAXBK/MAXIM
C     30     NC IS NONCONVERGENCE FLAG = 0 CONVERGED =1 DIDNT
      IF (IPR.GE.5) WRITE (6,220)
      NC=0
      FAIL=.FALSE.
      OLDNRM=PLINFY
      NCALLS=NCALLS+1
      IF (LV.NE.0) JR=IBV(LV)
      IF (IPR3.GE.4) WRITE (6,150) (X(I),I=1,NPMP1)
      DO 120 III = 1,ITLIM
          ITER=III-1
          IF (LV.EQ.0) GO TO 20
C
C     MODIFY SUPERBASICS IF IN BACKUP PHASE
C
          DO 10 I=1,NSUPER
              J=INBV(I)
10            X(J)=XSTAT(I)+XSTEP*D(I)
C
C     EVALUATE CONSTRAINTS
C
20        CALL GCOMP (G,X)
      NFTN = NFTN + 1
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
C
      IF ( .NOT. PENLTY ) GO TO 25
      IF ( VIOL ) GO TO 140
      CONST = 1.0D0
      IF (MAXIM) CONST = -CONST
      USEOBJ = G(NOBJ)
      G(NOBJ) = G(NOBJ) + RPEN*CONST*PENTRM
25    CONTINUE
C***********************************************************************
C
      IF (MAXIM) G(NOBJ)=-G(NOBJ)
C
C     COMPUTE NORM OF EQUATION ERROR
C
          IF (LV.EQ.0) GO TO 30
          IF (LV.GT.NB) X(JR)=G(JR-N)
          XNUM=XB-X(JR)
30        CONTINUE
          ERNORM=0.0D0
          IF (LV.NE.0) ERNORM=DABS(XNUM)
          IF (NB.EQ.0) GO TO 50
          DO 40 I=1,NB
              J=IBC(I)
              ERR=G(J)-X(N+J)
              IF (DABS(ERR).GT.ERNORM) ERNORM=DABS(ERR)
              ERROR(I)=ERR
40        CONTINUE
          IF (IPR3.GE.4) WRITE (6,180) (ERROR(I),I=1,NB)
C
C     TEST FOR CONVERGENCE
C
50        CONTINUE
          IF (IPR3.GE.4) WRITE (6,170) XNUM,ERNORM
          IF (ERNORM.LT.EPNEWT) GO TO 140
C
C     TESTS FOR NONCONVERGENCE
C
          IF (ERNORM.GT.OLDNRM.OR.ITER.EQ.ITLIM) GO TO 130
          OLDNRM=ERNORM
          NIT=NIT+1
C
C     PROCEED WITH NEWTON
C
          IF (LV.EQ.0) GO TO 90
C
C     COMPUTATIONS FOR BACKUP PHASE
C
          DENOM=-CORNER
          IF (NB.EQ.0) GO TO 70
          DO 60 I=1,NB
              RR=ROWB(I)
              XNUM=XNUM+RR*ERROR(I)
              DENOM=DENOM+RR*COLB(I)
60        CONTINUE
70        CONTINUE
          IF (IPR3.GE.4) WRITE (6,190) XNUM,DENOM
          IF (DABS(DENOM).LT.TOLZ) GO TO 130
          DELTA=XNUM/DENOM
          XSTEP=XSTEP-DELTA
          IF (IPR3.GE.4) WRITE (6,200) DELTA,XSTEP
          IF (NB.EQ.0) GO TO 120
          DO 80 I=1,NB
80            ERROR(I)=ERROR(I)-DELTA*COLB(I)
          IF (IPR3.GE.4) WRITE (6,180) (ERROR(I),I=1,NB)
C
C     COMPUTE NEWTON CORRECTION
C
90        DO 110 I=1,NB
              DEL=0.0D0
              DO 100 J=1,NB
100               DEL=DEL+BINV(I,J)*ERROR(J)
              K=IBV(I)
110           X(K)=X(K)-DEL
          IF (IPR3.GE.4) WRITE (6,150) (X(I),I=1,NPMP1)
120   CONTINUE
C
C     FAILURE
C
130   FAIL=.TRUE.
      NNFAIL=NNFAIL+1
140   CONTINUE
      IF (IPR.GE.5) WRITE (6,230)
      RETURN
C
150   FORMAT (5X,3HX =,/(1X,10E13.6))
160   FORMAT (21H NEWTON ITERATIONS = ,I5,2X,7HGNORM =,E12.5)
170   FORMAT (7H XNUM =,E15.8,9H ERNORM =,E13.6)
180   FORMAT (15H ERROR ARRAY =   /(1X,8E15.8))
190   FORMAT (7H XNUM =,E15.8,8H DENOM =,E15.8)
200   FORMAT (8H DELTA =,E15.8,8H XSTEP =,E15.8)
210   FORMAT (14H NEWTON FAILED   )
220   FORMAT (16H NEWTON ENTERED   )
230   FORMAT (18H NEWTON COMPLETED   )
C
C     END OF NEWTON
C
      END
      SUBROUTINE QUAD(IBV,X,V,XB1,XB2,XB3)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      DIMENSION IBV(M),X(NPMP1),V(NBMAX),XB1(NBMAX),XB2(NBMAX),
     1 XB3(NBMAX)
      COMMON/QUADBK/A1,A2,A3,ICON,IQUAD
      COMMON/BESTBK/STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND
      IF (ICON.EQ.2) GO TO 20
      AA=(STEP/A2)**2
      DO 10 I = 1,NB
          II=IBV(I)
10        X(II)=XB1(I)+STEP*V(I)+AA*(XB2(I)-XB1(I)-A2*V(I))
      GO TO 40
20    CONTINUE
      T2=A2-A1
      T3=A3-A1
      T32=A3-A2
      TA=STEP-A1
      AA=TA*TA
      W1=1.0D0-(TA*(T2+T3)-AA)/(T2*T3)
      W2=(TA*T3-AA)/(T2*T32)
      W3=(AA-TA*T2)/(T3*T32)
      DO 30 I = 1,NB
          II=IBV(I)
30        X(II)=W1*XB1(I)+W2*XB2(I)+W3*XB3(I)
40    RETURN
C
C     END OF QUAD
C
      END
      SUBROUTINE TANG(BINV,GRAD,IBC,INBV,V,D,RHS)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 GRAD
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      DIMENSION IBC(NBMAX),INBV(N),V(NBMAX),D(N)
      DIMENSION RHS(NBMAX)
      DIMENSION GRAD(MP1,N),BINV(NBMAX,NNBMAX)
C         **************************************************************
C     COMPUTES TANGENT VECTOR V = - BINV * JACOBIAN * DIRECTION
C         BINV IS THE BASIS INVERSE
C         JACOBIAN HAS I,J ELEMENT = PARTIAL I'TH BINDING CONSTRAINT
C              WITH RESPECT TO J'TH NONBASIC VARIABLE
C         **************************************************************
10    FORMAT (19H TANGENT VECTOR IS /(1X,10E13.6))
      DO 40 I=1,NB
          TRHS=0.0D0
          II=IBC(I)
          DO 30 J = 1,NSUPER
              JJ=INBV(J)
              IF (JJ.GT.N) GO TO 20
              TRHS=TRHS+GRAD(II,JJ)*D(J)
              GO TO 30
20            IF (JJ-N.EQ.II) TRHS=TRHS-D(J)
30        CONTINUE
          RHS(I)=TRHS
40    CONTINUE
      DO 60 I=1,NB
          TMP=0.0D0
          DO 50 J=1,NB
50            TMP=TMP+BINV(I,J)*RHS(J)
          V(I)=-TMP
60    CONTINUE
      RETURN
C
C     END OF TANG
C
      END
      SUBROUTINE CG(GRADFP,GRADF,D,MSGCG,INBV,IUB)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL MOVE,RESTRT,DROP,VARMET,CONJGR
      LOGICAL UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      DIMENSION GRADFP(N), GRADF(N), D(N)
      DIMENSION INBV(N), IUB(N)
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON /EPSCOM/ EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      COMMON /TOLS/ EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON /LOGBLK/ MOVE,RESTRT,DROP,VARMET,CONJGR
      COMMON /CGBK/ MODCG
      COMMON/BESTBK/STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON/MISC/MAXR,NSEAR,JP,LV,JQQ
      COMMON /SRCHLG/ UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES
      DATA INITCG/1/
C
C     CONJUGATE GRADIENT METHOD ON SUPERBASICS
C
      IF (IPR.GE.5) WRITE (6,330)
      IF (INITCG.NE.0.OR.NSEAR.EQ.0) MSGCG=1
      INITCG=0
      IF (MSGCG.GT.0) WRITE (6,220)
      IF ((.NOT.RESTRT).AND.(ITNCG.LE.NSUPER)) GO TO 20
C
C     RESTART
C
9     CONTINUE
      ITNCG=0
      CGBETA=0.0D0
      DO 10 I = 1,NSUPER
10        D(I)=-GRADF(I)
      IF (MSGCG.EQ.0) GO TO 210
      GO TO (11,12,13,14,15), MODCG
11    WRITE (6,230)
      GO TO 210
12    WRITE (6,240)
      GO TO 210
13    WRITE (6,250)
      GO TO 210
14    WRITE (6,260)
      GO TO 210
15    WRITE (6,270)
      GO TO 210
20    GO TO (30,60,90,110,120), MODCG
C
C     FLETCHER-REEVES
C
30    IF (MSGCG.GT.0) WRITE (6,230)
      G1N=0.0D0
      DO 40 I=1,NSUPER
40        G1N=G1N+GRADFP(I)**2
      IF (G1N.LE.EPS) GO TO 9
      G2N=0.0D0
      DO 50 I=1,NSUPER
50        G2N=G2N+GRADF(I)**2
      CGBETA=G2N/G1N
      IF (IPR.LT.5) GO TO 190
      WRITE (6,290) CGBETA
      WRITE (6,300) (GRADF(I),I=1,NSUPER)
      WRITE(6,310) (GRADFP(I),I=1,NSUPER)
      GO TO 190
C
C     POLAK-RIBIERE
C
60    IF (MSGCG.GT.0) WRITE (6,240)
      G1N=0.0D0
      DO 70 I=1,NSUPER
70        G1N=G1N+GRADFP(I)**2
      IF (G1N.LE.EPS) GO TO 9
      GTY=0.0D0
      DO 80 J=1,NSUPER
          G1=GRADFP(J)
          G2=GRADF(J)
          GTY=GTY+G2*(G2-G1)
80    CONTINUE
      CGBETA=GTY/G1N
      GO TO 190
C
C     PERRY
C
90    IF (MSGCG.GT.0) WRITE (6,250)
      GYD=0.0D0
      YTD=0.0D0
      DO 100 J=1,NSUPER
          YJ=GRADF(J)-GRADFP(J)
          GYD=GYD+GRADF(J)*(YJ-STEP*D(J))
          YTD=YTD+YJ*D(J)
100   CONTINUE
      IF (DABS(YTD).LE.EPS) GO TO 9
      CGBETA=GYD/YTD
      GO TO 190
C
C     DFP
C
110   IF (MSGCG.GT.0) WRITE (6,260)
      GO TO 130
C
C     COMPLEMENTARY DFP
C
120   IF (MSGCG.GT.0) WRITE (6,270)
130   GTY=0.0D0
      YTD=0.0D0
      YTY=0.0D0
      DO 140 J=1,NSUPER
          G1=GRADFP(J)
          G2=GRADF(J)
          YJ=G2-G1
          GTY=GTY+G2*YJ
          YTD=YTD+YJ*D(J)
          YTY=YTY+YJ**2
140   CONTINUE
      IF (DABS(YTD).LE.EPS) GO TO 9
      IF (DABS(YTY).LE.EPS) GO TO 9
      GTD=0.0D0
      DO 150 I=1,NSUPER
150       GTD=GTD+GRADF(I)*D(I)
      IF (MODCG.EQ.5) GO TO 160
C
      CGBETA=-STEP*GTD/YTD
      GAMMA=GTY/YTY
      GO TO 170
C
160   RHO=STEP+YTY/YTD
      CGBETA=(GTY-RHO*GTD)/YTD
      GAMMA=GTD/YTD
C
170   ITNCG=ITNCG+1
      IF (IPR.GE.5) WRITE (6,280) CGBETA,GAMMA,YTD,YTY,GTD,GTY
      DO 180 J=1,NSUPER
180       D(J)=-GRADF(J)+CGBETA*D(J)+GAMMA*(GRADF(J)-GRADFP(J))
      GO TO 210
C
C     SET UP NEXT CG-TYPE SEARCH DIRECTION
C
190   ITNCG=ITNCG+1
      DO 200 J=1,NSUPER
200       D(J)=-GRADF(J)+CGBETA*D(J)
210   MSGCG=0
      IF (IPR.GE.5) WRITE (6,320)(D(I),I=1,NSUPER)
      IF (IPR.GE.5) WRITE (6,340)
      RETURN
C
C
C
220   FORMAT (73H HESSIAN IS TOO LARGE FOR VARIABLE METRIC--SWITCH TO CO
     1NJUGATE GRADIENTS   )
230   FORMAT (40H FLETCHER-REEVES DIRECTION WILL BE USED   )
240   FORMAT(38H POLAK-RIBIERE DIRECTION WILL BE USED   )
250   FORMAT(41H PERRY (MARCH 76) DIRECTION WILL BE USED   )
260   FORMAT(28H DFP DIRECTION WILL BE USED   )
270   FORMAT (42H COMPLEMENTARY DFP DIRECTION WILL BE USED   )
280   FORMAT (9H CGBETA =,E14.7,8H GAMMA =,E14.7,6H YTD =,E14.7,
     1 6H YTY =,E14.7,6H GTD =,E14.7,6H GTY =,E14.7)
290   FORMAT (9H CGBETA =,E15.8)
300   FORMAT(8H GRADF =,8E15.8/(1X,8E15.8))
310   FORMAT(8H GRADFP=,8E15.8/(1X,8E15.8))
320   FORMAT (4H D =,8E15.8/(1X,8E15.8))
330   FORMAT (12H CG ENTERED   )
340   FORMAT (14H CG COMPLETED   )
C
C  END OF CG
C
      END
      SUBROUTINE COMDFP
     1 (R,Y,GG,NY,THETA,YTP,GTP,KADD,KSUB)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/TOLS/EPS,PLINFY,PLZERO,TOLX,TOLZ
      DIMENSION R(NRTOT)
      DIMENSION Y(N),GG(N)
C
C  MODIFY R USING COMPLEMENTARY DFP FORMULA
C
      IF (DABS(THETA).LE.EPS) RETURN
      IF (DABS(YTP).LT.EPS) RETURN
      IF (DABS(GTP).LT.EPS) GTP=EPS
      C1=1.0D0/(THETA*YTP)
      C2=1.0D0/GTP
      CALL R1MOD(R,Y,NY,C1,KADD)
      CALL R1MOD(R,GG,NY,C2,KSUB)
      RETURN
C
C     END OF COMDFP
C
      END
      SUBROUTINE RESETR(R,COND)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      LOGICAL MOVE,RESTRT,DROP,VARMET,CONJGR
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON/MISC/MAXR,NSEAR,JP,LV,JQQ
      COMMON /LOGBLK/ MOVE,RESTRT,DROP,VARMET,CONJGR
      DIMENSION R(NRTOT)
C
C  RESET THE CHOLESKY FACTOR OF THE HESSIAN
C
      IF (NSUPER.EQ.0) RETURN
      IF (IPR.GE.4) WRITE (6,100)
      COND=1.0D0
      NCOL=NSUPER
      IF (MAXR.LT.NCOL) NCOL=MAXR
      K=NCOL*(NCOL+1)/2
      DO 10 I = 1,K
10        R(I)=0.0D0
      K=0
      DO 20 I = 1,NCOL
          K=K+I
20        R(K)=1.0D0
C
      RESTRT=.TRUE.
      RETURN
100   FORMAT (40H CHOLESKY FACTOR OF HESSIAN RESET TO I.   )
C
C  END OF RESETR
C
      END
      SUBROUTINE RTRSOL
     1 (R,Y,NY)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      DIMENSION R(NRTOT)
      DIMENSION Y(N)
C
C  SOLVE R'R*Y = Y (R TRIANGULAR, STORED BY COLUMNS)
C
      Y(1) = Y(1)/R(1)
      K=1
      IF (NY.LE.1) GO TO 50
      DO 20 I = 2,NY
          S=Y(I)
          I1=I-1
          DO 10 J = 1,I1
              K=K+1
              S=S-R(K)*Y(J)
10        CONTINUE
          K=K+1
20    Y(I)=S/R(K)
C
50    DO 80 II = 1,NY
          I=NY+1-II
          T=Y(I)/R(K)
          Y(I)=T
          IF (I.LE.1) GO TO 80
          K=K-I
          I1=I-1
          DO 60 J = 1,I1
              Y(J)=Y(J) - R(K+J)*T
60        CONTINUE
80    CONTINUE
      RETURN
C
C  END OF RTRSOL
C
      END
      SUBROUTINE R1ADD
     1 (R,Y,NY)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/EPSCOM/EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      DIMENSION R(NRTOT)
      DIMENSION Y(N)
C
C  MODIFY R SUCH THAT R'R := R_R + YY' WHERE R IS UPPER-TRIANGULAR,
C  STORED BY COLUMNS.
C
      K1=1
      DO 100 K = 1,NY
      K0=K1
          T1=R(K1)
          T2=Y(K)
          D=DSQRT(T1*T1+T2*T2)
          CS=T1/D
          SN=T2/D
          J1=K1+K
          K1=J1+1
          IF (DABS(SN).LE.EPS0) GO TO 100
          R(K0)=D
          L=K+1
          IF (L.GT.NY) GO TO 100
          DO 50 J = L,NY
		      T1=R(J1)
              T2=Y(J)
              R(J1)=CS*T1 + SN*T2
              Y(J)= SN*T1 - CS*T2
          	  J1=J1+J
50        CONTINUE
100   CONTINUE
C
900   RETURN
C
C     END OF R1ADD
C	 
      END
      SUBROUTINE R1MOD
     1 (R,Y,NY,SIGMA,K)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL POSDEF
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/EPSCOM/EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      DIMENSION R(NRTOT)
      DIMENSION Y(N)
C
C  MODIFY R SUCH THAT R'R := R'R + SIGMA*YY'
C
      TOL=EPS0
      S=0.0D0
      T=DSQRT(DABS(SIGMA))
      DO 100 I = 1,NY
          S=Y(I)**2 + S
          Y(I)=Y(I)*T
100   CONTINUE
      S=SIGMA*S
      IF (DABS(S).LE.TOL) RETURN
      IF (S.LE.TOL) GO TO 200
      CALL R1ADD(R,Y,NY)
      RETURN
C
200   CALL R1SUB(R,Y,NY,TOL,POSDEF)
      IF (.NOT.POSDEF) K=-K
      RETURN
C
C  END OF R1MOD
C
      END
      SUBROUTINE R1SUB
     1 (R,Y,NY,TOL,POSDEF)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL POSDEF
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/EPSCOM/EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      DIMENSION R(NRTOT)
      DIMENSION Y(N)
C
C  MODIFY R SUCH THAT R'R := R'R - YY' WHERE R IS UPPER-TRIANGULAR,
C  STORED BY COLUMNS.
C   SEE SAUNDERS, STANFORD TECHNICAL REPORT STAN-CS-72-252, CHAPTER 7.
C
C  FIRST SOLVE R'P = Y
C
      T=Y(1)/R(1)
      D=T*T
      Y(1)=T
      K=1
      IF (NY.LE.1) GO TO 50
      DO 20 I = 2,NY
          S=Y(I)
          I1=I-1
          DO 10 J = 1,I1
              K=K+1
              S=S-R(K)*Y(J)
10        CONTINUE
          K=K+1
          T=S/R(K)
          D=T*T+D
          Y(I)=T
20    CONTINUE
C
C  SEE IF NEW R WILL BE POSITIVE DEFINITE.
C
50    D=1.0D0-D
      POSDEF=(D.GT.EPS0)
      IF (.NOT.POSDEF) RETURN
      S=DSQRT(D)
C
C  PERFORM BACKWARD SWEEP OF PLANE ROTATIONS.
C
      DO 80 II=1,NY
          I=NY+1-II
          U=S
          T=Y(I)
          D=T*T + D
          S=DSQRT(D)
          CS=U/S
          SN=T/S
          Y(I)=0.0D0
          L=K
          K=K-I
          IF (DABS(SN).LE.EPS0) GO TO 80
          DO 60 J = I,NY
              T1=Y(J)
              T2=R(L)
              Y(J)=CS*T1 + SN*T2
              R(L)=SN*T1 - CS*T2
              L=L+J
60        CONTINUE
80    CONTINUE
      RETURN
C
C     END OF R1SUB
C
      END
      SUBROUTINE PARSH(G,X,MP1,N,GRAD)
      WRITE (6,10)
10    FORMAT(60H0AN ATTEMPT IS BEING MADE TO USE THE DUMMY PARSH SUBROUT
     1INE   )
      STOP
      END
      SUBROUTINE PARSHC(X,G,GRAD,IFIK,GG,GBEST)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 GRAD
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      DIMENSION X(NPMP1),GBEST(MP1),GRAD(MP1,N),IFIK(N)
      DIMENSION GG(MP1),G(MP1)
      COMMON/BESTBK/STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      LOGICAL FAIL1
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
      COMMON /MNGRG/ EPSTOP,LIMSER,NSTOP,IERR,IPN4,IPN5,IPN6
      COMMON /MISC/ MAXR,NSEAR,JP,LV,JRR
      COMMON /USOBJ/ USEOBJ,USBEST,PENBST
      FAIL1 = .FALSE.
C***********************************************************************
C
C
C     THIS SUBROUTINE COMPUTES FINITE DIFFERENCE DERIVATIVES BY
C     CENTRAL DIFFERENCING
C
      PSTEP = 1.0D-4
      DXMIN = 0.1*PSTEP
      DO 200 I=1,N
          IF (IFIK(I).NE.0) GO TO 160
          DX=DABS(X(I))*PSTEP
          IF (DX.LT.DXMIN) DX=DXMIN
          TS=X(I)
          X(I) = X(I)+ DX
          CALL GCOMP(GG,X)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( .NOT. PENLTY ) GO TO 120
      IF ( .NOT. VIOL ) GO TO 120
      FAIL1 = .TRUE.
      DO 110 L = 1,MP1
         GG(L) = G(L)
110   CONTINUE
      GG(NOBJ) = USEOBJ
120   CONTINUE
C***********************************************************************
C
          X(I)=TS-DX
          CALL GCOMP(GBEST,X)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( .NOT. PENLTY ) GO TO 145
      IF ( .NOT. VIOL ) GO TO 140
      IF ( FAIL1 ) GO TO 210
      DO 130 L = 1,MP1
      GBEST(L) = G(L)
130   CONTINUE
      GBEST(NOBJ) = USEOBJ
      WRITE(6,320) DX,I
      DX = DX/2.0D0
140   CONTINUE
      IF ( FAIL1 ) WRITE (6,310) DX,I
      IF ( FAIL1 ) DX = DX/2.0D0
145   CONTINUE
C***********************************************************************
C
          DO 150 L=1,MP1
150           GRAD(L,I) = (GG(L)-GBEST(L))/(2.0D0*DX)
          X(I) = TS
          GO TO 200
160       DO 170 L = 1,MP1
170           GRAD(L,I) = 0.0D0
 200  CONTINUE
      RETURN
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
210   CONTINUE
C
C   SET LIMSER TO FORCE CODE TO STOP AND PRINT OUTPUT
C
      LIMSER = NSEAR
      WRITE (6,330) I,TS,DX
      WRITE ( 6,400)
      X(I) = TS
      RETURN
C
310   FORMAT(2X,45HDURING CENTRAL DIFFERENCING A HARD CONSTRAINT ,
     1 34H WAS VIOLATED USING AN INCREASE OF ,E14.5/2X,
     2 11HIN VARIABLE ,I5,28HFORWARD DIFFERENCING USED ON  ,
     3 14H THIS VARIABLE )
320   FORMAT(2X,45HDURING CENTRAL DIFFERENCING A HARD CONSTRAINT ,
     1 34H WAS VIOLATED USING A  DECREASE OF ,E14.5/2X,
     2 11HIN VARIABLE ,I5,28HFORWARD DIFFERENCING USED ON  ,
     3 14H THIS VARIABLE )
330   FORMAT(2X,45HDURING CENTRAL DIFFERENCING A HARD CONSTRAINT ,
     1 56HWAS VIOLATED BY PERTURBING A VARIABLE IN BOTH DIRECTIONS ,
     2 13H FOR VARIABLE ,I5 / 2X , 10HWITH X(I)= ,E14.5,8H AND DX= ,
     3 E14.5 )
400   FORMAT(2X,45HIN SUBROUTINE PARSHC, RESETTING SEARCH LIMIT,,1X,
     1 41HTO FORCE PROGRAM TO STOP AND PRINT OUTPUT )
C***********************************************************************
C
C
C     END OF PARSHC
C
      END
      SUBROUTINE PARSHF(X,G,GRAD,IFIK,GG,UB)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 GRAD
      LOGICAL MAXIM
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
      COMMON /USOBJ/ USEOBJ,USBEST,PENBST
      COMMON /MNGRG/ EPSTOP,LIMSER,NSTOP,IERR,IPN4,IPN5,IPN6
      COMMON /MISC/ MAXR,NSEAR,JP,LV,JRR
C***********************************************************************
C
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      DIMENSION X(NPMP1),G(MP1),GRAD(MP1,N),IFIK(N)
      DIMENSION UB(NPMP1)
      DIMENSION GG(MP1)
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON/BESTBK/STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON/MAXBK/MAXIM
C
C     THIS SUBROUTINE COMPUTES FINITE DIFFERENCE DERIVATIVES BY
C     FORWARD DIFFERENCING
C
      TMP=G(NOBJ)
      IF (MAXIM) TMP=-TMP
      IF (NINF.NE.0) TMP=TRUOBJ
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( PENLTY ) TMP = USEOBJ
C***********************************************************************
C
      PSTEP = 1.0D-4
      DXMIN = 0.1*PSTEP
      DO 200 I=1,N
          IF (IFIK(I).NE.0) GO TO 160
          DX=DABS(X(I))*PSTEP
          IF (DX.LT.DXMIN) DX=DXMIN
          TS=X(I)
          IF (UB(I)-TS.LT.DX) DX=-DX
          X(I) = X(I)+ DX
          CALL GCOMP(GG,X)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( .NOT. PENLTY ) GO TO 140
      IF ( .NOT. VIOL ) GO TO 140
      X(I) = TS-DX
      CALL GCOMP(GG,X)
      IF ( .NOT. VIOL ) GO TO 140
      WRITE (6,300) I,TS,DX
      WRITE(6,400)
C   SET LIMSER TO FORCE CODE TO STOP AND PRINT OUTPUT
      LIMSER = NSEAR
      X(I) = TS
      RETURN
140   CONTINUE
C***********************************************************************
C
          DO 150 L=1,MP1
 150          GRAD(L,I) = (GG(L)-G(L))/DX
          GRAD(NOBJ,I) = (GG(NOBJ)-TMP)/DX
          X(I) = TS
          GO TO 200
160       DO 170 L = 1,MP1
170           GRAD(L,I) = 0.0D0
 200  CONTINUE
      RETURN
C
C***********************************************************************
300   FORMAT(2X,45HDURING FORWARD DIFFERENCING A HARD CONSTRAINT ,
     1 56HWAS VIOLATED BY PERTURBING A VARIABLE IN BOTH DIRECTIONS ,
     2 13H FOR VARIABLE ,I5 / 2X , 10HWITH X(I)= ,E14.5,8H AND DX= ,
     3 E14.5 )
400   FORMAT(2X,45HIN SUBROUTINE PARSHF, RESETTING SEARCH LIMIT,,1X,
     1 41HTO FORCE PROGRAM TO STOP AND PRINT OUTPUT )
C***********************************************************************
C
C
C  END OF PARSHF
C
      END
      SUBROUTINE REDGRA
     1 (REDGR,BINV,GRAD,INBC,IBV,U,INBV,IBC,G,ALB,UB)
C         **************************************************************
C     COMPUTES THE REDUCED GRADIENT  = REDGRA(.)
C         **************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 GRAD
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
C***********************************************************************
C
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON/PH1BK/PHMULT,PH1EPS,INITPH
      DIMENSION REDGR(N),INBC(M),IBV(M),U(NBMAX),INBV(N),IBC
     1 (NBMAX),G(MP1),ALB(NPMP1),UB(NPMP1)
      DIMENSION GRAD(MP1,N),BINV(NBMAX,NNBMAX)
10    FORMAT (38H0REDGRA CAN'T FIND IBC INDEX TO MATCH ,I5)
      IF (NINF.EQ.0) GO TO 70
C     IF PHASE ONE THEN
C     THIS BLOCK COMPUTES GRADIENT OF SUM OF INFEASIBILITIES, STORES IT
C     IN ROW M+1 OF GRAD.
      DO 60 J=1,N
      TMP=0.0D0
      DO 40 I = 1,NNBC
      K=INBC(I)
      IF (G(K).LT.ALB(N+K)-EPNEWT) GO TO 20
      IF (G(K).GT.UB(N+K)+EPNEWT) GO TO 30
      GO TO 40
20    TMP=TMP-GRAD(K,J)
      GO TO 40
30    TMP=TMP+GRAD(K,J)
40    CONTINUE
      GRAD(NOBJ,J) = PHMULT*GRAD(NOBJ,J) + TMP
60    CONTINUE
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( .NOT. PENLTY ) GO TO 75
C
C   ADD IN GRADIENT OF PENALTY TERM IN PHASE 1
      DO 72 J = 1,N
          GRAD(NOBJ,J) = GRAD(NOBJ,J) + RPEN*GRDPEN(J)
72    CONTINUE
75    CONTINUE
C***********************************************************************
C
70    CONTINUE
      IF (NB.EQ.0) GO TO 100
C     COMPUTES LAGRANGE MULTIPLIERS U
      DO 90 I=1,NB
          UU=0.0D0
          DO 80 J=1,NB
              JJ=IBV(J)
              IF (JJ.GT.N) GO TO 80
              UU=UU+GRAD(NOBJ,JJ)*BINV(J,I)
80        CONTINUE
          U(I)=UU
90    CONTINUE
100   CONTINUE
C     NOW COMPUTE REDGRA USING MULTIPLIERS.
      DO 160 I=1,N
          II=INBV(I)
          IF (II.GT.N) GO TO 130
          TEMP=0.0D0
          IF (NB.EQ.0) GO TO 120
          DO 110 J=1,NB
              JJ=IBC(J)
              TEMP=TEMP+U(J)*GRAD(JJ,II)
110       CONTINUE
120       CONTINUE
          REDGR(I)=GRAD(NOBJ,II) - TEMP
          GO TO 160
C     HANDLES SLACK VARIABLES.
130       I2=II-N
          DO 140 J=1,NB
              ISV=J
              IF (IBC(J).EQ.I2) GO TO 150
140       CONTINUE
          GO TO 170
150       REDGR(I)=U(ISV)
160   CONTINUE
      RETURN
170   WRITE (6,10) I2
      STOP
C
C     END OF REDGRA
C
      END
      SUBROUTINE REDOBJ (BINV,GRAD,X,G,IBV,V,XB1,XB2,XB3,INBV,ALB,UB,XST
     1AT,D,ISTAT,INBC,GBEST,XBEST,IBC,ROWB,COLB,ROW,DBND)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES,MAXIM
      REAL *4GRAD
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
      COMMON /USOBJ/ USEOBJ,USBEST,PENBST
C***********************************************************************
C
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      DIMENSION X(NPMP1), G(MP1), IBV(M), V(NBMAX), XB1(NBMAX),
     1 XB2(NBMAX), XB3(NBMAX), INBV(N), ALB(NPMP1), UB(NPMP1), XSTAT(N),
     2 D(N),ISTAT(MP1),INBC(M),GBEST(MP1),XBEST(MPNBMX),IBC(NBMAX),
     3 ROWB(NBMAX),COLB(NBMAX)
      DIMENSION ROW(NBMAX)
      DIMENSION DBND(NPNBMX)
      DIMENSION GRAD(MP1,N),BINV(NBMAX,NNBMAX)
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON /COUNTS/ NFTN,NGRAD,NMINV,NNFAIL,NCALLS,NIT,NBS,NSTEPC,NDUB
      COMMON /SRCHLG/ UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES,UNCONP
      COMMON /REDNEW/ CORNER,XB,XSTEP
      COMMON/MISC/MAXR,NSEAR,JP,LV,JQQ
      COMMON /QUADBK/ A1,A2,A3,ICON,IQUAD
      COMMON/TOLS/EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON/BESTBK/STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON/REDPH/TRUBST
      COMMON/MAXBK/MAXIM
      COMMON/NEWSRC/ITER
      COMMON/REDSER/NINFB
      IF (IPR.GE.5) WRITE (6,540)
      FAIL=.FALSE.
      XB=0.0D0
      LV=0
      LV1=0
      MXSTEP=(STEP.GE.STEPMX)
      IF (MXSTEP) STEP=STEPMX
      DO 10 I=1,NSUPER
          J=INBV(I)
10        X(J)=XSTAT(I)+STEP*D(I)
      IF (NB.NE.0) GO TO 20
C     10     NO BINDING CONSTRAINTS
C     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I
      CALL GCOMP (G,X)
C     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( .NOT. PENLTY ) GO TO 15
      IF ( VIOL .AND. IPR3 .GE. 1 ) WRITE(6,600) STEP
      IF ( VIOL ) GO TO 440
      CONST = 1.0D0
      IF (MAXIM) CONST = -CONST
      USEOBJ = G(NOBJ)
      G(NOBJ) = G(NOBJ) + RPEN*CONST*PENTRM
15    CONTINUE
C***********************************************************************
C
      IF (MAXIM) G(NOBJ)=-G(NOBJ)
      NFTN=NFTN+1
      GO TO 60
20    CONTINUE
      IF (ICON.EQ.0) GO TO 30
C     OBTAIN INITIAL ESTIMATE OF BASIC VARIABLES BY QUADRATIC EXTRAPOLAT
C     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I
      CALL QUAD (IBV,X,V,XB1,XB2,XB3)
C     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I
      GO TO 50
30    CONTINUE
C
C     LINEAR EXTRAPOLATION
C
      DO 40 I=1,NB
          K=IBV(I)
          X(K)=XBEST(I)+(STEP-STPBST)*V(I)
40    CONTINUE
50    CONTINUE
C     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I
      CALL NEWTON(BINV,X,INBV,G,XSTAT,D,IBC,IBV,ROWB,COLB,INBC,ROW,
     1 DBND)
C     =I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I=I
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( .NOT. PENLTY ) GO TO 55
      IF ( VIOL .AND. IPR3 .GE. 1 ) WRITE(6,600) STEP
      IF ( VIOL ) GO TO 440
55    CONTINUE
C***********************************************************************
C
      IF (FAIL) GO TO 330
C
C     COMPUTE SLACKS FOR NONBINDING CONSTRAINTS AND MAKE THEM
C     TEMPORARILY BASIC
C
60    KNB=NB
      IF (NNBC.EQ.0) GO TO 80
      DO 70 I=1,NNBC
          J=INBC(I)
          X(N+J)=G(J)
          KNB=NB+I
70        IBV(KNB)=N+J
80    IF (KNB.EQ.0) GO TO 370
      IF (LV.EQ.0) GO TO 90
      IF (XSTEP.GE.STEP.OR.XSTEP.LE.STPBST) GO TO 320
      STEP=XSTEP
90    IF (NINF.EQ.0) GO TO 120
C
C     COME HERE IF INFEASIBLE
C
C
C                      -BEGIN CHECKING FOR CHANGES IN STATUS OF BASIC
C     VARIABLES
C
      LV=0
      THMIN=PLINFY
      DO 110 I=1,KNB
          J=IBV(I)
          XJ=X(J)
          XBESTJ=XBEST(I)
          DENOM=XJ-XBESTJ
          IF (DABS(DENOM).LE.EPS) GO TO 110
          TB=ALB(J)
          T=TB-EPNEWT
          IF ((XBESTJ.GE.T).AND.(XJ.LT.T)) GO TO 100
          TB=UB(J)
          T=TB+EPNEWT
          IF (XBESTJ.LE.T.AND.XJ.GT.T) GO TO 100
          GO TO 110
100       THET=(T -XBESTJ)/DENOM
          IF (THET.GE.THMIN) GO TO 110
          THMIN=THET
          LV=I
110   CONTINUE
C
C      TEST IF BACKUP PHASE NEEDED
       IF ( LV .NE. 0 ) GO TO 160
C      NO BACKUP PHASE NEEDED  CHECK IF FEASIBLE
      CALL PH1OBJ (INBC,G,UB,ALB)
      IF (NINF.EQ.0) GO TO 340
C    STILL INFEASIBLE
      IF (IPR.GE.2) WRITE (6,500) STEP,G(NOBJ),ITER
       GO TO 350
C
C     WE WERE FEASIBLE BEFORE NEWTON.  CHECK BOUNDS ON BASICS TO SEE IF
C     WE ARE STILL FEASIBLE.
C
120   LV=0
      THMIN=PLINFY
      DO 150 I=1,KNB
          J=IBV(I)
          XJ=X(J)
          XBESTJ=XBEST(I)
          DENOM=XJ-XBESTJ
          IF (DABS(DENOM).LE.EPS) GO TO 150
          IF (ALB(J)-XJ.LE.EPNEWT) GO TO 130
          T=ALB(J)-XBESTJ
          GO TO 140
130       IF (XJ-UB(J).LE.EPNEWT) GO TO 150
          T=UB(J)-XBESTJ
140       THET=T/DENOM
          IF (THET.GE.THMIN) GO TO 150
          THMIN=THET
          LV=I
150   CONTINUE
160   IF (IPR.GE.5) WRITE (6,490) LV,THMIN
      IF (THMIN.LT.0.0D0) GO TO 320
      IF (NINF.EQ.0.AND.MAXIM) G(NOBJ)=-G(NOBJ)
      IF (IPR3.GE.1) WRITE(6,500) STEP,G(NOBJ),ITER
      IF (NINF.EQ.0.AND.MAXIM) G(NOBJ)=-G(NOBJ)
      IF (LV.EQ.0) GO TO 350
      IF (IPR3.LT.1) GO TO 167
      IF (LV.GT.NB) GO TO 163
      WRITE (6,510) IBV(LV)
      GO TO 167
163   I=IBV(LV)-N
      WRITE (6,520) I
167   CONTINUE
C         **************************************************************
C     BOUND VIOLATED--START BACK UP PHASE
C         **************************************************************
C
C     SET XB = BOUND NEAREST X(LV)
C
      JR=IBV(LV)
      B1=ALB(JR)
      B2=UB(JR)
      D1=X(JR)-B1
      D2=X(JR)-B2
      XB=B1
      IF (DABS(D1).GT.DABS(D2)) XB=B2
      XSTEP=STPBST+THMIN*(STEP-STPBST)
      IF (IPR.GE.5) WRITE (6,470) G(NOBJ),LV,XSTEP
      IF (NB.EQ.0) GO TO 210
C
C     COMPUTE COLB=B2*D, COLUMN PART OF JACOBIAN BORDER
C
      DO 190 I=1,NB
          SUM=0.0D0
          K=IBC(I)
          DO 180 J=1,NSUPER
              KK=INBV(J)
              IF (KK.LE.N) GO TO 170
              TS=0.0D0
              IF (K.EQ.KK-N) TS=-1.0D0
              GO TO 180
170           TS=GRAD(K,KK)
180           SUM=SUM+TS*D(J)
190       COLB(I)=SUM
      IF (IPR.GE.5) WRITE (6,460) (COLB(I),I=1,NB)
C
C     DO ROW BORDER AND CORNER ELEMENT CALCULATIONS
C
      IF (LV.GT.NB) GO TO 210
C
C     CASE OF BASIC VARIABLE VIOLATING BOUND
C
      DO 200 I=1,NB
200       ROW(I)=0.0D0
      ROW(LV)=1.0D0
      CORNER=0.0D0
      GO TO 270
C
C     CASE OF CONSTRAINT VIOLATING BOUND
C
210   J=IBV(LV)
      L=J-N
      IF (NB.EQ.0) GO TO 240
      DO 230 I=1,NB
          K=IBV(I)
          IF (K.GT.N) GO TO 220
          ROW(I)=GRAD(L,K)
          GO TO 230
220       ROW(I)=0.0D0
230   CONTINUE
240   CORNER=0.0D0
      DO 260 I=1,NSUPER
          J=INBV(I)
          IF (J.LE.N) GO TO 250
          TS=0.0D0
          GO TO 260
250       TS=GRAD(L,J)
260       CORNER=CORNER+TS*D(I)
      IF (NB.EQ.0) GO TO 310
270   DO 290 I=1,NB
          TMP=0.0D0
          DO 280 J=1,NB
280           TMP=TMP+ROW(J)*BINV(J,I)
290       ROWB(I)=TMP
      IF (IPR.GE.5) WRITE (6,450) (ROWB(I),I=1,NB),CORNER
C
C     COMPUTE ESTIMATES OF BASIC VARIABLES
C
      DO 300 I=1,NB
          J=IBV(I)
          XBI=XBEST(I)
300       X(J)=XBI+(X(J)-XBI)*THMIN
310   LV1=LV
      GO TO 50
C
C     BACKUP FAILED
C
320   CONTINUE
      IF (IPR.GE.5) WRITE (6,480)
C
C     NEWTON CALL FAILED
C
330   CONTINUE
      FAIL=.TRUE.
      IF (IPR3.GE.1) WRITE (6,530) ITER
      GO TO 440
C
C     JUST BECAME FEASIBLE
C
340   CONTINUE
      JSTFES=.TRUE.
      G(NOBJ)=TRUOBJ
      GO TO 440
C
C     NORMAL TERMINATION
C
350   CONTINUE
      IF (LV.GT.0.OR.NNBC.EQ.0) GO TO 370
C
C     CHECK FOR NEW BINDING CONSTRAINTS
C
      DO 360 I=1,NNBC
          J=INBC(I)
          GJ=G(J)
          NJ=N+J
          IF (DABS(GJ-ALB(NJ)).LT.EPNEWT.OR.DABS(GJ-UB(NJ)).LT.EPNEWT) L
     *    V1=N
     1    B+I
360   CONTINUE
C
C     UPDATE BEST VALUES
C
370   CONTINUE
      IF (G(NOBJ).GT.OBJBST) GO TO 430
      IF (NB.EQ.0) GO TO 390
      DO 380 I=1,NB
          J=IBV(I)
380       XBEST(I)=X(J)
390   IF (NNBC.EQ.0) GO TO 410
      DO 400 I=1,NNBC
          K=INBC(I)
400       XBEST(NB+I)=G(K)
410   DO 420 I=1,MP1
420       GBEST(I)=G(I)
      STPBST=STEP
      OBJBST=G(NOBJ)
      TRUBST=TRUOBJ
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( PENLTY ) USBEST = USEOBJ
      IF ( PENLTY ) PENBST = PENTRM
C
600   FORMAT(25H HARD CONSTRAINT VIOLATED ,10X,6HSTEP = ,E15.6 )
C***********************************************************************
C
      NINFB=NINF
430   LV=LV1
440   IF (IPR.GE.5) WRITE (6,550)
      RETURN
C
450   FORMAT (33H ROW BORDER FOLLOWED BY CORNER =   /(1X,8E15.8))
460   FORMAT(17H COLUMN BORDER =   /(1X,8E15.8))
470   FORMAT(27H REDOBJ BACKING UP.  OBJ =   ,E15.8,5H LV =,I5,
     1 8H XSTEP =,E15.8)
480   FORMAT(51H REDOBJ CANNOT BACKUP TO FIRST VIOLATED CONSTRAINT   )
490   FORMAT(5H LV =,I5,8H THMIN =  ,E13.6)
500   FORMAT (30X,6HSTEP =,1PE15.6,2X,5HOBJ =,1PE15.6,2X,
     1 14HNEWTON ITERS =,I4)
510   FORMAT (30X,16HBASIC VARIABLE #,I4,15H VIOLATED BOUND)
520   FORMAT (30X,12HCONSTRAINT #,I4,15H VIOLATED BOUND)
530   FORMAT (30X,38HNEWTON FAILED TO CONVERGE.  NO.ITER. =,I4)
540   FORMAT (16H REDOBJ ENTERED   )
550   FORMAT (18H REDOBJ COMPLETED   )
C
C  END OF REDOBJ
C
      END
      SUBROUTINE PH1OBJ(INBC,G,UB,ALB)
C         **************************************************************
C     THIS SUBROUTINE COMPUTES PHASE 1 OBJECTIVE = SUM OF ABSOLUTE
C     VALUES OF CONSTRAINT VIOLATIONS AND STORES AS G(M+1). TRUE
C     OBJECTIVE SAVED AS TRUOBJ
C         **************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL MAXIM
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
C***********************************************************************
C
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON/BESTBK/STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      LOGICAL MOVE,RESTRT,DROP,VARMET,CONJGR
      COMMON /LOGBLK/ MOVE,RESTRT,DROP,VARMET,CONJGR
      COMMON/PH1BK/PHMULT,PH1EPS,INITPH
      COMMON/REDPH/TRUBST
      COMMON/MAXBK/MAXIM
      DIMENSION INBC(M),G(MP1),UB(NPMP1)
      DIMENSION ALB(NPMP1)
      DATA SINF0,TRUE0/0.0D0,0.0D0/
      FACTOR=100.0D0
      NINF=0
      SINF=0.0D0
      IF (NNBC.EQ.0) GO TO 300
      DO 200 I = 1,NNBC
          J=INBC(I)
          GJ=G(J)
          T=ALB(N+J)-GJ
          IF (T.LE.EPNEWT) GO TO 100
          NINF=NINF+1
          SINF=SINF+T
          GO TO 200
100       T=GJ-UB(N+J)
          IF (T.LE.EPNEWT) GO TO 200
          NINF=NINF+1
          SINF=SINF+T
200   CONTINUE
300   TRUOBJ=G(NOBJ)
      IF (NINF.EQ.0) RETURN
      T=TRUOBJ
      IF (MAXIM) T=-T
C
C     IF IN MIDDLE OF SEARCH, DO NOT RECOMPUTE PHMULT
C
      IF (INITPH.EQ.0) GO TO 400
C
C     IF STARTING PHASE 1, COMPUTE PHMULT FROM SCRATCH
C
      IF (INITPH.EQ.2) GO TO 350
      TRUOBJ=TRUBST
      T=TRUOBJ
      IF (MAXIM) T=-T
C
C     RIGHT AFTER CONSBS CALL, SEE IF SINF HAS CHANGED ENOUGH TO RESTART
C
      IF (PHMULT.EQ.0.0D0) GO TO 350
      IF (SINF0.GT.FACTOR*SINF) GO TO 350
      IF (TRUE0*TRUOBJ.LT.0.0) GO TO 350
      IF (DABS(TRUE0).LT.FACTOR*DABS(TRUOBJ).AND.DABS(TRUOBJ).LT.
     1 FACTOR*DABS(TRUE0)) GO TO 400
350   CONTINUE
      PHMULT=0.0D0
      SINF0=SINF
      TRUE0=TRUOBJ
      RESTRT=.TRUE.
      IF (DABS(TRUOBJ).LT.0.01D0) GO TO 400
      PHMULT=DABS(PH1EPS*SINF/TRUOBJ)
      IF (MAXIM) PHMULT=-PHMULT
      IF (IPR3.GE.2) WRITE (6,500) PHMULT,SINF,T
400   IF (NINF.NE.0) G(NOBJ)=SINF + PHMULT*TRUOBJ
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF (NINF .NE. 0 .AND. PENLTY) G(NOBJ)=G(NOBJ)+RPEN*PENTRM
C***********************************************************************
C
      IF (PH1EPS.NE.0.0D0.AND.IPR3.GE.2) WRITE (6,510) SINF,T
      RETURN
500   FORMAT (13H NEW PHMULT = ,E12.5,2X,25HSUM OF INFEASIBILITIES =   ,
     1 E15.8,2X,8HTRUOBJ =,E15.8)
510   FORMAT (7H SINF =,E15.8,9H TRUOBJ =,E15.8)
C
C     END OF PH1OBJ
C
      END
      SUBROUTINE SEARCH (BINV,GRAD,D,X,G,IBV,V,XB1,XB2,XB3,INBV,ALB,UB,X
     1STAT,ISTAT,INBC,GBEST,XBEST,IBC,ROWB,COLB,ROW,DBND)
C
C     THIS SUBROUTINE PERFORMS THE ONE-DIMENSIONAL MINIMIZATION
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      LOGICAL PENFLG
      LOGICAL  VIOL,PENLTY
      DIMENSION GRDPEN(300)
      COMMON /SUMT/ PENTRM,RPEN,GRDPEN,REDFAC,IPNCNT,NPEN,NHARD
      COMMON  /SUMTL/ PENLTY,VIOL
      COMMON /USOBJ/ USEOBJ,USBEST,PENBST
C***********************************************************************
C
      LOGICAL UNCONP,FAIL,FAILP,JSTFES,MXSTEP,UNBD,UNCON,SUCCES
      LOGICAL MAXIM
      REAL *4GRAD
      INTEGER *4M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON /DIMEN/ M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      DIMENSION BINV(NBMAX,NNBMAX), GRAD(MP1,N), D(N), X(NPMP1), G(MP1),
     1 IBV(M), V(NBMAX), XB1(NBMAX), XB2(NBMAX), XB3(NBMAX), INBV(N)
     2, ALB(NPMP1), UB(NPMP1), XSTAT(N), GBEST(MP1), XBEST(MPNBMX), ISTA
     3T(MP1), IBC(NBMAX), ROWB(NBMAX), COLB(NBMAX)
      DIMENSION INBC(M),ROW(NBMAX)
      DIMENSION DBND(NPNBMX)
      COMMON/NINTBK/NB,NNBC,NOBJ,NINF,NSUPER,IPR3,NCAND,IPR
      COMMON /LIMITS/ EPBOUN,EPNEWT,EPSPIV,ITLIM
      COMMON /COUNTS/ NFTN,NGRAD,NMINV,NNFAIL,NCALLS,NIT,NBS,NSTEPC,NDUB
      COMMON /NEWSRC/ ITER
      COMMON /EPSCOM/ EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      COMMON /TOLS/ EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON /SRCHLG/ UNCON,FAIL,JSTFES,MXSTEP,UNBD,SUCCES,UNCONP
      COMMON /BESTBK/ STPBST,OBJBST,STEP,STEPMX,TRUOBJ
      COMMON /MISC/ MAXR,NSEAR,JP,LV,JQQ
      COMMON /QUADBK/ A1,A2,A3,ICON,IQUAD
      COMMON/MAXBK/MAXIM
      COMMON/REDSER/NINFB
      COMMON /INITBK/  INIT,LASTCL
C         **************************************************************
C     INITIALIZATION FOR QUADRATIC EXTRAPOLATION SCHEME
C         **************************************************************
      IF (IPR.GE.5) WRITE (6,470)
      ICON=0
      A1=0.0D0
      IF (NB.EQ.0) GO TO 20
      DO 10 I=1,NB
          II=IBV(I)
10        XB1(I)=X(II)
20    CONTINUE
C
C     INITIALIZATION
C
      ITER=0
      LASTCL = 1
      FAIL=.FALSE.
      UNCON=.TRUE.
      SUCCES=.TRUE.
      MXSTEP=.FALSE.
      JSTFES=.FALSE.
      STEP=0.0D0
      LMQUIT=ITLIM/2
      IF (IQUAD.EQ.1) LMQUIT=4
      MAXCUT=7
      MAXDUB=30
      SA=0.0D0
      SB=STPBST
      FA=G(NOBJ)
      FTEMP=-PLINFY
      STPBST=0.0D0
      OBJBST=G(NOBJ)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      PENFLG = .FALSE.
      DD = 0.0D0
C***********************************************************************
C
      IF (NB.EQ.0) GO TO 40
      DO 30 I=1,NB
          J=IBV(I)
30        XBEST(I)=X(J)
40    CONTINUE
      DO 50 I=1,MP1
50        GBEST(I)=G(I)
      IF (NNBC.EQ.0) GO TO 70
      DO 60 I=1,NNBC
          J=INBC(I)
60        XBEST(NB+I)=G(J)
70    CONTINUE
      NINFB=NINF
      SMAX=0.0D0
      DO 80 I=1,NSUPER
          TS=DABS(D(I))
          IF (SMAX.LT.TS) SMAX=TS
80    CONTINUE
C
C     COMPUTE STEPMX
C
      STEPMX=PLINFY
      DO 90 I=1,NSUPER
          SI=D(I)
          IF (DABS(SI).LT.TOLZ) GO TO 90
          J=INBV(I)
          BT=UB(J)
          IF (SI.LT.0.0D0) BT=ALB(J)
          IF(DABS(BT) .GT. 1.0D20) BT = DSIGN(1.0D20,BT)
          TS=(BT-X(J))/SI
          IF (TS.GT.STEPMX) GO TO 90
          JQQ=I
          STEPMX=TS
90    CONTINUE
C
C     DETERMINE INITIAL STEP SIZE
C
      TS=EPS3/SMAX
      IF (SB.LT.TS) SB=TS
      PCTCHG=0.05D0
      TS=STEPMX
      DO 120 I=1,NSUPER
          J=INBV(I)
          XI=X(J)
          XSTAT(I)=XI
          XI=DABS(XI)
          SI=DABS(D(I))
          IF (SI.LT.TOLZ) GO TO 120
          IF (XI.LT.1.0D0) GO TO 100
          TS2=PCTCHG*XI/SI
          GO TO 110
100       TS2=PCTCHG/SI
110       IF (TS.GT.TS2) TS=TS2
120   CONTINUE
C
C     SET SB TO TS UNLESS PREVIOUS SEARCH WAS UNCONSTRAINED AND
C     SB IS SMALLER THAN TS.
C
      IF (.NOT.UNCONP.OR.SB.GT.TS) SB=TS
      IF (SB.GT.STEPMX) SB=STEPMX
      TMP=G(NOBJ)
      IF (MAXIM.AND.NINF.EQ.0) TMP=-TMP
      IF (IPR.GE.5) WRITE (6,310) TMP,SB,STEPMX
C
C     THIS LOOP COMPUTES FB AND CUTS BACK STEPSIZE IF FB > FA.
C
      DO 140 NCUT=1,MAXCUT
          FAILP=FAIL
          STEP=SB
          CALL REDOBJ (BINV,GRAD,X,G,IBV,V,XB1,XB2,XB3,INBV,ALB,UB,XSTAT
     1    ,D,ISTAT,INBC,GBEST,XBEST,IBC,ROWB,COLB,ROW,DBND)
          IF (FAIL) GO TO 125
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( .NOT. PENLTY ) GO TO 122
      IF ( VIOL ) GO TO 125
122   CONTINUE
C***********************************************************************
C
          FB=G(NOBJ)
          TMP=FB
          IF (MAXIM.AND.NINF.EQ.0) TMP=-TMP
          IF (JSTFES) GO TO 200
          IF (FB.LE.FA+EPS) GO TO 150
          FTEMP=FB
          SC=STEP
          GO TO 130
125       NSTEPC=NSTEPC+1
C
C     REDUCE STEPSIZE
C
130       SB=STEP/(2**NCUT)
140   CONTINUE
      GO TO 240
C
C     STEP REDUCTION PHASE COMPLETED--HAVE FOUND A BETTER POINT.
C     CONSIDER ALL POSSIBLE CASES.
C
150   CONTINUE
      IF (LV.NE.0) GO TO 210
      IF (MXSTEP) GO TO 220
      IF (FAILP) GO TO 190
      FC=FTEMP
C
C     BEGIN QUADRATIC INTERPOLATION BLOCK
C
      IF (IQUAD.EQ.0.OR.NB.EQ.0) GO TO 2160
      ICON=1
      A2=SB
      DO 2150 I = 1,NB
          II=IBV(I)
2150      XB2(I)=X(II)
      NEXT=3
2160   CONTINUE
C
C     INTERPOLATE IF A BRACKET HAS BEEN FOUND.
C
          IF (FC.GT.FB+EPS) GO TO 170
C
C     STEP DOUBLING PHASE
C
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      SC = SB + SB
C***********************************************************************
C
      DO 160 NDUB=1,MAXDUB
C
C     QUIT SEARCH IF NEWTON FAILURE ANTICIPATED.
C
          IF (ITER.GE.LMQUIT) GO TO 230
C***********************************************************************
C
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
C     THE NEXT CARD HAS BEEN COMMENTD OUT
C     SC = SB + SB
C
C***********************************************************************
          STEP=SC
          CALL REDOBJ (BINV,GRAD,X,G,IBV,V,XB1,XB2,XB3,INBV,ALB,UB,XSTAT
     1    ,D,ISTAT,INBC,GBEST,XBEST,IBC,ROWB,COLB,ROW,DBND)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
C
      IF ( .NOT. PENLTY ) GO TO 1300
      IF ( VIOL ) GO TO 1000
1300  CONTINUE
C***********************************************************************
C
          IF (FAIL) GO TO 180
          SC=STEP
          FC=G(NOBJ)
          TMP=FC
          IF (MAXIM.AND.NINF.EQ.0) TMP=-TMP
          IF (JSTFES) GO TO 200
C
C     BEGIN QUADRATIC INTERPOLATION BLOCK
C
      IF (IQUAD.EQ.0.OR.NB.EQ.0) GO TO 2260
      ICON=2
      GO TO (2200,2220,2240),NEXT
2200  DO 2210 I=1,NB
          II=IBV(I)
2210      XB1(I)=X(II)
      A1=SC
      NEXT=2
      GO TO 2260
2220  DO 2230 I=1,NB
          II=IBV(I)
2230      XB2(I)=X(II)
      A2=SC
      NEXT=3
      GO TO 2260
2240  DO 2250 I=1,NB
          II=IBV(I)
2250      XB3(I)=X(II)
      A3=SC
      NEXT=1
2260  CONTINUE
C
C     INTERPOLATE IF A BRACKET HAS BEEN FOUND.
C
          IF (FC.GT.FB+EPS) GO TO 170
      IF (NINF .GT. 0 .AND. FC .GE. (FB - EPS0)  ) GO TO 170
          IF (LV.NE.0) GO TO 210
          IF (MXSTEP) GO TO 220
C
C     MOVE 3 POINT PATTERN ONE STEP AHEAD.
C
          FA=FB
          SA=SB
          FB=FC
          SB=SC
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      SC = 2.0D0*SB - DD*SA
      GO TO 1200
1000  CONTINUE
C    CUT BACK STEP BY TWO THIRDS IF VIOLATION
      PENFLG = .TRUE.
      SC = SB + (SC - SB ) / 3.0D0
      DD = 1.0D0
      IF (DABS((SC - SB ) / SC ) .GE. 1.0D-07 ) GO TO 1200
      IF ( IPR .GE. 3 ) WRITE(6,1400)
      GO TO 270
1200  CONTINUE
C***********************************************************************
C
160   CONTINUE
      NDUB=MAXDUB
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( .NOT. PENLTY ) GO TO 250
      IF ( PENFLG ) GO TO 270
C***********************************************************************
C
      GO TO 250
C
C     INTERPOLATION PHASE
C
170   CONTINUE
      T2=SB-SA
      T3=SC-SA
      F2=(FB-FA)*T3
      F3=(FC-FA)*T2
      IF (DABS(F2-F3).LT.PLZERO) GO TO 260
C
C     SD IS MINIMUM POINT FOR QUADRATIC FIT.
C
      SD=SA+0.5D0*(T3*F2-T2*F3)/(F2-F3)
      IF (SD.LE.SA.OR.SD.GE.SC) GO TO 260
      IF (IPR3.LT.2) GO TO 178
      WRITE (6,330)
      TMP=FA
      TMP2=FB
      TMP3=FC
      IF (.NOT.MAXIM.OR.NINF.NE.0) GO TO 175
      TMP=-TMP
      TMP2=-TMP2
      TMP3=-TMP3
175   WRITE (6,320) SA,SB,SC,TMP,TMP2,TMP3,SD
178   CONTINUE
      IF (IPR3.GE.1) WRITE (6,450)
C
C     COMPUTE OBJECTIVE AT SD POINT.
C
      STEP=SD
      CALL REDOBJ (BINV,GRAD,X,G,IBV,V,XB1,XB2,XB3,INBV,ALB,UB,XSTAT,D,I
     1STAT,INBC,GBEST,XBEST,IBC,ROWB,COLB,ROW,DBND)
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( .NOT. PENLTY ) GO TO 1500
      IF ( VIOL ) LASTCL = 0
      IF ( VIOL ) GO TO 270
1500  CONTINUE
C***********************************************************************
C
      IF (FAIL) GO TO 180
      FD=G(NOBJ)
      IF (JSTFES) GO TO 200
      IF (LV.NE.0.AND.FD.LT.FB) GO TO 210
      IF ( FD .GT. FB ) LASTCL = 0
      GO TO 270
180   CONTINUE
C
C     QUIT BECAUSE NEWTON FAILED AND BETTER POINT FOUND.
C
      IF (IPR3.GE.1) WRITE (6,420)
      LASTCL = 0
      GO TO 270
190   CONTINUE
      IF (IPR3.GE.1) WRITE (6,380)
      GO TO 270
200   CONTINUE
C
C     HAVE JUST BECOME FEASIBLE
C
      IF (IPR3.GE.1) WRITE (6,460)
      RETURN
210   CONTINUE
C
C     BASIC VARIABLE HIT BOUND
C
      UNCON=.FALSE.
      NBS=NBS+1
      IF (IPR3.LT.1) GO TO 270
      I=IBV(LV)
      K=I-N
      IF (I.LE.N) WRITE (6,360) I
      IF (I.GT.N) WRITE (6,440) K
      GO TO 270
220   CONTINUE
C
C     SUPERBASIC VARIABLE HIT BOUND
C
      UNCON=.FALSE.
      IF (IPR3.GE.1) WRITE (6,370) INBV(JQQ)
      GO TO 270
230   CONTINUE
C
C     NEWTON TOOK TOO LONG
C
      IF (IPR3.GE.1) WRITE (6,390)
      GO TO 270
240   CONTINUE
C
C
      SUCCES=.FALSE.
      LASTCL = 0
      IF (IPR3.GE.1) WRITE (6,400) MAXCUT
      GO TO 270
250   CONTINUE
C
C     STEP SIZE DOUBLED NDUB TIMES
C
      UNBD=.TRUE.
      GO TO 270
260   CONTINUE
C
C     QUADRATIC INTERPOLATION OUT OF RANGE
C
      IF (IPR3.GE.1) WRITE (6,410)
      GO TO 270
270   CONTINUE
      UNCONP=UNCON
C
C     PICK UP BEST POINT ENCOUNTERED AND RETURN.
C
      STEP=STPBST
      NINF=NINFB
C
C***********************************************************************
C     CODE FOR INTERIOR PENALTY OPTION  -  HARD BOUNDS
      IF ( PENLTY ) USEOBJ = USBEST
      IF ( PENLTY ) PENTRM = PENBST
1400  FORMAT(52H0 C AND B POINTS ARE TOO CLOSE TOGETHER. STOP SEARCH   )
C***********************************************************************
C
      DO 280 I = 1,NSUPER
          J=INBV(I)
280       X(J)=XSTAT(I)+STEP*D(I)
       IF ( NB .EQ. 0 )  GO TO 302
      DO 300 I=1,NB
          J=IBV(I)
300       X(J)=XBEST(I)
  302   CONTINUE
      DO 290 I=1,MP1
          TS=GBEST(I)
          G(I)=TS
          X(N+I)=TS
290   CONTINUE
      IF (IPR.GE.5) WRITE (6,480)
      RETURN
C
310   FORMAT (12H0OBJECTIVE =,D15.7,5X,14HINITIAL STEP =,D15.7,5X,19HMAX
     1 STEP TO BOUND = ,D15.7)
320   FORMAT (3E13.6,11X,3E14.7,11X,E14.7)
330   FORMAT (7X,1HA,12X,1HB,12X,1HC,22X,2HFA,13X,2HFB,13X,2HFC,24X,1HD)
360   FORMAT(15X,16HBASIC VARIABLE #,I4,10H HIT BOUND )
370   FORMAT (15X,21HSUPERBASIC VARIABLE #,I4,10H HIT BOUND)
380   FORMAT(15X,22HEARLIER NEWTON FAILURE)
390   FORMAT (15X,26HANTICIPATED NEWTON FAILURE)
400   FORMAT (15X,30HNO OBJECTIVE IMPROVEMENT AFTER,I4,
     1 20H STEPSIZE REDUCTIONS)
410   FORMAT (35H QUADRATIC INTERPOLATION ABANDONED   )
420   FORMAT (15X,14HNEWTON FAILURE)
430   FORMAT (5H FD =,E15.8,5H FQ =,E15.8,8H FQVAR =,E12.5,7H NCUT =,I4)
440   FORMAT (15X,12HCONSTRAINT #,I4,14H NOW AT BOUND )
450   FORMAT (15X,23HQUADRATIC INTERPOLATION)
460   FORMAT (15X,74HALL VIOLATED CONSTRAINTS SATISFIED.  NOW BEGIN TO O
     1PTIMIZE TRUE OBJECTIVE   )
470   FORMAT (16H SEARCH ENTERED   )
480   FORMAT (18H SEARCH COMPLETED   )
C
C     END OF SEARCH
C
      END
      SUBROUTINE INITLZ
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TOLS/EPS,PLINFY,PLZERO,TOLX,TOLZ
      COMMON/EPSCOM/EPS0,EPS1,EPS2,EPS3,EPS4,EPS5
      COMMON/WORDBK/IWORDR,IWORDI
C
C  SET MACHINE DEPENDENT PARAMETERS
C
C  EPS IS MACHINE ACCURACY FOR VARIABLES STARTING WITH A-H,O-Z.
C  PLINFY IS 'PLUS INFINITY' FOR GIVEN MACHINE.
C  PLZERO IS POSITIVE VALUE NEAREST ZERO WITHOUT BEING ZERO ON MACHINE.
C
      EPS=1.0D-15
      PLINFY=1.0D30
      PLZERO=1.0D-30
C
C  SET MACHINE PRECISION CONSTANTS
C
      EPS0=EPS**(4.0/5.0)
      EPS1=EPS**(2.0/3.0)
      EPS2=EPS**(1.0/2.0)
      EPS3=EPS**(1.0/3.0)
      EPS4=EPS**(1.0/4.0)
      EPS5=EPS**(1.0/5.0)
C
C  SET TOLERANCES
C
      TOLZ=EPS0
      TOLX=DMAX1(EPS2,1.0D-6)
C
C  SET WORD LENGTHS FOR THE PARTICULAR MACHINE
C  IWORDI IS THE NUMBER OF INTEGER VARIABLES STARTING I-N CONTAINED
C  IN ONE WORD OF THE ZZ ARRAY.
C  IWORDR IS THE NUMBER OF REAL VARIABLES (GRAD ARRAY) CONTAINED IN
C  ONE WORD OF THE ZZ ARRAY.
C  FOR THE UNIVAC 1108 THERE ARE NO HALF-WORD INTEGERS.
C
      IWORDI=2
      IWORDR=2
C
C  END OF INITLZ
C
      RETURN
      END
      SUBROUTINE SETUP(NCORE)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 LASTZ
      INTEGER*4 NCORE
      INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/GOBK/NVARS,NROWS,MAXR,MAXB
      COMMON/DYNAM1/KX,KG,KALB,KUB,KICAND,KISTAT,KIFIX
      COMMON/DYNAM2/KU,KGRADF,KINBC,KIBC,KIBV,KINBV,KIUB
      COMMON/DYNAM3/KR,KV,KD,KGBEST,KXBEST,KXB1,KXB2,KXB3
      COMMON/DYNAM4/KDBND,KCNORM,KXSTAT,KGG,KRR,KY,KGRADP
      COMMON/DYNAM5/KROWB,KCOLB,KBINV,KGRAD,KICOLS,KINORM
       COMMON/DYNAM6/KX0,KG0,KCON,KVAR,KINBVP
      COMMON/DIMEN/M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
      COMMON/WORDBK/IWORDR,IWORDI
      COMMON/SETIN/LASTZ
C
C  CHECK FOR PUNCHING IN WRONG COLUMNS
C
      IF (NVARS.GT.0.AND.NVARS.LT.1000) GO TO 10
      WRITE (6,1060) NVARS,NROWS,MAXR,MAXB
      WRITE (6,1070)
      STOP
10    IF (NROWS.GT.0.AND.NROWS.LT.1000) GO TO 20
      WRITE (6,1060) NVARS,NROWS,MAXR,MAXB
      WRITE (6,1080)
      STOP
C
C  COMPUTE VARIABLES USED IN DIMENSION STATEMENTS
C
20    N=NVARS
      MP1=NROWS
      M=MP1-1
      IF (M.EQ.0) M=1
      NPMP1=N+MP1
      NBMAX=MAXB
      IF (NBMAX.LE.0.OR.NBMAX.GT.M) NBMAX=M
      IF (MAXR.LE.0.OR.MAXR.GT.N) MAXR=N
      NNBMAX=N
      IF (NBMAX.GT.N) NNBMAX=NBMAX
      MPNBMX=MP1
      NPNBMX=N+NBMAX
      NRTOT=MAXR*(MAXR+1)/2
C
C  COMPUTE INCREMENTS TO INTEGER AND REAL VARIABLE POSITIONS IN Z ARRAY
C
      NK=N/IWORDI+1
      MK=M/IWORDI + 1
      MP1K=MP1/IWORDI+1
      NPMP1K=NPMP1/IWORDI+1
      NBMAXK=NBMAX/IWORDI+1
      NNBK=NNBMAX/IWORDI+1
      NPNBK=NPNBMX/IWORDI + 1
      KGRADK=(MP1*N)/IWORDR + 1
C
C  COMPUTE INITIAL POSITIONS OF ALL ARRAYS IN Z VECTOR
C
      KX=1
      KG=KX+NPMP1
      KALB=KG+MP1
      KUB=KALB+NPMP1
      KR=KUB+NPMP1
      KGRADF=KR+NRTOT
      KV=KGRADF+N
      KD=KV+NBMAX
      KU=KD+N
      KGBEST=KU+NBMAX
      KXBEST=KGBEST+MP1
      KXB1=KXBEST+MP1
      KXB2=KXB1+NBMAX
      KXB3=KXB2+NBMAX
      KDBND=KXB3+NBMAX
      KCNORM=KDBND+NPNBMX
      KXSTAT=KCNORM+NPNBMX
      KGG=KXSTAT+N
      KRR=KGG+MP1
      KY=KRR+NBMAX
      KGRADP=KY+N
      KROWB=KGRADP+N
      KCOLB=KROWB+NBMAX
      KCON=KCOLB+NBMAX
      KVAR=KCON+MP1
      KX0=KVAR+N
      KG0 = KX0 + N
      KBINV = KG0 + MP1
      KGRAD=KBINV+NBMAX*NNBMAX
      KISTAT=KGRAD+KGRADK
      KINBV=KISTAT+MP1K
      KIUB = KINBV + NK
      KINBC=KIUB+NK
      KIBC=KINBC+MP1K
      KIBV=KIBC+NBMAXK
      KIFIX=KIBV+MK
      KICOLS=KIFIX+NK
      KINORM=KICOLS+NPNBK
      KINBVP=KINORM+NPNBK
      KICAND=KINBVP+NK
      KNEXT=KICAND+NPNBK
      LASTZ=KNEXT-1
      IF (LASTZ.LE.NCORE) RETURN
      WRITE (6,1000) LASTZ,NCORE
      WRITE (6,1060) N,MP1,MAXR,NBMAX
      STOP
1000  FORMAT(29H0ACTUAL LENGTH OF ZZ ARRAY =   ,I6,
     1 18H EXCEEDED NCORE =   ,I6)
1060  FORMAT(24H0NUMBER OF VARIABLES IS   ,I5/
     1 19H0NUMBER OF ROWS IS   ,I5/9H0MAXR IS   ,I5/
     2 35H0CEILING ON BINDING CONSTRAINTS IS   ,I5)
1070  FORMAT(26H0IMPROPER VALUE OF NVARS   )
1080  FORMAT(26H0IMPROPER VALUE OF NROWS   )
C
C  END OF SETUP
C
      END
