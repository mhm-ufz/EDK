!*****************************************************************************
!    Main_opti.for
!    AUTHOR
!         Luis E. Samaniego-Eguiguren, IREUS, 28.05.1999
!    DESCRIPTION
!         Initialization of the Nonlinear Optimization Subroutine GRG2 
!         The function to be optimized is suplied in subroutine GCOMP
!    DEFINITION OF INPUT/OUTPUT-CHANNELS
!         SCREEN:            *
!
!         SCREEN OUTPUT:     *
!         OUTPUT FILE:       6
!******************************************************************************
subroutine OPTI
  use VarFit  
  use mo_kind, only: i4, dp
  ! inputs for GRG2
  IMPLICIT  DOUBLE PRECISION(A-H,O-Z), INTEGER(I,J,L,M,N)
  INTEGER*4 NCORE,NNVARS,NFUN,MAXBAS,MAXHES,LASTZ
  INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
  LOGICAL   MAXIM,INPRNT,OTPRNT
  DIMENSION BLVAR(100),BUVAR(100),BLCON(10),BUCON(10)
  DIMENSION RAMCON(10),RAMVAR(10),INBIND(10),Z(5000)
  DIMENSION NONBAS(10),REDGR(10),DEFAUL(20),TITLE(19)
  DIMENSION XX(100),FCNS(10),RMULTS(10)
  DATA      BIG/1.D31/
  ! Initialization of GRG2
  !
  ! Minimize the Obj. Function
  !
  DEFAUL(19) = 1.0
  !
  ! Set problem size and objective
  !
  NCORE  = 5000	   !DIMENSION OF THE {Z} ARRAY         
  NNVARS = nParam	 !NUMBER OF VARIABLES
  NFUN   = 1			 !NUMBER OF FUNCTIONS INCLUDING OBJECTIVE
  MAXBAS = 1			 !UPPER LIMIT ON NUMBER OF BINDINGCONSTRAINTS: {NFUN}
  MAXHES = nParam  !MAXIMUM ALLOWABLE SIZE OF	HESSIAN	         : {NVARS}
	NNOBJ  = 1			 !INDEX OF COMPONENT OF VECTOR {G} CORRESPONDING TO
                   !OBJECTIVE FUNCTION
  DO I=1,18
    DEFAUL(I) = 1.0_dp
  end do
  !
  ! Set tolerances
  !
	DEFAUL(1) = -1.0
	FPNEWT    = 1.0D-7
	DEFAUL(2) = -1.0
	FPINIT    = 1.0D-6
	DEFAUL(3) = -1.0
	FPSTOP    = 1.0D-6
	DEFAUL(6) = -1.0
	NNSTOP    = 10
  !
  ! Quadratic extrapolation
  !
  DEFAUL(15) = -1.0_dp
	IIQUAD     = 1
  !
  ! Central difference approximation
  !
  DEFAUL(15) = -1.0_dp
  LDERIV     = 1
  !
  ! Set variable bounds
  !
  do i=1,nParam
    BLVAR(i) = 1.0d-2
    BUVAR(i) = 1.0_dp        
  end do

  ! lower boundary of variogram parameters
  ! BLVAR(1) = 1.6d-1    ! for pre
  ! BLVAR(2) = 6.00d-1
  ! BLVAR(3) = 1.0d-2

  ! upper boundary of variogram parameters
  ! BUVAR(1) = 2.0d-1    ! for pre   
  ! BUVAR(2) = 9.00d-1       
  ! BUVAR(3) = 3.0_dp       

  !
  ! Set fuctions bounds
  !
  BLCON(1) = 0._dp
  BUCON(1) = BIG
  !    
  ! Set initial values
  !
  !scaling range
  beta(3)=beta(3)/gmax(1)
  ! assigning
  do i=1,nParam
    XX(i) = beta(i)
  end do

  !
  ! Other stuff
  !
  ! INPRNT = .FALSE.
  !	OTPRNT = .FALSE.
  INPRNT = .TRUE.
  OTPRNT = .TRUE.
  !
  ! write report USE:  OPEN(UNIT = 6,FILE = 'Report_OPT.sol') and do not CLOSE 6
  !
  ! Call Nonlinear Optimization subroutine
  !
  CALL GRGSUB(INPRNT,OTPRNT,NCORE,NNVARS,NFUN,MAXBAS,&
      MAXHES,NNOBJ,TTITLE,BLVAR,BUVAR,BLCON,BUCON,DEFAUL,FPNEWT,FPINIT,&
      FPSTOP,FPSPIV,PPH1EP,NNSTOP,IITLIM,LLMSER,IIPR,IIPN4,IIPN5,&
      IIPN6,IIPER,IIDUMP,IIQUAD,LDERIV,MMODCG,&
      RAMCON,RAMVAR,XX,FCNS,INBIND,RMULTS,NONBAS,REDGR,&
      NBIND,NNONB,INFORM,Z)
  !
  ! delete report (partial results)
  !
  close (6, STATUS='DELETE')
  !
  ! keep parameters
  do i=1,nParam
    beta(i) = XX(i)
  end do
  ! estimate statistics
  call stats
  !
  ! scale up distance h
  where (gamma(:,1) > 0._dp) gamma(:,1) = gamma(:,1) * gmax(1)

  ! scale back parameters (only for range a)
  beta(3)=beta(3)*gmax(1)

  !print*, 'opt co =', beta(1)
  !print*, 'opt c  =', beta(2)
  !print*, 'opt a  =', beta(3)

end subroutine
