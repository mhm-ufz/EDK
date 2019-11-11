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
  use mo_obj_func, only: obj_func
  use mo_nelmin, only: nelmin

  ! parameters for Nelder-Mead algorithm
  real(dp) :: pstart(3) ! Starting point for the iteration.
  real(dp) :: pmin(3) ! optimized parameter set
  real(dp) :: prange(3, 2) ! Range of parameters (upper and lower bound).
  real(dp) :: varmin ! the terminating limit for the variance of the function values. varmin>0 is required
  real(dp) :: step(3) ! determines the size and shape of the initial simplex. The relative magnitudes of its elements should reflect the units of the variables. size(step)=size(start)
  integer(i4) :: konvge ! the convergence check is carried out every konvge iterations
  integer(i4) :: maxeval ! the maximum number of function evaluations. default: 1000
  real(dp) :: funcmin
  integer(i4) :: neval ! the number of function evaluations used.
  integer(i4) :: numrestart ! the number of restarts.
  integer(i4) :: ierror ! error indicator.
                        ! 0: no errors detected.
                        ! 1: varmin or konvge have an illegal value.
                        ! 2: iteration terminated because maxeval was exceeded without convergence.
  real(dp), allocatable :: history(:)
  
  
  ! ! inputs for GRG2
  ! IMPLICIT  DOUBLE PRECISION(A-H,O-Z), INTEGER(I,J,L,M,N)
  ! INTEGER*4 NCORE,NNVARS,NFUN,MAXBAS,MAXHES,LASTZ
  ! INTEGER*4 M,N,MP1,NPMP1,NBMAX,NNBMAX,NPNBMX,MPNBMX,NRTOT
  ! LOGICAL   MAXIM,INPRNT,OTPRNT
  ! DIMENSION BLVAR(100),BUVAR(100),BLCON(10),BUCON(10)
  ! DIMENSION RAMCON(10),RAMVAR(10),INBIND(10),Z(5000)
  ! DIMENSION NONBAS(10),REDGR(10),DEFAUL(20),TITLE(19)
  ! DIMENSION XX(100),FCNS(10),RMULTS(10)
  ! DATA      BIG/1.D31/

  ! Initialization of Nelder-Mead
  pstart = (/0.05, 100., 2./) ! Starting point for the iteration.
  prange(:, 1) = (/0., 0., 0./) ! Range of parameters (lower bound).
  prange(:, 2) = (/2., 300., 5./) ! Range of parameters (upper bound).
  varmin = 0.01 ! the terminating limit for the variance of the function values. varmin>0 is required
  step = (/0.001, 0.05, 0.01/) ! determines the size and shape of the initial simplex. The relative magnitudes of its elements should reflect the units of the variables. size(step)=size(start)
  konvge = 10 ! the convergence check is carried out every konvge iterations
  maxeval = 2000 ! the maximum number of function evaluations. default: 1000

  ! Call Nelder-Mead optimizer to reduce GCOMP
  pmin = nelmin(obj_func, pstart, varmin, step, konvge, maxeval, &
                funcmin, neval, numrestart, ierror, history)
  ! pmin = nelminrange(obj_func, pstart, prange, varmin, step, konvge, maxeval, &
  !                    funcmin, neval, numrestart, ierror)

  ! scale up distance h
  where (gamma(:,1) > 0._dp) gamma(:,1) = gamma(:,1) * gmax(1)

  ! scale back parameters (only for range a)
  pmin(2)=pmin(2)*gmax(1)


  print *, neval, ' of ', maxeval
  print *, "funcmin: ", funcmin
  print *, "p_obj:   ", pmin
  print *, 'error: ', ierror
  print *, 'varmin: ', varmin
  print *, 'history: ', history(1), history(size(history))
  
  ! stop 'Testing Nelder Mead algorithm'

  ! ! estimate statistics
  ! call stats
  ! !
  ! ! scale up distance h
  ! where (gamma(:,1) > 0._dp) gamma(:,1) = gamma(:,1) * gmax(1)

  ! ! scale back parameters (only for range a)
  ! beta(3)=beta(3)*gmax(1)

  ! !print*, 'opt co =', beta(1)
  ! !print*, 'opt c  =', beta(2)
  ! !print*, 'opt a  =', beta(3)

end subroutine OPTI
