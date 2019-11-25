!******************************************************************************
! Function to be minimised for the nelder mead algorithm
!******************************************************************************
module mo_obj_func
  implicit none
  PRIVATE
  PUBLIC :: obj_func
contains
  function obj_func(p)
    use mo_kind, only     : i4, dp
    use varfit, only      : nbins, gamma, nh
    use mo_setVario, only : tVar
    implicit none
    real(dp), dimension(:), intent(in) :: p
    real(dp) :: obj_func
    real(dp) :: gcal
    integer(i4) :: k 
    !
    obj_func = 0.0_dp
    !
    do k=1,nbins
      if (gamma(k,1) > 0.0_dp) then
        ! gcal = tvar (gamma(k,1), xu(1), xu(2), xu(3))
        gcal = tvar(gamma(k,1), p(1), p(2), p(3))
        !
        ! estimator l1 weighted
        if (gamma(k,1) <= 1.0_dp ) obj_func = obj_func + dabs(gamma(k,2) - gcal) * dble(Nh(k))
      end if
    end do
  end function obj_func
end module mo_obj_func
