  ! modulo com a declaração da subrotina que calcula a jacobiana do
  !  sistema

module jacob
  use tipos, only : dp
  implicit none

  private
  public :: jacobian

contains
!!$
!!$  subroutine jacobian( fex, n, x, t, jac )
!!$    external fex
!!$
!!$    integer, intent(in) :: n
!!$    real(dp), dimension(n), intent(in) :: x
!!$    real(dp), intent(in) :: t
!!$    real(dp), dimension(n,n), intent(out) :: jac
!!$
!!$    real(dp), parameter :: h = 1.0d-3
!!$    real(dp), dimension(n) :: xdot, xdot2, x2
!!$    real(dp) :: aux
!!$    integer :: i
!!$
!!$    x2 = x
!!$    call fex(n,t,x2,xdot)
!!$
!!$    do i = 1, n
!!$       aux = x2(i)
!!$       x2(i) = x2(i) + h
!!$       call fex(n,t,x2,xdot2)
!!$       jac(:,i) = ( xdot2 - xdot ) / h
!!$       x2(i) = aux
!!$    end do
!!$
!!$    return
!!$  end subroutine jacobian

  
  subroutine jacobian( fex, n, x, t, jac )
    external fex

    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: x
    real(dp), intent(in) :: t
    real(dp), dimension(n,n), intent(out) :: jac

    real(dp), parameter :: h = 1.0d-3
    real(dp), dimension(n), save :: xdot, xdot2, x2
    real(dp) :: aux
    integer :: i

    x2 = x
    call fex(n,t,x2,xdot)

    do i = 1, n
       aux = x2(i)
       x2(i) = x2(i) + h
       call fex(n,t,x2,xdot2)
       jac(:,i) = ( xdot2 - xdot ) / h
       x2(i) = aux
    end do

    return
  end subroutine jacobian
end module jacob

    
