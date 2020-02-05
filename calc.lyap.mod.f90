! modulo que contém as subrotinas que são utilizadas para calcular os
!  expoentes de lyapunov a tempo finito de um sistema arbitrário

module feltf
  use tipos, only : dp
  implicit none

  procedure(), pointer :: fex
  integer :: nn
  integer :: neq, liw, lrw
  real(dp), dimension(:), allocatable :: xx, rwork
  real(dp), dimension(:), allocatable :: y_jac, ydotp
  integer, dimension(:), allocatable :: iwork

  real(dp), dimension(:,:), allocatable :: m,jac

  save
  private
  public :: calcula_lyap, inicializa_lyap, finaliza_lyap, gs_ortho

contains

  ! ****************************************************************
  ! Apenas aloca os vetores necessários para o cálculo dos expoentes
  !  de lyapunov
  subroutine inicializa_lyap(f,n)
    external f
    integer, intent(in) :: n

    nn = n
    fex => f

    neq = n + n * n
    liw = 20
    lrw = 20 + neq * 16

    allocate(xx(neq), rwork(lrw), iwork(liw))
    allocate(y_jac(nn),ydotp(nn))
    allocate(m(nn,nn), jac(nn,nn) )
    return
  end subroutine inicializa_lyap

  ! ****************************************************************
  ! Apenas dealloca os vetores allocados
  subroutine finaliza_lyap
    deallocate(xx, rwork, iwork)
    deallocate(y_jac,ydotp)
    deallocate(m, jac )
    return
  end subroutine finaliza_lyap

  ! ****************************************************************
  ! Subrotina que calcula o valor dos expoentes de lyapunov
  !  propriamente ditos

  subroutine calcula_lyap( n, x, t, tlyap, tstep, lyap, itol,&
       & itask, iopt, jt, mf, rtol, atol, icall, vectors )

    integer, intent(in), target :: n
    real(dp), dimension(n), intent(inout) :: x
    real(dp), dimension(n), intent(out) :: lyap
    real(dp), dimension(n,n), intent(out), optional :: vectors
    real(dp), intent(inout) :: t
    real(dp), intent(in) :: tlyap, rtol, atol, tstep
    integer(dp) :: itol, itask, iopt, jt, mf
    integer, intent(inout) :: icall

    integer :: istate, i, j
    real(dp) :: tout, tcount

    istate = 1
    icall = icall + 1

    ! optional argument
    iopt = 1
    iwork(5) = 4

    if( icall == 1 ) then

       ! Cond. Inicial para a parte nao linear
       xx(1:nn) = x

       ! cond. inicial para a parte linear
       xx(nn+1:neq) = 0.0_dp
       ident : do i = 1, n
          xx((nn+1)*i) = 1.0_dp
       end do ident
    end if

    tcount = 0.0

    do while ( tcount < tlyap )
       tout = t + tstep
       tcount = tcount + tstep
       call dlsode(flyap,neq,xx,t,tout,itol,rtol,atol,itask,istate,iopt&
            &,rwork,lrw,iwork,liw,jt,mf)
    end do

    call gs_ortho(n,xx,lyap)
    lyap = lyap / tlyap

    if( present(vectors) ) then
       do i = 1, n
          do j = 1, n
             vectors(i,j) = xx( i * n + j )
          end do
       end do
    end if

    x = xx(1:nn)

    return
  end subroutine calcula_lyap

  ! ****************************************************************
  ! Subrotina que retorna o valor de f como função de x para o caso
  !  da matriz de monodromia

  subroutine flyap(nin,t,y,ydot)
    implicit none
    integer, intent(in) :: nin
    real(dp), dimension(nin), intent(inout) :: y, ydot
    real(dp), intent(in) :: t

    integer :: i, j

    ! Calcula a derivada da parte não linear
    call fex(nn,t,y(1:nn),ydot(1:nn))

    ! Calculo das derivadas da matriz de monodromia
    forall( i = 1 : nn, j = 1 : nn )
       m(i,j) = y(i*nn+j)
    end forall

    call jacobian(y(1:nn),ydot(1:nn),t)

    m = matmul(jac,m)

    forall( i = 1 : nn, j = 1 : nn )
       ydot(i*nn+j) = m(i,j)
    end forall

    return
  end subroutine flyap

  ! ****************************************************************
  ! ****************************************************************
  ! Subrotina que calcula a jacobiana
  subroutine jacobian(yin,ydotin,t)
    real(dp), dimension(nn), intent(in) :: yin, ydotin
    real(dp), intent(in) :: t
    real(dp), parameter :: h = 1.0d-6

    real(dp) :: aux
    integer :: i

    y_jac = yin

    do i = 1, nn
       aux = y_jac(i)
       y_jac(i) = aux + h
       call fex(nn,t,y_jac,ydotp)
       jac(:,i) = ( ydotp - ydotin ) / h
       y_jac(i) = aux
    end do
    return
  end subroutine jacobian

  ! ****************************************************************
  ! Subrotina que ortonormaliza as equações e retorna os expoentes de
  !  lyapunov
  subroutine gs_ortho(n,y,cum)
    use tipos, only: dp
    implicit none
    integer :: n,j,k,l
    real(dp), dimension(n) :: cum,znorm,gsc
    real(dp), dimension(n*n+n) :: y

    ! construct a new orthonormal basis by gram-schimidt
    ! normalize first vector
    znorm(1)=0.
    do j=1,n
       znorm(1)=znorm(1)+y(n*j+1)**2
    end do
    znorm(1)=sqrt(znorm(1))
    do j=1,n
       y(n*j+1)=y(n*j+1)/znorm(1)
    end do

    ! generate new orthonormal set
    ! make j-1 gsr coefficients
    do j=2,n
       do k=1,j-1
          gsc(k)=0.
          do l=1,n
             gsc(k)=gsc(k)+y(n*l+j)*y(n*l+k)
          end do
       end do

       ! construct a new vector
       do k=1,n
          do l=1,j-1
             y(n*k+j)=y(n*k+j)-gsc(l)*y(n*k+l)
          end do
       end do

       ! calculate the vector's norm
       znorm(j)=0.
       do k=1,n
          znorm(j)=znorm(j)+y(n*k+j)**2
       end do
       znorm(j)=sqrt(znorm(j))

       ! normalize the new vector
       do k=1,n
          y(n*k+j)=y(n*k+j)/znorm(j)
       end do
    end do

    ! update running vector magnitudes
    do k=1,n
       cum(k)=log(znorm(k))
    end do

    return
  end subroutine gs_ortho



end module feltf


