! modulo que contém as subrotinas que são utilizadas para calcular os
!  expoentes de lyapunov a tempo finito de um sistema arbitrário

module clv
  use tipos, only : dp
  implicit none

  procedure(), pointer :: fex
  integer :: nn
  integer :: neq, liw, lrw
  real(dp), dimension(:), allocatable, save :: xx, rwork, u
  integer, dimension(:), allocatable, save :: iwork

  real(dp), dimension(:,:), allocatable :: m,jac,r,cm,q

  save
  private
  public :: calcula_lyap, inicializa_lyap, finaliza_lyap,usang

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
    
    allocate(xx(neq), rwork(lrw), iwork(liw),u(nn))
    allocate(m(nn,nn), jac(nn,nn), r(nn,nn), cm(nn,nn), q(nn,nn) )
    return
  end subroutine inicializa_lyap

  ! ****************************************************************
  ! Apenas dealloca os vetores allocados
  subroutine finaliza_lyap
    deallocate(xx, rwork, iwork, u)
    deallocate(m, jac, r , q)
    return
  end subroutine finaliza_lyap
  
  ! ****************************************************************
  ! Subrotina que calcula o valor dos expoentes de lyapunov
  !  propriamente ditos
  
  subroutine calcula_lyap( n, x, t, tlyap, norto, tstep, lyap, itol,&
       & itask, iopt, jt, mf, rtol, atol, icall, qmatrix, rout )

    use matqrinv, only : qr
    
    integer, intent(in), target :: n
    integer, intent(in) :: norto
    real(dp), dimension(n), intent(inout) :: x
    real(dp), dimension(n), intent(out) :: lyap
    real(dp), dimension(n,n), intent(out), optional :: rout, qmatrix
    real(dp), intent(inout) :: t
    real(dp), intent(in) :: tlyap, rtol, atol, tstep
    integer(dp), intent(in) :: itol, itask, iopt, jt, mf
    integer, intent(inout) :: icall
    
    integer :: istate, i, j, iorto
    real(dp) :: tout, tcount, torto

    lyap = 0.0_dp
    
    if ( norto /= 0 ) then
       torto = tlyap / real(norto)
    else
       write(0,"(' NORTO MUST BE DIFFERENT FROM ZERO ')")
       call exit(2)
    end if

    icall = icall + 1
    
    testa_icall : if( icall == 1 ) then
       ! Cond. Inicial para a parte nao linear
       xx(1:n) = x
       
       ! cond. inicial para a parte linear       
       xx(n+1:neq) = 0.0_dp
       identidade : do i = 1, n
          xx((n+1)*i) = 1.0_dp
       end do identidade
    end if testa_icall
    
    lyap_loop : do iorto = 1, norto
       
       istate = 1
       tcount = 0.0
       
       orto_loop : do while ( tcount <= torto ) 
          tout = t + tstep
          tcount = tcount + tstep
          
          call dlsode(flyap,neq,xx,t,tout,itol,rtol,atol,itask,istate,iopt&
               &,rwork,lrw,iwork,liw,jt,mf)
       end do orto_loop
       
       do i = 1, nn
          do j = 1, nn
             m(i,j) = xx(i*nn+j)
          end do
       end do
       
       call qr(m,q,r)
       
       forall(i=1:n)
          lyap(i) = lyap(i) + log(r(i,i))
       end forall

     do i = 1, nn
          do j = 1, nn
             xx(i*nn+j) = q(i,j)
          end do
       end do
       
    end do lyap_loop
    
    lyap = lyap / tlyap 
    
    if( present(rout) ) then
       rout = r
    end if

    if( present(qmatrix) ) then
       qmatrix = q
    end if

    x = xx(1:n)
    
    return
  end subroutine calcula_lyap

  ! ****************************************************************
  ! Subrotina que retorna o valor de f como função de x para o caso
  !  da matriz de monodromia

  subroutine flyap(n,t,x,xdot)
    use jacob, only: jacobian
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(n), intent(inout) :: x, xdot
    real(dp), intent(in) :: t
    integer, save :: count = 0

    integer :: i, j
    
    count = count + 1
    
    ! Calcula a derivada da parte não linear
    call fex(nn,t,x(1:nn),xdot(1:nn))
    
    ! Calculo das derivadas da matriz de monodromia
    do i = 1, nn
       do j = 1, nn
          m(i,j) = x(i*nn+j)
       end do
    end do

    call jacobian(fex,nn,x(1:nn),t,jac)

    m = matmul(jac,m)
    
    do i = 1, nn
       do j = 1, nn
          xdot(i*nn+j) = m(i,j)
       end do
    end do

    return
  end subroutine flyap

  ! ****************************************************************
  ! Função que simplesmente normaliza um vetor dado x
  
  function normaliza(x) result(y)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: y
    y = x / sqrt(dot_product(x,x))
    return
  end function normaliza

  ! ****************************************************************
  ! Função que retorna uma matriz simétrica, cuja diagonal principal
  !  é nula, e o elemento i,j contém o ângulo entre os vetores de
  !   lyapunov i e j
  
  function usang(q,rold) result(theta)
    use matqrinv, only : inv
    real(dp), dimension(:,:), intent(in) :: q
    real(dp), dimension(:,:,:), intent(in) :: rold
    real(dp), dimension(size(q,1),size(q,1)) :: theta

    integer :: j, k, i, n, nr, istatus, ir
    real(dp) :: covcos
    
    n = size(q,1)
    nr = size(rold,3)
    
    linha : do j = 1, n
       coluna : do k = 1, n
          u = 0.0_dp
          cria_u : do i = 1, k
             u = u + q(:,i)
          end do cria_u
          cm(j,k)=dot_product(q(:,j),u)
       end do coluna
    end do linha

    ! normaliza
    do j = 1, n
       cm(:,j) = normaliza(cm(:,j))
    end do
    
    ! itera pra trás
    do ir = nr, 1, -1
       cm = matmul(inv(rold(:,:,ir),istatus),cm)
       do j = 1, n
          cm(:,j) = normaliza(cm(:,j))
       end do
    end do
    
    ! calcula o ângulo
    do i = 1, n-1
       do j = i+1,n
          covcos = dot_product(cm(:,i),cm(:,j))
          covcos = max( -1.0 , min( 1.0, covcos ) )
          theta(i,j) = acos(covcos)
          theta(j,i) = theta(i,j)
       end do
    end do
    
    do i = 1, n
       theta(i,i) = 0.0_dp
    end do
    return
    
  end function usang
  
  
end module clv

  
