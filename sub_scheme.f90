!----------------------------------------------------------------------------------
Subroutine init_scheme
use var
implicit none
integer:: i,j,alpha_scd

	!Coefficients of Classic and Optimized Schemes
	allocate (C_cls(-alpha:alpha),C_opt(-alpha:alpha),C_trunc(-alpha-tau:alpha+tau))
	
	!Coefficients of the classical finite-difference scheme
	!for solving n-th derivatives with a stencil of [-alpha:alpha]
	call cal_cof_cls(n,alpha,C_cls)
	
	!initialize C_opt
	C_opt=C_cls
	
	!Coefficients of the classical finite-difference scheme
	!for solving n-th derivatives with a larger stencil of 
	![-alpha-tau:alpha+tau]
	call cal_cof_cls(n,alpha+tau,C_trunc)

	!Calculate coefficients of truncation error(approximated by
	!classical scheme with larger stencils)
	C_trunc(-alpha:alpha)=C_trunc(-alpha:alpha)-C_cls(-alpha:alpha)

	!Coefficients of standard second order central format
	allocate (der_scd(-alpha:alpha,0:2*alpha))
	der_scd(:,:)=0.q0
	do i=0,2*alpha
		alpha_scd=(i+mod(i,2))/2
		call cal_cof_cls(i,alpha_scd,der_scd(-alpha_scd:alpha_scd,i))
	enddo 

	!Discrete Taylor Expansion Coefficients of 
	!the Classical Finite Difference Scheme of the Nth Order Derivative
	allocate (tlc_cls(0:2*P))
	!call a subroutine to solve for these coefficients
	call cal_cof_tlc(n,P,tlc_cls)
	
End
!----------------------------------------------------------------------------------
!Calculate the coefficients of the classical finite difference scheme
subroutine cal_cof_cls(n,alpha,C_cls)
implicit none
Real(Kind=16):: fact
integer :: n,alpha
integer(kind=8):: i,j !Large number calculations require double precision
real(Kind=16):: C_cls(-alpha:alpha)
!矩阵和解向量
real(Kind=16):: T(0:2*alpha,-alpha:alpha),IM(0:2*alpha)	

	do i=0,2*alpha 
	do j=-alpha,alpha 
		!The i-th order Taylor expansion coefficient of u(0) at point j
		T(i,j)=j**i/(1.Q0*fact(i))	
	enddo
	enddo
	IM(:)=0.Q0
	IM(n)=1.q0
	call Guass_solver(size(C_cls),1,T(0:2*alpha,-alpha:alpha),IM(0:2*alpha),C_cls(-alpha:alpha))
	!~ 	  print*, "C_cls:","n=",n,"alpha=", alpha 
	!~ 	  do i=-alpha,alpha
	!~ 	      Write(*,"(30f20.15)") C_cls(i)
	!~ 	  end do

end
!-------------------------------------------------------
! Calculate derivatives
subroutine Cal_der(nx,KL,KR,alpha,u,der,C)
implicit none
integer:: i,j
integer:: nx,KL,KR,alpha
real(Kind=16):: u(1+KL:nx+KR),der(1:nx),C(-alpha:alpha)
	do i=1,nx
		der(i)=sum(C(-alpha:alpha)*u(i-alpha:i+alpha))
	enddo
end
!----------------------------------------------------------------------------------
!Calculate discrete Taylor expansion coefficients
subroutine cal_cof_tlc(n,P,tlc)
implicit none
Real(Kind=16):: fact
integer :: n,P,ii
!Large number calculations must use double precision
integer(kind=8):: i,j
real(Kind=16):: T(0:2*P,-P:P),IM(0:2*P)
!Temporary classical finite difference coefficients
real(Kind=16):: tmp_cls(-P:P,1:P)
real(Kind=16):: tlc(0:2*P)

	do i=0,2*P 
	do j=-P,P 
		!The i-th order Taylor expansion coefficient of u(0) at point j
		T(i,j)=j**i/(1.Q0*fact(i))	
	enddo
	enddo

	IM(:)=0.q0
	IM(n)=1.q0

	tmp_cls(-P:P,1:P)=0.q0
	tlc(0:2*P)=0.q0
	tlc(n)=1.q0
	do i=(n+mod(n,2))/2,P-1 
		call Guass_solver(2*i+1,1,T(0:2*i,-i:i),IM(0:2*i),tmp_cls(-i:i,i))
		tlc(2*i+2-mod(n,2))=dot_product(T(2*i+2-mod(n,2),-i:i),tmp_cls(-i:i,i))
	enddo

	!~ 	  print*, "tlc"
	!~ 	  do i=0,2*(P-1)
	!~ 	      Write(*,"(30f20.15)") tlc(i)
	!~ 	  end do
	!~ stop

end
!-------------------------------------------------------
!  Gaussian Elimination Solver
subroutine Guass_solver(n,m,A_in,B_in,x)
implicit none
integer i,j,k,n,m
!增广矩阵和解向量
real(Kind=16):: A_in(n,n),B_in(n,m),Arr(n,n+m),x(n,m)	
real(Kind=16):: a,b(m),error
	Arr(1:n,1:n)=A_in(1:n,1:n)
	Arr(1:n,n+1:n+m)=B_in(1:n,1:m)

	do k=1,n-1 
		do i=k+1,n 
			a=Arr(i,k)/Arr(k,k)
			do j=k,n+m 
				Arr(i,j)=Arr(i,j)-a*Arr(k,j)
			enddo
		enddo
	enddo

	x(n,1:m)=Arr(n,n+1:n+m)/Arr(n,n)

	do i=n-1,1,-1 
		b(:)=0
		do j=n,i+1,-1 
			b(:)=b(:)+Arr(i,j)*x(j,:)
		enddo
		x(i,1:m)=(Arr(i,n+1:n+m)-b(1:m))/Arr(i,i)
	enddo

	error=maxval(abs(matmul(A_in(1:n,1:n),x(1:n,1:m))-B_in(1:n,1:m)))
	if(error .gt.1.q-20)then
		print*, "hahaha",error
		stop
	endif

end
!----------------------------------------------------------------------------------
!factorial function
Real(Kind=16) function fact(n)
integer(kind=8), intent(in) :: n
integer(kind=8) :: i
integer(Kind=8) :: f

	if (n < 0) error stop 'factorial is singular for negative integers'
	f = 1
	do i = 2, n
		f = f * i
	enddo
	fact=Real(f,kind=16)
end function fact
