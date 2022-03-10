!-------------------------------------------------------
! Calculate the derivative of odd order,
! using "Building Blocks" acceleration
subroutine Cal_der_odd_BB
use var
implicit none
integer:: i,j
			u_hlf(1,0)=0.5d0*(u_itg(1,0)+u_itg(nx,0))
		do i=2,nx
			u_hlf(i,0)=0.5d0*(u_itg(i,0)+u_itg(i-1,0))
		enddo
		do i=1,nx-1
			u_itg(i,1)=u_hlf(i+1,0)-u_hlf(i,0)
		enddo
			u_itg(nx,1)=u_hlf(1,0)-u_hlf(nx,0)
	do j=2,P-2,2
			u_hlf(1,j)=u_itg(1,j-1)-u_itg(nx,j-1)
		do i=2,nx
			u_hlf(i,j)=u_itg(i,j-1)-u_itg(i-1,j-1)
		enddo
		do i=1,nx-1
			u_itg(i,j+1)=u_hlf(i+1,j)-u_hlf(i,j)
		enddo
			u_itg(nx,j+1)=u_hlf(1,j)-u_hlf(nx,j)
	enddo
end
!-------------------------------------------------------
! Calculate the derivative of even order,
! using "Building Blocks" acceleration
subroutine Cal_der_even_BB
use var
implicit none
integer:: i,j
	do j=1,P-1,2
			u_hlf(1,j)=u_itg(1,j-1)-u_itg(nx,j-1)
		do i=2,nx
			u_hlf(i,j)=u_itg(i,j-1)-u_itg(i-1,j-1)
		enddo
		do i=1,nx-1
			u_itg(i,j+1)=u_hlf(i+1,j)-u_hlf(i,j)
		enddo
			u_itg(nx,j+1)=u_hlf(1,j)-u_hlf(nx,j)
	enddo
end
!-------------------------------------------------------
! Calculate the statistics of odd order,
! using "Building Blocks" acceleration
subroutine Cal_sta_odd_BB
use var
implicit none
integer:: i,j,k,it,jt
real(kind=16):: diagonal(2*(alpha+1-dof)-1:2*(alpha+tau)-1)
	do i=2*(alpha+1-dof)-1,2*(alpha+tau)-3,2 !1，3，5，7，，，
		diagonal(i)=sum(u_itg(1:nx,i)*u_itg(1:nx,i))
		diagonal(i+1)=-sum(u_hlf(1:nx,i+1)*u_hlf(1:nx,i+1))
	enddo		
		i=2*(alpha+tau)-1 
		diagonal(i)=sum(u_itg(1:nx,i)*u_itg(1:nx,i))
	do i=1,dof
		it=2*(alpha+i-dof)-1
		do j=1,dof
			jt=2*(alpha+j-dof)-1
			k=(it+jt)/2
			A_matrix(i,j)=diagonal(k)
		enddo
		B_vector(i)=0.q0
		do j=1,2*tau
			jt=2*(alpha+j)-1
			k=(it+jt)/2
			B_vector(i)=B_vector(i)-tlc_cls(jt)*diagonal(k)
		enddo
	enddo
end
!-------------------------------------------------------
! Calculate the statistics of even order,
! using "Building Blocks" acceleration
subroutine Cal_sta_even_BB
use var
implicit none
integer:: i,j,k,it,jt
real(kind=16):: diagonal(2*(alpha+1-dof):2*(alpha+tau))
	do i=2*(alpha+1-dof),2*(alpha+tau)-2,2 !2，4，6，8，，，
		diagonal(i)=sum(u_itg(1:nx,i)*u_itg(1:nx,i))
		diagonal(i+1)=-sum(u_hlf(1:nx,i+1)*u_hlf(1:nx,i+1))
	enddo		
		i=2*(alpha+tau) 
		diagonal(i)=sum(u_itg(1:nx,i)*u_itg(1:nx,i))
	do i=1,dof
		it=2*(alpha+i-dof)
		do j=1,dof
			jt=2*(alpha+j-dof)
			k=(it+jt)/2
			A_matrix(i,j)=diagonal(k)
		enddo
		B_vector(i)=0.q0
		do j=1,2*tau
			jt=2*(alpha+j)
			k=(it+jt)/2
			B_vector(i)=B_vector(i)-tlc_cls(jt)*diagonal(k)
		enddo
	enddo
end
!-------------------------------------------------------
! Approximate vector B using "Building Blocks" acceleration
! Available only when n is odd
subroutine Cal_approx_BB
use var
implicit none
integer:: i,j,k,it,jt
real*16:: tmp
	do i=1,dof
		it=2*(alpha+i-dof)-1
		jt=2*alpha+4*tau+1
		k=(it+jt)/2
		if(mod(k-it,2)==1)then
			tmp=-sum(u_itg(1:nx,k)*u_itg(1:nx,k))
		else
			tmp=sum(u_hlf(1:nx,k)*u_hlf(1:nx,k))		
		endif
		B_vector(i)=B_vector(i)-tlc_cls(jt)*tmp
	enddo
end
