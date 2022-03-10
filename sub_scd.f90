!-------------------------------------------------------
! Calculate the derivative of odd order
subroutine Cal_der_odd_scd
use var
implicit none
integer:: i,j,k
	do j=2*(alpha+1-dof)-1,2*alpha-1,2
		do i=1,nx
!~ 		u_itg(i,j)=0.d0
!~ 			do k=1,(j+1)/2
!~ 				u_itg(i,j)=u_itg(i,j)+der_scd(k,j)*(u(i+k)-u(i-k))
!~ 			enddo
			u_itg(i,j)=sum(der_scd(-(j+1)/2:(j+1)/2,j)*u(i-(j+1)/2:i+(j+1)/2))
		enddo
	enddo
end
!-------------------------------------------------------
! Calculate the derivative of even order
subroutine Cal_der_even_scd
use var
implicit none
integer:: i,j,k
	do j=2*(alpha+1-dof),2*alpha,2
		do i=1,nx
!~ 		u_itg(i,j)=der_scd(0,j)*u(i)
!~ 			do k=1,j/2
!~ 				u_itg(i,j)=u_itg(i,j)+der_scd(k,j)*(u(i+k)+u(i-k))
!~ 			enddo
			u_itg(i,j)=sum(der_scd(-j/2:j/2,j)*u(i-j/2:i+j/2))
		enddo
	enddo
end
!-------------------------------------------------------
! Calculate the derivative of odd order
subroutine Cal_sta_odd_scd
use var
implicit none
integer:: i,j,k,it,jt
	do i=1,dof
		it=2*(alpha+i-dof)-1
		do j=1,dof
			jt=2*(alpha+j-dof)-1
			A_matrix(i,j)=sum(u_itg(1:nx,it)*u_itg(1:nx,jt))
		enddo
		B_vector(i)=sum(u_itg(1:nx,it)*u_trunc(1:nx))
	enddo
end
!-------------------------------------------------------
! Calculate the derivative of even order
subroutine Cal_sta_even_scd
use var
implicit none
integer:: i,j,k,it,jt
	do i=1,dof
		it=2*(alpha+i-dof)
		do j=1,dof
			jt=2*(alpha+j-dof)
			A_matrix(i,j)=sum(u_itg(1:nx,it)*u_itg(1:nx,jt))
		enddo
		B_vector(i)=sum(u_itg(1:nx,it)*u_trunc(1:nx))
	enddo
end
