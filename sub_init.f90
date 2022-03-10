!---------------------------------------------------------------------------------
subroutine init
use var
implicit none
integer:: i

	n=2; KL=-10; KR=10
	nx=2**8
	Building_Blocks_acceleration=1; Building_Blocks_approximation=0
	alpha=3; tau=4; P=2*(alpha+tau); dof=3

	dx=1.q0/nx

	allocate(x(1+KL:nx+KR),u(1+KL:nx+KR),der_u_n(1:nx),exact(1+KL:nx+KR),&
			 u_itg(1:nx,0:P),u_hlf(1:nx,0:P-1),&
			 u_trunc(1:nx),A_matrix(dof,dof),B_vector(dof),e_vector(dof))
	   
	do i=1,nx
	 x(i)=(i-1.q0/2.q0)*dx
	enddo
	   
	u(:)=0.d0
	u_itg(:,:)=0.q0

	do  i=1+KL,nx+KR    
		u(i)=exp(-300.q0*(x(i)-0.5q0)**2)
		select case(n)
		case(1)
			exact(i)=(-600*(-0.5 + x(i)))/Exp(300*(-0.5 + x(i))**2)
		case(2)
			exact(i)=(600*(149 + 600*(-1 + x(i))*x(i)))/Exp(75*(1 - 2*x(i))**2)
		case default
			exact(i)=(-600*(-0.5 + x(i)))/Exp(300*(-0.5 + x(i))**2)
		end select
	enddo 

	u_itg(1:nx,0)=u(1:nx)

end subroutine init
