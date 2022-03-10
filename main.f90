!------------------------------------------------------------------------------
! This code is an implementation of 
! "On Finite Difference Schemes Approximated minimized L2 norm error 
! and "Building Blocks" acceleration"
! by Wei HuangSheng, Huang Zhu and Xi Guang.
!==============================================================================
module var    
	implicit none
!------------------------------------------------------------------------------
! Global
	integer:: n !optimize and compute the nth derivative
	integer:: nx
	integer:: KL,KR
    real(Kind=16):: dx
	real*16,dimension(:),Pointer:: x,u,der_u_n,exact
!------------------------------------------------------------------------------
! Finite difference scheme approximated minimized L2 norm error
	integer:: alpha,tau,P,dof
	integer:: Building_Blocks_acceleration	!0：no; 1:yes
	integer:: Building_Blocks_approximation  !For odd derivatives, 0：no; 1:yes
	real*16,dimension(:,:),Pointer:: u_itg,u_hlf
	real*16,dimension(:,:),Pointer:: A_matrix
	real*16,dimension(:),Pointer:: B_vector,e_vector
	real*16,dimension(:),Pointer:: C_cls,C_opt
!------------------------------------------------------------------------------
! Use Building Blocks
	real*16,dimension(:),Pointer:: tlc_cls
!------------------------------------------------------------------------------
! If not use Building Blocks, then need
	real*16,dimension(:,:),Pointer:: der_scd
	real*16,dimension(:),Pointer:: C_trunc
	real*16,dimension(:),Pointer:: u_trunc
end module var
!------------------------------------------------------------------------------
program main
use var
implicit none
integer:: i,j,k,times
real*16:: time_start,time_end,cost
real*16:: e_L2_cls,e_L2_opt
real*16:: opt_cost

	times=2**16
	opt_cost=0.q0

	print*, "Initializing parameters"
		call init
	print*, "Initializing classicl and standard second order accurate schemes"
		call init_scheme
	print*, "End initialization"
	print*, "****************************************************************"
	print*, "Whether to use Building_Blocks_acceleration", Building_Blocks_acceleration
	print*, "Whether to use Building_Blocks_approximation", Building_Blocks_approximation
	print*, "n=",n,"alpha=",alpha,"tau=",tau,"dof=",dof
	print*, "nx=",nx,"The number of repeated calculations is",times
	print*, "****************************************************************"
	print*, "Calculating derivatives"
		call cpu_time(time_start)
			do i=1,times
				call Cal_derivatives
			enddo
		call cpu_time(time_end)
		cost=time_end-time_start
		opt_cost=opt_cost+cost
	print*, "End Calculation"
	print*, "Cal_derivatives cost:", cost
	print*, "****************************************************************"
	print*, "Calculating statistics"
		call cpu_time(time_start)
			do i=1,times
				call Cal_statistics
			enddo
		call cpu_time(time_end)
		cost=time_end-time_start
		opt_cost=opt_cost+cost
	print*, "End Calculation"
	print*, "Cal_statistics cost:", cost
	print*, "****************************************************************"
	print*, "Optimizing"
		call cpu_time(time_start)
			do i=1,times
				call Cal_optimization
			enddo
		call cpu_time(time_end)
		cost=time_end-time_start
		opt_cost=opt_cost+cost
	print*, "End Optimization"
	print*, "Cal_optimization cost:", cost
	print*, "****************************************************************"
	print*, "calculating derivatives using classical scheme"
		call cpu_time(time_start)
			do i=1,times
				call Cal_der(nx,KL,KR,alpha,u(1+KL:nx+KR),der_u_n(1:nx),C_cls(-alpha:alpha))
			enddo
			der_u_n(1:nx)=der_u_n(1:nx)/(dx**n)
		call cpu_time(time_end)
		cost=time_end-time_start
		e_L2_cls=sqrt(sum( (der_u_n(1:nx)-exact(1:nx))**2 )/sum(exact(1:nx)**2) )
	print*, "End Calculation"
	print*, "Cal_der cost:", cost
	print*, "****************************************************************"
	print*, "calculating derivatives using optimized scheme"
		call cpu_time(time_start)
			do i=1,times
				call Cal_der(nx,KL,KR,alpha,u(1+KL:nx+KR),der_u_n(1:nx),C_opt(-alpha:alpha))
			enddo
			der_u_n(1:nx)=der_u_n(1:nx)/(dx**n)
		call cpu_time(time_end)
		cost=time_end-time_start
		e_L2_opt=sqrt(sum( (der_u_n(1:nx)-exact(1:nx))**2 )/sum(exact(1:nx)**2) )
	print*, "End Calculation"
	print*, "Cal_der cost:", cost
	print*, "****************************************************************"
	print*, "A_matrix ="
		do i=1,dof
			Write(*,"(20f20.15)") A_matrix(i,1:dof)
		enddo
	print*, "B_vector ="
			Write(*,"(20f20.15)") B_vector(1:dof)
	print*, "e_vector ="
			Write(*,"(20f20.15)") e_vector(1:dof)
	print*, "****************************************************************"
	print*, "Classical finite difference scheme is"
			Write(*,"(20f20.15)") C_cls(-alpha:alpha)
	print*, "Its relative error is"
	print*, e_L2_cls	
	print*, "****************************************************************"
	print*, "Optimized finite difference scheme is"
			Write(*,"(20f20.15)") C_opt(-alpha:alpha)
	print*, "Its relative error is"
	print*, e_L2_opt	
	print*, "****************************************************************"
	print*, "Calculation error decreased by", (e_L2_cls-e_L2_opt)/e_L2_cls
	print*, "Calculations increased by", opt_cost/cost
	print*, "****************************************************************"
end

!   -----------------------
include "sub_init.f90"
include "sub_scheme.f90"
include "sub_buildingblocks.f90"
include "sub_scd.f90"
include "sub_opt.f90"



