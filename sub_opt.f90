!---------------------------------------------------------------------------------
subroutine cal_derivatives
use var
implicit none

if (Building_Blocks_acceleration==1)then
	if(mod(n,2)==1)then
		call Cal_der_odd_BB
		if (Building_Blocks_approximation==1)then
			call Cal_der_even_BB
		endif
	else
		call Cal_der_even_BB
	endif
else
	if(mod(n,2)==1)then
		call Cal_der_odd_scd
	else
		call Cal_der_even_scd
	endif
	call Cal_der(nx,KL,KR,alpha+tau,u(1+KL:nx+KR),u_trunc(1:nx),C_trunc(-alpha-tau:alpha+tau))
endif

end
!---------------------------------------------------------------------------------
subroutine cal_statistics
use var
implicit none

if (Building_Blocks_acceleration==1)then
	if(mod(n,2)==1)then
		call Cal_sta_odd_BB
		if (Building_Blocks_approximation==1)then
			call Cal_approx_BB
		endif
	else
		call Cal_sta_even_BB
	endif
else
	if(mod(n,2)==1)then
		call Cal_sta_odd_scd
	else
		call Cal_sta_even_scd
	endif
endif

end
!---------------------------------------------------------------------------------
subroutine cal_optimization
use var
implicit none
integer:: i,it

	call Guass_solver(size(e_vector),1,A_matrix,B_vector,e_vector) 
	C_opt(-alpha:alpha)=C_cls(-alpha:alpha)
	do i=1,dof
		it=2*(alpha+i-dof)-mod(n,2)
		C_opt(-alpha:alpha)=C_opt(-alpha:alpha)+e_vector(i)*der_scd(-alpha:alpha,it)
	enddo

end
