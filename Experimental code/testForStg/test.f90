program test
	integer::k1,k2,counter=0;
	integer::nx=8
	do k1=0,nx
		do k2=0,nx
			if (k1==0 .and. k2==nx/2) then
				write(*,*)"0,pi ",counter
			else if (k1==nx/2 .and. k2==nx/2) then
				write(*,*)"pi,pi ",counter
			else if (k1==nx/2 .and. k2==0) then
				write(*,*)"pi,0 ",counter
			else if (k1==nx .and. k2==nx) then
				write(*,*)"2pi, 2pi", counter
			else if (k1==nx .and. k2==0) then
				write(*,*)"2pi, 0", counter
			else if (k1==0 .and. k2==nx) then
				write(*,*)"0, 2pi", counter
			endif
			counter=counter+1;
		enddo
	enddo
	stop
end program test