module rfft
	use params
      use constants
      use telescope1

contains

	subroutine fourt(data,nn,ndim,isign,iform,work)

! Wrapper for the SCSL library fft (2D,single precision,complex)
! with a 'fourt' interface.
! The resulting data array is the same as fourt would return (same scale..)
! NB - the SCSL routine is more accurate (smaller remainder after f+inv)
! NB2 - would be more efficient to allocate table automatically 

  	double complex  data(*)
  	double precision work(*)
  	double precision scale
  	double precision,allocatable :: table(:)
  	integer nn(2)
  	integer ndim,isign,iform,init,n
  	integer isys(2)

  	n=max(nn(1),nn(2))
  	allocate(table(2*nn(1)+2*nn(2)+512))
  	scale=1.0
  	isys(1)=1
  	init=0
	
!	initialize table
!  	call ZZFFT2D(init,nn(1),nn(2),scale,data,nn(1),data,nn(1),table,work,isys)
	
!	compute
!  	call ZZFFT2D(isign,nn(1),nn(2),scale,data,nn(1),data,nn(1),table,work,isys)
	
  
	end subroutine fourt


end module rfft
