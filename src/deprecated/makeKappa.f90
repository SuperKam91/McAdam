program makeKappa
      
      use fits_utils
      use RandomNS

	integer i,j,k
      double precision kappa(256,256)
      character*100 file
      
     	open(unit=23,file='/local/cosmos-med/tws29/McAdam_v3/data/GL/sim/090707/gcat13_kappa.cat',status='unknown')
            
	do i=1,256
      	do j=1,256
      		read(23,*),k,k,kappa(i,j)
		enddo
	enddo
      
      close(23)
      
      call initRandomNS(1)
      
      do i=1,256
      	do j=1,256
            	kappa(i,j)=kappa(i,j)+0.256*Gaussian1NS(0)
		enddo
	enddo
      
      file='gcat13_n.fits'
      call writefits(kappa,256,256,7.03125d0,file)
      
end program makeKappa
      
      
      
      
