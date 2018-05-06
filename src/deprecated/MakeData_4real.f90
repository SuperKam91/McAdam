module MakeData1
	use params
      use constants
      use telescope1
      use params_MakeData
      use lensing
      use ReadWrite
      use utilities
      use RandomNS
      use GasModels
      use MassModels
      use fits_utils
      use massfunction
      use CheckPars1

contains

!=======================================================================

	subroutine MakeData
	
	implicit none
	
	double precision obsarea
	double precision ehat(Ngals),pphihat(Ngals)
      double precision e1hat(Ngals),e2hat(Ngals)
	double precision e(Ngals),pphi(Ngals),mode
      integer eflag(Ngals)
	double precision noise1(Ngals),noise2(Ngals)
	double precision ehats(nxpix,nypix),pphihats(nxpix,nypix)
      double precision e1hats(nxpix,nypix),e2hats(nxpix,nypix)
	double precision es(nxpix,nypix),pphis(nxpix,nypix)
      integer eflags(nxpix,nypix)
	double precision noise1s(nxpix,nypix),noise2s(nxpix,nypix)
      double precision junk(Ngals),z_s(Ngals)

	integer seed,ndata,nx2,ny2,idummy

	double precision xmin,xmax,ymin,ymax,pixel
	double precision massmap(nxpix,nypix)

	double precision dummy,dummy2
	integer i,j,k,m,flag,ic,nbad,mAMI,mVSA,mRT,mSEA,mCBI
	character*100 string,catfile
	character*100 visfile(Nvisfiles),visfile2(Nvisfiles)
	character*1 char
	
	double complex vis_buffer
	double precision weight_buffer,uv_buffer(2),rms_buffer
	integer baseline(large)
	double precision ra,dec,ra2,dec2
	double precision unitnoise(2*large),scalednoise(2*large)
	double precision row(2*NLCM),sum
	double complex vis_dump(large)
	double precision weight_dump(large),uv_dump(2,large)
	double precision szsky(nx,ny),urv
      
      double precision area,error,d1,xrange,yrange
      integer i1,j1

!-----------------------------------------------------------------------

! Explanation:

	write(*,'(a)')' In simulation mode - generate mock data sets, then stop.'
	write(*,'(a)')

      if(ModelClass==2)then
      	write(*,'(a,a,/,a,a)')' The data will be created using ',Sobsfile,' and ',covmatfile
      	write(*,*)
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! Gravitational Lensing data
! **************************

	if (GL==1.and.ModelClass==1) then

		write(*,'(a)') ' Gravitational lensing data:'
		write(*,'(a)') ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^'

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! Set up output filename:

		string = GLroot
		do ic=1,100
        		read(string(ic:ic),'(a1)') char
        		if (char==' ') goto 5
      	enddo
 5    	write(catfile(1:ic-1),'(a)') string(1:ic-1)
		do j=ic,100
        		write(catfile(j:j),'(a1)') ' '
      	enddo
	      
		if (nreals==1) then
        		write(catfile(ic:ic+3),'(a4)') '.cat'
		else
        		write(catfile(ic:ic+7),'(a8)') '_000.cat'
      	endif
	 
	
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! For each noise realisation,
            
		do k=1,Nreals

! Initialise Random Number generator:
	  
			seed = 13*k
			call initRandomNS(1)
	!i = RandomInit(state,seed)

! Either use existing catalogue for galaxy positions and redshifts, or
! generate random galaxy positions within observing area. 

			if (GL==1 .and. readCatalog) then
        			call readcat(Ngals,x,y,junk,junk,junk,junk,junk,i1,GLobsfile)
        
        			if (k==1) then
          				write(*,'(a,i5,a)') ' Observed ',Ngals,' galaxy images.'
        			endif
        			ndata = Ngals
                        
                        !mass map setup
                        xmin=minval(x)
                        xmax=maxval(x)
                        ymin=minval(y)
                        ymax=maxval(y)
                        !range
                        xrange=xmax-xmin
                        yrange=ymax-ymin
                        
      		else if(.not.readCatalog) then
        			if(survey) then
                        	ndata=nxpix*nypix
				else
                        	!calculate the total no. of galaxies
                        	if(k==1) then
                  			ndata=ndensity*xrange_s*yrange_s
                        		if (ndata>Ngals) then
            					write(*,*) 'Error: in MakeData.inc Ngals = ',ndata
            					write(*,*) '       in McAdam.inc Ngals = ',Ngals
            					write(*,*) ' --> Change McAdam.inc to avoid overflow.'
            					write(*,*)
            					stop
          					else if (ndata<Ngals) then
            					write(*,*) 'Warning: in MakeData.f Ngals = ',ndata
            					write(*,*) '         in McAdam.inc Ngals = ',Ngals
            					write(*,*)
          					endif
          					write(*,'(a,i4,a,i5,a)') ' Observed ',nint(xrange_s*yrange_s), &
          					' square arcmin of sky (',ndata,' galaxy images)'
          					write(*,'(a,f4.2,a,f8.1)') ' Mean background galaxy redshift = ', &
          					zsource
					end if
                        	!generate random galaxy positions
                        	do i=1,ndata
                  			urv=ranmarNS(0)
                              	x(i)=(-xrange_s/2.+xrange_s*urv)*60.
                  			urv=ranmarNS(0)
                              	y(i)=(-yrange_s/2.+yrange_s*urv)*60.
					end do
                              
                              !mass map setup
                              xmin=-xrange/2
                              xmax=xrange/2
                              ymin=-yrange/2
                              ymax=yrange/2
                              !ranges
                              xrange=xrange_s
                              yrange=yrange_s
				endif
			endif
                  	
			nbad=0
            	if(.not.survey) then
				! Allocate intrinsic shapes and critical density
      			do i=1,ndata
         				e1hat(i)=Gaussian1NS(0)*sigmaint
         				e2hat(i)=Gaussian1NS(0)*sigmaint
	  				
                              ehat(i) = sqrt(e1hat(i)*e1hat(i)+e2hat(i)*e2hat(i))
	  				pphihat(i) = 0.5*atan2(e2hat(i),e1hat(i))
	 				eflag(i) = 0
	   
!      				deal with bad intrinsic galaxies (hack)
	 				if (ehat(i)>=1.) then
	   					eflag(i) = 1
	   					nbad = nbad + 1
	 				endif
	  
!      				zero variables	  
	 				kappa(i) = 0.
	 				gamma1(i) = 0.
	 				gamma2(i) = 0.
      			end do
 			else
				! Allocate intrinsic shapes and critical density 
      			do i=1,nxpix
                        	do j=1,nypix
         					e1hats(i,j)=Gaussian1NS(0)*sigmaint
         					e2hats(i,j)=Gaussian1NS(0)*sigmaint
	  					
                                    ehats(i,j)=sqrt(e1hats(i,j)**2.+e2hats(i,j)**2.)
	  					pphihats(i,j)=0.5*atan2(e2hats(i,j),e1hats(i,j))
	 				
                              	eflags(i,j)=0
	   
!      					deal with bad intrinsic galaxies (hack)
	 					if(ehats(i,j)>=1.) then
	   						eflags(i,j)=1
	   						nbad=nbad+1
	 					endif
	  
!      					zero variables	  
	 					kappas(i,j)=0.
	 					gamma1s(i,j)=0.
	 					gamma2s(i,j)=0.
					enddo
      			end do
			endif
                  
                  !open the file in which the true parameters will be written
                  open(unit=23,file=trim(n_root)//'true.param',status='unknown')
                  
                  !allocate the memory for mass function lookup table
                  if(genTrue .and. MassModel==1 .and. NFWstyle(1)==2 .and.  &
                  (z_PriorType==8 .and. Mass_PriorType(1,2)==8)) then
                  	allocate(lookM(n,2),lookZ(n,n,2))
            		call makeMZlookup
			endif
                  
			do j=1,natoms

!       			Copy parameters:
				
                        !sample the true parameters from the prior?
				if(genTrue) then
                        	NFWStyle(j)=NFWStyle(1)
                              !rescale the actual parameters
                        	do m=1,NPars+1
                              	urv=ranmarNS(0)
                                    if(m<=NPars) then
                              		call rescale(urv,m)
						else
                                    	call rescale(urv,NAtoms*NPars+1)
						endif
					enddo
                              GeoPars(j,1:NGeoPars)=GeoPars(1,1:NGeoPars)
      				MassPars(j,1:NMassPars)=MassPars(1,1:NMassPars)
				else   
	  				GeoPars(j,1) = x0_s(j)
	  				GeoPars(j,2) = y0_s(j)
	  
!       				Sort out elliptical geometry:
        				if (GeoModel==2) then
	    					GeoPars(j,3) = theta_s(j)
	    					GeoPars(j,4) = axratio_s(j)
	  				endif
	    
	  				if (MassModel==1) then
	    					if (NFWstyle(j)==1) then
	      					MassPars(j,1) = rs_s(j)
	      					MassPars(j,2) = ps_s(j)
	    					elseif (NFWstyle(j)==2) then
	      					MassPars(j,1) = c200_s(j) 
	      					MassPars(j,2) = M200_s(j)
	    					elseif (NFWstyle(j)==3) then
	      					MassPars(j,1) = M200_s(j) 
	    					elseif (NFWstyle(j)==4) then
	      					MassPars(j,1) = thetaE_s(j) 
	      					MassPars(j,2) = M200_s(j) 
	    					endif  
	  				elseif (MassModel==2) then 
	    					MassPars(j,1) = M200_s(j)
	  				endif
                              z=zcluster_s(j)
				endif
                        if (GeoModel==2) call EllGeometry(j)
                        
                        !now write the true parameters to the file
				write(23,*)GeoPars(j,1:NGeoPars)
      			if(Mass) write(23,*)MassPars(j,1:NMassPars)
				if(varyzs==1) then
                             	write(23,*)z,zs
				else
                             	write(23,*)z
				endif
                        write(23,*)
                        write(23,*)
                        
                        D=calcAngDiamDis(0.d0,z)
                        if(varyzs==0) zs=zsource
                        zsmin=zsource
                        zsmax=zsource
                        z_s=zsource
!       			Now increase the shear and convergence at each galaxy position:
	  			call LensFields(j,flag,.false.)
        			if(flag==1) stop 'Flag = 1 in LensFields'
			enddo
                  close(23)
                  if(genTrue .and. MassModel==1 .and. NFWstyle(1)==2 .and.  &
                  (z_PriorType==8 .and. Mass_PriorType(1,2)==8)) then
                  	deallocate(lookM,lookZ)
			endif

			if(survey) then
                  	!work out the gaussian noise per pixel
                        
                        !fist calculate the no. of galaxies per pixel
                        !calculating the are of each pixel in square arcmin
                        area=pix_sizex*pix_sizey/3600.
                        d1=area*ndensity
                        error=sqrt(sigmaint**2.+sigmaobs**2.)/sqrt(d1)
                        
                        
! 				Now construct reduced shear, transform galaxy elipticities and deal 
! 				with critical regions
				do i=1,nxpix
					do j=1,nypix
                              	!deal with bad intrinsic galaxies (hack)
        					if (eflags(i,j)==1) then
          						e1s(i,j)=0.
          						e1errs(i,j)=0.
          						e2s(i,j)=0.
          						e2errs(i,j)=0.
                                          g1s(i,j)=0.
                                          g2s(i,j)=0.
                                  	  	cycle
        					endif
	  
	 					g1s(i,j) = gamma1s(i,j)/(1.-kappas(i,j))	   
	 					g2s(i,j) = gamma2s(i,j)/(1.-kappas(i,j))

		 				call lensdble(e1hats(i,j),e2hats(i,j),g1s(i,j),g2s(i,j),e1s(i,j),e2s(i,j))
	 
	      				!alter data inside strong lensed region

	 					mode=e1s(i,j)*e1s(i,j)+e2s(i,j)*e2s(i,j)
       					if (mode>=1.) then
	   						e1s(i,j) = e1s(i,j)/mode
	   						e2s(i,j) = e2s(i,j)/mode
	 					endif
	 
	 					es(i,j) = sqrt(e1s(i,j)*e1s(i,j)+e2s(i,j)*e2s(i,j))
	 					pphis(i,j) = 0.5*atan2(e2s(i,j),e1s(i,j))
                              
                              	!Add experimental noise:
                              	noise1s(i,j)=Gaussian1NS(0)*sigmaobs
         					noise2s(i,j)=Gaussian1NS(0)*sigmaobs
	  
        					e1s(i,j)=e1s(i,j)+noise1s(i,j)
        					e2s(i,j)=e2s(i,j)+noise2s(i,j)
	  
	  					es(i,j)=sqrt(e1s(i,j)*e1s(i,j)+e2s(i,j)*e2s(i,j))
	  					pphis(i,j)=0.5*atan2(e2s(i,j),e1s(i,j))

	   					e1errs(i,j)=error
	   					e2errs(i,j)=error
					enddo
				enddo
                  
                  else
! 				Now construct reduced shear, transform galaxy elipticities and deal 
! 				with critical regions
				dummy = sqrt(sigmaint*sigmaint + sigmaobs*sigmaobs)
				do i=1,Ndata
                              
                              !deal with bad intrinsic galaxies (hack)
        				if (eflag(i)==1) then
          					e1(i) = 0.
          					e1err(i) = 0.
          					e2(i) = 0.
          					e2err(i) = 0.
						g1(i)=0.
						g2(i)=0.
                                    cycle
        				endif
	  
	 				g1(i) = gamma1(i)/(1.-kappa(i))	   
	 				g2(i) = gamma2(i)/(1.-kappa(i))

		 			call lensdble(e1hat(i),e2hat(i),g1(i),g2(i),e1(i),e2(i))
	 
!      				alter data inside strong lensed region

	 				mode = e1(i)*e1(i) + e2(i)*e2(i)
       				if (mode >= 1.0) then
	   					e1(i) = e1(i)/mode
	   					e2(i) = e2(i)/mode
	 				endif
	 
	 				e(i) = sqrt(e1(i)*e1(i)+e2(i)*e2(i))
	 				pphi(i) = 0.5*atan2(e2(i),e1(i))
                              
                              !Add experimental noise:
                              noise1(i)=Gaussian1NS(0)*sigmaobs
         				noise2(i)=Gaussian1NS(0)*sigmaobs
	  
        				e1(i) = e1(i) + noise1(i)
        				e2(i) = e2(i) + noise2(i)
	  
	  				e(i) = sqrt(e1(i)*e1(i)+e2(i)*e2(i))
	  				pphi(i) = 0.5*atan2(e2(i),e1(i))

	   				e1err(i) = dummy
	   				e2err(i) = dummy
				enddo
			endif

			if(nreals==1) write(*,*) '  No. of galaxies in catalogue = ',ndata 
			if(nreals==1) write(*,*) '  (Including ',nbad,' bad images)' 
			if(nreals==1) write(*,*)  

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! 			Write out new catalogue(s)
! 	
			if (nreals==1) then
	  
        			write(*,*) '  Writing data to file: ',catfile(1:ic+7)
        			write(*,*) '   No. of galaxies in catalogue = ',ndata
                        
                        if(survey) then
                        	call writecat_survey(nxpix,nypix,e1s,e2s,error,zsource,catfile)
                        else
        				write(*,*) '   (Including ',nbad,' bad images)'
        				write(*,*)
                              
        				call writecat(ndata,x,y,e1,e2,e1err,e2err,z_s,catfile)	
				endif
	
			else
	  
        			if (k<10) then
	    				write(catfile(ic+3:ic+3),'(i1)') k  
	  			elseif (k<100) then
	    				write(catfile(ic+2:ic+3),'(i2)') k  
	  			else 
	    				write(catfile(ic+1:ic+3),'(i3)') k  
        			endif 
	   
        			write(*,*) '  Writing data to file: ',catfile(1:ic+7)
       	 		write(*,*) '   No. of galaxies in catalogue = ',ndata
        			write(*,*) '   (Including ',nbad,' bad images)'
        			write(*,*)
	  			call writecat(ndata,x,y,e1,e2,e1err,e2err,z_s,catfile)	
	
! 	  			if (mod(k,10)==0) write(*,*) k 
	    	
			endif
		
		enddo

		write(catfile(ic:ic+5),'(a6)') '_g.cat'
		do j=ic+10,100
      		write(catfile(j:j),'(a1)') ' '
      	enddo
      	write(*,*) '  Writing true g to file: ',catfile(1:ic+9)
            if(survey) then
            	call writecat_survey(nxpix,nypix,g1s,g2s,error,zsource,catfile)
                  
                  write(catfile(ic:ic+9),'(a10)') '_kappa.cat'
      		write(*,*) '  Writing true map to file: ',catfile(1:ic+10) 
            	call writekappa_survey(nxpix,nypix,kappas,catfile)
                  
			!Write out in fits format
                  if(nxpix==nypix .and. pix_sizex==pix_sizey) then
	      		write(catfile(ic:ic+10),'(a11)') '_kappa.fits'
      			write(*,*) '  Writing true map to file: ',catfile(1:ic+10) 
				call writefits(kappas,nxpix,nypix,pix_sizex,catfile)
      			write(*,*)
			else
                  	write(*,*)"can not generate the fits file for true kappa since &
                        either the no. of pixels or the pixel sizes in x & y are not the &
                        same"
			endif
		endif
		write(*,*)  
	
! 		Finished creating simulated data files.  

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		
            if(.not.survey) then
! 			Generate true mass map:

			if (nxpix*nypix>Ngals) then
	  			write(*,*) '  Warning: No. of pixels = ',nxpix*nypix
	  			write(*,*) '           Length of x,y arrays = ',Ngals
	  			write(*,*) '    --> Cannot make true map.'
	  			goto 40
			endif

! 			pixel size:
			pixel=max(xrange/nxpix,yrange/nypix)

! 			Redefine x and y arrays to be the centres of the pixels:

			do j=1,nypix
	  			do i=1,nxpix
	    				k = (j-1)*nxpix + i
                       	 	x(k)=(2.*i-1.-nxpix)*pixel/2.
                       	 	y(k)=(2.*j-1.-nypix)*pixel/2.
!	    				x(k) = -float(nxpix/2+1)*pixel + i*pixel
!	    				y(k) = -float(nypix/2+1)*pixel + j*pixel
!	    				if (x(k)==0. .and. y(k)==0.) then
!	      				x(k) = 0.
!	      				y(k) = pixel/2.
!	    				endif
	    				kappa(k) = 0.
	  			enddo
			enddo
		
			do j=1,natoms
	
!       			Copy parameters:
                        
				z=zcluster_s(j)
				D=calcAngDiamDis(0.d0,z)
				zs=zsource
				zsmin=zsource
				zsmax=zsource
				z_s=zsource

	  			GeoPars(j,1) = x0_s(j)
	  			GeoPars(j,2) = y0_s(j)
	  
!       			Sort out elliptical geometry:
        			if (GeoModel==2) then
	   		 		GeoPars(j,3) = theta_s(j)
	    				GeoPars(j,4) = axratio_s(j)
          				call EllGeometry(j)
	  			endif
	    
	  			if (MassModel==1) then
	    				if (NFWstyle(j)==1) then
	      				MassPars(j,1) = rs_s(j)
		      			MassPars(j,2) = ps_s(j)
		    			elseif (NFWstyle(j)==2) then
	      				MassPars(j,1) = c200_s(j) 
	      				MassPars(j,2) = M200_s(j)
		    			elseif (NFWstyle(j)==3) then
		      			MassPars(j,1) = M200_s(j) 
	    				elseif (NFWstyle(j)==4) then
	      				MassPars(j,1) = thetaE_s(j) 
		      			MassPars(j,2) = M200_s(j) 
		    			endif  
	  			elseif (MassModel==1) then 
		    			MassPars(j,1) = M200_s(j)
		  		endif  
	
!     	  		Now calculate the convergence at each pixel position:
	 	
	  			call LensFields(j,flag,.false.)
	 	
			enddo
	
! 			Chop out the first part of the kappa vector into the mass array:	
	
 			do j=1,nypix
	  			do i=1,nxpix
	    				k=(j-1)*nxpix+i
		    			massmap(i,j) = kappa(k)
		    			if (massmap(i,j)==0.) write(*,*) i,j,k,massmap(i,j)
	 			enddo
			enddo
            
            	!write kappa in ascii catalogue
            	write(catfile(ic:ic+9),'(a10)') '_kappa.cat'
      		write(*,*) '  Writing true map to file: ',catfile(1:ic+10) 
            	call writekappa(nxpix*nypix,x,y,kappa,catfile)
		
! 			Write out in fits format
      	
	      	write(catfile(ic:ic+10),'(a11)') '_kappa.fits'
      		write(*,*) '  Writing true map to file: ',catfile(1:ic+10) 
			call writefits(massmap,nxpix,nypix,pixel,catfile)
      		write(*,*)
		endif  

! 		End of gravitational lensing data simulation.	
      	write(*,*) 'End of gravitational lensing data simulation.'	
	      write(*,*)  
      	write(*,*)  

	endif

!-----------------------------------------------------------------------


! Sunyaev-Zel'dovich effect data and Cosmic strings data
! ******************************************************

 40	if (SZ==1.and.ModelClass==1) then

	   write(*,'(a)') ' SZ visibility data:'
	   write(*,'(a)') ' ^^^^^^^^^^^^^^^^^^^'
	   write(*,'(a)') ' '
	   write(*,'(a)') ' Can only write out .vis files at SLAC...'
	   
!FF	   SZscale = y2Jy(cell)
	   
	elseif (ModelClass==2) then

	   write(*,'(a)') ' Cosmic String visibility data:'
	   write(*,'(a)') ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'

	   SZscale = 0.001
	endif     

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! Loop over pointings:

        if (SZ==1.or.ModelClass==2) then

	   mVSA = 0
	   mSEA = 0
	   mAMI = 0
	   mRT = 0
	   mCBI = 0
	   DO m=1,Nvisfiles

! First get parameters of each observation:

	  call phitsread_open(Sobsfile(m))
	  pb_x(m) = CRVAL5
	  pb_y(m) = CRVAL6
	  ra = map_x*deg2rad
	  dec = map_y*deg2rad
	  ra2 = pb_x(m)*deg2rad
	  dec2 = pb_y(m)*deg2rad
	  call calc_ra_dec_off(ra,dec,ra2,dec2,pb_x_offset(m),pb_y_offset(m)) 	  
! NB. RA-Dec is left handed, need -ve sign here. 
!         pb_x_offset(m) = pb_x_offset(m)*rad2sec
	  pb_x_offset(m) = -pb_x_offset(m)*rad2sec
	  pb_y_offset(m) = pb_y_offset(m)*rad2sec

	  ph_x(m) = OBSRA
	  ph_y(m) = OBSDEC
	  ra = map_x*deg2rad
	  dec = map_y*deg2rad
	  ra2 = ph_x(m)*deg2rad
	  dec2 = ph_y(m)*deg2rad
	  call calc_ra_dec_off(ra,dec,ra2,dec2,ph_x_offset(m),ph_y_offset(m)) 	  
! NB. RA-Dec is left handed, need -ve sign here. 
!         ph_x_offset(m) = ph_x_offset(m)*rad2sec
	  ph_x_offset(m) = -ph_x_offset(m)*rad2sec
	  ph_y_offset(m) = ph_y_offset(m)*rad2sec
	  telescope(m) = TELESCOP
	  write(*,*)    
	  write(*,*) 'Pointing number ',m     
	  write(*,*) 'Template filename: ',Sobsfile(m)(1:70)     
	  write(*,*) 'Telescope: ',telescope(m)(1:3)     
	  write(*,*) 'Pointing centre: ',pb_x(m),pb_y(m)     
	  write(*,*) 'Phase centre: ',ph_x(m),ph_y(m)     
	  write(*,*) 'Relative pointing centre: ',pb_x_offset(m),pb_y_offset(m)     
	  write(*,*) 'Relative phase centre: ',ph_x_offset(m),ph_y_offset(m)     
      
! Set up output filename:

	  string = Sroot
	  do ic=1,100
	     read(string(ic:ic),'(a1)') char
	     if (char==' ') goto 45
	  enddo
 45	  write(visfile(m)(1:ic-1),'(a)') string(1:ic-1)
	  do j=ic,100
	     write(visfile(m)(j:j),'(a1)') ' '
	  enddo
	  visfile2(m) = visfile(m)

	  if (nreals==1) then
	     if (telescope(m)(1:3)=='VSA') then
		write(visfile(m)(ic:ic+6),'(a7)') '_VSAe_00'
		write(visfile2(m)(ic:ic+6),'(a7)') '_VSAe_00'
		mVSA = mVSA + 1
		if (mVSA<10) then
		   write(visfile(m)(ic+6:ic+6),'(i1)') mVSA
		   write(visfile2(m)(ic+6:ic+6),'(i1)') mVSA
		elseif (mVSA<100) then
		   write(visfile(m)(ic+5:ic+6),'(i2)') mVSA
		   write(visfile2(m)(ic+5:ic+6),'(i2)') mVSA
		else
		   stop 'Too many pointings in mosaic! Stopping'
		endif
!       write(visfile(m)(ic+7:ic+15),'(a9)') '_vis.fits'
!       write(visfile2(m)(ic+7:ic+21),'(a15)') '_noCMB_vis.fits'
		write(visfile(m)(ic+7:ic+14),'(a8)') '_vis.vis'
		write(visfile2(m)(ic+7:ic+20),'(a14)') '_noCMB_vis.vis'

	     elseif (telescope(m)(1:3)=='VSAs') then
		write(visfile(m)(ic:ic+6),'(a7)') '_VSAs_00'
		write(visfile2(m)(ic:ic+6),'(a7)') '_VSAs_00'
		mSEA = mSEA + 1
		if (mSEA<10) then
		   write(visfile(m)(ic+6:ic+6),'(i1)') mSEA
		   write(visfile2(m)(ic+6:ic+6),'(i1)') mSEA
		elseif (mSEA<100) then
		   write(visfile(m)(ic+5:ic+6),'(i2)') mSEA
		   write(visfile2(m)(ic+5:ic+6),'(i2)') mSEA
		else
		   stop 'Too many pointings in mosaic! Stopping'
		endif
!       write(visfile(m)(ic+7:ic+15),'(a9)') '_vis.fits'
!       write(visfile2(m)(ic+7:ic+21),'(a15)') '_noCMB_vis.fits'
		
	     elseif (telescope(m)(1:3)=='AMI') then
		write(visfile(m)(ic:ic+6),'(a7)') '_AMI_00'
		write(visfile2(m)(ic:ic+6),'(a7)') '_AMI_00'
		mAMI = mAMI + 1
		if (mAMI<10) then
		   write(visfile(m)(ic+6:ic+6),'(i1)') mAMI
		   write(visfile2(m)(ic+6:ic+6),'(i1)') mAMI
		elseif (mAMI<100) then
		   write(visfile(m)(ic+5:ic+6),'(i2)') mAMI
		   write(visfile2(m)(ic+5:ic+6),'(i2)') mAMI
		else
		   stop 'Too many pointings in mosaic! Stopping'
		endif
!       write(visfile(m)(ic+7:ic+15),'(a9)') '_vis.fits'
!       write(visfile2(m)(ic+7:ic+21),'(a15)') '_noCMB_vis.fits'
		write(visfile(m)(ic+7:ic+14),'(a8)') '_vis.vis'
		write(visfile2(m)(ic+7:ic+20),'(a14)') '_noCMB_vis.vis'

	     elseif (telescope(m)(1:2)=='RT') then
		write(visfile(m)(ic:ic+5),'(a6)') '_RT_00'
		write(visfile2(m)(ic:ic+5),'(a6)') '_RT_00'
		mRT = mRT + 1
		if (mRT<10) then
		   write(visfile(m)(ic+5:ic+5),'(i1)') mRT
		   write(visfile2(m)(ic+5:ic+5),'(i1)') mRT
		elseif (mRT<100) then
		   write(visfile(m)(ic+4:ic+5),'(i2)') mRT
		   write(visfile2(m)(ic+4:ic+5),'(i2)') mRT
		else
		   stop 'Too many pointings in mosaic! Stopping'
		endif
!       write(visfile(m)(ic+6:ic+14),'(a9)') '_vis.fits'
!       write(visfile2(m)(ic+6:ic+20),'(a15)') '_noCMB_vis.fits'
		write(visfile(m)(ic+6:ic+13),'(a8)') '_vis.vis'
		write(visfile2(m)(ic+6:ic+19),'(a14)') '_noCMB_vis.vis'

	     elseif (telescope(m)(1:3)=='CBI') then
		write(visfile(m)(ic:ic+6),'(a7)') '_CBI_00'
		write(visfile2(m)(ic:ic+6),'(a7)') '_CBI_00'
		mCBI = mCBI + 1
		if (mCBI<10) then
		   write(visfile(m)(ic+6:ic+6),'(i1)') mCBI
		   write(visfile2(m)(ic+6:ic+6),'(i1)') mCBI
		elseif (mCBI<100) then
		   write(visfile(m)(ic+5:ic+6),'(i2)') mCBI
		   write(visfile2(m)(ic+5:ic+6),'(i2)') mCBI
		else
		   stop 'Too many pointings in mosaic! Stopping'
		endif
!       write(visfile(m)(ic+7:ic+15),'(a9)') '_vis.fits'
!       write(visfile2(m)(ic+7:ic+21),'(a15)') '_noCMB_vis.fits'
		write(visfile(m)(ic+7:ic+14),'(a8)') '_vis.vis'
		write(visfile2(m)(ic+7:ic+20),'(a14)') '_noCMB_vis.vis'

	     else
		write(*,*) 'Unknown SZ telescope:'
		write(*,*) '  Template file = ',Sobsfile(m)(1:70)   
		write(*,*) '      telescope = ',telescope(m)(1:3)   
		stop 
	     endif
	     
	  elseif (nreals > 1) then
          
          write(visfile(m)(ic:ic+3),'(a4)') '_000'
          write(visfile2(m)(ic:ic+3),'(a4)') '_000'
          
          if (telescope(m)(1:3)=='VSA') then
            write(visfile(m)(ic+4:ic+10),'(a7)') '_VSA_00'
            write(visfile2(m)(ic+4:ic+10),'(a7)') '_VSA_00'
            mVSA = mVSA + 1
            if (mVSA<10) then
              write(visfile(m)(ic+10:ic+10),'(i1)') mVSA
              write(visfile2(m)(ic+10:ic+10),'(i1)') mVSA
            elseif (mVSA<100) then
              write(visfile(m)(ic+9:ic+10),'(i2)') mVSA
              write(visfile2(m)(ic+9:ic+10),'(i2)') mVSA
            else
              stop 'Too many pointings in mosaic! Stopping'
            endif
!             write(visfile(m)(ic+11:ic+19),'(a9)') '_vis.fits'
!             write(visfile2(m)(ic+11:ic+25),'(a15)') '_noCMB_vis.fits'
            write(visfile(m)(ic+11:ic+18),'(a8)') '_vis.vis'
            write(visfile2(m)(ic+11:ic+24),'(a14)') '_noCMB_vis.vis'

	 elseif (telescope(m)(1:3)=='VSAs') then
            write(visfile(m)(ic+4:ic+10),'(a7)') '_VSAs_00'
            write(visfile2(m)(ic+4:ic+10),'(a7)') '_VSAs_00'
            mSEA = mSEA + 1
            if (mSEA<10) then
	       write(visfile(m)(ic+10:ic+10),'(i1)') mSEA
	       write(visfile2(m)(ic+10:ic+10),'(i1)') mSEA
            elseif (mSEA<100) then
	       write(visfile(m)(ic+9:ic+10),'(i2)') mSEA
	       write(visfile2(m)(ic+9:ic+10),'(i2)') mSEA
            else
	       stop 'Too many pointings in mosaic! Stopping'
            endif
!       write(visfile(m)(ic+11:ic+19),'(a9)') '_vis.fits'
!       write(visfile2(m)(ic+11:ic+25),'(a15)') '_noCMB_vis.fits'
            write(visfile(m)(ic+11:ic+18),'(a8)') '_vis.vis'
            write(visfile2(m)(ic+11:ic+24),'(a14)') '_noCMB_vis.vis'	  

	 elseif (telescope(m)(1:3)=='AMI') then
            write(visfile(m)(ic+4:ic+10),'(a7)') '_AMI_00'
            write(visfile2(m)(ic+4:ic+10),'(a7)') '_AMI_00'
            mAMI = mAMI + 1
            if (mAMI<10) then
	       write(visfile(m)(ic+10:ic+10),'(i1)') mAMI
	       write(visfile2(m)(ic+10:ic+10),'(i1)') mAMI
            elseif (mAMI<100) then
	       write(visfile(m)(ic+9:ic+10),'(i2)') mAMI
	       write(visfile2(m)(ic+9:ic+10),'(i2)') mAMI
            else
	       stop 'Too many pointings in mosaic! Stopping'
            endif
!       write(visfile(m)(ic+11:ic+19),'(a9)') '_vis.fits'
!       write(visfile2(m)(ic+11:ic+25),'(a15)') '_noCMB_vis.fits'
            write(visfile(m)(ic+11:ic+18),'(a8)') '_vis.vis'
            write(visfile2(m)(ic+11:ic+24),'(a14)') '_noCMB_vis.vis'

	 elseif (telescope(m)(1:2)=='RT') then
            write(visfile(m)(ic+4:ic+9),'(a6)') '_RT_00'
            write(visfile2(m)(ic+4:ic+9),'(a6)') '_RT_00'
            mRT = mRT + 1
            if (mRT<10) then
	       write(visfile(m)(ic+9:ic+9),'(i1)') mRT
	       write(visfile2(m)(ic+9:ic+9),'(i1)') mRT
            elseif (mRT<100) then
	       write(visfile(m)(ic+8:ic+9),'(i2)') mRT
	       write(visfile2(m)(ic+8:ic+9),'(i2)') mRT
            else
	       stop 'Too many pointings in mosaic! Stopping'
            endif
!       write(visfile(m)(ic+10:ic+18),'(a9)') '_vis.fits'
!       write(visfile2(m)(ic+10:ic+24),'(a15)') '_noCMB_vis.fits'
            write(visfile(m)(ic+10:ic+17),'(a8)') '_vis.vis'
            write(visfile2(m)(ic+10:ic+23),'(a14)') '_noCMB_vis.vis'
	    
	 elseif (telescope(m)(1:3)=='CBI') then
            write(visfile(m)(ic+4:ic+10),'(a7)') '_CBI_00'
            write(visfile2(m)(ic+4:ic+10),'(a7)') '_CBI_00'
            mCBI = mCBI + 1
            if (mCBI<10) then
	       write(visfile(m)(ic+10:ic+10),'(i1)') mCBI
	       write(visfile2(m)(ic+10:ic+10),'(i1)') mCBI
            elseif (mCBI<100) then
	       write(visfile(m)(ic+9:ic+10),'(i2)') mCBI
	       write(visfile2(m)(ic+9:ic+10),'(i2)') mCBI
            else
	       stop 'Too many pointings in mosaic! Stopping'
            endif
!       write(visfile(m)(ic+11:ic+19),'(a9)') '_vis.fits'
!       write(visfile2(m)(ic+11:ic+25),'(a15)') '_noCMB_vis.fits'
            write(visfile(m)(ic+11:ic+18),'(a8)') '_vis.vis'
            write(visfile2(m)(ic+11:ic+24),'(a14)') '_noCMB_vis.vis'

          else
            write(*,*) 'Unknown SZ telescope:'
            write(*,*) '  Template file = ',Sobsfile(m)(1:70)   
            write(*,*) '      telescope = ',telescope(m)(1:3)   
            stop 
          endif
        endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Read in uv points (Nb read in thermal noise values too - output file
! contains 'unknown' CMB noise):

        if (GCOUNT>large) then
          write(*,*) 'Error: in MakeData.f Nvis = ',GCOUNT
          write(*,*) '       in McAdam.inc Nmax = ',large
          write(*,*) ' --> Change McAdam.inc to avoid overflow.'
          write(*,*)
          stop
        elseif (GCOUNT<Nvis(m)) then
          write(*,*) 'Error: in MakeData.f Nvis = ',GCOUNT
          write(*,*) '       in McAdam.inc Nvis = ',Nvis(m)
          write(*,*) ' --> Change McAdam.inc.'
          write(*,*)
          stop
        endif

          write(*,'(a,i5,a,e9.2,a)') ' Observed ',GCOUNT, &
          ' complex visibilities, at a frequency of ',CRVAL4,' Hz.'
          write(*,*)

        do i=1,GCOUNT
          call phitsread_one(i,vis_buffer,uv_buffer,weight_buffer,rms_buffer,j)
          u(m,i) = 1.*uv_buffer(1)
          v(m,i) = 1.*uv_buffer(2)
          visrms(m,i) = 1./sqrt(weight_buffer)
          viswt(m,i) = 1.*weight_buffer
          baseline(i) = j
        enddo

        call phitsread_close

	ENDDO
      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! Make the y map - first initialise to zero:
	  
	do i=1,nx
	  do j=1,ny
	    ymap(i,j) = 0.
	  enddo
	enddo

! Make sure SZscale parameter is set correctly - in inference mode this
! is calculated in Initialise, which is called AFTER the vis file has
! been read in. See start of SZ/strings simulation.
	  	  
! Then add atoms
	  
	do j=1,Natoms

	   if (ModelClass==1) then
!       Copy parameters:

	      GeoPars(j,1) = x0_s(j)
	      GeoPars(j,2) = y0_s(j)
	  
!       Sort out elliptical geometry:
	      if (GeoModel==2) then
		 GeoPars(j,3) = theta_s(j)
		 GeoPars(j,4) = axratio_s(j)
		 call EllGeometry(j)
	      endif
	    

	      if (GasModel > 1) then
!       Potential parameters:
		 if (MassModel==1) then
		    if (NFWstyle(j)==1) then
		       MassPars(j,1) = rs_s(j)
		       MassPars(j,2) = ps_s(j)
		    elseif (NFWstyle(j)==2) then
		       MassPars(j,1) = c200_s(j)
		       MassPars(j,2) = M200_s(j)
		    elseif (NFWstyle(j)==2) then  
		       MassPars(j,1) = M200_s(j)
		    else
		       stop 'Only 3 NFW styles implemented.'
		    endif 
		 else
		    stop 'Only NFW HSE model implemented.'
		 endif
	      endif

!       Gas parameters:
	      if (GasModel==0 .or. GasModel==1) then
			GasPars(j,1) = rc_s(j)
		 	GasPars(j,2) = beta_s(j)
		 	GasPars(j,3) = Mgas_s(j)
			TPars(j,1) = Tgas_s(j)
	      elseif (GasModel==2) then
		 	if (TModel==1) then
				GasPars(j,1) = Mgas_s(j)
		    		TPars(j,1) = Tgas_s(j)
		 	elseif (TModel==2) then
		    		GasPars(j,1) = Mgas_s(j)
		    		TPars(j,1) = Tgas_s(j)
		    		TPars(j,2) = Tgamma_s(j)
		 	endif	
	      elseif (GasModel >= 4) then 
		 	stop 'Only 3 gas models implemented.'
	      endif  

!       Now fill up the relevant gas distribution arrays:

	  idummy = 0
	  call MakeGasDistributions(j,idummy)
	  if (idummy.ne.0) stop ' T goes negative with these parameters.'
	  

	elseif (ModelClass==2) then
	   
	   StringGeoPars(1) = x0_s(j)
	   StringGeoPars(2) = y0_s(j)
	   StringPars(1) = flux_s(j)
	   StringPars(2) = angle_s(j)
	   StringPars(3) = size_s(j)
	   
	endif

!       Now increase the y parameter at each grid position:

	  call MakeYMap(j)

	enddo
	
!	write(*,*) 'central y value is ',ymap(nx/2+1,ny/2+1)
	
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!     Calculate visibilities at each uv position - these are stored in
!     the pvisr and pvisi arrays:

	call MakeVisibilities
	
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!     Add point sources in the u-v plane:

	if (SourceSubtract==1) then
	
	do i=1,NSrc
	  SrcPars(i,1) = srcx_s(i)
	  SrcPars(i,2) = srcy_s(i)
	  SrcPars(i,3) = srcf_s(i)
	enddo
	
	call AddSources
	
	endif
	
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! For each noise realisation,
            
	do k=1,Nreals

! and each pointing,

      do m=1,Nvisfiles

! Initialise Random Number generator:
	  
	seed = 13*k
	!i = RandomInit(state,seed)
      call initRandom()

      if (IncludeCMB(m)==1) then

! Approximation: skip this part for Ryle data and only produce noCMB
! data.

        if (telescope(m)(1:2)=='RT') goto 48

! Add noise, using covariance matrix generated by writechol, which
! includes both thermal noise and primordial CMB:

        !call gauss(0.0,1.0,2*large,unitnoise)
        do i=1,2*large
        	unitnoise(i)=Gaussian1NS(0)
        end do

! Open L matrix data file

        open(unit=9,file=covmatfile(m),form='unformatted',status='old')
        read(9) idummy
        if (idummy.ne.2*NLCM) then
          write(*,*)  'Array size mismatch - expecting ',2*NLCM
          write(*,*)  '                      received  ',idummy
          write(*,*)  '-> Check sizes of template files and .LCM file.'
          stop
        endif

! Read in matrix one row at a time and compute y = L*x

        do i=1,2*NLCM
          read(9) (row(j),j=1,i)
          sum = 0.
          do j=1,i
            sum=sum+row(j)*unitnoise(j)
          enddo
          scalednoise(i)=sum
        enddo
        close(9)

! Add noise to predicted visibilities to get observed visibilities:

        do i=1,NLCM
          visr(m,i)=pvisr(m,i)+scalednoise(i)
        enddo
        do i=NLCM+1,2*NLCM
          j =i-NLCM
          visi(m,j)=pvisi(m,j)+scalednoise(i)
        enddo

! Package up buffers ready for writing to fits file:

        nbad = 0
        do i=1,Nvis(m)
          if (viswt(m,i) .ne. 0.) then
            nbad = nbad + 1
            uv_dump(1,nbad) = 1.*u(m,i)
            uv_dump(2,nbad) = 1.*v(m,i)
            vis_dump(nbad) = dcmplx(1.*visr(m,i),1.*visi(m,i))
            weight_dump(nbad) = 1.*viswt(m,i)
!       write(*,*) i,u(i),v(i),visr(i),visi(i),viswt2(i)
          endif
        enddo

! And set new keywords!

        GCOUNT = nbad
        TELESCOP = telescope(m)
	  CRVAL5 = pb_x(m)
	  CRVAL6 = pb_y(m)
        OBSRA = ph_x(m) 
        OBSDEC = ph_y(m)

        nbad = Nvis(m) - nbad

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Write out data to files

        if (nreals==1) then

         write(*,*) 'Writing data to file: ',visfile(m)(1:ic+25)
!          call write_phits(vis_dump,uv_dump,weight_dump,baseline,visfile(m))
         call write_vis(vis_dump,uv_dump,weight_dump,baseline,visfile(m))

        else

          if (k<10) then
            write(visfile(m)(ic+3:ic+3),'(i1)') k
          elseif (k<100) then
            write(visfile(m)(ic+2:ic+3),'(i2)') k
          else
            write(visfile(m)(ic+1:ic+3),'(i3)') k
          endif

          write(*,*) 'Writing data to file: ',visfile(m)(1:ic+25)
!           call write_phits(vis_dump,uv_dump,weight_dump,baseline,visfile(m))
         call write_vis(vis_dump,uv_dump,weight_dump,baseline,visfile(m))

!           if (mod(k,10)==0) write(*,*) k

        endif
	
	endif
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Extra bit - make a vis.fits file with no added CMB (but thermal noise
! realisation is also different...)
	    
 48	do i=1,2*large
	   unitnoise(i)=Gaussian1NS(0)
      end do	
      !call gauss(0.d0,1.d0,2*large,unitnoise)

	do i=1,Nvis(m)
	  visr(m,i) = pvisr(m,i) + visrms(m,i)*unitnoise(i)
	enddo
	do i=Nvis(m)+1,2*Nvis(m)
	  j = i - Nvis(m)  
	  visi(m,j) = pvisi(m,j) + visrms(m,j)*unitnoise(i)
	enddo

	nbad = 0
	do i=1,Nvis(m)
	  if (viswt(m,i).ne.0.0) then 	  
	    nbad = nbad + 1
	    uv_dump(1,nbad) = 1.0*u(m,i)
	    uv_dump(2,nbad) = 1.0*v(m,i)
	    vis_dump(nbad) = dcmplx(1.0*visr(m,i),1.0*visi(m,i))
	    weight_dump(nbad) = 1.0*viswt(m,i)
	  endif
	enddo
	
	GCOUNT = nbad
      TELESCOP = telescope(m)
      CRVAL5 = pb_x(m)
      CRVAL6 = pb_y(m)
      OBSRA = ph_x(m)
      OBSDEC = ph_y(m)
        
	nbad = Nvis(m) - nbad

	if (nreals==1) then
        write(*,*) 'Writing data (noCMB) to file: ',visfile2(m)(1:ic+25)
        write(*,*)
!         call write_phits(vis_dump,uv_dump,weight_dump,baseline,visfile2(m))     
        call write_vis(vis_dump,uv_dump,weight_dump,baseline,visfile2(m))     
	else
        if (k<10) then
	    write(visfile2(m)(ic+3:ic+3),'(i1)') k  
	  elseif (k<100) then
	    write(visfile2(m)(ic+2:ic+3),'(i2)') k  
	  else 
	    write(visfile2(m)(ic+1:ic+3),'(i3)') k  
        endif 
        write(*,*) ' Writing data (noCMB) to file: ',visfile2(m)(1:ic+25)
        write(*,*)
!         call write_phits(vis_dump,uv_dump,weight_dump,baseline,visfile2(m))     	    	
        call write_vis(vis_dump,uv_dump,weight_dump,baseline,visfile2(m))     	    	
	endif
	
! ! TESTING:
! ! write out noise-free visibilities for mapping checks:
! 	nbad = 0
! 	do i=1,Nvis(m)
! 	  if (viswt(m,i).ne.0.0) then 	  
! 	    nbad = nbad + 1
! 	    uv_dump(1,nbad) = 1.0*u(m,i)
! 	    uv_dump(2,nbad) = 1.0*v(m,i)
! 	    vis_dump(nbad) = dcmplx(1.0*pvisr(m,i),1.0*pvisi(m,i))
! 	    weight_dump(nbad) = 1.0*viswt(m,i)
! 	  endif
! 	enddo
! 	GCOUNT = nbad
! 	nbad = Nvis(m) - nbad
!       call write_phits(vis_dump,uv_dump,weight_dump,baseline,visfile(m))     
! 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     	
	enddo
	enddo

! Finished creating simulated data files.  

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
! Write out noise-free map, in units of y parameter, and Jy:	
	
	string = Sroot
	do ic=1,100
        read(string(ic:ic),'(a1)') char
        if (char==' ') goto 50
      enddo
 50   write(visfile(1)(1:ic-1),'(a)') string(1:ic-1)
	do j=ic,100
        write(visfile(1)(j:j),'(a1)') ' '
      enddo
	      
      write(visfile(1)(ic:ic+11),'(a12)') '_true_y.fits'	
	call writefits(ymap,nx,ny,cell,visfile(1))  
	
      sum = 0.0
      do i=1,nx
         do j=1,ny
            sum = sum + ymap(i,j)
         end do
      end do
      sum = sum/float(nx*ny)	
	do i=1,nx
         do j=1,ny
!FF            szsky(i,j) = (ymap(i,j) - sum)*SZscale
         end do
      end do
	
      write(visfile(1)(ic:ic+11),'(a12)') '_true_S.fits'	
	call writefits(szsky,nx,ny,cell,visfile(1))  
	
! Save arrays too:

	if (SZ==1.and.GasModel > 0.and.ModelClass==1) then
	          
	  write(visfile(1)(ic:ic+11),'(a12)') '_rhogas.txt '	
	  open (unit=23,file=visfile(1),status='unknown')
	  do i=1,n
	    write(23,*) r(i),Rhogas(i)
	  enddo
	  close(23)
	  
	  write(visfile(1)(ic:ic+11),'(a12)') '_T.txt      '	
	  open (unit=23,file=visfile(1),status='unknown')
	  do i=1,n
	    write(23,*) r(i),T(i)
	  enddo
	  close(23)
	  
	  write(visfile(1)(ic:ic+11),'(a12)') '_Pgas.txt   '	
	  open (unit=23,file=visfile(1),status='unknown')
	  do i=1,n
	    write(23,*) r(i),Pgas(i)
	  enddo
	  close(23)
	  
	  write(visfile(1)(ic:ic+11),'(a12)') '_y.txt      '	
	  open (unit=23,file=visfile(1),status='unknown')
	  do i=1,n
	    write(23,*) r(i),Yarray(i)
	  enddo
	  close(23)
	endif	
	
! End of SZ data simulation.	
	endif
!-----------------------------------------------------------------------
	
      write(*,*)
      write(*,*)
      
	return
	end subroutine MakeData	
  
!======================================================================

	subroutine observe(obsarea,ndata)
	
	implicit none

	integer nx1,ny1,ndata
	parameter(nx1=128,ny1=128)
	double precision pixel,obsregion(nx1,ny1),obsarea
      double precision urv
	
	double precision xmin,xmax,ymin,ymax,tr(6),xx,yy

	integer i,j,k

!-----------------------------------------------------------------------

! read in a flag map (first realisation only) and calculate observing 
! area in sq. arcmin

	if (obsarea==0.) then	
	  
	  call readfits(obsregion,nx1,ny1,pixel,GLobsfile)	  
	  	  
	  do i=1,nx1
	    do j=1,ny1
	      if (obsregion(i,j)==1.0) then
	        obsarea = obsarea + pixel*pixel
            endif
	    enddo
	  enddo  
	
        obsarea = obsarea/3600.0
	  ndata = obsarea*ndensity
	
	endif  
	
! Set up geometry
	
	i = nx1/2+1
      j = ny1/2+1
	tr(1)=-float(i)*pixel
      tr(2)=pixel
      tr(3)=0.
      tr(4)=-float(j)*pixel
      tr(5)=0.
      tr(6)=pixel
      xmin=-(float(i-1)+0.5)*pixel
      xmax=(float(i-1)-0.5)*pixel
      ymin=-(float(j-1)+0.5)*pixel
      ymax=(float(j-1)-0.5)*pixel

! Now generate randomimage positions within the obsregion map
 	
	do k=1,ndata
     	  urv=ranmarNS(0)
10	  xx=xmin+(xmax-xmin)*urv
        urv=ranmarNS(0)
	  yy=ymin+(ymax-ymin)*urv
	  i = nint((xx - tr(1))/tr(2))
	  j = nint((yy - tr(4))/tr(6))
	  if (obsregion(i,j)==0.) then
	    goto 10
	  else
	    x(k) = 1.*xx	  
	    y(k) = 1.*yy
	  endif  	  
	enddo
	
	return
	end subroutine observe
  
!======================================================================


      subroutine write_vis(vis_buffer,uv_buffer,weight_buffer,baselines,outfile)     	    	

      implicit none
	
      character*100 outfile
      double complex vis_buffer(large)
      double precision weight_buffer(large),uv_buffer(2,large)
      integer baselines(large)

	integer i,n
      double precision u,v,re,im,rms
	
!-----------------------------------------------------------------------

!  Open the output file

      open(unit=43,file=outfile,status='unknown')

!  Write visibilities for each baseline as a single group of data

      do i = 1,GCOUNT
        n = baselines(i)
        u = uv_buffer(1,i)
        v = uv_buffer(2,i)
        re = dble(vis_buffer(i))
        im = aimag(vis_buffer(i))
        rms = 1.0/sqrt(weight_buffer(i))
	  write(43,*) n,u,v,re,im,rms		
      enddo

!  Close the file

      close(43)

      end subroutine write_vis
        
!=======================================================================
end module MakeData1
