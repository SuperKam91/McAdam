module GasModels
	use params
      	use constants
      	use telescope1
	use utilities
      	use matrix_utils
      	use fftw_inc
      	!use rfft
      	use cosmology
      
      	implicit none
      	double precision lmin,lmax,gd1,gd2

contains

!=======================================================================

	subroutine MakeGasDistributions(k,flag)

! Make 1-dimensional distributions needed for modelling cluster gas in
! the kth atom.
! Almost certainly needs tidying up, including further modularisation.
! Models to be included should be a HSE profile based on beta model for
! gas.

	implicit none
			
      integer i,j,k,flag
	
	double precision A1,Gm,rc,beta,yc,rhocritz,index,index_g
      parameter(A1=1.424d-19,Gm=3.12d-14)
	
      double precision Mgas200,rs,ps,cc,Rhogas0,logRhogas0,T0,logT0,rr,rlimit1
      double precision prefactor,prefactor2
      double precision thetaE
	double precision eps
	parameter(eps=1d-4)
	
!-----------------------------------------------------------------------

      rhocritz=rhocritofz(z(i))
      call lookUp1D(Dn,z(i),lookD(:,1),lookD(:,2),D)
      
      Phi = 0.d0
	Rhogas = 0.d0
	T = 0.d0
	logPhi = 0.d0
	logRhogas = 0.d0
	logT = 0.d0
	
!-----------------------------------
! Isothermal beta models:
!-----------------------------------

! Hobson and McLachlan beta model, normalised to central y parameter:
! Parameters of kth atom are 
!    GasPars(k,1)=apparent core radius (arcmin)
!    GasPars(k,2)=central y parameter\
! Beta is fixed at 2/3, and profile is "truncated" at rt=3*rc
! Note no need for D or z to take sensible values-the r array is
! interpreted as being in arcminutes not Mpc. This is ok even at z=0.05,
! where rmax=100 corresponds to 4 Mpc or so. Also note that this profile
! falls to 0.1% of its maximum at about 18 core radii.

      if(GasModel==0) then
		rc=GasPars(k,1)
	  	yc=GasPars(k,2)
	   
        	do i=1,n
	    		Yarray(i)=BetaModelHM(r(i),rc,yc)
	    		if(Yarray(i)<0.0) then
	      		write(*,*) 'about to take log of negative number: '
	      		write(*,*) '  i,Yarray(i)=',i,Yarray(i)
	      		write(*,*) '  yc,rc=',yc,rc
	      		write(*,*) '  GasPars=',GasPars
	    		endif
	    		logYarray(i)=phlog10(Yarray(i))
	  	enddo
	  
	  	flag=0
	  	return
	    	
!-----------------------------------

! Isothermal beta model, normalised to gas mass within 3D radius rmass if 
! Mass=0, or within 3D radius r200 if Mass=1; projection is numerical.
! Note tying of beta model core radii together by using CPLstyle-1:
         
      elseif(GasModel==1) then

 	  	if(Mass==1.and.MassModel==3.and.CPLstyle(k)==-1) then 
	    		rc=MassPars(k,1)/1000.d0
        	else
          		rc=GasPars(k,1)/1000.d0
	  	endif
        
        	beta=GasPars(k,2) 
		  	    
        	if(Mass==1) then
	   		if(MassModel==1) then
!        			NFW profile for potential:
          			if(NFWstyle(k)==1) then 
	      			stop 'Cannot calculate r200 in this NFW style.'    
	    			elseif(NFWstyle(k)==2) then
	      			M200=MassPars(k,2)
	      			r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
	    			elseif(NFWstyle(k)==3) then
	      			M200=MassPars(k,1)
	      			r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
	    			endif
         		elseif(MassModel==2) then
!        			SIS profile for potential:
          			if(SISstyle(k)==1) then 
	      			M200=MassPars(k,1)
	      			r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
          			elseif(SISstyle(k)==2) then 
	      			r200=sqrt(MassPars(k,1)*sec2rad*D*SigcritE/2.0)
	      			r200=r200*sqrt(3.0/(200.0*rhocritz*pi))
	    			endif
         		elseif(MassModel==3.and.CPLstyle(k)==-1) then
 	     			M200=MassPars(k,2)
           			r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
         		else  
	    			stop 'Can only do joint fit with NFW, SIS and CIS potentials!'
	   		endif

	  	else
          		if(Betastyle(k)==1) then 
!           		Unknown potential, normalise to fixed radius:
	      		r200=rmass

	    		elseif(Betastyle(k)==2.or.Betastyle(k)==3) then      
!           		HSE potential, compute r200:
		 		r200=max((9.d0/(4.d0*pi*Gm*200.d0*rhocritz))*beta*TPars(k,1)-rc*rc,0.d0)
             		r200=sqrt(r200)
             
            		if(Tmodel==2) then 
              			index_g=TPars(k,2)
              			rr=2.0/(3.0*beta*(index_g-1.0)+2.0)
              			r200=(r200*TPars(k,2))**rr
				endif
          		else  
	      		stop 'Unrecognised Beta style.'
	    		endif
	  
          		M200=4.0*pi*200.0*rhocritz*r200*r200*r200/3.0
                
        	endif
		if( M200 <= 0d0 .or. r200 <= 0d0 .or. beta <= 0d0 .or. rc <= 0d0 ) then
            		flag = 1
                  	return
		endif

! 		Gas density profile normalisation: by mass within fixed radius rmass 
! 		(Betastyle=1), by mass within r200 from HSE (Betastyle=2), or by gas 
! 		fraction within within r200 from HSE:    
        	if(Betastyle(k)==1) then
          		Mgas200=GasPars(k,3)
        	elseif(Betastyle(k)==2) then
          		Mgas200=GasPars(k,3)
        	elseif(Betastyle(k)==3) then
          		Mgas200=GasPars(k,3)*M200
        	endif
		
		!null run
		if( Mgas200 == 0d0 ) then
			flag = 2
			return
		endif
		
         
!       	Integrate profile to ensure model has correct gas mass
	  	do i=1,n
	    		Rhogas(i)=BetaModel3D(r(i),rc,beta,1.d0)
	    		logRhogas(i)=phlog10(Rhogas(i))
	  	enddo
        	
            rr=r200
        	Rhogas0=Mgas200/GasMass(rr,rc,beta)
        	Rhogas_central=Rhogas0
	  	logRhogas0=phlog10(Rhogas0)
	  
!       	Tabulate T profile and compute pressure:        
	  	do i=1,n
          		if(Tmodel==1) then
            		T(i)=TPars(k,1)
          		elseif(Tmodel==2) then
	      		index_g=TPars(k,2)
            		T(i)=TPars(k,1)*(Rhogas(i))**(index_g-1)
          		else
	      		stop 'Unknown temperature model.'
          		endif
          		
                  prefactor=A1*T(i)
	    		Rhogas(i)=Rhogas0*Rhogas(i)
	    		logRhogas(i)=logRhogas0+logRhogas(i)
	    		Pgas(i)=prefactor*Rhogas(i)
	    		logPgas(i)=phlog10(Pgas(i))
	  	enddo
	
!-----------------------------------
! Hydrostatic equilibrium models, numerical projection:
!-----------------------------------

! 	Isothermal temperature profile:

      elseif(GasModel==2.and.TModel==1) then
	
!     	First calculate gas density and temperature:	
	
        	T(1:n)=TPars(k,1)
	    
        	if(MassModel==1) then

!       		NFW profile for potential:
          		if(NFWstyle(k)==1) then 
	      		rs=MassPars(k,1) 
	      		ps=MassPars(k,2) 
	    		elseif(NFWstyle(k)==2) then
	      		cc=MassPars(k,1)
	      		M200=MassPars(k,2)
	      		r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
	      		rs=r200/cc
	      		ps=M200/(4.0*pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
	    		elseif(NFWstyle(k)==3) then

				!M200 now contains M100, so that Komatsu and Seljak's relation can be used!!
            		M200=MassPars(k,2)
            		cc=6.0*(M200/1e14)**(-0.2)
            		r200=(M200/(100.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
            		rs=r200/cc
            		ps=M200/(4.0*Pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
          		elseif(NFWstyle(k)==4) then
            		thetaE=MassPars(k,1)
            		M200=MassPars(k,2)
            		rE=ThetaE*sec2rad*D
            		scE=sigcritE
            		r200=(M200/(200.0*(4.0*Pi/3.0)*rhocritz))**(1.0/3.0)
            		rs=FINDR(rsfunc,1.0d-4,1.0d4,1.0d-4)
            		cc=r200/rs
            		ps=M200/(4.0*Pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
          		else
            		STOP 'Unrecognised NFWstyle.'
          		endif
                  
                  !sanity check
            	if( rs <= 0d0 .or. ps <= 0d0 ) then
            		flag = 1
                  	return
			endif
           	   
	    		prefactor=-4.d0*Pi*Gmu*ps*rs*rs/TPars(k,1)
	    
	    		do i=1,n
	      		rr=r(i)/rs
	      		Rhogas(i)=1.0*dexp(prefactor*((rr-log(1.0+rr))/rr))
	      		logRhogas(i)=phlog10(Rhogas(i))
	    		enddo

!         		Now normalise to Mgas200, and calculate pressure:
	    		rr=r200
          		Rhogas0=GasMass(rr,rc,beta)
          		Rhogas_central=Rhogas0

 	    		if(Rhogas0 > 1d-20) then
	      		Rhogas0=GasPars(k,1)/GasMass(rr,rc,beta)
            		Rhogas_central=Rhogas0
	      		logRhogas0=phlog10(Rhogas0)
	      		prefactor=A1*TPars(k,1)
	      		do i=1,n
	        			Rhogas(i)=Rhogas0*Rhogas(i)
	        			logRhogas(i)=logRhogas(i)+logRhogas0
	        			Pgas(i)=1.0*prefactor*Rhogas(i)
	        			logPgas(i)=phlog10(Pgas(i))
	      		enddo
	    		else
	      		Rhogas0=0d0
            		Rhogas_central=Rhogas0
	      		logRhogas0=-20d0
	      		Rhogas(1:n)=0d0
	        		logRhogas(1:n)=-45d0
	        		Pgas(1:n)=0d0
	        		logPgas(1:n)=-45d0
	        		Yarray(1:n)=0d0
	        		logYarray(1:n)=-45d0
				return
	    		endif
		else
	    		stop 'Cannot do non-NFW potentials!'
	  	endif
	   	   	
!-----------------------------------

! 	Polytropic temperature profile:
         
	elseif(GasModel==2.and.TModel==2) then
		  
!     	First calculate potential, then the temperature, then the density:	
	
        	if(MassModel==1) then

!       		NFW profile for potential:
          		if(NFWstyle(k)==1) then 
	      		rs=MassPars(k,1) 
	      		ps=MassPars(k,2) 
	    		elseif(NFWstyle(k)==2) then
	      		cc=MassPars(k,1)
	      		M200=MassPars(k,2)
	      		r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
	      		rs=r200/cc
	      		ps=M200/(4.0*pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
	    		elseif(NFWstyle(k)==3) then
! 				M200 now contains M100, so that Komatsu and Seljak's relation can be used!!
            		M200=MassPars(k,2)
            		cc=6.0*(M200/1e14)**(-0.2)
            		r200=(M200/(100.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
            		rs=r200/cc
            		ps=M200/(4.0*Pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
          		elseif(NFWstyle(k)==4) then
            		thetaE=MassPars(k,1)
            		M200=MassPars(k,2)
            		rE=ThetaE*sec2rad*D
            		scE=sigcritE
            		r200=(M200/(200.0*(4.0*Pi/3.0)*rhocritz))**(1.0/3.0)
            		rs=FINDR(rsfunc,1.0d-4,1.0d4,1.0d-4)
            		cc=r200/rs
            		ps=M200/(4.0*Pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
          		else
            		STOP 'Unrecognised NFWstyle.'
          		endif
         
	    		if(TPars(k,2)==1.0d0) then
	    			prefactor=-4.d0*Pi*Gmu*ps*rs*rs/TPars(k,1)
				T0=TPars(k,1)
				do i=1,n
	      			T(i)=TPars(k,1)
					rr=r(i)/rs
	      			Rhogas(i)=1.0d0*exp(prefactor*((rr-log(1.0d0+rr))/rr))
	      			logRhogas(i)=phlog10(Rhogas(i))
	    			enddo
	    		else
	    			prefactor=4.d0*Pi*Gmu*ps*rs*rs*(TPars(k,2)-1.0d0)/TPars(k,2)
	    			T0=1.d0*TPars(k,1)
	    			logT0=phlog10(T0)
	    			prefactor2=(1.d0/(TPars(k,2)-1.0d0))
	    			do i=1,n
	      			rr=r(i)/rs
	      			T(i)=TPars(k,1)-1.d0*prefactor*(((rr-log(1.d0+rr))/rr))
	      			logT(i)=phlog10(T(i))
	      			logRhogas(i)=prefactor2*(logT(i)-logT0)
	      			Rhogas(i)=10.d0**(logRhogas(i))
	    			enddo
	    		endif
                  
                  !sanity check
            	if( rs <= 0d0 .or. ps <= 0d0 ) then
            		flag = 1
                  	return
			endif

!         		Now normalise to Mgas200, and calculate pressure:

	    		rr=1.0*r200
          		Rhogas0=GasMass(rr,rc,beta)
          		Rhogas_central=Rhogas0
 	    		if(Rhogas0 > 1d-20) then
	      		Rhogas0=GasPars(k,1)/GasMass(rr,rc,beta)
            		Rhogas_central=Rhogas0
	      		logRhogas0=phlog10(Rhogas0)
	      		prefactor=A1
	      		do i=1,n
	        			Rhogas(i)=Rhogas0*Rhogas(i)
	        			logRhogas(i)=logRhogas(i)+logRhogas0
	        			Pgas(i)=1.0*prefactor*T(i)*Rhogas(i)
	        			logPgas(i)=phlog10(Pgas(i))
	      		enddo
	    		else
	      		Rhogas0=0.0
            		Rhogas_central=Rhogas0
	      		logRhogas0=-20d0
	      		Rhogas(1:n)=0d0
	        		logRhogas(1:n)=-45d0
	        		Pgas(1:n)=0d0
	        		logPgas(1:n)=-45d0
	        		Yarray(1:n)=0d0
	        		logYarray(1:n)=-45d0
				return
	    		endif
	  	else
	    		stop 'Cannot do non-NFW potentials!'
	  	endif
	   	   	
!-----------------------------------
	
      endif
     
	
! 	Finally project pressure to generate y arrays:

	if(GasModel /= 0) then
	  	do i=1,n
	    		uu=r(i)
          		rlimit1=rlimit
          
	    		call qtrap(PgasIntegrand,-rlimit1,rlimit1,eps,Yarray(i))
	    		if(Yarray(i) > 1d10) then
	      		write(*,*) 'i,r(i),rlimit1,Yarray(i)=',i,r(i),rlimit1,Yarray(i)
	      		write(*,*) 'QTRAP error. Printing pressure/density arrays:'
	      		do j=1,n
	          			write(*,*)j,r(j),Pgas(j),logPgas(j),logRhogas(j),Rhogas(j)
		  		enddo
	      		write(*,*) 'rhogas0=',rhogas0
	      		write(*,*) 'prefactor=',prefactor
	      		write(*,*) 'r200,rc,rlimit=',r200,rc,rlimit
	      		write(*,*) 'GasPars(1)=',GasPars(k,1)
	      		write(*,*) 'TPars(1)=',TPars(k,1)
	      		write(*,*) 'MassPars(1)=',MassPars(k,1)
	      		write(*,*) 'MassPars(2)=',MassPars(k,2)
            		pause
	    		endif
	    		logYarray(i)=phlog10(Yarray(i))
	  	enddo 
	endif   
	flag=0
      
      !store derived parameters
      aux(k,1)=rc/(sec2rad*D) !r_{core} in arcsec
      aux(k,2)=rhogas_central !central gas density in [h Msun Mpc^{-3}]
      aux(k,3)=rhogas_central/(1.14d0*m_p*(Mpc2m**3)) !central electron number density in [m^{-3}]
      aux(k,4)=totMass(1500d0,r200,rc,rhocritz,aux(k,22))
      aux(k,5)=totMass(1000d0,r200,rc,rhocritz,aux(k,23))
      aux(k,6)=totMass(500d0,r200,rc,rhocritz,aux(k,24))
      aux(k,7)=totMass(200d0,r200,rc,rhocritz,aux(k,25))
      aux(k,8)=totMass(178d0,r200,rc,rhocritz,aux(k,26))
      aux(k,9)=totMass(150d0,r200,rc,rhocritz,aux(k,27))
      call getGassMass(rc,beta,rhogas_central,r200,aux(k,10:14))
      aux(k,15)=aux(k,14)
      aux(k,14)=aux(k,13)
      aux(k,13)=Mgas200
      do i=16,21
      	if(aux(k,i-6)==0d0 .and. aux(k,i-12)==0d0) then
            	aux(k,i)=0d0
		else
            	aux(k,i)=aux(k,i-6)/aux(k,i-12)
		endif
      	if(.not.(aux(k,i)/h>=0d0 .and. aux(k,i)/h<=1d0)) then
            	flag=1
                  return
		endif
	enddo
      aux(k,28)=rhocritz
      
 999	return 
	end subroutine MakeGasDistributions

!=======================================================================

	subroutine MakeYMap(k)

! Generate a map of comptonisation parameter y for a given set of
! cluster parameters corresponding to the kth atom

	implicit none
	
	integer i,j,k	
	double precision x0,y0,rc,beta,yc
	double precision xx,yy,dx(2),rr,rold,angfactor
	double precision theta,polar,flux,angle,size
	double precision A1,rp1,fmin,fmax
	parameter(A1=1.424d-19,rp1=1.0)
	
!-----------------------------------------------------------------------

	x0=GeoPars(k,1)
	y0=GeoPars(k,2)

!------------------------------------------------------------------------
      
      if(GasModel==0) then
        	angfactor=1.0/60.0
      else
        	angfactor=sec2rad*D
      endif
      
!-----------------------------------

! Spherical geometry:
	
	if(GeoModel==1) then
	   	do i=1,nx
			do j=1,ny
		 		xx=trans(1)+i*trans(2)
		 		yy=trans(4)+j*trans(6)
		 		dx(1)=(xx-x0)
		 		dx(2)=(yy-y0)
		 		rr=sqrt(dx(1)*dx(1)+dx(2)*dx(2))*angfactor
		 		ymap(i,j)=ymap(i,j)+yfunc(rr)
	      	enddo
	   	enddo

!-----------------------------------

! Elliptical geometry: 

	elseif(GeoModel==2) then
         	   
	   	do i=1,nx
	      	do j=1,ny
		 		xx=trans(1)+i*trans(2)
		 		yy=trans(4)+j*trans(6)
		 		dx(1)=(xx-x0)
		 		dx(2)=(yy-y0)
		 		call QuadForm(2,dx,Q,dx,rr)
		 		rr=sqrt(rr)*angfactor
		 		ymap(i,j)=ymap(i,j)+yfunc(rr)
	      	enddo
	   	enddo
	endif
      
	end subroutine MakeYMap
		
!=======================================================================

	subroutine MakeVisibilities

! Convert a map of comptonisation parameter y to microwave sky
! brightness, fourier transform and sample at the relevant uv points.
! Might well be worth investigating more efficient fourier transform
! routines, if dataset is small.

	implicit none
	
	integer i,j,m
	double precision szsky(nx,ny),re1(nx,ny),im(nx,ny)
	double precision xx,yy,r2,sum,beam,taper,fu(large),fv(large),rr,ii
	double precision ll,mm,theta      
      double complex wkspce(nx*ny)
      double precision work(2*nx*ny)
      integer nn(2),k,nd
      integer job,iform	
	double precision zmin,zmax		
      character(len=100) file
!-----------------------------------------------------------------------
			
! Make sure map has zero mean, and convert to Janskys:
      
      sum=0.0
      do i=1,nx
         do j=1,ny
            sum=sum+ymap(i,j)
         end do
      end do
      sum=sum/float(nx*ny)	
	do i=1,nx
         do j=1,ny
            !szsky(i,j)=(ymap(i,j)-sum)*SZscale
            !!FF Hack
            !szsky(i,j)=(ymap(i,j)-sum)
            szsky(i,j)=(ymap(i,j))
         end do
      end do

! Loop over data files:

      do m=1,Nvisfiles
	  if(verbose) write(*,*) 'for pointing ',m,'...'

!  Apply primary beam 

        if(telescope(m)(1:3)=='VSA') then
          !beam=pb_vsa_ext
          beam=pb_sig(m)
          !lmax=lmax_vsa
          !lmin=lmin_vsa
	  elseif(telescope(m)(1:4)=='VSAs') then
	    !beam=pb_vsa_sea
          beam=pb_sig(m)
          !lmax=lmax_vsas
          !lmin=lmin_vsas
        elseif(telescope(m)(1:3)=='AM3') then
          beam=pb_ami(1)
          lmax=lmax_ami(1)
          lmin=lmin_ami(1)
        elseif(telescope(m)(1:3)=='AM4') then
          beam=pb_ami(2)
          lmax=lmax_ami(2)
          lmin=lmin_ami(2)
        elseif(telescope(m)(1:3)=='AM5') then
          beam=pb_ami(3)
          lmax=lmax_ami(3)
          lmin=lmin_ami(3)
        elseif(telescope(m)(1:3)=='AM6') then
          beam=pb_ami(4)
          lmax=lmax_ami(4)
          lmin=lmin_ami(4)
        elseif(telescope(m)(1:3)=='AM7') then
          beam=pb_ami(5)
          lmax=lmax_ami(5)
          lmin=lmin_ami(5)
        elseif(telescope(m)(1:3)=='AM8') then
          beam=pb_ami(6)
          lmax=lmax_ami(6)
          lmin=lmin_ami(6)
        elseif(telescope(m)(1:2)=='RT') then
          !beam=pb_ryl
          beam=pb_sig(m)
          !lmax=lmax_ryl
          !lmin=lmin_ryl
        elseif(telescope(m)(1:3)=='CBI') then
	    !beam=pb_cbi
          beam=pb_sig(m)
          !lmax=lmax_cbi
          !lmin=lmin_cbi
        endif
        if(verbose) then
        write(*,*) ' telescope=',telescope(m)(1:3)
        write(*,*) ' SZscale=',SZscale(m)
        write(*,*) ' beam size=',beam
        write(*,*) ' map centre (deg)=',map_x,map_y
        write(*,*) ' beam centre (deg)=',pb_x(m),pb_y(m)
        write(*,*) ' phase centre (deg)=',ph_x(m),ph_y(m)
        write(*,*) ' beam offset (sec)=',pb_x_offset(m),pb_y_offset(m)
        write(*,*) ' phase offset (sec)= ',ph_x_offset(m),ph_y_offset(m)
        write(*,*) ' model offset (sec)=',GeoPars(1,1),GeoPars(1,2)
        endif
        beam=2.0*beam*beam
!write(file,'(i,a)') m,'.ymap'
!open(2,file=file,form='formatted',status='replace')

     
        do j=1,ny
           do i=1,nx
           
!write(2,'(2i6,E20.6)')i,j,szsky(i,j)*SZscale(m)/ &
!	((cell*sec2rad)*(cell*sec2rad))* &
!	((clight/nu(1))**2.)*(10.**-26)/(kboltzmann*2.)

             xx=trans(1)+i*trans(2)-pb_x_offset(m)
             yy=trans(4)+j*trans(6)-pb_y_offset(m)
             r2=xx*xx+yy*yy
             taper=exp(-r2/beam)

             re1(i,j)=taper*szsky(i,j)*SZscale(m)
             im(i,j)=0.0

             if(verbose) then
             if(i==nx/2.and.j==ny/2) then
               write(*,*) 'szsky central pixel:',szsky(nx/2,ny/2)*SZscale(m)
               write(*,*) 'taper at this point=',taper
               write(*,*) 're map central pixel:',re1(nx/2,ny/2)
             endif
             endif

           enddo
        enddo
! close(2)
        
! FFT into aperture plane:
      
! Old code (profile) used NR FFT, painfully slow.
!
!        call xzapit(re,nx,ny)
!         nn(1)=nx
!         nn(2)=ny
!         k=1
!         do j=1,ny
!            do i=1,nx
!                work(k)=re(i,j)
!                k=k+1
!                work(k)=im(i,j)
!                k=k+1
!            enddo
!         enddo
!         call fourn2(work,nn,2,1)
!         k=1
!         do j=1,ny
!            do i=1,nx
!               re(i,j)=work(k)
!               k=k+1
!               im(i,j)=work(k)
!               k=k+1
!            end do
!         end do
!         call xzapit(re,nx,ny)
!         call xzapit(im,nx,ny)
 	
! New FFT code: Hobson's rfft.f code used, but only if the FFTW library
! is not available.

!      if(FFTW) then

        do j=1,ny
          do i=1,nx
            arr(i,j)=dcmplx(re1(i,j),im(i,j))
          end do
        end do
        call cheqboard(arr,nx,ny)
        call dfftw_execute(fftwplan)
        call cheqboard(arr,nx,ny)
        do j=1,ny
          do i=1,nx
            re1(i,j)=dble(arr(i,j))
            im(i,j)=aimag(arr(i,j))
          end do
        end do
!
!      else
!
!        nd=2
!        nn(1)=nx
!        nn(2)=ny
!        k=0
!        do j=1,ny
!          do i=1,nx
!            k=k+1
!            wkspce(k)=dcmplx(re1(i,j),im(i,j))
!          end do
!        end do
!        job=+1                   ! forward transform (-ve exponential)
!        iform=0                  ! data are real
!        call cheqboard(wkspce,nx,ny)
!        call fourt(wkspce,nn,nd,job,iform,work)
!        call cheqboard(wkspce,nx,ny)
!        k=0
!        do j=1,ny
!          do i=1,nx
!            k=k+1
!            re1(i,j)=dble(wkspce(k))
!            im(i,j)=aimag(wkspce(k))
!          end do
!        end do

!      endif
      
!  Now sample aperture plane at given u-v points, and rotate to phase
!  centre:
	
        ll=sin(-ph_x_offset(m)*sec2rad)*cos(-ph_y_offset(m)*sec2rad)
        mm=sin(-ph_y_offset(m)*sec2rad)

	
!      write(file,'(i,a)') m,'.vis'
!	open(2,file=file,form='formatted',status='replace')

        do k=1,Nvis(m)

          if(SZeflag(m,k)==0) then

            fu(k)=-1.0*u(m,k)
            fv(k)=1.0*v(m,k)
            
            call extract_visibility(re1,im,nx,ny,fu(k),fv(k),cell,rr,ii)
		
            theta=1.0*TwoPi*(-u(m,k)*ll+v(m,k)*mm)
            
!           Phase rotation by theta:
!           (vr+i*vi)*(cos+i*sin)=(vr*cos-vi*sin)+i*(vr*sin+vi*cos)

            pvisr(m,k)=dble(rr)*dble(cos(theta))-dble(ii)*dble(sin(theta))
            pvisi(m,k)=dble(rr)*dble(sin(theta))+dble(ii)*dble(cos(theta))
            
!            write(2,'(i4,2F10.3,2E14.4,F4.1)')k,u(m,k),v(m,k),pvisr(m,k),pvisi(m,k),0.0

          else
            pvisr(m,k)=0.d0
            pvisi(m,k)=0.d0
          endif

        enddo
!      close(2)
	
      enddo
!      stop
      
	return
	end subroutine MakeVisibilities
	
!=======================================================================

	subroutine AddSources

! Add a number of point sources to the predicted visibilities.

	implicit none
	
	integer i,j,m
	double precision ll,mm,ff,r2,beam,taper,ra,dec,sra,sdec,ra_off,dec_off,flux
	double precision theta,lambda
	
!-----------------------------------------------------------------------

! Loop over pointings

      do m=1,Nvisfiles

       if(telescope(m)(1:3)=='VSA') then
         !beam=pb_vsa_ext
         beam=pb_sig(m)
       elseif(telescope(m)(1:4)=='VSAs') then
         !beam=pb_vsa_sea
         beam=pb_sig(m)
       elseif(telescope(m)(1:4)=='AM3') then
         beam=pb_ami(1)
       elseif(telescope(m)(1:4)=='AM4') then
         beam=pb_ami(2)
       elseif(telescope(m)(1:4)=='AM5') then
         beam=pb_ami(3)
       elseif(telescope(m)(1:4)=='AM6') then
         beam=pb_ami(4)
       elseif(telescope(m)(1:4)=='AM7') then
         beam=pb_ami(5)
       elseif(telescope(m)(1:4)=='AM8') then
         beam=pb_ami(6)
       elseif(telescope(m)(1:2)=='RT') then
         !beam=pb_ryl
         beam=pb_sig(m)
       elseif(telescope(m)(1:3)=='CBI') then
	   !beam=pb_cbi
         beam=pb_sig(m)
	 else
         beam=pb_sig(m)
       endif
       beam=2.d0*beam*beam

       !lambda=clight/CRVAL4
       lambda=clight/nu(m)

! Loop over sources, and visibilities, calculating amplitude and phase
! factors at each uv point and applying these to the predicted data:
! Cannot use small angle approximation for VSA fields!

       do i=1,NSrc

!       Apply primary beam to source flux:

         ra=pb_x(m)*deg2rad
         dec=pb_y(m)*deg2rad
         sra=SrcPars(i,1)*deg2rad
         sdec=SrcPars(i,2)*deg2rad

         call calc_sepn(ra,dec,sra,sdec,r2)
         r2=r2*rad2sec
         r2=r2*r2
         taper=exp(-r2/beam)
         !ff=SrcPars(i,3)*taper
         !calculate the source flux at the channel frequency using spectral index
         flux=flux0(i)*((nu(m)/nu0(i))**-SrcPars(i,4))
         ff=flux*taper

!        Find direction cosines of source position vector
!        For small angles, ll \approx x/radians etc:
!        ll=SrcPars(i,1)*sec2rad
!        mm=SrcPars(i,2)*sec2rad
!        Otherwise, do offsets properly and use full formula:
!        (Note that it is relative to the phase centre we calculate)

         ra=ph_x(m)*deg2rad
         dec=ph_y(m)*deg2rad
         call calc_ra_dec_off(ra,dec,sra,sdec,ra_off,dec_off)
         ra_off=-1.0*ra_off

         ll=sin(ra_off)*cos(dec_off)
         mm=sin(dec_off)

!        Finally add the point source flux directly to all visibilities

         do j=1,NVis(m)
           theta=1.0*TwoPi*(-u(m,j)*ll+v(m,j)*mm)
           pvisr(m,j)=pvisr(m,j)+ff*cos(theta)
           pvisi(m,j)=pvisi(m,j)+ff*sin(theta)
         enddo

       enddo

      enddo 
       
	return
	end subroutine AddSources
	
!-----------------------------------------------------------------------

      subroutine extract_visibility(re_map,im_map,nx,ny,u,v,cell,re1,im)

! De-grids a visibility at specified u,v from an unconvolved gridded 
! aperture

      implicit none

	integer nx,ny
      double precision re_map(nx,ny),im_map(nx,ny)
      double precision u,v,cell,re1,im,uvdist,l
      double precision sec2rad
      parameter(sec2rad=4.8481368d-6)
      double precision centre,x,y,uvcell

      uvcell=1.0/(nx*cell*sec2rad)
      centre=(nx/2)+1
      x=(u/uvcell)+centre
      y=(v/uvcell)+centre

      if((x<0).or.(y<0).or.(x.gt.nx).or.(y.gt.ny)) then
         re1=0.0
         im=0.0
         write(*,*) 'extract_visibility: aperture plane overshoot.'
      else
         uvdist=sqrt(u*u+v*v)
         l=2.0*pi*uvdist
         if((l>lmax .or. l<lmin) .and. .false.) then
         	re1=0.
            im=0.
	   else
         	call interp2d(re_map,nx,ny,x,y,re1)
         	call interp2d(im_map,nx,ny,x,y,im)
  	   endif
      end if

      return
	end subroutine extract_visibility

!======================================================================


	function BetaModel2D(r,rc,beta,y0)
	
	implicit none
	
	double precision r,rc,beta,y0
	double precision BetaModel2D
	
	BetaModel2D=y0/((1.0+(r/rc)*(r/rc))**(3.0*beta/2.0-0.5))

	return
	end function BetaModel2D
      
!=======================================================================

	function BetaModelHM(r,rc,y0)
	
	implicit none
	
	double precision r,rc,y0,rt,bracket
	double precision BetaModelHM
	
	rt=3.0*rc
      bracket=1.0/sqrt(rc*rc+r*r)	
      bracket=bracket-1.0/sqrt(rt*rt+r*r)	
      BetaModelHM=bracket*y0*rc*rt/(rt-rc)
       
	return
	end function BetaModelHM
      
!=======================================================================

	function BetaModel3D(r,rc,beta,rhogas0)
	
	implicit none
	
	double precision r,rc,beta,rhogas0
	double precision BetaModel3D
	
	BetaModel3D=rhogas0/((1.0+(r/rc)*(r/rc))**(3.0*beta/2.0))
 
	return
	end function BetaModel3D
	
!=======================================================================

      function BetaMass(rc,beta,y0,rlimit1)
      
      implicit none
      
      double precision rc,beta,y0,rlimit1,result,BetaMass
      double precision eps
      parameter(eps=1d-4)
      
      a=rc
      b=beta
      c=y0
      call qtrap(Mgas_Beta_Integrand,0.d0,rlimit1,eps,result)
      
      BetaMass=result
      
      return
      end function BetaMass

!=======================================================================

      function GasMass(rlimit1,rc,beta)
      
      implicit none
      
      double precision rlimit1,GasMass,result,rc,beta
	double precision eps
	parameter(eps=1d-4)
	
      gd1=rc
      gd2=beta

      !call qtrap(Mgas_Integrand,0.d0,rlimit1,eps,result)
      call qtrap(betaProfileInt,0.d0,rlimit1,eps,result)
      GasMass=result
      
      return
      end function GasMass

!=======================================================================

      function Mgas_Beta_Integrand(r)
      
      implicit none
      
      double precision r,Mgas_Beta_Integrand
      double precision TwoPi
      parameter(TwoPi=6.283185307)
      
      Mgas_Beta_Integrand=TwoPi*r*c/((1.0+(r/a)*(r/a))**(3.0*b/2.0-0.5))
      
      return
      end function Mgas_Beta_Integrand

!=======================================================================

      function Mgas_Integrand(rr)
      
      implicit none
      
      double precision rr,Mgas_Integrand
      double precision t3
      
	if(rr<rmin) then 
	  t3=Rhogasfunc(rmin)
        Mgas_Integrand=4.0*Pi*rmin*rmin*t3
	else
	  t3=Rhogasfunc(rr)
        Mgas_Integrand=4.0*Pi*rr*rr*t3
	endif
	
      return
      end function Mgas_Integrand

!======================================================================

	function yfunc(rr)
	
	implicit none
	
	double precision rr,yfunc,result
	
	if(rr<rmin) then 
	  call interp1d(logYarray,logr,n,phlog10(rmin),result)
      else if(rr>rlimit) then
        yfunc=0.
        return
	else
	  call interp1d(logYarray,logr,n,phlog10(rr),result)
	endif
	yfunc=10.0**result
	
	return
	end function yfunc

!======================================================================

	function Rhogasfunc(rr)
	
	implicit none
	
	double precision rr,Rhogasfunc,result
	
	if(rr<rmin) then 
	  call interp1d(logRhogas,logr,n,phlog10(rmin),result)
	else if(rr>rlimit) then
        Rhogasfunc=0.
        return
	else
	  call interp1d(logRhogas,logr,n,phlog10(rr),result)
	endif
	Rhogasfunc=10.0**result
	
	return
	end function Rhogasfunc

!======================================================================

	function Pgasfunc(rr)
	
	implicit none
	
	double precision rr,Pgasfunc,result
	
	if(rr<rmin) then 
	  call interp1d(logPgas,logr,n,phlog10(rmin),result)
	else if(rr>rlimit) then
        Pgasfunc=0.
        return
      else
	  call interp1d(logPgas,logr,n,phlog10(rr),result)
	endif
	
	Pgasfunc=10.0**result
	
	return
	end function Pgasfunc

!======================================================================

	function PgasIntegrand(zz)
	
	implicit none
	
	double precision PgasIntegrand,zz,rr

	rr=uu*uu+zz*zz
	if(rr<rmin*rmin) then 
	  rr=rmin
	  PgasIntegrand=Pgasfunc(rr)
	else if(rr>rlimit*rlimit) then
        PgasIntegrand=0.
	else
	  rr=sqrt(rr)
	  PgasIntegrand=Pgasfunc(rr) 
	endif
	
	return
	end function PgasIntegrand

!======================================================================

	double precision function totMass(x,r200,rc,rhocrit,r)
	
	implicit none
	
	double precision x,r200,rc,rhocrit
      double precision r

	r=rx(x,r200,rc)
      totMass=(4d0*pi/3d0)*(r**3)*x*rhocrit
      
	end function totMass

!======================================================================

	double precision function rx(x,r200,rc)
	
	implicit none
	
	double precision x,r200,rc

	rx=sqrt(max(200d0*(r200**2+rc**2)/x-rc**2,0.d0))
      
	end function rx

!======================================================================

	subroutine getGassMass(rc,beta,rho0,r200,Mgas)
	
	implicit none
	
	double precision rc,beta,r200,rho0
      double precision Mgas(5) !Mgas value at r1500, r1000, r150, r178 & r500
      double precision rhi,rlo
	double precision eps
	parameter(eps=1d-4)
	
      gd1=rc
      gd2=beta
      
      !r1500
      rlo=0.d0
      rhi=rx(1500d0,r200,rc)
	call qtrap(betaProfileInt,rlo,rhi,eps,Mgas(1))
      Mgas(1)=Mgas(1)*rho0
      
      !r1000
      rlo=rhi
      rhi=rx(1000d0,r200,rc)
	call qtrap(betaProfileInt,rlo,rhi,eps,Mgas(2))
      Mgas(2)=Mgas(2)*rho0+Mgas(1)
      
      !r500
      rlo=rhi
      rhi=rx(500d0,r200,rc)
	call qtrap(betaProfileInt,rlo,rhi,eps,Mgas(3))
      Mgas(3)=Mgas(3)*rho0+Mgas(2)
      
      !r178
      rlo=rhi
      rhi=rx(178d0,r200,rc)
	call qtrap(betaProfileInt,rlo,rhi,eps,Mgas(4))
      Mgas(4)=Mgas(4)*rho0+Mgas(3)
      
      !r150
      rlo=rhi
      rhi=rx(150d0,r200,rc)
	call qtrap(betaProfileInt,rlo,rhi,eps,Mgas(5))
      Mgas(5)=Mgas(5)*rho0+Mgas(4)
      
	end subroutine getGassMass

!======================================================================

	double precision function betaProfileInt(r)
	
	implicit none
	
	double precision r
	
      if(r<rmin) then
		betaProfileInt=4d0*pi*(rmin**2)/((1d0+((rmin/gd1)**2))**(1.5d0*gd2))
	else
		betaProfileInt=4d0*pi*(r**2)/((1d0+((r/gd1)**2))**(1.5d0*gd2))
      endif
      
	end function betaProfileInt

!======================================================================

end module GasModels
