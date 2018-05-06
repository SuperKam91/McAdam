c=====================================================================
c
c                            McAdam
c
c  Bayesian parameterised analysis of astronomical data sets
c
c  Author:	Phil Marshall
c  Editors:	Nutan Rajguru, Joern Geisbuesch, Marko Velic
c  Uses:	BayeSys3 (Skilling 2003)
c 		Jelly (MacLachlan 2003)
c
c  History: 
c	8/2002	Original version, weak lensing + SZ effect
c	1/2003	Multiple atoms combined into single parameter set
c	3/2003	Multiple data files added
c	4/2003	Group release, CVS access
c	3/2004	FFTW added 
c	6/2004	More SZ beta model styles, convergence GL constraint 
c
c=======================================================================

      program McAdam 	
       
 	implicit none
 
	include 'McAdam.inc'
	
	integer index,Iseed,i
      integer JellyInit
	real*8 Evidence

c-----------------------------------------------------------------------
	
	if (bayesys) then 
	  index = 1  
	else
	  index = JellyInit()
	endif
	  
	if (bayesys) call Welcome
		
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	if (simulate) then 
	  call Initialise(index)	
	  call MakeData
	  stop
	endif

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	call ReadInData
c       if (HelpObj.or.HelpSrc) call HelperReadin
      
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      
	do i=1,Nruns
	
	if (index.ne.0) write(*,*) 
	if (index.ne.0.and.Nruns.gt.1) write(*,*) 'Starting run no. ',i	
	if (index.ne.0.and.Nruns.gt.1) write(*,*) 
	
	call Initialise(index)	
		
	if (plot.and.index.ne.0.and.i.eq.1) call PlotSetup
		 
	if (bayesys) then
	  Iseed = -1
	  call Bayes3_Init(Nwalks,127)
	  call Bayes3_Run(Rate,Iseed,Ndim,Minatoms,Maxatoms,
     &                  Alpha,Evidence)
      else
	  Iseed = -1
	  call JellyRun(Rate,Iseed,Ndim,Evidence,NWalks,Method)
	endif
      call CloseFiles
	 
	if (index.ne.0) write(*,*) 
	if (index.ne.0) write(*,*) 
	if (index.ne.0.and.Nruns.gt.1) write(*,*) '********************'
	
	enddo
		
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	if (bayesys) then
	  write(*,*) 'BayeSys3 finished normally.'
	  write(*,*) ' (Took ',Nsamples,' samples after burn-in.'
	  write(*,*) '  During burn-in, ',Hoops-Nsamples,' were accepted.'
	  write(*,*) '  Total no. of calls to UserBuild = ',Shots
	  write(*,*) '  Total no. of unphysical samples = ',Airballs
	  write(*,*) '  Total no. of rejected samples = ',
     &                Shots - Airballs - Hoops
	  write(*,*) '  Total no. of accepted samples = ',Hoops,')'
	else  
	  call JellyFinish
	endif
	if (plot) call pgend
	if (index.ne.0) then
	  write(*,*) 
	endif  

      end

c=======================================================================

      subroutine FUserbuild(Cool,Nsystem,Lhood,Nstore,idummy,ndummy,
     &                      Cube,retcode)
      
	implicit none

        include 'McAdam.inc'
	
	integer Nsystem,Nstore,idummy,ndummy,retcode,i,j,Uflag,flag
        real*8 Cool,Lhood,Cube(Ndim)
      
	real*8 GLLhood,SZLhood,dLhood

c-----------------------------------------------------------------------
	
	Uflag = 1
	flag = 0
	shots = shots + 1
	
	GLLhood = 0.d0
	SZLhood = 0.d0
	
c c Update helper functions as necessary:	
c       
c       if ((HelpObj.or.HelpSrc).and.HelpUpdate) then
c         call HelperInit(Cool)
c 	  HelpUpdate = .false.
c       endif  
      
c Additional prior constraints: check parameters and return both an
c additional contribution to the likelihood and a flag value (for
c truncated distributions)

	dLhood = 0.0d0
      call CheckPars(Cube,dLhood,flag)
      if (flag.eq.1) goto 999

c Skip L calculation if sampling prior:

      if (SamplePrior) goto 999

c Make predicted data:

	if (GL.eq.1.and.ModelClass.eq.1) call PredictShear(Cube,flag)
	if (SZ.eq.1.or.ModelClass.eq.2) call PredictVisibilities(Cube,flag)
	if (flag.eq.1) goto 999 

c Calculate Likelihood	

	if (GL.eq.1.and.ModelClass.eq.1) call ShearLhood(GLLhood)
	if (SZ.eq.1.or.ModelClass.eq.2) call VisibilitiesLhood(SZLhood)

	Lhood = GLLhood + SZLhood + dLhood

      if (verbose) write(*,'(a,e10.3,2x,i5,2x,a,e13.6,a)') 
     &'UserBuild: ',Cool,Shots,'Lhood = ',Lhood,'                      '

c Return
 999	if (flag.eq.0) then
	  retcode = 1
	else
	  retcode = 0
	  Lhood = B_A_D
	endif  

      return
      end

c=======================================================================

      subroutine FUsermonitor(Cool,Nsystem,Lhood,Nstore,idummy,ndummy,
     &                        Cube,Walk,retcode,Evidence)
      
	implicit none

      include 'McAdam.inc'
	
      integer Nsystem,Nstore,idummy,ndummy,Walk,retcode,i,j,Uflag
      real*8 Cool,Lhood,Cube(Ndim),Evidence,diff

      integer Nrem,tot_dim
      real rejects
	
c-----------------------------------------------------------------------

	if (verbose) write(*,*) 'UserMonitor called:'
	Uflag = 2

        if(ModelClass.eq.1)then
        tot_dim = NDim
        elseif(ModelClass.eq.2)then
        tot_dim = NDim2
        endif   

c Count successful samples:

	hoops = hoops + 1
	
c Examine Cooling progress:

	if (maximise) then
	
        Nrem = Nsystem +1
	  diff = Lhood - OldLhood
	  if (diff.gt.0.0) then
	    OldLhood = Lhood
	    if (diff.le.Aim) then
	      write(*,*) 'stopping criterion reached!'
		Nrem = 0
	      call OpenFiles 
	      do j=1,Ndim
  	        call Rescale(Cube(j),j)
		enddo  
	      call WritePars(Lhood,Evidence)	    
	    endif  
	  endif 
	  if (Cool.ge.1.d0) Uflag=1 
	
      else
      
	  if (Cool.ge.1.d0) Cool = 1.d0
	
      endif

c c Debugging - do 20 samples then stop burnin
c       if (hoops.ge.2) Cool = 1.d0

c If doing ordinary sampling, then get the process stage right:
c If sampling from prior, start writing out a.s.a.p.	
	
	if (SamplePrior.and.hoops.eq.Nwalks+1) then
	  Cool = 1.d0
	endif
	
	if (Cool.lt.1.d0.and.Burnin.eq.1) then
c       [During burn-in]	

	  Nrem = Nsystem*Nwalks
	
	elseif (Cool.eq.1.d0.and.Burnin.eq.1) then
c       [Burn-in has just finished]	
	
	  Burnin = 0
	  NS_afterburn = 0
	  Nrem = Nsamples - NS_afterburn
        write(*,*) 
	  write(*,*) 
	  write(*,*) 'Burn-in complete, now exploring posterior.'
	  call OpenFiles 
	
	elseif (Cool.eq.1.d0.and.Burnin.eq.0) then 
c       [After burn-in]	
	
	  NS_afterburn = NS_afterburn + 1
	  Nrem = Nsamples - NS_afterburn 
	  
	  do j=1,Ndim
	    call Rescale(Cube(j),j)
	  enddo
	  call WritePars(Lhood,Evidence)	    
        
	endif

c Log progress
	
	if (maximise) then
	
        write(*,'(a,e10.3,2x,i5,2x,a,e13.6,2x,a,e15.8)')
     &  'Progress: ',Cool,Nrem,'Lhood = ',Lhood,'Diff. = ',diff

      else

        write(*,'(a,e10.3,2x,i5,2x,a,e13.6,2x,a,e15.8)')
     &  'Progress: ',Cool,Nrem,'Lhood = ',Lhood,'Evidence = ',Evidence
c      &  'Progress: ',Cool,Nrem,'Lhood = ',Lhood,'Lphi = ',Lphi

	endif

	call rewindline

c c Calculate mean log lhood and helper functions
c 
c 	if (HelpObj.or.HelpSrc) then
c         do j=1,Ndim
c 	    call Rescale(Cube(j),j)
c 	  enddo
c 	endif
c       
c       if (NS_afterburn.le.Nwalks+1) then
c         if (Walk.eq.1) then
c           Lbar = Lhood
c           Lphibar = Lphi
c         elseif (Walk.lt.NWalks) then  
c           Lbar = (Lbar*float(Walk-1) + Lhood)/float(Walk)
c           Lphibar = (Lphibar*float(Walk-1) + Lphi)/float(Walk)
c         elseif (Walk.eq.NWalks) then
c           Lbar = (Lbar*float(Walk-1) + Lhood)/float(Walk)
c           Lphibar = (Lphibar*float(Walk-1) + Lphi)/float(Walk)
c           call RecordProgress(Cool)
c         else
c         endif  
c 	endif

c Plot progress

      if (plot) then
	  do j=1,Ndim
	    call Rescale(Cube(j),j)
	  enddo
	  call PlotPoints(Uflag,Lhood)
	endif
	
c c Change helper function next time UserBuild is called - only set this
c c flag in last call to UserMonitor:
c 
c       if ((HelpObj.or.HelpSrc).and.Walk.eq.NWalks) HelpUpdate = .true.

c Exit Bayesys3 if ensemble is large enough, or if desired maximisation 
c accuracy has been reached:

c c Debugging - do 40 samples then stop sampling
c       if (hoops.ge.4) Nrem = 0
      
	if (Nrem.eq.0.and.Burnin.ne.1) then
	  retcode = -1
	else
	  retcode = 0
	endif  
	  
c c Force exit after first sample accepted (debugging)
c 	if (Nrem.gt.20) retcode = -1 
 	  
	return
      end

c=======================================================================

      subroutine Jellybuild(Cool,Nsystem,Lhood,Nstore,ndummy,
     &                      Cube,retcode)
      
	implicit none

      include 'McAdam.inc'
	
	integer Nsystem,Nstore,idummy,ndummy,retcode,i,j,Uflag,flag
      real*8 Cool,Lhood,Cube(Ndim)

c-----------------------------------------------------------------------

	idummy = 1
	call FUserbuild(Cool,Nsystem,Lhood,Nstore,idummy,ndummy,
     &                      Cube,retcode)

      return
      end

c=======================================================================

      subroutine Jellymonitor(Cool,Nsystem,Lhood,Nstore,ndummy,
     &                        Cube,Walk,retcode,Evidence)
      
	implicit none

      include 'McAdam.inc'
	
      integer Nsystem,Nstore,idummy,ndummy,Walk,retcode
      real*8 Cool,Lhood,Cube(Ndim),Evidence
	
c-----------------------------------------------------------------------
      
      idummy = 1
      call FUsermonitor(Cool,Nsystem,Lhood,Nstore,idummy,ndummy,
     &                  Cube,Walk,retcode,Evidence)
 	
      return
      end

c=======================================================================

	subroutine PredictShear(Cube,flag)
	
	implicit none

      include 'McAdam.inc'
	
	integer i,j,k,flag
      real*8 Cube(Ndim)
	
	real*8 gsq,sigma2,lookupSigCrit
	
c-----------------------------------------------------------------------

c First zero the working arrays:

      do k=1,Ngals
	  gamma1(k) = 0.d0
	  gamma2(k) = 0.d0
	  kappa(k) = 0.d0
	enddo

c All parameters are contained in one "atom" in the
c Bayesys3 sense of the word. Atom is now taken to mean "cluster
c component". This should all be taken care of in the Rescale
c function, which is responsible for mapping the unit hypercube into
c the defined parameter set. This now happens once only, all the atoms
c being parameterised in the same function call.
c
c Note, Cube is now always a 1-dimensional array.
c
c Rescale parameters:

	do j=1,Ndim
	  call Rescale(Cube(j),j)
	enddo

c If we are fitting atoms,

      if (Atoms.eq.1) then

c For each atom,

c 	if (verbose) write(*,*) 'Calculating shears:'
	do i=1,Natoms

c       Sort out elliptical geometry:
        if (GeoModel.eq.2) then
          call EllGeometry(i)
	  endif
	    
c       Now increase the shear and convergence at each galaxy position:

	  call LensFields(i,flag)

	enddo

c 	if (verbose) write(*,*) 'Calculating reduced shears...'
	
      do k=1,Ngals

	  if (GLeflag(k).eq.0) then
	  
c       Now calculate reduced shear fields

	  g1(k) = gamma1(k)/(1.0-kappa(k))
	  g2(k) = gamma2(k)/(1.0-kappa(k))

c       Correct for strong lensing region
	
	  gsq = g1(k)*g1(k) + g2(k)*g2(k)
	  if (gsq.gt.1.0) then
	    g1(k) = g1(k)/gsq
	    g2(k) = g2(k)/gsq
	    gsq = 1.d0/gsq
        endif

c       Correct intrinsic errors for effect of lensing - does not seem
c       to work properly!

	  if (errcorrect) then
	    sigma2 = e1err(k)*e1err(k) - sigmaobs*sigmaobs
	    sigma2 = sigma2*(1.d0-gsq)*(1.d0-gsq)
	    wt1(k) = 1.d0/(sigma2 + sigmaobs*sigmaobs)
	    sigma2 = e2err(k)*e2err(k) - sigmaobs*sigmaobs
	    sigma2 = sigma2*(1.d0-gsq)*(1.d0-gsq)
	    wt2(k) = 1.d0/(sigma2 + sigmaobs*sigmaobs)	
	  endif
	  
	  endif
	  
	enddo	

	endif
      
c 	if (verbose) write(*,*) 'Finished.'

      return
      end

c=======================================================================

	subroutine PredictVisibilities(Cube,flag)
	
	implicit none

        include 'McAdam.inc'
	
	integer i,j,m,k,flag,tot_dim
        real*8 Cube(*)
	
c-----------------------------------------------------------------------

        if(ModelClass.eq.1)then
        tot_dim = NDim
        elseif(ModelClass.eq.2)then
        tot_dim = NDim2
        endif

c-----------------------------------------------------------------------

c If we are fitting atoms, first set y map to zero:

      if (Atoms.eq.1) then
        do i=1,nx
          do j=1,ny
            ymap(i,j) = 0.0
          enddo
        enddo
      endif
      	
c     Rescale parameters:

	do j=1,tot_dim
	  call Rescale(Cube(j),j)
	enddo

c If we are fitting atoms,

      if (Atoms.eq.1) then

c For each atom,

        do i=1,Natoms

c       Sort out elliptical geometry:

          if (GeoModel.eq.2.and.ModelClass.eq.1) then
            call EllGeometry(i)
          endif

c       Generate working arrays (r,T,rhogas etc):

          if (verbose) write(*,*) 'Making gas distribution ',i,'...'
          call MakeGasDistributions(i,flag)
          if (verbose) write(*,*) '...done'
          if (flag.eq.1) goto 999

c       Generate Comptonisation parameter map:

          if (verbose) write(*,*) 'Updating ymap...'
          call MakeYMap(i)
          if (verbose) write(*,*) '...done'

        enddo

c     Scale, Fourier transform and sample to get visibilities:

        if (verbose) write(*,*) 'Making visibilities...'
        call MakeVisibilities
        if (verbose) write(*,*) '...done'

c Otherwise, zero all visibility arrays:

      else

        if (verbose) write(*,*) 'Making visibilities...'
        do m=1,Nvisfiles
          if (verbose) write(*,*) 'for pointing ',m,'...'
          do k=1,Nvis(m)
            if (SZeflag(m,k).eq.0) then
              pvisr(m,k) = 0.d0
              pvisi(m,k) = 0.d0
            endif
          enddo
        enddo

      endif
      
c     Add point sources if required:

	if (Nuisance.eq.1.and.SourceSubtract.eq.1) call AddSources
      
 999	return
      end

c=======================================================================

	subroutine ShearLhood(GLLhood)
	
	implicit none

      include 'McAdam.inc'
	include 'constants.inc'
	
	integer k
      real*8 GLLhood,Chisq
	
	Chisq = 0.d0
	if (errcorrect) then
	  GLLhood0 = 0.d0
	endif	     
	
	do k=1,Ngals
	  if (GLeflag(k).eq.0) then
	    Chisq = Chisq + (e1(k)-g1(k))*wt1(k)*(e1(k)-g1(k))
	    Chisq = Chisq + (e2(k)-g2(k))*wt2(k)*(e2(k)-g2(k))
	    if (errcorrect) then
	      GLLhood0 = GLLhood0 + 0.5d0*dlog(wt1(k))
	      GLLhood0 = GLLhood0 + 0.5d0*dlog(wt2(k))
	    endif	     
	  endif	     
	enddo
		
	if (errcorrect) then
	  GLLhood0 = GLLhood0 - ((2*Ngals)/2)*dlog(TwoPi)
	endif	     
	
	GLLhood = GLLhood0 - Chisq/2.0d0	  
		
	return
	end
	
c=======================================================================
c Note: when move was made to multiple covariance matrices, the memory
c usage grew considerably, leading (I think) to a bug in this part of
c the code where the b vector was being calculated but then replaced
c with INFs and NANs when compiled with g77. Using the -fno-automatic
c seems to have fixed this. 

	subroutine VisibilitiesLhood(SZLhood)
	
	implicit none

        include 'McAdam.inc'
	
	integer i,j,k,m
        real*8 SZLhood,Chisq
	real*8 a(2*large),b(2*large),row(2*large),sum

c-----------------------------------------------------------------------
	
	Chisq = 0.d0

c Approximation: CMB=0 for Ryle Telescope dishes. Only apply
c full covariance matrix if IncludeCMB(m)=1, 
c otherwise just use noise values in the data file.
	
c Loop over data files:      
     
      do m=1,Nvisfiles

	  if (verbose) write(*,*) 'Calculating Lhood for file ',m,'...'
     
        if (IncludeCMB(m).eq.1) then

c         Compute chi-squared as the dot product b^T b, where
c         Lb = a and a is the vector of residuals:
	  	            
          do i=1,Nvis(m)
            a(i) = visr(m,i)-pvisr(m,i)
          enddo
          do j=Nvis(m)+1,2*Nvis(m)
            i = j - Nvis(m)
            a(j) = visi(m,i)-pvisi(m,i)
          enddo
	  
          k = 1
          do i=1,2*Nvis(m)

            do j=1,i
              row(j) = LCM(m,k)
              k = k + 1
            enddo

            if (i.eq.1) then
              b(1) = a(1)/row(1)
            else
              sum = 0.d0
              do j=1,i-1
                sum = sum + row(j)*b(j)
              enddo
              b(i) = (a(i)-sum)/row(i)
            endif
          
          enddo

          do i=1,2*Nvis(m)
            Chisq = Chisq + b(i)*b(i)
          enddo
	
	  else
	
c 	  Compute chisquared for diagonal thermal noise only covariance
c       matrix, ie the usual independent Gaussian errors routine:
	
 10       do k=1,Nvis(m)
            if (SZeflag(m,k).eq.0) then
              Chisq = Chisq + (visr(m,k)-pvisr(m,k))*viswt(m,k)*
     &                        (visr(m,k)-pvisr(m,k))
              Chisq = Chisq + (visi(m,k)-pvisi(m,k))*viswt(m,k)*
     &                        (visi(m,k)-pvisi(m,k))
            endif
          enddo

        endif
	
      enddo
      	
	SZLhood = SZLhood0 - Chisq/2.0d0
		
	return
	end
	
c=======================================================================

	subroutine Rescale(value,ival)
	
	implicit none

        include 'McAdam.inc'
	
	integer ival,i,ii,j,ih,jh,tot_pars
	logical nuis
	real*8 value,Prior

c----------------------------------------------------------------------

        if(ModelClass.eq.1)then
        tot_pars = NPars
        elseif(ModelClass.eq.2)then
        tot_pars = NPars2
        endif	

c----------------------------------------------------------------------

c First determine whether 1D coordinate ival refers to an atomic or 
c a nuisance parameter:	
	
      if (ival.gt.Atoms*NAtoms*tot_pars) then
        nuis = .true.
      else
        nuis = .false.
      endif
         
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	if (.not.nuis.and.ModelClass.eq.1) then
		
c     Rescale the parameters of the atom, labelled j = 1..NAtoms;
c     i runs from 1 to NPars:

           j = (ival-1)/NPars + 1
           i = ival - (j-1)*NPars
	
           if (i.le.NGeoPars) then	

              ii = i	 

c 	  if (HelpObj.and.ii.eq.1) then
c c         Convert first dimension sample to xy position:
c           call HelperGetij(value,ih,jh)
c           call HelperGetlogphi(ih,jh)
c           call HelperGetxy(ih,jh,GeoPars(j,1),GeoPars(j,2))
c 	  elseif (HelpObj.and.ii.eq.2) then 
c c         Ignore second dimension sample:
c           goto 999
c         else 
	    
              GeoPars(j,ii) = Prior(Geo_PriorType(j,ii),value,
     &             Geo_Prior(j,ii,1),Geo_Prior(j,ii,2))
c 	  endif
              
           elseif (Mass.eq.1.and.
     &           i.le.(NGeoPars+NMassPars*Mass)) then
	  
              ii = i - NGeoPars	 

              MassPars(j,ii) = Prior(Mass_PriorType(j,ii),value,
     &             Mass_Prior(j,ii,1),Mass_Prior(j,ii,2))

           elseif (Gas.eq.1.and.
     &             i.le.(NGeoPars+NMassPars*Mass+NGasPars*Gas)) then

              ii = i - NGeoPars - NMassPars*Mass	 

              GasPars(j,ii) = Prior(Gas_PriorType(j,ii),value,
     &                     Gas_Prior(j,ii,1),Gas_Prior(j,ii,2))

           elseif (Gas.eq.1.and.
     &             i.le.(NGeoPars+NMassPars*Mass+
     &             (NGasPars+NTPars)*Gas)) then

              ii = i - NGeoPars - NMassPars*Mass - NGasPars*Gas	 

              TPars(j,ii) = Prior(T_PriorType(j,ii),value,
     &             T_Prior(j,ii,1),T_Prior(j,ii,2))

           endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 
        elseif (.not.nuis.and.ModelClass.eq.2) then

c     Rescale the parameters of the atom, labelled j = 1..NAtoms;
c     i runs from 1 to NPars:

	j = (ival-1)/tot_pars + 1
	i = ival - (j-1)*tot_pars

        if (i.le.NStringGeoPars) then	

           ii = i	 

           StringGeoPars(ii) = Prior(StringGeo_PriorType(ii),value,
     &          StringGeo_Prior(ii,1),StringGeo_Prior(ii,2))
	        
        elseif (i.le.(NStringGeoPars+NStringPars)) then	

           ii = i - NGeoPars 

           StringPars(ii) = Prior(String_PriorType(ii),value,
     &          String_Prior(ii,1),String_Prior(ii,2))
    
	endif

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      elseif (nuis) then
		
c       Rescale the nuisance parameters:

         i = ival - Atoms*NAtoms*NPars
         j = (i-1)/3 + 1
        
         if (SourceSubtract.eq.1.and.(i/3).le.NSrc) then	
	    
	    ii = i - (j-1)*3 

c           if (HelpSrc.and.ii.eq.1) then
c c         Convert first dimension sample to xy position:
c             call HelperGetij(value,ih,jh)
c             call HelperGetlogphi(ih,jh)
c             call HelperGetxy(ih,jh,SrcPars(j,1),SrcPars(j,2))
c           elseif (HelpSrc.and.ii.eq.2) then
c c         Ignore second dimension sample:
c             goto 999
c           else

            SrcPars(j,ii) = Prior(Src_PriorType(j,ii),value,
     &                        Src_Prior(j,ii,1),Src_Prior(j,ii,2))	
c           endif
        
         elseif (SourceSubtract.eq.1.and.Varyzs.eq.1)then

	    ii = i - NSrc*3	 
	    zs = Prior(zs_PriorType,value,zs_Prior(1),zs_Prior(2))	
            
         elseif (SourceSubtract.eq.0.and.Varyzs.eq.1)then

	    ii = i 
	    zs = Prior(zs_PriorType,value,zs_Prior(1),zs_Prior(2))	

         endif

      endif

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

 999  return
      end

c=======================================================================

	subroutine Initialise(index)
	
	implicit none

        include 'McAdam.inc'
	include 'fftw.inc'
	include 'constants.inc'
	
	integer i,j,k,m,index,iend,NSZcovmats,Wastage
	real*8 row(2*large),SZLhood,GLLhood,Lhood
	real*8 y2Jy
	real zsmin,zsmax,makesigcrit
        character*100 string
        character*25 ModClassName
        character*25 GeoModelName,MassModelName,GasModelName,TModelName
        character*15 engine

c-----------------------------------------------------------------------
        
c Write out to screen status of parameters


        if(ModelClass.eq.1)then
           ModClassName = 'Clusters'
 
           if(GeoModel.eq.1)then
              GeoModelName = 'Circular symmetry'
           elseif(GeoModel.eq.2)then
              GeoModelName = 'Elliptical symmetry'
           endif

           if(MassModel.eq.1)then
              MassModelName = 'NFW'
           elseif(MassModel.eq.2)then
              MassModelName = 'SIS'
           elseif(MassModel.eq.3)then
              MassModelName = 'Cored Power Law'
           endif
           
           if(GasModel.eq.0)then
              GasModelName = 'H+MO2 survey beta model'
           elseif(GasModel.eq.1)then
              GasModelName = 'Beta Model'
           elseif(GasModel.eq.2)then
              GasModelName = 'Hydrostatic equilibrium'
           endif
                
           if (TModel.eq.1) then
              TModelName = 'Isothermal'
           elseif (TModel.eq.2) then
              TModelName = 'Polytropic'
           endif

        elseif(ModelClass.eq.2)then
        ModClassName = 'Cosmic Strings'
        endif   

        if(bayesys)then
           engine='bayesys'
        else
           engine='jelly'
        endif

        if(simulate)then
        write(*,'(a,I1,1x,a)')
     &  ' In simulation mode. Generating mock data sets with ', NAtoms,
     &  ModClassName
        write(*,*)
       elseif(.not.simulate)then
        write(*,'(a,I1,1x,a)') 
     &  ' In inference mode. Recovering parameters of ', NAtoms, 
     &   ModClassName
        write(*,*)  

        if(ModelClass.eq.1)then
        write(*,'(a27,a)')
     &  ' The geometrical model is ',GeoModelName

              if(Mass.eq.1)then
              write(*,'(a,a)')
     &        ' The mass model is ', MassModelName        
              else
              write(*,'(a)') 'Not fitting mass content'   
              endif
        
              if(Gas.eq.1)then
              write(*,'(a18,a)')
     &        ' The gas model is ', GasModelName
              else
              write(*,'(a30)') ' Not fitting gas content'
              endif

        write(*,'(a,a)')
     &  ' The temperature model is ', TModelName

      endif

      endif

        if(Nuisance.eq.1)then 
        write(*,'(a,i1,1x,a)') 
     &  'The dataset contains ', NSrc, 'point sources.'
        else
        write(*,'(a)') ' The data contains no point sources'   
        endif

        if(ModelClass.eq.1.and..not.simulate)then

c           do i=1,Nvisfiles
c           write(*,'(/,a,a)') 
c     &     ' The dataset used is called ', visdatafile(i)
c           enddo
c           write(*,'(a,a)')
c     &     ' The covariance matrix is ', covmatfile
c               
c           if(GL.eq.1)then
c           write(*,'(a,a)')
c     &     ' The lensing datafile is ', GLdatafile
c           endif

        elseif(ModelClass.eq.2.and..not.simulate)then
           
           do i=1,Nvisfiles
           write(*,'(/,a,a)')
     &     ' The dataset used is called ', visdatafile(i)
           enddo
           write(*,'(a,a)')
     &     ' The covariance matrix is ', covmatfile

        endif

        if(.not.simulate)then
         write(*,'(/,a20,a7,a14,i2,a8,i5,a,/,a,f6.3,/)') 
     &  'The MCMC engine is ',engine,'. It is using ', Nwalks, 
     &  ' chains, ', Nsamples,' samples and has a burn in rate of ', 
     &    Rate
         endif
c-----------------------------------------------------------------------	
c MCMC and general stuff:
	
	if (MinAtoms.ne.MaxAtoms) then
	  write(*,*) '  Cannot handle variable atom numbers... stopping'
	  stop
	endif

	if (maximise) then
	  Burnin = 2
	  OldLhood = -1.d32
	else  
	  Burnin = 1
	endif
	
	shots = 0	
	hoops = 0	
	bricks = 0	
	airballs = 0	

	call Prior2Ranges('Geometry')
	
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

c Gravitational lensing stuff:
	
	if (Mass.eq.1.and.ModelClass.eq.1) then
	  call Prior2Ranges('Mass')
	endif
	  
	if (GL.eq.1) then
	
	  if (Varyzs.eq.1) then
	    call Prior2Ranges('SigCrit')
	    zsmin = z
	    zsmax = 6.0
	    do i=1,n
	      logzs(i) = log10(zsmin) + 
     &		(log10(zsmax)-log10(zsmin))*(i-1) / float(n-1)
            zs = 10.d0**logzs(i)
	      logSigmaCrit(i) = log10(makesigcrit(1.d0*z,zs))
c 		write(*,*) i,z,zs,logzs(i),makesigcrit(z,zs),
c      &                 10.0**(logSigmaCrit(i)),logSigmaCrit(i)
	    enddo
c 	    pause	
	  endif
	  
	  do i=1,Ngals
          gamma1(i) = 0.d0
          gamma2(i) = 0.0d0
          kappa(i) = 0.d0
          g1(i) = 0.d0
          g2(i) = 0.d0          
        enddo

	  GLLhood0 = 0.d0
	  if (simulate) goto 5
	  
	  do k=1,Ngals
	    if (GLeflag(k).eq.0) then
	      GLLhood0 = GLLhood0 - dlog(e1err(k))
	      GLLhood0 = GLLhood0 - dlog(e2err(k))
	    endif
	  enddo
      
	  GLLhood0 = GLLhood0 - ((2*Ngals)/2)*dlog(TwoPi)

	  if (index.ne.0) write(*,*) '       GLLhood0 = ',GLLhood0	 
        if (index.ne.0) write(*,*)
	  
	endif
	
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

c SZ stuff:
	
 5	if (Gas.eq.1.and.ModelClass.eq.1) then
           call Prior2Ranges('Gas')
           call Prior2Ranges('Temperature')
        elseif (ModelClass.eq.2)  then
           call Prior2Ranges('StringObject')
	endif

	if (SourceSubtract.eq.1) then
	  call Prior2Ranges('Source')
	endif
	  
	if (SZ.eq.1.or.ModelClass.eq.2) then
	
	  do m=1,NVisfiles 
          if (Nvis(m).gt.large) then
            write(*,*) 'Error: not enough memory allocated for SZ data;'
            write(*,*) '  Nvis(',m,') = ',Nvis(m)
            write(*,*) '        large = ',large
            write(*,*) 'Increase value of large in .inc file and '
            write(*,*) '  recompile, but note that total memory'
            write(*,*) '  allocated is ~2*Nvisfiles*large^2*8 bytes.'
            stop
          endif
        enddo  

	  SZLhood0 = 0.d0
	  
        do m=1,Nvisfiles
          do i=1,Nvis(m)
            pvisr(m,i) = 0.0d0
            pvisi(m,i) = 0.0d0
          enddo
        enddo
        do i=1,nx
          do j=1,ny
            ymap(i,j) = 0.0
          enddo
        enddo

	  if (simulate) goto 10

c       Loop over all data files:

        NSZcovmats = 0
        Wastage = 0
        do m=1,Nvisfiles

          if (IncludeCMB(m).eq.1) then
	            
c         Read in Cholesky decomposition of covariance matrices, one row at
c         a time, and compute the determinant by the product of the
c         diagonal values. Also store the matrix for future use:
	  
	    open(unit=9,file=covmatfile(m),form='unformatted',status='old')
	    read(9) i
          if (i.ne.2*Nvis(m)) then
            write(*,*)  'Array size mismatch - expecting ',2*Nvis(m)
            write(*,*)  '                      received  ',i
            write(*,*)  'Possible causes: '
            write(*,*)  '  1) vis.fits / .LCM file mismatch'
            write(*,*)  '  2) byte-swapped .LCM file'
            stop
          endif
          k = 1
          do i=1,2*Nvis(m)
            read(9) (row(j),j=1,i)
            SZLhood0 = SZLhood0 + dlog(row(i))
            do j=1,i
              LCM(m,k) = row(j)
              k = k + 1
c             write(*,*) m,i,j,row(j)
c             read(*,*)
            enddo
          enddo
          close(9)

          if (bayesys) then
            string = covmatfile(m)
            do iend=100,1,-1
              if (string(iend:iend).ne.' ') goto 6
            enddo
 6          i = nint(Nvis(m)*(2*Nvis(m)+1)*8.0/1048576.0)
            write(*,'(a,i3,a)') ' Read in covariance matrix (',i,' Mb)'
            write(*,*) ' from file ',string(1:iend)
          endif

          NSZcovmats = NSZcovmats + 1
          j = nint(large*(2*large+1)*8.0/1048576.0)
          Wastage = Wastage + (j-i)

          SZLhood0 =  -((2*Nvis(m))/2)*dlog(TwoPi) - 0.5d0*SZLhood0

	    else

c         Use noise values in visibility data file in diagonal likelihood
	  
 7          do k=1,Nvis(m)
              if (SZeflag(m,k).eq.0) then
                SZLhood0 = SZLhood0 - dlog(visrms(m,k))
              endif
            enddo

            j = nint(large*(2*large+1)*8.0/1048576.0)
            Wastage = Wastage + j

            SZLhood0 =  SZLhood0 - ((2*Nvis(m))/2)*dlog(TwoPi)
      
          endif
          
       enddo   
       
       if (NSZcovmats.gt.0) then
         i = Nvisfiles*nint(large*(2*large+1)*8.0/1048576.0)
         write(*,'(a,i3,a)') ' Memory allocated = ',i,' Mb'
         write(*,'(a,i3,a)') '  (',Wastage,' Mb left unused).'
         write(*,*)
       endif   

	 if (index.ne.0) write(*,*) '       SZLhood0 = ',SZLhood0	 

c Gridding setup:
	
 10    trans(1)=-float(nx/2+1)*cell
       trans(2)=cell
       trans(3)=0.
       trans(4)=-float(ny/2+1)*cell
       trans(5)=0.
       trans(6)=cell
	
       uvcell = 1.0/(nx*cell*sec2rad)
       uvtrans(1)=-float(nx/2+1)*uvcell
       uvtrans(2)=uvcell
       uvtrans(3)=0.
       uvtrans(4)=-float(ny/2+1)*uvcell
       uvtrans(5)=0.
       uvtrans(6)=uvcell

c      Compute SZ conversion factor:

       if(ModelClass.eq.1)then
       SZscale = y2Jy(cell)
       elseif(ModelClass.eq.2)then
       SZscale = 0.001
       endif 

c      Set FFTW to work finding the most efficient transform...

        if (FFTW) then
          call sfftw_plan_dft_2d(fftwplan,nx,ny,arr,arr,
     &                         FFTW_BACKWARD,FFTW_MEASURE)
        endif
	
      endif

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      
c Initialise working arrays:

	if (Gas.eq.1.and.ModelClass.eq.1) then
	
	  do i=1,n
	    logr(i) = log10(rmin) + (log10(rmax)-log10(rmin))*(i-1) / 
     &                                                float(n-1)
	    r(i) = 10.0**logr(i)
	    Phi(i) = 0.0
	    Rhogas(i) = 0.0
	    T(i) = 0.0
	    logPhi(i) = 0.0
	    logRhogas(i) = 0.0
	    logT(i) = 0.0
	  enddo
	
	endif
	
	if (.not.simulate) then
	  if (index.ne.0) 
     &    write(*,*) '         Lhood0 = ', SZLhood0 + GLLhood0	 
	  GLLhood = 0.d0
	  SZLhood = 0.d0
        if (GL.eq.1) call ShearLhood(GLLhood)
	  if (SZ.eq.1) call VisibilitiesLhood(SZLhood)
	  NullEv = GLLhood + SZLhood
	  if (index.ne.0) write(*,*) '  Null Evidence = ',NullEv	 
	endif

c       call RecordProgress(0.0d0)
	
	return
	end

c=======================================================================

	subroutine Prior2Ranges(ParType)

	implicit none
	
	include 'McAdam.inc'
	
	character*(*) ParType
	integer i,j
	real Range
	external Range
	
c Loop over NAtoms, setting plotting ranges:	
	
	do j=1,NAtoms
	
	if (ParType.eq.'Geometry') then
	  do i=1,NGeoPars
	    Geo_Ranges(j,i,1) = Range(1,Geo_PriorType(j,i),
     &                              Geo_Prior(j,i,1),Geo_Prior(j,i,2))
	    Geo_Ranges(j,i,2) = Range(2,Geo_PriorType(j,i),
     &                              Geo_Prior(j,i,1),Geo_Prior(j,i,2))
	  enddo
	elseif (ParType.eq.'Mass') then
	  do i=1,NMassPars
	    Mass_Ranges(j,i,1) = Range(1,Mass_PriorType(j,i),
     &                            Mass_Prior(j,i,1),Mass_Prior(j,i,2))
	    Mass_Ranges(j,i,2) = Range(2,Mass_PriorType(j,i),
     &                            Mass_Prior(j,i,1),Mass_Prior(j,i,2))
	  enddo
	elseif (ParType.eq.'Gas') then
	  do i=1,NGasPars
	    Gas_Ranges(j,i,1) = Range(1,Gas_PriorType(j,i),
     &                              Gas_Prior(j,i,1),Gas_Prior(j,i,2))
	    Gas_Ranges(j,i,2) = Range(2,Gas_PriorType(j,i),
     &                              Gas_Prior(j,i,1),Gas_Prior(j,i,2))
	  enddo
	elseif (ParType.eq.'Temperature') then
           do i=1,NTPars
              T_Ranges(j,i,1) = Range(1,T_PriorType(j,i),
     &             T_Prior(j,i,1),T_Prior(j,i,2))
              T_Ranges(j,i,2) = Range(2,T_PriorType(j,i),
     &             T_Prior(j,i,1),T_Prior(j,i,2))
           enddo

        elseif (ParType.eq.'StringGeometry') then
           do i=1,NStringGeoPars
              StringGeo_Ranges(i,1) = Range(1,StringGeo_PriorType(i),
     &             StringGeo_Prior(i,1),StringGeo_Prior(i,2))
              StringGeo_Ranges(i,2) = Range(2,StringGeo_PriorType(i),
     &             StringGeo_Prior(i,1),StringGeo_Prior(i,2))
           enddo
	elseif (ParType.eq.'StringObject') then
           do i=1,NStringPars
              String_Ranges(i,1) = Range(1,String_PriorType(i),
     &             String_Prior(i,1),String_Prior(i,2))
              String_Ranges(i,2) = Range(2,String_PriorType(i),
     &             String_Prior(i,1),String_Prior(i,2))
           enddo

	endif
	
      enddo
	
c Now deal with nuisance parameters:	
	
	if (ParType.eq.'Source') then
	
        do j=1,NSrc
	    do i=1,3
            Src_Ranges(j,i,1) = Range(1,Src_PriorType(j,i),
     &                                Src_Prior(j,i,1),Src_Prior(j,i,2))
	      Src_Ranges(j,i,2) = Range(2,Src_PriorType(j,i),
     &                                Src_Prior(j,i,1),Src_Prior(j,i,2))
	    enddo
        enddo
	
      elseif (ParType.eq.'SigCrit') then
	  zs_Ranges(1) = Range(1,zs_PriorType,zs_Prior(1),zs_Prior(2))
	  zs_Ranges(2) = Range(2,zs_PriorType,zs_Prior(1),zs_Prior(2))
	endif
	
	return
	end

c=======================================================================
c Return the lower (lim=1) or upper (lim=2) plotting limit corresponding
c to the prior denoted by flag with characteristic values x1 and x2.

	function Range(lim,flag,x1,x2)

	implicit none
	
	integer lim,flag
	real*8 x1,x2,x(2)
	real Range
	integer m
	parameter(m=5)
	
	if (lim.ne.1.and.lim.ne.2) stop 'Illegal range limits requested.'
	x(1) = x1
	x(2) = x2

c Delta function prior:
	if (flag.eq.0) then
	  if (lim.eq.1) then
	    Range = 0.9*x(1) - 1.0
	  elseif (lim.eq.2) then
	    Range = 1.1*x(1) + 1.0
	  endif	    
c Uniform prior:
	elseif (flag.eq.1) then
	  Range = x(lim)
c Uniform prior in log:
	elseif (flag.eq.2) then 
	  Range = x(lim)
c Gaussian prior - plot between +/- m sigma:
	elseif (flag.eq.3) then 
	  if (lim.eq.1) then
	    Range = x(1) - m*x(2)
	  elseif (lim.eq.2) then
	    Range = x(1) + m*x(2)
	  endif	    
c LogNormal prior - plot between max/10 and max*4:
	elseif (flag.eq.4) then 
	  if (lim.eq.1) then
	    Range = x(1)/10.0
	  elseif (lim.eq.2) then
	    Range = x(1)*4.0
	  endif	    
c Sinusoidal prior:
	elseif (flag.eq.5) then 
	  Range = x(lim)
c Cauchy prior - plot between +/- m sigma:
	elseif (flag.eq.6) then 
	  if (lim.eq.1) then
	    Range = x(1) - m*x(2)
	  elseif (lim.eq.2) then
	    Range = x(1) + m*x(2)
	  endif	    
	endif
		
	return
	end

c=======================================================================

	subroutine Welcome
	
	write(*,'(a)')
     &'****************************************************************'
	write(*,'(a)')
	write(*,'(a)')
     &'                           McAdam'
	write(*,'(a)')
	write(*,'(a)')
     &'     Bayesian parameterised fitting of astronomical datasets.'
	write(*,'(a)')
	write(*,'(a)')
     &'               P.J. Marshall et al (April 2003)'

	write(*,'(a)')
	write(*,'(a)')
     &'****************************************************************'
	write(*,'(a)')
	
      return
      end

c=======================================================================
