module mass_stats
	use params
      use randomNS
      use lensing
      use MassModels
      use fits_utils
      use ReadWrite
	
      
contains

	subroutine readPostFiles(nPar,nCls,nullZ,baseroot,locZ,nParCls,params,bestfit)
      
      	!input variable
            integer nCls !total number of modes found
            double precision nullZ !log of null evidence
            integer npar !total no. of parameters to read
            character*100 baseroot
            
            !output variable
            double precision locZ(:) !difference of the local evidence from the null evidence
            integer nParCls(:) !total no. of posterior samples per cluster
            double precision params(:,:,:) !cumulative prob, log-like & posterior parameters
            double precision bestfit(:,:,:) !mean, max-like & map params of each cluster
            
            !work variables
      	character(len=32) fmt
            integer idummy(nCls)
            double precision rdummy(nCls),d1,maxl,maxp,sigma(nCls,nPar),temp(nPar*2)
            integer i,j,k,indx(1)
      
      
      	!read in the local evidences
      	open(11,file=TRIM(baseroot)//'p.dat',form='formatted',status='old')
      	write(fmt,'(a,i4,a,i4,a)')  '(',nCls,'i4,',2*nCls+2,'E20.12)' 
     		read(11,fmt) idummy(1:nCls),rdummy(1),locZ(1:nCls),rdummy(1),rdummy(1:nCls)
            write(fmt,'(a,i4,a)')  '(',nCls,'i6)'
            read(11,fmt)nParCls(1:nCls)
            write(fmt,'(a,i4,a)')  '(',nCls,'E20.12)'
            read(11,fmt)locZ(1:nCls)
      	close(11)
            
            !now read the separated posterior samples file
            open(5,file=TRIM(baseroot)//'post_separate_comm.dat',form='formatted',status='old')
      	write(fmt,'(a,i2.2,a)')  '(',nPar+2,'E20.12)'
      	read(5,*)
      	read(5,*)
            do i=1,nCls
            	locZ(i)=locZ(i)-nullZ
                  
                  maxp=0.
                  maxl=-huge(1.d0)
                  bestfit(i,1,1:nPar)=0.
                  sigma(i,1:nPar)=0.
                  
            	do j=1,nParCls(i)
                  	read(5,fmt)params(i,j,1:nPar+2)
                        !mean
                        bestfit(i,1,1:nPar)=bestfit(i,1,1:nPar)+params(i,j,1)*params(i,j,3:nPar+2)
                        !sigma
                       	sigma(i,1:nPar)=sigma(i,1:nPar)+params(i,j,1)*params(i,j,3:nPar+2)**2
                        !max-like
                        if(params(i,j,2)>maxl) then
                        	maxl=params(i,j,2)
                              bestfit(i,2,1:nPar)=params(i,j,3:nPar+2)
				endif
                        !map
                        if(params(i,j,1)>maxp) then
                        	maxp=params(i,j,1)
                              bestfit(i,3,1:nPar)=params(i,j,3:nPar+2)
				endif
                        !convert the probabilities into cumulative probabilities
                       	if(j>1) params(i,j,1)=params(i,j,1)+params(i,j-1,1)
			enddo
                  !normalize the mean & sigma
                  bestfit(i,1,1:nPar)=bestfit(i,1,1:nPar)/params(i,nParCls(i),1)
                  sigma(i,1:nPar)=sigma(i,1:nPar)/params(i,nParCls(i),1)
                  sigma(i,1:nPar)=sqrt(sigma(i,1:nPar)-bestfit(i,1,1:nPar)**2)
                  
                  !nomralize the cumulative probabilities
                  params(i,1:nParCls(i),1)=params(i,1:nParCls(i),1)/params(i,nParCls(i),1)
                  
                  !read the blank lines
                  read(5,*)
                  read(5,*)
		enddo
            close(5)
            
            maxp=1.d99
            !now write the cluster parameters in the descending order of evidence
           	open(5,file=TRIM(baseroot)//'clusters.cat',form='formatted',status='replace')
            write(fmt,'(a,i4,a)')'(i4,',(nPar+1)*2,'E14.6)'
            do i=1,nCls
            	!find the cluster with max evidence
                  indx=maxloc(locZ(1:nCls),MASK=locZ<maxp)
                  maxp=locZ(indx(1))
                  do j=1,nPar
                  	temp((j-1)*2+1)=bestfit(indx(1),1,j)
                        temp(j*2)=sigma(indx(1),j)
			enddo
                  write(5,fmt)i,temp,locZ(indx(1)),rdummy(indx(1))
		enddo
            close(5)
            
	end subroutine readPostFiles

!-----------------------------------------------------------------------

	subroutine calcStats(nPar,nCls,locZ,nParCls,params,bestfit,Ztol,n,infile,outfile,fitsFile)
      	!input variables
            integer npar !total number of parameters
            integer nCls !total no. of modes
            double precision locZ(:) !difference between local logZ & log-nullZ
            integer nParCls(:) !total no. of posterior samples per mode
            double precision params(:,:,:) !cumulative prob, log-like & posterior parameters
            double precision bestfit(:,:,:) !mean, max-like & map params of each cluster
            double precision Ztol !what difference in Ztol is to be tolerated
      	integer n !total number of monte carlo realizations
            character*100 infile !file containing true mass map
            character*100 outfile !output file
            character*100 fitsFile !output kappa map file
            
            !work variable
            integer i,j,k,m,flag
            double precision kappar(nxpix,nypix) !true kappa map
            double precision kappam(nxpix,nypix) !inference kappa map
            double precision kappabf(nxpix,nypix) !best kappa map
            double precision meankappa(nxpix,nypix) !mean kappa map
            double precision urv,bestStat
            double precision meanr,meanm,meanc,sdr,sdm,sdc,mean_err,sd_err
            
            !read the true catalogue first
            !call readInCat(256,256,infile,kappar)
            
            !initialize the random no. generator
            call initRandomNS(1)
            
            NFWStyle(2:nCls)=NFWStyle(1)
            
            mean_err=0.
            sd_err=0.
            bestStat=1.d99
            meankappa=0.
            
           	open(unit=23,file='./data/GL/sim/starck/kappa.dat',form='formatted',status='unknown')
            do i=256,1,-1
                  read(23,*)kappar(1:256,i)
		enddo
            close(23)
            
            open(unit=23,file=outfile,form='formatted',status='unknown')
            do i=1,n+3
			!zero variables
			kappas=0.
			gamma1s=0.
			gamma2s=0.
                  
            	do j=1,nCls
                  	if(locZ(j)<Ztol) cycle
                        
                        if(i>3) then
                        	!pick a random sample from the posterior
            			urv=ranmarNS(0)
                        	do k=1,nParCls(j)
                        		if(params(j,k,1)>=urv) then
                                    	GeoPars(j,1)=params(j,k,3)
	  						GeoPars(j,2)=params(j,k,4)
	  
       						!Sort out elliptical geometry:
        						if(GeoModel==2) then
	    							GeoPars(j,3)=params(j,k,5)
	    							GeoPars(j,4)=params(j,k,6)
          							call EllGeometry(j)
                                          	m=2
							else
                                    		m=0
	  						endif
	    
	  						if (MassModel==1) then
	    							if (NFWstyle(j)==1) then
	      							MassPars(j,1)=params(j,k,5+m)
	      							MassPars(j,2)=params(j,k,6+m)
                                                	z=params(j,k,7+m)
	    							elseif (NFWstyle(j)==2) then
	      							MassPars(j,1)=params(j,k,5+m)
		      						MassPars(j,2)=params(j,k,6+m)
      	                                          z=params(j,k,7+m)
	    							elseif (NFWstyle(j)==3) then
	      							MassPars(j,1)=params(j,k,5+m)
                        	                        z=params(j,k,6+m)
	    							elseif (NFWstyle(j)==4) then
	      							MassPars(j,1)=params(j,k,5+m)
	      							MassPars(j,2)=params(j,k,6+m)
                                                	z=params(j,k,7+m)
		    						endif  
		  					elseif (MassModel==2) then 
	    							MassPars(j,1)=params(j,k,5+m)
                  	                        z=params(j,k,6+m)
	  						endif
                                    	exit
						endif
					enddo
				else
                        	GeoPars(j,1)=bestfit(j,i,1)
	  				GeoPars(j,2)=bestfit(j,i,2)
	  
       				!Sort out elliptical geometry:
        				if(GeoModel==2) then
	    					GeoPars(j,3)=bestfit(j,i,3)
	    					GeoPars(j,4)=bestfit(j,i,4)
          					call EllGeometry(j)
						m=2
					else
						m=0
	 				endif
	   
 					if (MassModel==1) then
	    					if (NFWstyle(j)==1) then
	      					MassPars(j,1)=bestfit(j,i,3+m)
	      					MassPars(j,2)=bestfit(j,i,4+m)
							z=bestfit(j,i,7+m)
	  					elseif (NFWstyle(j)==2) then
	      					MassPars(j,1)=bestfit(j,i,3+m)
		      				MassPars(j,2)=bestfit(j,i,4+m)
							z=bestfit(j,i,5+m)
						elseif (NFWstyle(j)==3) then
	      					MassPars(j,1)=bestfit(j,i,5+m)
							z=bestfit(j,i,6+m)
	    					elseif (NFWstyle(j)==4) then
	      					MassPars(j,1)=bestfit(j,i,3+m)
	      					MassPars(j,2)=bestfit(j,i,4+m)
							z=bestfit(j,i,5+m)
						endif  
					elseif (MassModel==2) then 
	    					MassPars(j,1)=bestfit(j,i,3+m)
						z=bestfit(j,i,4+m)
					endif
				endif
                  	
				D=calcAngDiamDis(0.d0,z)
                        if(z_s(1)==0.) then
                        	zs=z_ss(1,1)
				else
                        	zs=z_s(1)
				endif
                       	zsmin=zs
                      	zsmax=zs
				Vary_zs=0
                        
                        !Now increase the shear and convergence at each galaxy position:
                        call LensFields(j,flag,.false.)
        			if(flag==1) stop 'Flag = 1 in LensFields'
			enddo
                  
                  !error statistic calculation
                  meanr=0.
                  meanm=0.
                  meanc=0.
                  sdr=0.
                  sdm=0.
                  sdc=0.
                  do j=1,256
      			do k=1,256
                        	!real kappa mean & variance
                        	meanr=meanr+kappar(j,k)
                              sdr=sdr+kappar(j,k)**2.
                              !model kappa mean & variance
                        	meanm=meanm+kappas(j,k)
                              sdm=sdm+kappas(j,k)**2.
                              !kappar error mean & variance
                  		meanc=meanc+kappar(j,k)-kappas(j,k)
                  		sdc=sdc+(kappar(j,k)-kappas(j,k))**2.
				enddo
			enddo
      		meanr=meanr/(256*256)
      		sdr=sdr/(256*256)-meanr**2.
      		meanm=meanm/(256*256)
      		sdm=sdm/(256*256)-meanm**2.
      		meanc=meanc/(256*256)
      		sdc=sdc/(256*256)-meanc**2.
                  
                  if(i<=3) then
                  	if(i==1) then
                        	write(*,*)"error statistic for parameter means",sqrt(sdc)/sqrt(sdr)
				elseif(i==2) then
                        	write(*,*)"error statistic for max likelihood parameter",sqrt(sdc)/sqrt(sdr)
				elseif(i==3) then
                        	write(*,*)"error statistic for map parameter",sqrt(sdc)/sqrt(sdr)
				endif
                  else
                  	meankappa=meankappa+kappas
                  	mean_err=mean_err+sqrt(sdc)/sqrt(sdr)
                        sd_err=sd_err+(sqrt(sdc)/sqrt(sdr))**2.
      			write(23,*)sqrt(sdc)/sqrt(sdr)
			endif
                  
                  !find the parameters for which the error statistic is the best
                  if(sqrt(sdc)/sqrt(sdr)<bestStat) then
                  	bestStat=sqrt(sdc)/sqrt(sdr)
                        kappabf=kappas
			endif
		enddo
            
            mean_err=mean_err/n
            sd_err=sqrt(sd_err/n-(mean_err)**2.)
            meankappa=meankappa/n
            
            !find the no. of cluster included
            j=0
            do i=1,nCls
            	if(locZ(i)>=Ztol) j=j+1
		enddo
            write(*,*)"posterior error statis:",mean_err,"+/-",sd_err
            write(*,*)j,"modes were included in the analysis"
            close(23)
            
            write(*,*)"best error statistic",bestStat
            
            
            !Write out in fits format
		if(nxpix==nypix .and. pix_sizex==pix_sizey) then
	      	write(*,*) 'Writing true map to file: ',fitsFile
			call writefits(meankappa,nxpix,nypix,pix_sizex,fitsFile)
      		write(*,*)
		else
                  write(*,*)"can not generate the fits file for true kappa since &
				either the no. of pixels or the pixel sizes in x & y are not the same"
		endif
            call killRandomNS()
                        
	end subroutine calcStats

!-----------------------------------------------------------------------

	!reads in the true shear catalogue
	subroutine readInCat(nx,ny,filename,kappar)
      
            !input/out variables
      	integer nx,ny !no. of x & y pixels
            character*100 filename !file containting the catalogue
            double precision kappar(nx,ny)
            
            integer i,j
      
      	open(unit=23,file=filename,form='formatted',status='unknown')
            do
            	read(23,*)i,j,d1
			if(j>ny) cycle
                  if(i>nx) exit
			kappar(i,j)=d1
		enddo
            close(23)
            
	end subroutine readInCat

!-----------------------------------------------------------------------

	subroutine makeCatalogue
      
      	integer i,j,k
      	character*100 kappafile,gammafile,gfile
            double precision kappa(300,300),gamma1(300,300),gamma2(300,300)
            double precision g1(300,300),g2(300,300),gsq,err,zs,d1
            
           	kappafile='./data/GL/sim/starck/300/kappa.tab'
           	gammafile='./data/GL/sim/starck/300/g.tab'
           	gfile='./data/GL/sim/starck/300/g300.cat'
            
            zs=1.
            err=0.256
            
            open(unit=23,file=gammafile,form='formatted',status='unknown')
            do
            	read(23,*)i,j,k,d1
			if(j>300) cycle
                  if(i>300) exit
			if(k==1) then
                  	gamma1(i,j)=d1
			else
            		gamma2(i,j)=d1
			endif
		enddo
            close(23)
            
            open(unit=23,file=kappafile,form='formatted',status='unknown')
            do
            	read(23,*)i,j,d1
			if(j>300) cycle
                  if(i>300) exit
                  kappa(i,j)=d1
		enddo
            close(23)
            
            do i=1,300
            	do j=1,300
                  	!Now calculate reduced shear fields
				g1(i,j)=gamma1(i,j)/(1.-kappa(i,j))
				g2(i,j)=gamma2(i,j)/(1.-kappa(i,j))

				!Correct for strong lensing region
				gsq=g1(i,j)**2+g2(i,j)**2
				if (gsq>1.) then
					g1(i,j)=g1(i,j)/gsq
					g2(i,j)=g2(i,j)/gsq
					gsq=1./gsq
				endif
			enddo
		enddo
            call writecat_survey(300,300,g1,g2,err,zs,gfile)
	end subroutine makeCatalogue
                  
                  
end module mass_stats
