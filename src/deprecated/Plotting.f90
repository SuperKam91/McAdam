module plotting1
	use params

contains

!=======================================================================

	subroutine PlotSetup

	implicit none
	
	integer pgopen
	integer i,j,ix,iy,k,Nvistotal
	double precision xmin,xmax,ymin,ymax,aspect
	character*7 xchars,ychars
	
      Nvistotal=0
      do i=1,Nvisfiles
        Nvistotal=Nvistotal+Nvis(i)
      enddo
      
! If we are fitting atoms,

      if (Atoms.eq.1) then

! First plot atom parameters:

        do i=1,NAtoms

          dev(i)=pgopen('/xs')
          call pgslct(dev(i))

          ix=max(NGeoPlots,NMassPlots*Mass)
          ix=max(ix,(NGasPlots+NTplots)*Gas)
          ix=max(ix,NSrcPlots*SourceSubtract)
          iy=1+Mass+Gas
          aspect=float(iy)/float(ix)
          call pgpap(0.0,aspect)
          call pgsch(1.5)
          call pgsubp(ix,iy)

          k=1
          do j=1,NGeoPlots
            call pgpanl(j,k)
            call pgvport(0.15,0.85,0.15,0.85)
            xmin=Geo_Ranges(i,Geo_xplot(j),1)
            xmax=Geo_Ranges(i,Geo_xplot(j),2)
            ymin=Geo_Ranges(i,Geo_yplot(j),1)
            ymax=Geo_Ranges(i,Geo_yplot(j),2)
            call pgwindow(xmin,xmax,ymin,ymax)
            call pgbox('BCTSNI',0.0,0,'BCTSNI',0.0,0)
            call pglab(Geo_Labels(Geo_xplot(j)),Geo_Labels(Geo_yplot(j)), &
            'Geometry Parameters')
          enddo

          if (Mass.eq.1) then
          k=k+1
          do j=1,NMassPlots
            call pgpanl(j,k)
            call pgvport(0.15,0.85,0.15,0.85)
            xmin=Mass_Ranges(i,Mass_xplot(j),1)
            xmax=Mass_Ranges(i,Mass_xplot(j),2)
            if (Mass_xplot(j).ne.Mass_yplot(j)) then
              ymin=Mass_Ranges(i,Mass_yplot(j),1)
              ymax=Mass_Ranges(i,Mass_yplot(j),2)
            else
              ymin=GL*(GLLhood0-2*0.5*2*Ngals)+SZ*(SZLhood0-2*0.5*2*Nvistotal)
              ymax=GL*GLLhood0+SZ*SZLhood0
            endif
            if (Mass_PriorType(i,Mass_xplot(j)).eq.2) then
              xmin=log10(xmin)
              xmax=log10(xmax)
              xchars='BCTSNIL'
            else
              xchars='BCTSNI'
            endif
            call pgwindow(xmin,xmax,ymin,ymax)
            call pgbox(xchars,0.0,0,'BCTSNI',0.0,0)
            if (Mass_xplot(j).ne.Mass_yplot(j)) then
              call pglab(Mass_Labels(Mass_xplot(j)),Mass_Labels(Mass_yplot(j)), &
              'Mass Parameters')
            else
              call pglab(Mass_Labels(Mass_xplot(j)),'Likelihood','Mass Parameter')
            endif
          enddo
          endif

          if (Gas.eq.1) then

          k=k+1
          do j=1,NGasPlots
            call pgpanl(j,k)
            call pgvport(0.15,0.85,0.15,0.85)
            xmin=Gas_Ranges(i,Gas_xplot(j),1)
            xmax=Gas_Ranges(i,Gas_xplot(j),2)
            if (Gas_xplot(j).ne.Gas_yplot(j)) then
              ymin=Gas_Ranges(i,Gas_yplot(j),1)
              ymax=Gas_Ranges(i,Gas_yplot(j),2)
            else
              ymin=GL*(GLLhood0-2*0.5*2*Ngals)+SZ*(SZLhood0-2*0.5*2*Nvistotal)
              ymax=GL*GLLhood0+SZ*SZLhood0
            endif
            call pgwindow(xmin,xmax,ymin,ymax)
            call pgbox('BCTSNI',0.0,0,'BCTSNI',0.0,0)
            if (Gas_xplot(j).ne.Gas_yplot(j)) then
              call pglab(Gas_Labels(Gas_xplot(j)),Gas_Labels(Gas_yplot(j)), &
              'Gas Parameters')
            else
              call pglab(Gas_Labels(Gas_xplot(j)),'Likelihood','Gas Parameter')
            endif
          enddo

          do j=NGasPlots+1,NGasPlots+NTPlots
            call pgpanl(j,k)
            call pgvport(0.15,0.85,0.15,0.85)
            xmin=T_Ranges(i,T_xplot(j-NGasPlots),1)
            xmax=T_Ranges(i,T_xplot(j-NGasPlots),2)
            if (T_xplot(j-NGasPlots).ne.T_yplot(j-NGasPlots)) then
              ymin=T_ranges(i,T_yplot(j-NGasPlots),1)
              ymax=T_ranges(i,T_yplot(j-NGasPlots),2)
            else
              ymin=NullEv-GL*2*0.5*2*Ngals-SZ*2*0.5*2*Nvistotal
              ymax=NullEv+GL*2*0.5*2*Ngals+SZ*2*0.5*2*Nvistotal
            endif
            call pgwindow(xmin,xmax,ymin,ymax)
            call pgbox('BCTSNI',0.0,0,'BCTSNI',0.0,0)
            if (T_xplot(j-NGasPlots).ne.T_yplot(j-NGasPlots)) then
              call pglab(T_labels(T_xplot(j-NGasPlots)),T_labels(T_yplot(j-NGasPlots)), &
              'Gas Parameters')
            else
              call pglab(T_labels(T_xplot(j-NGasPlots)),'Likelihood','Gas Parameter')
            endif
          enddo

          endif

        enddo

	endif
      
! Now do nuisance parameters:
	  	
	if (Nuisance.eq.1) then
	
	  i=Atoms*NAtoms+1
	  dev(i)=pgopen('/xs')
	  call pgslct(dev(i))
		  	  
	  ix=NSrcPlots+NzsPlots
! 	  ix=1
	  iy=1
	  aspect=float(iy)/float(ix)	
	  call pgpap(0.0,aspect)
! 	  call pgsch(1.5)
	  call pgsubp(ix,iy)
	  
	  j=0
	  k=1
	
	  if (SourceSubtract.eq.1) then
	  
	    call pgpanl(1,k)
	    call pgvport(0.15,0.85,0.15,0.85)
          xmin=map_x
          xmax=map_x
          ymin=map_y
          ymax=map_y
	    do i=1,NSrc
            xmin=min(Src_Ranges(i,1,1),xmin)
	      xmax=max(Src_Ranges(i,1,2),xmax)
	      ymin=min(Src_Ranges(i,2,1),ymin)
	      ymax=max(Src_Ranges(i,2,2),ymax)
	    enddo
          xmin=map_x-1.1*max(abs(xmax-map_x),abs(map_x-xmin))
          xmax=map_x+1.1*max(abs(xmax-map_x),abs(map_x-xmin))
          ymin=map_y-1.1*max(abs(ymax-map_y),abs(map_y-ymin))
          ymax=map_y+1.1*max(abs(ymax-map_y),abs(map_y-ymin))
          call pgwindow(xmax,xmin,ymin,ymax)
          call pgbox('BCTSNI',0.0,0,'BCTSNI',0.0,0)
          call pglab(Src_Labels(1),Src_Labels(2),'Point Source Positions') 

          do j=2,NSrcPlots
            call pgpanl(j,k)
            call pgvport(0.15,0.85,0.15,0.85)
            i=j-1
            xmin=Src_Ranges(Src_xplot(i),3,1)
            xmax=Src_Ranges(Src_xplot(i),3,2)
            if (Src_xplot(i).ne.Src_yplot(i)) then
              ymin=Src_Ranges(Src_yplot(i),3,1)
              ymax=Src_Ranges(Src_yplot(i),3,2)
            else
              ymin=GL*(GLLhood0-2*0.5*2*Ngals)+SZ*(SZLhood0-2*0.5*2*Nvistotal)
              ymax=GL*GLLhood0+SZ*SZLhood0
            endif
            if (Src_PriorType(Src_xplot(i),3).eq.2) then
              xmin=log10(xmin)
              xmax=log10(xmax)
              xchars='BCTSNIL'
            else
              xchars='BCTSNI'
            endif
            if (Src_xplot(i).ne.Src_yplot(i).and.Src_PriorType(Src_yplot(i),3).eq.2) then
              ymin=log10(ymin)
              ymax=log10(ymax)
              ychars='BCTSNIL'
            else
              ychars='BCTSNI'
            endif
            call pgwindow(xmin,xmax,ymin,ymax)
            call pgbox(xchars,0.0,0,ychars,0.0,0)
            if (Src_xplot(i).ne.Src_yplot(i)) then
              call pglab(Src_Labels(3),Src_Labels(3),'Source Fluxes')
            else
              call pglab(Src_Labels(3),'Likelihood','Source Flux')
            endif
          enddo

	  endif
	
	  if (Mass.eq.1.and.Varyzs.eq.1) then
	  
	    j=j+1
	  
	    call pgpanl(j,k)
	    call pgvport(0.15,0.85,0.15,0.85)
!  	  write(*,*) 'zs_Ranges: ',zs_Ranges
	    xmin=zs_Ranges(1)
	    xmax=zs_Ranges(2)
	    ymin=GL*(GLLhood0-2*0.5*2*Ngals)+SZ*(SZLhood0-2*0.5*2*Nvistotal)
	    ymax=GL*GLLhood0+SZ*SZLhood0
	    if (zs_PriorType.eq.2) then
	      xmin=log10(xmin)
	      xmax=log10(xmax)
		xchars='BCTSNIL'
	    else	
		xchars='BCTSNI'
	    endif	
	    ychars='BCTSNI'
	    call pgwindow(xmin,xmax,ymin,ymax)
          call pgbox(xchars,0.0,0,ychars,0.0,0)
	    call pglab(zs_Label,'Likelihood','Source plane redshift')
	
	  endif
	
	
	endif
	  
	return
	end subroutine PlotSetup

!=======================================================================

	subroutine PlotPoints(Uflag,Lhood)

	implicit none
	
	double precision Lhood
	integer i,j,iy,Uflag,k,Nvistotal
	double precision xmin,xmax,ymin,ymax,xpt,ypt
	
	
      Nvistotal=0
      do i=1,Nvisfiles
        Nvistotal=Nvistotal+Nvis(i)
      enddo
      
! If we are fitting atoms,

      if (Atoms.eq.1) then

! First plot atom parameters:

        iy=2

        do i=1,NAtoms

          call pgslct(dev(i))

          if (Burnin.ge.1) then
            if (Uflag.eq.1) then
              call pgsci(7)
            elseif (Uflag.eq.2) then
              call pgsci(2)
            endif
          else
            if (Uflag.eq.1) then
              call pgsci(7)
            elseif (Uflag.eq.2) then
              call pgsci(5)
            endif
          endif

          k=1
          do j=1,NGeoPlots
            call pgpanl(j,k)
            call pgvport(0.15,0.85,0.15,0.85)
            xmin=Geo_Ranges(i,Geo_xplot(j),1)
            xmax=Geo_Ranges(i,Geo_xplot(j),2)
            ymin=Geo_Ranges(i,Geo_yplot(j),1)
            ymax=Geo_Ranges(i,Geo_yplot(j),2)
            call pgwindow(xmin,xmax,ymin,ymax)
            xpt=GeoPars(i,Geo_xplot(j))
            ypt=GeoPars(i,Geo_yplot(j))
            call pgpt1(xpt,ypt,17)
          enddo

          if (Mass.eq.1) then
          k=k+1
          do j=1,NMassPlots
            call pgpanl(j,k)
            call pgvport(0.15,0.85,0.15,0.85)
            xmin=Mass_Ranges(i,Mass_xplot(j),1)
            xmax=Mass_Ranges(i,Mass_xplot(j),2)
            if (Mass_xplot(j).ne.Mass_yplot(j)) then
              ymin=Mass_Ranges(i,Mass_yplot(j),1)
              ymax=Mass_Ranges(i,Mass_yplot(j),2)
              ypt=MassPars(i,Mass_yplot(j))
            else
              ymin=GL*(GLLhood0-2*0.5*2*Ngals)+SZ*(SZLhood0-2*0.5*2*Nvistotal)
              ymax=GL*GLLhood0+SZ*SZLhood0
              ypt=1.0*Lhood
            endif
            if (Mass_PriorType(i,Mass_xplot(j)).eq.2) then
              xmin=log10(xmin)
              xmax=log10(xmax)
              xpt=dlog10(MassPars(i,Mass_xplot(j)))
            else
              xpt=MassPars(i,Mass_xplot(j))
            endif
            call pgwindow(xmin,xmax,ymin,ymax)
            call pgpt1(xpt,ypt,17)
          enddo
          endif

          if (Gas.eq.1) then
          k=k+1
          do j=1,NGasPlots
            call pgpanl(j,k)
            call pgvport(0.15,0.85,0.15,0.85)
            xmin=Gas_Ranges(i,Gas_xplot(j),1)
            xmax=Gas_Ranges(i,Gas_xplot(j),2)
            if (Gas_xplot(j).ne.Gas_yplot(j)) then
              ymin=Gas_Ranges(i,Gas_yplot(j),1)
              ymax=Gas_Ranges(i,Gas_yplot(j),2)
            else
              ymin=GL*(GLLhood0-2*0.5*2*Ngals)+SZ*(SZLhood0-2*0.5*2*Nvistotal)
              ymax=GL*GLLhood0+SZ*SZLhood0
            endif
            call pgwindow(xmin,xmax,ymin,ymax)
            xpt=GasPars(i,Gas_xplot(j))
            if (Gas_xplot(j).ne.Gas_yplot(j)) then
              ypt=GasPars(i,Gas_yplot(j))
            else
              ypt=1.0*Lhood
            endif
            call pgpt1(xpt,ypt,17)
          enddo
          do j=NGasPlots+1,NGasPlots+NTPlots
            call pgpanl(j,k)
            call pgvport(0.15,0.85,0.15,0.85)
            xmin=T_Ranges(i,T_xplot(j-NGasPlots),1)
            xmax=T_Ranges(i,T_xplot(j-NGasPlots),2)
            if (T_xplot(j-NGasPlots).ne.T_yplot(j-NGasPlots)) then
              ymin=T_ranges(i,T_yplot(j-NGasPlots),1)
              ymax=T_ranges(i,T_yplot(j-NGasPlots),2)
            else
              ymin=NullEv-GL*2*0.5*2*Ngals-SZ*2*0.5*2*Nvistotal
              ymax=NullEv+GL*2*0.5*2*Ngals+SZ*2*0.5*2*Nvistotal
            endif
            call pgwindow(xmin,xmax,ymin,ymax)
            xpt=TPars(i,T_xplot(j-NGasPlots))
            if (T_xplot(j-NGasPlots).ne.T_yplot(j-NGasPlots)) then
              ypt=TPars(i,T_yplot(j-NGasPlots))
            else
              ypt=1.0*Lhood
            endif
            call pgpt1(xpt,ypt,17)
          enddo
          endif

        enddo
        
      endif  
	
! Now plot nuisance parameters:	  
	
	if (Nuisance.eq.1) then
	
	  i=Atoms*NAtoms+1
	  call pgslct(dev(i))
	  
	  if (Burnin.eq.1) then
	    if (Uflag.eq.1) then
	      call pgsci(7)
	    elseif (Uflag.eq.2) then
	      call pgsci(2)
	    endif	
	  else
	    if (Uflag.eq.1) then
	      call pgsci(7)
	    elseif (Uflag.eq.2) then
	      call pgsci(5)
	    endif	
	  endif  
	
	  k=1
	  j=0

	  if (SourceSubtract.eq.1) then
	  
	    call pgpanl(1,k)
	    call pgvport(0.15,0.85,0.15,0.85)
          xmin=map_x
          xmax=map_x
          ymin=map_y
          ymax=map_y
	    do i=1,NSrc
            xmin=min(Src_Ranges(i,1,1),xmin)
	      xmax=max(Src_Ranges(i,1,2),xmax)
	      ymin=min(Src_Ranges(i,2,1),ymin)
	      ymax=max(Src_Ranges(i,2,2),ymax)
	    enddo
          xmin=map_x-1.1*max(abs(xmax-map_x),abs(map_x-xmin))
          xmax=map_x+1.1*max(abs(xmax-map_x),abs(map_x-xmin))
          ymin=map_y-1.1*max(abs(ymax-map_y),abs(map_y-ymin))
          ymax=map_y+1.1*max(abs(ymax-map_y),abs(map_y-ymin))
          call pgwindow(xmax,xmin,ymin,ymax)
          do i=1,NSrc
            xpt=SrcPars(i,1)
            ypt=SrcPars(i,2)
            call pgpt1(xpt,ypt,1)
          enddo  

          
          do j=2,NSrcPlots
            call pgpanl(j,k)
            call pgvport(0.15,0.85,0.15,0.85)
            i=j-1
            xmin=Src_Ranges(Src_xplot(i),3,1)
            xmax=Src_Ranges(Src_xplot(i),3,2)
            if (Src_xplot(i).ne.Src_yplot(i)) then
              ymin=Src_Ranges(Src_yplot(i),3,1)
              ymax=Src_Ranges(Src_yplot(i),3,2)
            else
              ymin=GL*(GLLhood0-2*0.5*2*Ngals)+SZ*(SZLhood0-2*0.5*2*Nvistotal)
              ymax=GL*GLLhood0+SZ*SZLhood0
            endif
            if (Src_PriorType(Src_xplot(i),3).eq.2) then
              xmin=log10(xmin)
              xmax=log10(xmax)
              xpt=dlog10(SrcPars(Src_xplot(i),3))
            else
              xpt=SrcPars(Src_xplot(i),3)
            endif
            if (Src_xplot(i).ne.Src_yplot(i).and.Src_PriorType(Src_yplot(i),3).eq.2) then
              ymin=log10(ymin)
              ymax=log10(ymax)
              ypt=dlog10(SrcPars(Src_yplot(i),3))
            elseif (Src_xplot(i).ne.Src_yplot(i)) then
              ypt=SrcPars(Src_yplot(i),3)
            else
              ypt=1.0*Lhood
            endif
            call pgwindow(xmin,xmax,ymin,ymax)
            call pgpt1(xpt,ypt,17)
	    enddo

	  endif

	  if (Mass.eq.1.and.Varyzs.eq.1) then
	  
	    j=j+1
	  
	    call pgpanl(j,k)
	    call pgvport(0.15,0.85,0.15,0.85)
	    xmin=zs_Ranges(1)
	    xmax=zs_Ranges(2)
	    if (zs_PriorType.eq.2) then
	      xmin=log10(xmin)
	      xmax=log10(xmax)
	      xpt=dlog10(zs)
	    else	
	      xpt=1.0*zs
	    endif	
	    ymin=GL*(GLLhood0-2*0.5*2*Ngals)+SZ*(SZLhood0-2*0.5*2*Nvistotal)
	    ymax=GL*GLLhood0+SZ*SZLhood0
! 	    ymin=GL*(GLLhood0-5*0.5*2*Ngals)+SZ*(SZLhood0-5*0.5*2*Nvistotal)
! 	    ymax=GL*GLLhood0+SZ*SZLhood0
	    ypt=1.0*Lhood
	    call pgwindow(xmin,xmax,ymin,ymax)
! 	    write(*,*) 
! 	    write(*,*) 'Plotting zs=',zs,xmin,xmax,xpt
! 	    write(*,*) 
          call pgpt1(xpt,ypt,17)
	  endif
	
	endif
	
	return
	end subroutine PlotPoints

!=======================================================================

end module plotting1
