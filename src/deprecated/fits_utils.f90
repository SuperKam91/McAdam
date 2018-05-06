module fits_utils

contains
!=======================================================================
!
	subroutine readfits(array,nx,ny,cell,filename)
!
      implicit none
	
      integer  nx,ny
      double precision   array(nx,ny)
      double precision 	   array_sp(nx,ny)
      double precision     cell
      character filename*100
!
	integer unit,status,blocksize,readwrite,group
      integer maxdim,bitpix,naxis,naxes(2),pcount,gcount,decimals
	double precision nullval,skycell
	logical anyf
      logical simple,extend
	character*100 dummy
      integer i,j
	
!-----------------------------------------------------------------------

      status=0

      call ftgiou(unit,status)
	if (status.ne.0) goto 999
      
! Open the input file, with read-only access

      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
	if (status.ne.0) goto 999

! Initialize variables

      group = 1
      nullval = 0.0
	decimals = 10

! Read data from file:

      call ftghpr(unit,maxdim,simple,bitpix,naxis,naxes,pcount,gcount,extend,status)
	if (status.ne.0) goto 999

      if (naxis.ne.maxdim) then
!        write(*,*) 'WARNING: array dimension mismatch:'
!        write(*,*) 'maxdim is expected to be ',naxis,' but is ', &
!	   maxdim,'. Continuing with maxdim reset...'
        maxdim = naxis
	end if
	if (naxes(1).ne.nx.or.naxes(2).ne.ny) then
        write(*,*) 'ERROR: file sizes do not match in readfits:'
        write(*,*) 'data is expected to be ',nx,' by ',ny,' but is ', &
        naxes(1),' by ',naxes(2)
        write(*,*) 'File read in = ',filename
        write(*,*) 'Naxis,Naxes(1),Naxes(2) = ',Naxis,Naxes(1),Naxes(2)
        write(*,*) 'Bitpix = ',bitpix
        write(*,*) 'Status = ',status
        stop
      end if
		
      call ftgkyd(unit,'CDELT2',skycell,dummy,status)
	if (status.ne.0) goto 999

	cell = skycell*3600.0
	    
      call ftg2de(unit,group,nullval,nx,nx,ny,array_sp,anyf,status)
      
      do i=1,nx
      	do j=1,ny
      		array=dble(array_sp(i,j))
		end do
	end do
      
	if (status.ne.0) goto 999

! Close files and finish 
	
      call ftclos(unit,status)
	if (status.ne.0) goto 999
      call ftfiou(unit,status)
	if (status.ne.0) goto 999

! Check for any error, and if so print out error messages

 999  if (status .gt. 0) call printerror(status)
	
	return
	end subroutine readfits
!
!=======================================================================
!
	subroutine readfits_1(array,nx,ny,cell,filename)
!
      implicit none
	
      integer  nx,ny
      double precision     array(nx,ny)
      double precision     array_sp(nx,ny)
      double precision     cell
      character filename*100
!
	integer unit,status,blocksize,readwrite,group
      integer maxdim,bitpix,naxis,naxes(2),pcount,gcount,decimals
	double precision nullval
	logical anyf
      logical simple,extend
      integer i,j
!	
!-----------------------------------------------------------------------
!
      status=0
!
      call ftgiou(unit,status)
!
! Open the input file, with read-only access
!
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
!
!
! Initialize variables
!
      group = 1
      nullval = 0.0
	decimals = 10
!
! Read data from file:
!
!
      call ftghpr(unit,maxdim,simple,bitpix,naxis,naxes,pcount,gcount,extend,status)
      if (naxis.ne.maxdim) then
!        write(*,*) 'WARNING: array dimension mismatch:'
!        write(*,*) 'maxdim is expected to be ',naxis,' but is ', &
!	   maxdim,'. Continuing with maxdim reset...'
        maxdim = naxis
	end if
	if (naxes(1).ne.nx.or.naxes(2).ne.ny) then
        write(*,*) 'ERROR: file sizes do not match in readfits:'
        write(*,*) 'data is expected to be ',nx,' by ',ny,' but is ', &
        naxes(1),' by ',naxes(2)
        stop
      end if
!		
	cell = 1.0/float(nx)
!	    
      call ftg2de(unit,group,nullval,nx,nx,ny,array_sp,anyf,status)
      do i=1,nx
      	do j=1,ny
            	array(i,j)=dble(array_sp(i,j))
		end do
	end do
!
! Close files and finish 
!	
      call ftclos(unit,status)
      call ftfiou(unit,status)
!
! Check for any error, and if so print out error messages
!
      if (status .gt. 0) call printerror(status)
!	
	return
	end subroutine readfits_1
!
!=======================================================================
!
	subroutine readfitsn(nx,ny,filename)
!
      implicit none
	
      integer  nx,ny
      character filename*100
!
	integer unit,status,blocksize,readwrite
      integer maxdim,bitpix,naxis,naxes(2),pcount,gcount
      logical simple,extend
!	
!-----------------------------------------------------------------------
!
      status=0
!
      call ftgiou(unit,status)
!
! Open the input file, with read-only access
!
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
!
! Read data from file:
!
      call ftghpr(unit,maxdim,simple,bitpix,naxis,naxes,pcount,gcount,extend,status)
!
	nx = naxes(1)
	ny = naxes(2)		
!
! Close files and finish 
!	
      call ftclos(unit,status)
      call ftfiou(unit,status)
!
! Check for any error, and if so print out error messages
!
      if (status.gt.0) call printerror(status)
!	  	
	return
	end subroutine readfitsn
!
!=======================================================================
!
	subroutine readfitsint(array,nx,ny,cell,filename)
!
      implicit none
	
      integer  nx,ny
      integer  array(nx,ny)
      double precision  cell
      character filename*100
!
	integer unit,status,blocksize,readwrite,group,nullval
      integer maxdim,bitpix,naxis,naxes(2),pcount,gcount,decimals
	double precision skycell
	logical anyf
      logical simple,extend
	character*100 dummy
!	
!-----------------------------------------------------------------------
!
      status=0
!
      call ftgiou(unit,status)
!
! Open the input file, with read-only access
!
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
!
!
! Initialize variables
!
      group = 1
      nullval = 0
	decimals = 10
!
! Read data from file:
!
!
      call ftghpr(unit,maxdim,simple,bitpix,naxis,naxes,pcount,gcount,extend,status)
      if (naxis.ne.maxdim) then
        write(*,*) 'WARNING: array dimension mismatch:'
        write(*,*) 'maxdim is expected to be ',naxis,' but is ', &
        maxdim,'. Continuing with maxdim reset...'
        maxdim = naxis
	end if
	if (naxes(1).ne.nx.or.naxes(2).ne.ny) then
        write(*,*) 'ERROR: file sizes do not match in readfits:'
        write(*,*) 'data is expected to be ',nx,' by ',ny,' but is ', &
        naxes(1),' by ',naxes(2)
        stop
      end if
!		
      call ftgkyd(unit,'CDELT2',skycell,dummy,status)
!
	cell = skycell*3600.0
!	    
      call ftg2dj(unit,group,nullval,nx,nx,ny,array,anyf,status)
!
! Close files and finish 
!	
      call ftclos(unit,status)
      call ftfiou(unit,status)
!
! Check for any error, and if so print out error messages
!
      if (status .gt. 0) call printerror(status)
!	
	return
	end subroutine readfitsint
!
!=======================================================================
!
	subroutine writefits(array,nx,ny,cell,filename)
!
      implicit none
	
      integer nx,ny
      double precision array(nx,ny)
      double precision array_sp(nx,ny)
      double precision cell
      character filename*100
!
	integer unit,status,blocksize,group,fpixel,nelements
      integer bitpix,naxis,naxes(2),xrefpix,yrefpix,decimals
	double precision ra,dec,z,skycell
      logical simple,extend
      integer i,j
!	
!-----------------------------------------------------------------------
!
      status=0
      blocksize=1
!
! Create the input file:
!
      call deletefile(filename,status)
      call ftgiou(unit,status)
      call ftinit(unit,filename,blocksize,status)
!
! Initialize variables:
!
!     FITS image (nx x ny 4 byte floats)
!
      group = 1
      fpixel=1
      nelements=nx*ny
      simple=.true.
      bitpix=-32
      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      extend=.true.

!     Header (world coordinate system)

      ra = 0.0
      dec = 0.0
	skycell = dble(cell/3600.0)
      xrefpix = nx/2+1
      yrefpix = ny/2+1
      z = 0.0
      decimals = 10
!
!
! Write compulsory part of header:
!
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	
!     and then the optional keywords:

      call ftpkys(unit,'CTYPE1', 'RA---SIN','Axis 1', status)
      call ftpkys(unit,'CTYPE2', 'DEC--SIN','Axis 2', status)
      call ftpkyg(unit,'CRVAL1',ra,decimals,'Coordinate 1 at ref pixel',status)
      call ftpkyg(unit,'CRVAL2',dec,decimals,'Coordinate 2 at ref pixel',status)
      skycell = -1.0 * skycell
      call ftpkyg(unit,'CDELT1', skycell,decimals,'Coordinate 1/pixel', status)
      skycell = -1.0 * skycell
      call ftpkyg(unit,'CDELT2', skycell,decimals,'Coordinate 1/pixel', status)
      call ftpkyj(unit,'CRPIX1', xrefpix,'Reference pixel 1', status)
      call ftpkyj(unit,'CRPIX2', yrefpix,'Reference pixel 2', status)
      call ftpkyg(unit,'CROTA1',z,decimals,'Rotation 1', status)
      call ftpkyg(unit,'CROTA2',z,decimals,'Rotation 2', status)
      call ftpkyj(unit,'LTV1',xrefpix,'Physical to image',status)
      call ftpkyj(unit,'LTV2',yrefpix,'Physical to image',status)
      call ftpkyg(unit,'LTM1_1',1.0/cell,decimals,'Physical to image',status)
      call ftpkyg(unit,'LTM2_2',1.0/cell,decimals,'Physical to image',status)
!
! Write data to file:
	do i=1,nx
      	do j=1,ny
            	array_sp(i,j)=dble(array(i,j))
            end do
	end do
!
      call ftp2de(unit,group,nx,nx,ny,array_sp,status)
!
! Close files and finish 
!	
      call ftclos(unit,status)
      call ftfiou(unit,status)
!
! Check for any error, and if so print out error messages
!
      if (status.gt.0) call printerror(status)
!	
	return
	end subroutine writefits
!
!=======================================================================
!
	subroutine writefitsWCS(array,nx,ny,ra,dec,cell,R,filename)
!
      implicit none
	
      integer nx,ny
      double precision array(nx,ny)
      double precision array_sp(nx,ny)
      double precision cell,ra,dec,R(2,2)
      character filename*100
!
	integer unit,status,blocksize,group,fpixel,nelements
      integer bitpix,naxis,naxes(2),xrefpix,yrefpix,decimals
	double precision z,skycell
      logical simple,extend
      integer i,j
!	
!-----------------------------------------------------------------------
!
      status=0
      blocksize=1
!
! Create the input file:
!
      call deletefile(filename,status)
      call ftgiou(unit,status)
      call ftinit(unit,filename,blocksize,status)
!
! Initialize variables:
   
!     FITS image (nx x ny 4 byte floats) 

      group = 1
      fpixel=1
      nelements=nx*ny
      simple=.true.
      bitpix=-32
      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      extend=.true.

!     Header (world coordinate system)
      
      skycell = dble(cell/3600.0)
      xrefpix = nx/2+1
      yrefpix = ny/2+1
      z = 0.0
      decimals = 10
!
!
! Write compulsory part of header:
!
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	
!     and then the optional keywords:

      call ftpkyj(unit,'CRPIX1', xrefpix,'Reference pixel 1', status)
      call ftpkyg(unit,'CRVAL1',dble(ra),decimals,'Coordinate 1 at ref pixel',status)
      call ftpkys(unit,'CTYPE1', 'RA---TAN','Axis 1', status)
      call ftpkyg(unit,'CD1_1',dble(R(1,1)),decimals, &
      'Coordinate rotation matrix R(1,1) / degrees', status)
      call ftpkyg(unit,'CD2_1',dble(R(2,1)),decimals, &
      'Coordinate rotation matrix R(2,1) / degrees', status)
      
      
      call ftpkyj(unit,'CRPIX2', yrefpix,'Reference pixel 2', status)
      call ftpkyg(unit,'CRVAL2',dble(dec),decimals,'Coordinate 2 at ref pixel',status)
      call ftpkys(unit,'CTYPE2', 'DEC--TAN','Axis 2', status)
      call ftpkyg(unit,'CD1_2',dble(R(1,2)),decimals, &
      'Coordinate rotation matrix R(1,2) / degrees', status)
      call ftpkyg(unit,'CD2_2',dble(R(2,2)),decimals, &
      'Coordinate rotation matrix R(2,2) / degrees', status)

! Old fits wcs numbers:
!       call ftpkyg(unit,'CDELT1', skycell,decimals,'Coordinate 1/pixel', status)
!       skycell = skycell * -1.0
!       call ftpkyg(unit,'CDELT2', skycell,decimals,'Coordinate 1/pixel', status)
!       call ftpkyg(unit,'CROTA1',z,decimals,'Rotation 1', status)
!       call ftpkyg(unit,'CROTA2',z,decimals,'Rotation 2', status)

      call ftpkyj(unit,'LTV1',xrefpix,'Physical to image',status)
      call ftpkyj(unit,'LTV2',yrefpix,'Physical to image',status)
      call ftpkyg(unit,'LTM1_1',1.0/cell,decimals,'Physical to image',status)
      call ftpkyg(unit,'LTM2_2',1.0/cell,decimals,'Physical to image',status)
!
! Write data to file:
!
	do i=1,nx
      	do j=1,ny
      		array_sp=dble(array(i,j))
		end do
	end do
      call ftp2de(unit,group,nx,nx,ny,array_sp,status)
!
! Close files and finish 
!	
      call ftclos(unit,status)
      call ftfiou(unit,status)
!
! Check for any error, and if so print out error messages
!
      if (status.gt.0) call printerror(status)
!	
	return
	end subroutine writefitsWCS
!
!=======================================================================
!
	subroutine readfits3d(array,nx,ny,nz,cell,filename)
!
      implicit none
	
      integer  nx,ny,nz
      double precision     array(nx,ny,nz)
      double precision     array_sp(nx,ny,nz)
      double precision     cell
      character filename*100
!
	integer unit,status,blocksize,readwrite,group
      integer maxdim,bitpix,naxis,naxes(3),pcount,gcount,decimals
	double precision nullval,skycell
	logical anyf
      logical simple,extend
	character*100 dummy
      integer i,j,k
!	
!-----------------------------------------------------------------------
!
      status=0
!
      call ftgiou(unit,status)
!
! Open the input file, with read-only access
!
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
!
!
! Initialize variables
!
      group = 1
      nullval = 0.0
	decimals = 10
!
! Read data from file:
!
!
      call ftghpr(unit,maxdim,simple,bitpix,naxis,naxes,pcount,gcount,extend,status)
         
      if (naxis.ne.maxdim) then
!        write(*,*) 'WARNING: array dimension mismatch:'
!        write(*,*) 'maxdim is expected to be ',naxis,' but is ', &
!	   maxdim,'. Continuing with maxdim reset...'
        maxdim = naxis
	end if
	if (naxes(1).ne.nx.or.naxes(2).ne.ny.or.naxes(3).ne.nz) then
        write(*,*) 'ERROR: file sizes do not match in readfits:'
        write(*,*) 'data is expected to be ',nx,' by ',ny,' by ',nz, &
        ' but is ',naxes(1),' by ',naxes(2),' by ',naxes(3)
        stop
      end if
!		
      call ftgkyd(unit,'CDELT2',skycell,dummy,status)
!
	cell = skycell*3600.0
      
      call ftg3de(unit,group,nullval,nx,ny,nx,ny,nz,array_sp,anyf,status)
!	    
	do i=1,nx
      	do j=1,ny
            	do k=1,nz
                  	array(i,j,k)=dble(array_sp(i,j,k))
			end do
		end do
	end do
!
! Close files and finish 
!	
      call ftclos(unit,status)
      call ftfiou(unit,status)
!
! Check for any error, and if so print out error messages
!
      if (status .gt. 0) call printerror(status)
!	
	return
	end subroutine readfits3d
!
!=======================================================================
!
	subroutine readfits3dn(nx,ny,nz,filename)
!
      implicit none
	
      integer  nx,ny,nz
      character filename*100
!
	integer unit,status,blocksize,readwrite
      integer maxdim,bitpix,naxis,naxes(3),pcount,gcount
      logical simple,extend
!	
!-----------------------------------------------------------------------
!
      status=0
!
      call ftgiou(unit,status)
!
! Open the input file, with read-only access
!
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
!
! Read data from file:
!
      call ftghpr(unit,maxdim,simple,bitpix,naxis,naxes,pcount,gcount,extend,status)
!
	nx = naxes(1)
	ny = naxes(2)		
	nz = naxes(3)		
!
! Close files and finish 
!	
      call ftclos(unit,status)
      call ftfiou(unit,status)
!
! Check for any error, and if so print out error messages
!
      if (status.gt.0) call printerror(status)
!	  	
	return
	end subroutine readfits3dn
!
!=======================================================================
!
	subroutine writefits3d(array,nx,ny,nz,cell,filename)
!
      implicit none
	
      integer nx,ny,nz
      double precision array(nx,ny,nz)
      double precision array_sp(nx,ny,nz)
      double precision cell
      character filename*100
!
	integer unit,status,blocksize,group,fpixel,nelements
      integer bitpix,naxis,naxes(3),xrefpix,yrefpix,zrefpix,decimals
	double precision ra,dec,z,skycell
      logical simple,extend
      integer i,j,k
!	
!-----------------------------------------------------------------------
!
      status=0
      blocksize=1
!
! Create the input file:
!
      call deletefile(filename,status)
      call ftgiou(unit,status)
      call ftinit(unit,filename,blocksize,status)
!
! Initialize variables:
!
!     FITS image (nx x ny x nz 4 byte floats)
!
      group = 1
      fpixel = 1
      nelements = nx*ny*nz
      simple = .true.
      bitpix = -32
      naxis = 3
      naxes(1) = nx
      naxes(2) = ny
      naxes(3) = nz
      extend = .true.

!     Header (world coordinate system)

      ra = 0.0
      dec = 0.0
      z = 0.0
	skycell = dble(cell/3600.0)
      xrefpix = nx/2+1
      yrefpix = ny/2+1
      zrefpix = nz/2+1
      decimals = 10
!
!
! Write compulsory part of header:
!
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	
!     and then the optional keywords:

      call ftpkyg(unit,'CRVAL1',ra,decimals,'Coordinate 1 at ref pixel',status)
      call ftpkyg(unit,'CRVAL2',dec,decimals,'Coordinate 2 at ref pixel',status)
      call ftpkyg(unit,'CRVAL3',z,decimals,'Coordinate 3 at ref pixel',status)
      skycell = -1.0 * skycell
      call ftpkyg(unit,'CDELT1', skycell,decimals,'Coordinate 1/pixel', status)
      skycell = -1.0 * skycell
      call ftpkyg(unit,'CDELT2', skycell,decimals,'Coordinate 1/pixel', status)
      call ftpkyg(unit,'CDELT3', skycell,decimals,'Coordinate 1/pixel', status)
      call ftpkyj(unit,'CRPIX1', xrefpix,'Reference pixel 1', status)
      call ftpkyj(unit,'CRPIX2', yrefpix,'Reference pixel 2', status)
      call ftpkyj(unit,'CRPIX3', zrefpix,'Reference pixel 3', status)
      call ftpkyg(unit,'CROTA1',z,decimals,'Rotation 1', status)
      call ftpkyg(unit,'CROTA2',z,decimals,'Rotation 2', status)
      call ftpkyg(unit,'CROTA3',z,decimals,'Rotation 3', status)
      call ftpkyj(unit,'LTV1',xrefpix,'Physical to image',status)
      call ftpkyj(unit,'LTV2',yrefpix,'Physical to image',status)
      call ftpkyj(unit,'LTV3',zrefpix,'Physical to image',status)
      call ftpkyg(unit,'LTM1_1',1.0/cell,decimals,'Physical to image',status)
      call ftpkyg(unit,'LTM2_2',1.0/cell,decimals,'Physical to image',status)
      call ftpkyg(unit,'LTM3_3',1.0/cell,decimals,'Physical to image',status)
!
! Write data to file:
	do i=1,nx
      	do j=1,ny
            	do k=1,nz
                  	array_sp(i,j,k)=dble(array(i,j,k))
			end do
		end do
	end do
!
      call ftp3de(unit,group,nx,ny,nx,ny,nz,array_sp,status)
!
! Close files and finish 
!	
      call ftclos(unit,status)
      call ftfiou(unit,status)
!
! Check for any error, and if so print out error messages
!
      if (status.gt.0) call printerror(status)
!	
	return
	end subroutine writefits3d
!
!=======================================================================
!
      subroutine deletefile(filename,status)

!     A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

!     simply return if status is greater than zero
      if (status .gt. 0)return

!     Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

!     try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
!         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103)then
!         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
!	
      else
!         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

!     free the unit number for later reuse
      call ftfiou(unit, status)
      
	end subroutine deletefile
!	
!=======================================================================
!
      subroutine printerror(status)

!     Print out the FITSIO error messages to the user and stop

      integer status
      character errtext*30,errmessage*80
!
!-----------------------------------------------------------------------
!
!     check if status is OK (no error); if so, simply return
      if (status .le. 0)return

!     get the text string which describes the error
 1    call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

!     read and print out all the error messages on the FITSIO stack
 2    call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
      
	stop
	
	end subroutine printerror
!
! !=======================================================================
! !
! !  VSA Project - Reduce program, write visibilities in FITS format.
! !
! !  Writes selected samples from the reduced data arrays to a file using
! !  the FITS random groups format, suitable for import to AIPS, etc.
! !
! !  History:  12/7/99 - original version [DJT].
! !            1/11/00 - output of antenna tables added [DJT].
! !            13/6/01 - added output of report to log file [DJT].
! !
! ! Newer version using values stored in common block in the include file
! ! (pjm)
! 
!       subroutine write_phits(vis_buffer,uv_buffer,weight_buffer,baselines,outfile)
! 
!       implicit none
! 	
!       include 'telescope.inc'
! 
!       character*100 outfile
!       complex vis_buffer(maxvis)
!       real weight_buffer(maxvis),uv_buffer(2,maxvis)
!       integer baselines(maxvis)
! 
!       integer status
!       integer nbytes,i,iout
!       real*4 group(9)
! 	real freq,uscale,vscale
! 	
! !-----------------------------------------------------------------------
! 
!       if (status.ne.0) return
! 
!       call io_enqout(iout)
! 
! !  Open the output file
! 
!       call phits_open(outfile,'WRITE',status)
! 
! !  Write FITS header records
! 
!       call phits_wrhdr(status)
! 
! !  Write visibilities for each baseline as a single group of data
! 
! 	freq = CRVAL4
! 	uscale = PSCAL1
! 	vscale = PSCAL2
! 	
!       do i = 1,GCOUNT
!            
!         group(1) = uv_buffer(1,i)/(uscale*freq)
!         group(2) = uv_buffer(2,i)/(vscale*freq)
!         group(3) = 0.0
!         group(4) = 0.0
!         group(5) = 0.0
!         group(6) = float(baselines(i))
!         group(7) = dble(vis_buffer(i))/bscale
!         group(8) = aimag(vis_buffer(i))/bscale
!         group(9) = weight_buffer(i)/bscale
! 	  
!         nbytes = 9*4
!         call phits_write(group,nbytes,status)
! 		
!       enddo
! 
! !  End of data, fill current block with zeros
! 
!       group(1)=0
!       call phits_fill(group,status)
! 
! !  Close the file
! 
!       call phits_close(status)
! 
!       end
! 	
! !=======================================================================
! !
!       subroutine phits_open(file, access, status)
! 
! !  VSA Project - Reduce program, open file for UVFITS output
! !
! !  Open output file, set file and output buffer byte pointers.
! !
! !  History:  12/7/99 - original version [DJT].
! !            25/4/01 - give warning if cnnot open file [KG].
! !
! *-----------------------------------------------------------------------
! 
!        character file*(*), access*(*)
!        integer status
! 
!        include 'telescope.inc'
! 
!        integer bytep
! 
!        integer fseek,getfilep
!        external  fseek,getfilep
! 
!        if (status.ne.0) return
! 
! !  Open output file
! 
!        call io_nxtlun(ifile, status)
! 
!        if (access.eq.'WRITE') then
! 
!           open(ifile, file=file, status='UNKNOWN', iostat=status)
!           if (status.eq.0) then
!              bytep = 0
!              status = fseek(ifile,bytep,0)
!              filep = getfilep(ifile)
!           else 
!              write(*,*) 'Cannot open FITS file'
!           endif
! 
!        endif
! 
!        nblock = 0
!        ibuff = 0
! 
!        end
! !=======================================================================
! !
!       subroutine phits_wrhdr(status)
! !
! !  VSA Project - Reduce program, write FITS header.
! !
! !  Constructs entries in the FITS header and writes the complete header
! !  to disc.
! !
! !  History:  12/7/99 - original version, adapted from UV2FITS program [DJT].
! !            1/11/00 - IF coordinate suppressed [DJT].
! !
! *-
!       integer status
! 
!       include 'telescope.inc'
! 
!       character scratch*32
!       integer   i
!       integer   iscr
! 
!       integer   nline
!       parameter (nline=80)
!       character line*(nline)
!       real*4    rline(nline/4)
!       equivalence (line,rline)
! 
!       real*8    one, zero
!       data      one, zero / 1.d0, 0.d0 /
! 
!       character stokes(7)*4
!       data      stokes / 'BEAM', 'I', 'Q', 'U', 'V', 'I-Q', 'I+Q' /
! 
!       integer   chr_lenb
!       external  chr_lenb
! 
!       if (status.ne.0) return
! 
! !  Entry, open scratch file
! 
! !     scratch = 'scr.out'
! !     call io_opefil(iscr,scratch,'write',0,status)
!       call io_opescr(iscr,scratch,'write',0,status)
! 
!       naxes(1) = 0
!       naxes(2) = 3
!       naxes(3) = 1
!       naxes(4) = 1
!       naxes(5) = 1
!       naxes(6) = 1
!       naxes(7) = 1
!       write(iscr,10) bitpix,naxis,(naxes(i),i=1,naxis)
!  10   format(
!      :  'SIMPLE  = ',19x,'T',1x,'/ Standard FITS format'/
!      :  'BITPIX  = ',16x,i4,1x,'/ Bits per pixel'/
!      :  'NAXIS   = ',17x,i3,1x,'/ Number of axes'/
!      :  'NAXIS1  = ',16x,i4,1x,'/ No image: uv data'/
!      :  'NAXIS2  = ',16x,i4,1x,'/ Complex: cos, sin, weight'/
!      :  'NAXIS3  = ',16x,i4,1x,'/ Number of polarisations'/
!      :  'NAXIS4  = ',16x,i4,1x,'/ Number of frequencies'/
!      :  'NAXIS5  = ',16x,i4,1x,'/ RA'/
!      :  'NAXIS6  = ',16x,i4,1x,'/ Dec'/
!      :  'BLOCKED = ',19x,'T',1x,'/ Tape may be blocked')
! 
! !  Random groups format
! 
!       write(iscr,11) pcount,gcount
!  11   format(
!      :  'GROUPS  = ',19x,'T',1x,'/ Groups data structure'/
!      :  'PCOUNT  = ',16x,i4,1x,'/ Parameters per group'/
!      :  'GCOUNT  = ',14x,i6,1x,'/ Number of groups'/
!      :  'EXTEND  = ',19x,'T',1x,'/ Extension is antenna table')
! 
! !  Object, telescope, date of observation
! 
!       write(iscr,12) object(1:16),telescop(1:3),telescop(1:3)
!  12   format(
!      :  'OBJECT  = ',1h',a,1h',3x,'/ Source/field name'/
!      :  'TELESCOP= ',1h',a,5x,1h',11x,'/ Radio telescope'/
!      :  'INSTRUME= ',1h',a,5x,1h',11x,'/ Instrument')
! !     :  'DATE-OBS= ',1h',i4,2('-',i2.2),1h',9x,
! !     :  'DATE-OBS= ',1h',2(i2.2,'/'),i2.2,1h',8x,
! !     :                             '/ Start date of observation')
! 
! !  Bunit, bzero, bscale, blank, epoch
! 
!       write(iscr,13) bunit,bzero,bscale
!  13   format(
!      :  'BUNIT   = ',1h',a,1h',11x,'/ Units of data'/
!      :  'BZERO   = ',4x,1pe16.9,1x,'/ Data offset'/
!      :  'BSCALE  = ',4x,1pe16.9,1x,'/ Data = FITS*BSCALE + BZERO')
! 
! !      blank = zero
! !      if(blank.ne.0.d0) write(iscr,14) blank
! ! 14   format(
! !     :  'BLANK   = ',4x,1pe16.9,1x,'/ Undefined values on tape')
! 
! !      if (obs_epoch.eq.1) epoch=1950.0d0
! !      if (obs_epoch.eq.2) epoch=2000.0d0
!       write(iscr,15) epoch
!  15   format(
!      :  'EPOCH   = ',10x,f10.4,1x,'/ Equinox of RA, Dec')
! 
! !  Pointing centre
! 
!       write(iscr,16) obsra,obsdec
!  16   format(
!      :  'OBSRA   = ',4x,1pe16.9,1x,'/ Antenna pointing RA'/
!      :  'OBSDEC  = ',4x,1pe16.9,1x,'/ Antenna pointing DEC')
! 
! 
! !  Details of group data
! 
!       write(iscr,17) crval2,cdel2,crpix2,crota2
!  17   format(
!      :  'CTYPE2  = ',10h'COMPLEX ',11x,'/'/
!      :  'CRVAL2  = ',4x,1pe16.9,1x,'/'/
!      :  'CDELT2  = ',4x,1pe16.9,1x,'/'/
!      :  'CRPIX2  = ',4x,1pe16.9,1x,'/'/
!      :  'CROTA2  = ',4x,1pe16.9,1x,'/')
! 
!       write(iscr,18) stokes(1),crval3,cdel3,crpix3,crota3
!  18   format(
!      :  'CTYPE3  = ',10h'STOKES  ',11x,'/ Stokes parameter : ',A/
!      :  'CRVAL3  = ',4x,1pe16.9,1x,'/'/
!      :  'CDELT3  = ',4x,1pe16.9,1x,'/'/
!      :  'CRPIX3  = ',4x,1pe16.9,1x,'/'/
!      :  'CROTA3  = ',4x,1pe16.9,1x,'/')
! 
!       write(iscr,19) crval4,cdel4,crpix4,crota4
!  19   format(
!      :  'CTYPE4  = ',10h'FREQ    ',11x,'/ Frequency, Hertz'/
!      :  'CRVAL4  = ',4x,1pe16.9,1x,'/'/
!      :  'CDELT4  = ',4x,1pe16.9,1x,'/'/
!      :  'CRPIX4  = ',4x,1pe16.9,1x,'/'/
!      :  'CROTA4  = ',4x,1pe16.9,1x,'/')
! 
! !     write(iscr,195) one,one,one,zero
! c195  format(
! !    :  'CTYPE5  = ',10h'IF      ',11x,'/'/
! !    :  'CRVAL5  = ',4x,1pe16.9,1x,'/'/
! !    :  'CDELT5  = ',4x,1pe16.9,1x,'/'/
! !    :  'CRPIX5  = ',4x,1pe16.9,1x,'/'/
! !    :  'CROTA5  = ',4x,1pe16.9,1x,'/')
! 
!       write(iscr,20) crval5,cdel5,crpix5,crota5
!  20   format(
!      :  'CTYPE5  = ',10h'RA      ',11x,'/ Right Ascension, degrees'/
!      :  'CRVAL5  = ',4x,1pe16.9,1x,'/'/
!      :  'CDELT5  = ',4x,1pe16.9,1x,'/'/
!      :  'CRPIX5  = ',4x,1pe16.9,1x,'/'/
!      :  'CROTA5  = ',4x,1pe16.9,1x,'/')
! 
!       write(iscr,21) crval6,cdel6,crpix6,crota6
!  21   format(
!      :  'CTYPE6  = ',10h'DEC     ',11x,'/ Declination, degrees'/
!      :  'CRVAL6  = ',4x,1pe16.9,1x,'/'/
!      :  'CDELT6  = ',4x,1pe16.9,1x,'/'/
!      :  'CRPIX6  = ',4x,1pe16.9,1x,'/'/
!      :  'CROTA6  = ',4x,1pe16.9,1x,'/')
! 
! !  Details of group parameters
! 
!       write(iscr,22) pscal1,pzero1
!  22   format(
!      :  'PTYPE1  = ',10h'UU      ',11x,'/ U coordinate, seconds'/
!      :  'PSCAL1  = ',4x,1pe16.9,1x,'/ Real = FITS*PSCAL + PZERO'/
!      :  'PZERO1  = ',4x,1pe16.9,1x,'/')
! 
!       write(iscr,23) pscal2,pzero2
!  23   format(
!      :  'PTYPE2  = ',10h'VV      ',11x,'/ V coordinate, seconds'/
!      :  'PSCAL2  = ',4x,1pe16.9,1x,'/ Real = FITS*PSCAL + PZERO'/
!      :  'PZERO2  = ',4x,1pe16.9,1x,'/')
! 
!       write(iscr,24) pscal3,pzero3
!  24   format(
!      :  'PTYPE3  = ',10h'WW      ',11x,'/ W coordinate, seconds'/
!      :  'PSCAL3  = ',4x,1pe16.9,1x,'/ Real = FITS*PSCAL + PZERO'/
!      :  'PZERO3  = ',4x,1pe16.9,1x,'/')
! 
! !  Date parameter, as two entries:
! !   ptype1: pzero1 = intpt(jd at start-time of observations)
! !           pscal1 = 1/128 days
! !           value of type1 parameter = jd of sample to nearest
! !              integral value of 1/128 days
! !   ptype2: pzero2 = 0
! !           pscal2 = 1/6e6; precision is about 0.01 second
! 
! !      jdzero = int(obs_mjd_start+2400000.5d0)
! !      jdscal1 = 1.d0/128.d0
! !      jdscal2 = 1.d0/6.d6
! !      write(iscr,25) jdscal1,jdzero,jdscal2,zero
!       write(iscr,25) pscal4,pzero4,pscal5,pzero5
!  25   format(
!      :  'PTYPE4  = ',10h'DATE    ',11x,'/ Julian Date' /
!      :  'PSCAL4  = ',4x,1pe16.9,1x,'/'/
!      :  'PZERO4  = ',4x,1pe16.9,1x,'/'/
!      :  'PTYPE5  = ',10h'DATE    ',11x,'/ Julian Date'/
!      :  'PSCAL5  = ',4x,1pe16.9,1x,'/'/
!      :  'PZERO5  = ',4x,1pe16.9,1x,'/')
! 
!       write(iscr,26) pscal6,pzero6
!  26   format(
!      :  'PTYPE6  = ',10h'BASELINE',11x,'/ Ant1*256 + Ant2'/
!      :  'PSCAL6  = ',4x,1pe16.9,1x,'/'/
!      :  'PZERO6  = ',4x,1pe16.9,1x,'/')
! 
! !  End of generating header to scratch file
! 
!       write(iscr,29)
!  29   format('END')
!       endfile(iscr)
! 
! !  Rewind scratch file, copy records to FITS output file
! 
!       rewind(iscr)
!  30   read(iscr,'(A)',end=31) line
!       call phits_write(rline,nline,status)
!       goto 30
! 
! !  End of header, fill current block with blanks
! 
!  31   line=' '
!       call phits_fill(rline,status)
! 
!       close(iscr)
! 
!       end
! 
! !=======================================================================
! *+
!       subroutine phits_fill (data, status)
! !
! !  VSA Project - Reduce program, fill FITS block
! !
! !  Pads the current FITS buffer using input data, and writes to file.
! !
! !  History:  12/7/99 - original version [DJT].
! !
! *-
!        real*4   data
!        integer  status
! 
!        include 'telescope.inc'
! 
!        integer  i, block, bytes, ndata
! 
!        integer  io_wfile
!        external io_wfile
! 
!        if (status.ne.0) return
! 
!        if (ibuff.eq.0) return
! 
! !  Fill current buffer and write to file
! 
!        block = bsize/4
!        ndata = ((ibuff-1)/block+1)*block
! 
!        if (ibuff.lt.ndata) then
!           do i = ibuff+1, ndata
!              buffer(i) = data
!           enddo
!        endif
! 
!        bytes = ndata*4
!        status = io_wfile(filep, buffer, bytes)
!        if (status.ne.0) call io_wrerr(status,' on write')
!        nblock = nblock + ndata/block
!        ibuff = 0
! 
!        end
! 
! !=======================================================================
! !
!       subroutine phits_write(data, nbytes, status)
! !
! !  VSA Project - Reduce program, write data to FITS output file
! !
! !  Writes data to output buffer and thence to file when full.
! !
! !  History:  12/7/99 - original version [DJT].
! !
! *-----------------------------------------------------------------------
! 
!        real*4   data(*)
!        integer  nbytes, status
! 
!        include 'telescope.inc'
! 
!        integer  i, bytes, ndata
! 
!        integer  io_wfile
!        external io_wfile
! 
!        if (status.ne.0) return
! 
!        if (nbytes.le.0) return
! 
! !  Write data to output buffer, write to file when full
! 
!        ndata = (nbytes-1)/4+1
! 
!        do i = 1, ndata
!           ibuff = ibuff+1
!           buffer(ibuff) = data(i)
!           if (ibuff.ge.mbuff) then
!              bytes = mbuff*4
!              status = io_wfile(filep, buffer, bytes)
!              if (status.ne.0) call io_wrerr(status,' on write')
!              nblock = nblock+blocks
!              ibuff = 0
!           endif
!        enddo
! 
!        end
! 
! !=======================================================================
! !
!       subroutine phits_close(status)
! !
! !  VSA Project - Reduce program, close FITS output file
! !
! !  History:  12/7/99 - original version [DJT].
! !
! *-----------------------------------------------------------------------
!        integer  status
! 
!        include 'telescope.inc'
! 
!        character chstr*8
!        integer   iout, ls
! 
!        if (status.ne.0) return
! 
!        call io_enqout(iout)
! 
! !  Close output file
! 
!        call chr_chitoc(nblock,chstr,ls)
! !       write(iout,'(/X,A)')
! !     :            chstr(1:ls)//' FITS blocks of 2880 bytes written'
! 
!        close(ifile)
! 
!        end
! 	 
! !=======================================================================


end module fits_utils
