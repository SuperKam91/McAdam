module fits_utils_MRAO
	use telescope1

contains

!=======================================================================
! Cambridge - specific SZ fits file CREATION routines.
!=======================================================================
!
!  VSA Project - Reduce program, write visibilities in FITS format.
!
!  Writes selected samples from the reduced data arrays to a file using
!  the FITS random groups format, suitable for import to AIPS, etc.
!
!  History:  12/7/99 - original version [DJT].
!            1/11/00 - output of antenna tables added [DJT].
!            13/6/01 - added output of report to log file [DJT].
!
! Newer version using values stored in common block in the include file
! (pjm)

      subroutine write_phits(vis_buffer,uv_buffer,
     &                           weight_buffer,baselines,outfile)

      implicit none
	
      character*100 outfile
      double complex vis_buffer(maxvis)
      double precision weight_buffer(maxvis),uv_buffer(2,maxvis)
      integer baselines(maxvis)

      integer status
      integer nbytes,i,iout
      double precision group(9)
	double precision freq,uscale,vscale
	
!-----------------------------------------------------------------------

      if (status.ne.0) return

      call io_enqout(iout)

!  Open the output file

      call phits_open(outfile,'WRITE',status)

!  Write FITS header records

      call phits_wrhdr(status)

!  Write visibilities for each baseline as a single group of data

	freq = CRVAL4
	uscale = PSCAL1
	vscale = PSCAL2
	
      do i = 1,GCOUNT
           
        group(1) = uv_buffer(1,i)/(uscale*freq)
        group(2) = uv_buffer(2,i)/(vscale*freq)
        group(3) = 0.0
        group(4) = 0.0
        group(5) = 0.0
        group(6) = float(baselines(i))
        group(7) = dble(vis_buffer(i))/bscale
        group(8) = aimag(vis_buffer(i))/bscale
        group(9) = weight_buffer(i)/bscale
	  
        nbytes = 9*4
        call phits_write(group,nbytes,status)
		
      enddo

!  End of data, fill current block with zeros

      group(1)=0
      call phits_fill(group,status)

!  Close the file

      call phits_close(status)

      end subroutine write_phit
	
!=======================================================================
!
      subroutine phits_open(file, access, status)

!  VSA Project - Reduce program, open file for UVFITS output
!
!  Open output file, set file and output buffer byte pointers.
!
!  History:  12/7/99 - original version [DJT].
!            25/4/01 - give warning if cnnot open file [KG].
!
*-----------------------------------------------------------------------

       character file*(*), access*(*)
       integer status
       integer bytep
       integer fseek,getfilep
       external  fseek,getfilep

       if (status.ne.0) return

!  Open output file

       call io_nxtlun(ifile, status)

       if (access.eq.'WRITE') then

          open(ifile, file=file, status='UNKNOWN', iostat=status)
          if (status.eq.0) then
             bytep = 0
             status = fseek(ifile,bytep,0)
             filep = getfilep(ifile)
          else 
             write(*,*) 'Cannot open FITS file'
          endif

       endif

       nblock = 0
       ibuff = 0

       end subroutine phits_open
!=======================================================================
!
      subroutine phits_wrhdr(status)
!
!  VSA Project - Reduce program, write FITS header.
!
!  Constructs entries in the FITS header and writes the complete header
!  to disc.
!
!  History:  12/7/99 - original version, adapted from UV2FITS program [DJT].
!            1/11/00 - IF coordinate suppressed [DJT].
!
*-
      integer status
      character scratch*32
      integer   i
      integer   iscr

      integer   nline
      parameter (nline=80)
      character line*(nline)
      double precision    rline(nline/4)
      equivalence (line,rline)

      double precision    one, zero
      data      one, zero / 1.d0, 0.d0 /

      character stokes(7)*4
      data      stokes / 'BEAM', 'I', 'Q', 'U', 'V', 'I-Q', 'I+Q' /

      integer   chr_lenb
      external  chr_lenb

      if (status.ne.0) return

!  Entry, open scratch file

!     scratch = 'scr.out'
!     call io_opefil(iscr,scratch,'write',0,status)
      call io_opescr(iscr,scratch,'write',0,status)

      naxes(1) = 0
      naxes(2) = 3
      naxes(3) = 1
      naxes(4) = 1
      naxes(5) = 1
      naxes(6) = 1
      naxes(7) = 1
      write(iscr,10) bitpix,naxis,(naxes(i),i=1,naxis)
 10   format(
     :  'SIMPLE  = ',19x,'T',1x,'/ Standard FITS format'/
     :  'BITPIX  = ',16x,i4,1x,'/ Bits per pixel'/
     :  'NAXIS   = ',17x,i3,1x,'/ Number of axes'/
     :  'NAXIS1  = ',16x,i4,1x,'/ No image: uv data'/
     :  'NAXIS2  = ',16x,i4,1x,'/ Complex: cos, sin, weight'/
     :  'NAXIS3  = ',16x,i4,1x,'/ Number of polarisations'/
     :  'NAXIS4  = ',16x,i4,1x,'/ Number of frequencies'/
     :  'NAXIS5  = ',16x,i4,1x,'/ RA'/
     :  'NAXIS6  = ',16x,i4,1x,'/ Dec'/
     :  'BLOCKED = ',19x,'T',1x,'/ Tape may be blocked')

!  Random groups format

      write(iscr,11) pcount,gcount
 11   format(
     :  'GROUPS  = ',19x,'T',1x,'/ Groups data structure'/
     :  'PCOUNT  = ',16x,i4,1x,'/ Parameters per group'/
     :  'GCOUNT  = ',14x,i6,1x,'/ Number of groups'/
     :  'EXTEND  = ',19x,'T',1x,'/ Extension is antenna table')

!  Object, telescope, date of observation

      write(iscr,12) object(1:16),telescop(1:3),telescop(1:3)
 12   format(
     :  'OBJECT  = ',1h',a,1h',3x,'/ Source/field name'/
     :  'TELESCOP= ',1h',a,5x,1h',11x,'/ Radio telescope'/
     :  'INSTRUME= ',1h',a,5x,1h',11x,'/ Instrument')
!     :  'DATE-OBS= ',1h',i4,2('-',i2.2),1h',9x,
!     :  'DATE-OBS= ',1h',2(i2.2,'/'),i2.2,1h',8x,
!     :                             '/ Start date of observation')

!  Bunit, bzero, bscale, blank, epoch

      write(iscr,13) bunit,bzero,bscale
 13   format(
     :  'BUNIT   = ',1h',a,1h',11x,'/ Units of data'/
     :  'BZERO   = ',4x,1pe16.9,1x,'/ Data offset'/
     :  'BSCALE  = ',4x,1pe16.9,1x,'/ Data = FITS*BSCALE + BZERO')

!      blank = zero
!      if(blank.ne.0.d0) write(iscr,14) blank
! 14   format(
!     :  'BLANK   = ',4x,1pe16.9,1x,'/ Undefined values on tape')

!      if (obs_epoch.eq.1) epoch=1950.0d0
!      if (obs_epoch.eq.2) epoch=2000.0d0
      write(iscr,15) epoch
 15   format(
     :  'EPOCH   = ',10x,f10.4,1x,'/ Equinox of RA, Dec')

!  Pointing centre

      write(iscr,16) obsra,obsdec
 16   format(
     :  'OBSRA   = ',4x,1pe16.9,1x,'/ Antenna pointing RA'/
     :  'OBSDEC  = ',4x,1pe16.9,1x,'/ Antenna pointing DEC')


!  Details of group data

      write(iscr,17) crval2,cdel2,crpix2,crota2
 17   format(
     :  'CTYPE2  = ',10h'COMPLEX ',11x,'/'/
     :  'CRVAL2  = ',4x,1pe16.9,1x,'/'/
     :  'CDELT2  = ',4x,1pe16.9,1x,'/'/
     :  'CRPIX2  = ',4x,1pe16.9,1x,'/'/
     :  'CROTA2  = ',4x,1pe16.9,1x,'/')

      write(iscr,18) stokes(1),crval3,cdel3,crpix3,crota3
 18   format(
     :  'CTYPE3  = ',10h'STOKES  ',11x,'/ Stokes parameter : ',A/
     :  'CRVAL3  = ',4x,1pe16.9,1x,'/'/
     :  'CDELT3  = ',4x,1pe16.9,1x,'/'/
     :  'CRPIX3  = ',4x,1pe16.9,1x,'/'/
     :  'CROTA3  = ',4x,1pe16.9,1x,'/')

      write(iscr,19) crval4,cdel4,crpix4,crota4
 19   format(
     :  'CTYPE4  = ',10h'FREQ    ',11x,'/ Frequency, Hertz'/
     :  'CRVAL4  = ',4x,1pe16.9,1x,'/'/
     :  'CDELT4  = ',4x,1pe16.9,1x,'/'/
     :  'CRPIX4  = ',4x,1pe16.9,1x,'/'/
     :  'CROTA4  = ',4x,1pe16.9,1x,'/')

!     write(iscr,195) one,one,one,zero
c195  format(
!    :  'CTYPE5  = ',10h'IF      ',11x,'/'/
!    :  'CRVAL5  = ',4x,1pe16.9,1x,'/'/
!    :  'CDELT5  = ',4x,1pe16.9,1x,'/'/
!    :  'CRPIX5  = ',4x,1pe16.9,1x,'/'/
!    :  'CROTA5  = ',4x,1pe16.9,1x,'/')

      write(iscr,20) crval5,cdel5,crpix5,crota5
 20   format(
     :  'CTYPE5  = ',10h'RA      ',11x,'/ Right Ascension, degrees'/
     :  'CRVAL5  = ',4x,1pe16.9,1x,'/'/
     :  'CDELT5  = ',4x,1pe16.9,1x,'/'/
     :  'CRPIX5  = ',4x,1pe16.9,1x,'/'/
     :  'CROTA5  = ',4x,1pe16.9,1x,'/')

      write(iscr,21) crval6,cdel6,crpix6,crota6
 21   format(
     :  'CTYPE6  = ',10h'DEC     ',11x,'/ Declination, degrees'/
     :  'CRVAL6  = ',4x,1pe16.9,1x,'/'/
     :  'CDELT6  = ',4x,1pe16.9,1x,'/'/
     :  'CRPIX6  = ',4x,1pe16.9,1x,'/'/
     :  'CROTA6  = ',4x,1pe16.9,1x,'/')

!  Details of group parameters

      write(iscr,22) pscal1,pzero1
 22   format(
     :  'PTYPE1  = ',10h'UU      ',11x,'/ U coordinate, seconds'/
     :  'PSCAL1  = ',4x,1pe16.9,1x,'/ Real = FITS*PSCAL + PZERO'/
     :  'PZERO1  = ',4x,1pe16.9,1x,'/')

      write(iscr,23) pscal2,pzero2
 23   format(
     :  'PTYPE2  = ',10h'VV      ',11x,'/ V coordinate, seconds'/
     :  'PSCAL2  = ',4x,1pe16.9,1x,'/ Real = FITS*PSCAL + PZERO'/
     :  'PZERO2  = ',4x,1pe16.9,1x,'/')

      write(iscr,24) pscal3,pzero3
 24   format(
     :  'PTYPE3  = ',10h'WW      ',11x,'/ W coordinate, seconds'/
     :  'PSCAL3  = ',4x,1pe16.9,1x,'/ Real = FITS*PSCAL + PZERO'/
     :  'PZERO3  = ',4x,1pe16.9,1x,'/')

!  Date parameter, as two entries:
!   ptype1: pzero1 = intpt(jd at start-time of observations)
!           pscal1 = 1/128 days
!           value of type1 parameter = jd of sample to nearest
!              integral value of 1/128 days
!   ptype2: pzero2 = 0
!           pscal2 = 1/6e6; precision is about 0.01 second

!      jdzero = int(obs_mjd_start+2400000.5d0)
!      jdscal1 = 1.d0/128.d0
!      jdscal2 = 1.d0/6.d6
!      write(iscr,25) jdscal1,jdzero,jdscal2,zero
      write(iscr,25) pscal4,pzero4,pscal5,pzero5
 25   format(
     :  'PTYPE4  = ',10h'DATE    ',11x,'/ Julian Date' /
     :  'PSCAL4  = ',4x,1pe16.9,1x,'/'/
     :  'PZERO4  = ',4x,1pe16.9,1x,'/'/
     :  'PTYPE5  = ',10h'DATE    ',11x,'/ Julian Date'/
     :  'PSCAL5  = ',4x,1pe16.9,1x,'/'/
     :  'PZERO5  = ',4x,1pe16.9,1x,'/')

      write(iscr,26) pscal6,pzero6
 26   format(
     :  'PTYPE6  = ',10h'BASELINE',11x,'/ Ant1*256 + Ant2'/
     :  'PSCAL6  = ',4x,1pe16.9,1x,'/'/
     :  'PZERO6  = ',4x,1pe16.9,1x,'/')

!  End of generating header to scratch file

      write(iscr,29)
 29   format('END')
      endfile(iscr)

!  Rewind scratch file, copy records to FITS output file

      rewind(iscr)
 30   read(iscr,'(A)',end=31) line
      call phits_write(rline,nline,status)
      goto 30

!  End of header, fill current block with blanks

 31   line=' '
      call phits_fill(rline,status)

      close(iscr)

      end subroutine phits_wrhdr

!=======================================================================
*+
      subroutine phits_fill (data, status)
!
!  VSA Project - Reduce program, fill FITS block
!
!  Pads the current FITS buffer using input data, and writes to file.
!
!  History:  12/7/99 - original version [DJT].
!
*-
       double precision   data
       integer  status
       integer  i, block, bytes, ndata
       integer  io_wfile
       external io_wfile

       if (status.ne.0) return

       if (ibuff.eq.0) return

!  Fill current buffer and write to file

       block = bsize/4
       ndata = ((ibuff-1)/block+1)*block

       if (ibuff.lt.ndata) then
          do i = ibuff+1, ndata
             buffer(i) = data
          enddo
       endif

       bytes = ndata*4
       status = io_wfile(filep, buffer, bytes)
       if (status.ne.0) call io_wrerr(status,' on write')
       nblock = nblock + ndata/block
       ibuff = 0

       end subroutine phits_fill

!=======================================================================
!
      subroutine phits_write(data, nbytes, status)
!
!  VSA Project - Reduce program, write data to FITS output file
!
!  Writes data to output buffer and thence to file when full.
!
!  History:  12/7/99 - original version [DJT].
!
*-----------------------------------------------------------------------

       double precision   data(*)
       integer  nbytes, status
       integer  i, bytes, ndata
       integer  io_wfile
       external io_wfile

       if (status.ne.0) return

       if (nbytes.le.0) return

!  Write data to output buffer, write to file when full

       ndata = (nbytes-1)/4+1

       do i = 1, ndata
          ibuff = ibuff+1
          buffer(ibuff) = data(i)
          if (ibuff.ge.mbuff) then
             bytes = mbuff*4
             status = io_wfile(filep, buffer, bytes)
             if (status.ne.0) call io_wrerr(status,' on write')
             nblock = nblock+blocks
             ibuff = 0
          endif
       enddo

       end subroutine phits_write

!=======================================================================
!
      subroutine phits_close(status)
!
!  VSA Project - Reduce program, close FITS output file
!
!  History:  12/7/99 - original version [DJT].
!
*-----------------------------------------------------------------------
       integer  status
       character chstr*8
       integer   iout, ls

       if (status.ne.0) return

       call io_enqout(iout)

!  Close output file

       call chr_chitoc(nblock,chstr,ls)
!       write(iout,'(/X,A)')
!     :            chstr(1:ls)//' FITS blocks of 2880 bytes written'

       close(ifile)

       end subroutine phits_close
	 
!=======================================================================

end module fits_utils_MRAO
