module fourn1

contains

!----------------------------------------------------------------------
! 2-d FFT from Numerical Recipies. Faster than NAG.
      subroutine fourn(data,nn,ndim,isign)
      double precision wr,wi,wpr,wpi,wtemp,theta
      integer ndim, isign
      integer nn(ndim)
      double precision data(*)
      integer ntot,idim,nprev,n,nrem,ip1,ip2,ip3,i1,i2,i3,i2rev,i3rev,ibit,ifp1,ifp2,k1,k2
      double precision tempr,tempi
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          go to 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*dsin(0.5d0*theta)**2
          wpi=dsin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                if((data(k2).ne.0).or.(data(k2+1).ne.0)) then
                   tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                   tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                else
                   tempr = 0.0
                   tempi = 0.0
                end if
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        go to 2
        endif
        nprev=n*nprev
18    continue
      return
      end subroutine fourn
!----------------------------------------------------------------------

! 2-d FFT from Numerical Recipies. Faster than NAG.

      subroutine fourn2(data,nn,ndim,isign)

	implicit none

	double precision data(*)
	integer ndim,nn(ndim),isign

      double precision wr,wi,wpr,wpi,wtemp,theta
	double precision tempr,tempi
	integer ntot,idim,nprev,n,nrem,ip1,ip2,ip3,i2rev,i2,i1,i3,i3rev
	integer ibit,ifp1,ifp2,k1,k2

! Initialise variables:

	wr = 0.d0
	wi = 0.d0
	wpr = 0.d0
	wpi = 0.d0
	wtemp = 0.d0
	theta = 0.d0
	tempr = 0.0
	tempi = 0.0
	ntot = 0
	idim = 0
	nprev = 0
	n = 0
	nrem = 0
	ip1 = 0
	ip2 = 0
	ip3 = 0
	i2rev = 0
	i2 = 0
	i1 = 0
	i3 = 0
	i3rev = 0
	ibit = 0
	ifp1 = 0
	ifp2 = 0
	k1 = 0
	k2 = 0

! FFT algorithm:

      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          go to 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*dsin(0.5d0*theta)**2
          wpi=dsin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                if((data(k2).ne.0).or.(data(k2+1).ne.0)) then
                   tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                   tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                else
                   tempr = 0.0
                   tempi = 0.0
                end if
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        go to 2
        endif
        nprev=n*nprev
18    continue
      return
      end subroutine fourn2


end module fourn1
