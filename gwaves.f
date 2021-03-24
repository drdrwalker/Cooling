      subroutine gwaves(igrd)
c     calculate quantities for gravitational waves
 
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'
 
      real*8 tmp1(q), tmp2(q), rxyz(3), rkpm
      real*8 g11, g12, g13, g21, g22, g23, g31, g32, g33

      real*8 vxyz(qx,qy,qz,ngrd,3)
      equivalence ( vxyz(1,1,1,1,1), velx(1,1,1,1) ) ,
     &            ( vxyz(1,1,1,1,2), vely(1,1,1,1) ) ,
     &            ( vxyz(1,1,1,1,3), velz(1,1,1,1) )

      integer*4 igrd
c *** commented in to test influence on gw lum structure
c      call gpress  ( igrd )      ! commented out to save CPU time 
c *** --------------------------------------------------
      call dtgravty(igrd)
      itopgr = indxgr(lone,lone)
      cmgaz  = dble( 0/2) + half    !  cmgax,cmgay given in compct.cmn
  
      do j = 1, 3
      do i = 1, 3
         gdij3dt(i,j,igrd) = zero
      enddo
      enddo

      g11=zero
      g12=zero
      g13=zero
      g21=zero
      g22=zero
      g23=zero
      g31=zero
      g32=zero
      g33=zero

      d1x   = delx(igrd)
      d3x   = d1x**3

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i), 
C$OMP+            SHARED(nx,ny,nz,ro,po,d3x)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
c         po(i,j,k) = ro(i,j,k)   ! set in dtgravty !!
         ro(i,j,k) = two * d3x
c     one factor two because of symmetry about z-plane -> double volume
      enddo 
      enddo 
      enddo 
C$OMP END PARALLEL DO

c     zero area covered by finer grids
      ifgr = idxfgr(igrd)
      if( ifgr .ne. -lone ) call rozero(ifgr)

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+       PRIVATE(j,jp,jm,k,kp,km,i,ip,im,rxyz,divrhov,d1xd,d1yd,
C$OMP+              rkpm,rd1zd,rd1yd,rd1xd,mxyz,tmp1,tmp2,d1zd,sg,gadd),
C$OMP+       SHARED(nx,ny,nz,d1x,igrd,cmgax,cmgay,cmgaz,g13,g23,g31,g32,
C$OMP+              densty,vxyz,gpot,gtpot,press,ro,velx,vely,velz),
C$OMP+       REDUCTION(+:g11,g12,g21,g22,g33) 
      do j = 1, ny
        jp = min ( j+lone, ny)
        jm = max ( j-lone, lone )
        rxyz(2) = (dble(j)-cmgay) * d1x
      do k = 1, nz
        sg = +one
        if ( k .eq. lone )  sg = -one    ! reflection of velz at k=half
        kp = min ( k+lone, nz)
        km = max ( k-lone, lone )
        rxyz(3) = (dble(k)-cmgaz) * d1x
        rkpm = dble(kp-km)
        if ( k .eq. lone ) rkpm = two   ! reflection at k=0.5
      do i = 1, nx
        ip = min ( i+lone, nx)
        im = max ( i-lone, lone )
        rxyz(1) = (dble(i)-cmgax) * d1x

        d1zd = one / ( d1x * rkpm )
        d1yd = one / ( d1x * dble(jp-jm) )
        d1xd = one / ( d1x * dble(ip-im) )
 
        divrhov = 
     &    (   densty(ip,j,k,igrd)*velx(ip,j,k,igrd)
     &       -densty(im,j,k,igrd)*velx(im,j,k,igrd) )*d1xd
     &   +(   densty(i,jp,k,igrd)*vely(i,jp,k,igrd)
     &       -densty(i,jm,k,igrd)*vely(i,jm,k,igrd) )*d1yd
     &   +(   densty(i,j,kp,igrd)*velz(i,j,kp,igrd)
     &       -densty(i,j,km,igrd)*velz(i,j,km,igrd) )*d1zd

        rd1zd = ro(i,j,k) * d1zd
        rd1yd = ro(i,j,k) * d1yd
        rd1xd = ro(i,j,k) * d1xd
 
c *** 
        mxyz = lone    ! unrolled loop
 
        tmp1(i) =  - two * densty(i,j,k,igrd) * vxyz(i,j,k,igrd,mxyz)
     &                                + rxyz(mxyz) * divrhov
        tmp2(i) =  - densty(i,j,k,igrd) * rxyz(mxyz)

c  -- derivative in x direction --
        gadd = rd1xd* ( two * press(i,j,k,igrd)
     &               * (vxyz(ip,j,k,igrd,mxyz) - vxyz(im,j,k,igrd,mxyz))
c     &   - two * vxyz(i,j,k,igrd,mxyz)
c     &          * (press(ip,j,k,igrd) - press(im,j,k,igrd))
     &   + tmp1(i) * (  gpot(ip,j,k,igrd) -  gpot(im,j,k,igrd) )
     &   + tmp2(i) * ( gtpot(ip,j,k,igrd) - gtpot(im,j,k,igrd) ) )
        g11 = g11  +  gadd

c  -- derivative in y direction --

        gadd = rd1yd * (
     &   + two * press(i,j,k,igrd)
     &        * (vxyz(i,jp,k,igrd,mxyz) - vxyz(i,jm,k,igrd,mxyz))
c     &   - two * vxyz(i,j,k,igrd,mxyz)
c     &          * (press(i,jp,k,igrd) - press(i,jm,k,igrd))
     &   + tmp1(i) * (  gpot(i,jp,k,igrd) -  gpot(i,jm,k,igrd) )
     &   + tmp2(i) * ( gtpot(i,jp,k,igrd) - gtpot(i,jm,k,igrd) ) )
        g12 = g12  +  gadd

c  -- derivative in z direction --
 
        g13 = zero
cs       gdij3dt(mxyz,3,igrd) = gdij3dt(mxyz,3,igrd)  +  rd1zd * (
ccs    &   + two * press(i,j,k,igrd)
ccs    &          * (  vxyz(i,j,kp,igrd,mxyz) - sg*vxyz(i,j,km,igrd,mxyz))
cs    &   - two*vxyz(i,j,k,igrd,mxyz)*(press(i,j,kp,igrd)-press(i,j,km,igrd))
cs   &   + tmp1(i) * (  gpot(i,j,kp,igrd) -  gpot(i,j,km,igrd) )
cs   &   + tmp2(i) * ( gtpot(i,j,kp,igrd) - gtpot(i,j,km,igrd) ) )
 
c ***
        mxyz = ltwo    ! unrolled loop
 
        tmp1(i) =  - two * densty(i,j,k,igrd) * vxyz(i,j,k,igrd,mxyz)
     &                                + rxyz(mxyz) * divrhov
        tmp2(i) =  - densty(i,j,k,igrd) * rxyz(mxyz)

c  -- derivative in x direction --
 
        
        gadd = rd1xd * (
     &   + two * press(i,j,k,igrd)
     &        * (vxyz(ip,j,k,igrd,mxyz) - vxyz(im,j,k,igrd,mxyz))
c     &   - two * vxyz(i,j,k,igrd,mxyz)
c     &          * (press(ip,j,k,igrd) - press(im,j,k,igrd))
     &   + tmp1(i) * (  gpot(ip,j,k,igrd) -  gpot(im,j,k,igrd) )
     &   + tmp2(i) * ( gtpot(ip,j,k,igrd) - gtpot(im,j,k,igrd) ) )
        g21 = g21  +  gadd

c  -- derivative in y direction --
        
       gadd = rd1yd * (
     &   + two * press(i,j,k,igrd)
     &        * (vxyz(i,jp,k,igrd,mxyz) - vxyz(i,jm,k,igrd,mxyz))
c     &   - two * vxyz(i,j,k,igrd,mxyz)
c     &          * (press(i,jp,k,igrd) - press(i,jm,k,igrd))
     &   + tmp1(i) * (  gpot(i,jp,k,igrd) -  gpot(i,jm,k,igrd) )
     &   + tmp2(i) * ( gtpot(i,jp,k,igrd) - gtpot(i,jm,k,igrd) ) )
       g22 = g22  +  gadd

c  -- derivative in z direction --
 
        g23 = zero
cs      gdij3dt(mxyz,3,igrd) = gdij3dt(mxyz,3,igrd)  +  rd1zd * (
ccs   &   + two * press(i,j,k,igrd)
ccs   &        * (  vxyz(i,j,kp,igrd,mxyz) - sg*vxyz(i,j,km,igrd,mxyz))
cs    &   - two*vxyz(i,j,k,igrd,mxyz)*(press(i,j,kp,igrd)-press(i,j,km,igrd))
cs   &   + tmp1(i) * (  gpot(i,j,kp,igrd) -  gpot(i,j,km,igrd) )
cs   &   + tmp2(i) * ( gtpot(i,j,kp,igrd) - gtpot(i,j,km,igrd) ) )
 
c ***
        mxyz = lthree    ! unrolled loop
c       because of symmetry about z=0.5, two parts have to be added:
c       one with sgp = +1 and the other with sgp = -1 , with
c       vxyz -> sgp*vxyz  and   rxyz -> sgp*rxyz    => dij3dt = 0
  
        tmp1(i) =  - two * densty(i,j,k,igrd) * vxyz(i,j,k,igrd,mxyz)
     &                                + rxyz(mxyz) * divrhov
        tmp2(i) =  - densty(i,j,k,igrd) * rxyz(mxyz)

c  -- derivative in x direction --
 
        g31 = zero
cs      gdij3dt(mxyz,1,igrd) = gdij3dt(mxyz,1,igrd)  +  rd1xd * (
ccs   &   + two * press(i,j,k,igrd)
ccs           * (vxyz(ip,j,k,igrd,mxyz) - vxyz(im,j,k,igrd,mxyz))
cs     &   - two*vxyz(i,j,k,igrd,mxyz)*(press(ip,j,k,igrd)-press(im,j,k,igrd))
cs   &   + tmp1(i) * (  gpot(ip,j,k,igrd) -  gpot(im,j,k,igrd) )
cs   &   + tmp2(i) * ( gtpot(ip,j,k,igrd) - gtpot(im,j,k,igrd) ) )
 
c  -- derivative in y direction --
 
        g32 = zero
cs      gdij3dt(mxyz,2,igrd) = gdij3dt(mxyz,2,igrd)  +  rd1yd * (
ccs   &   + two * press(i,j,k,igrd)
ccs           * (vxyz(i,jp,k,igrd,mxyz) - vxyz(i,jm,k,igrd,mxyz))
cs     &   - two*vxyz(i,j,k,igrd,mxyz)*(press(i,jp,k,igrd)-press(i,jm,k,igrd))
cs   &   + tmp1(i) * (  gpot(i,jp,k,igrd) -  gpot(i,jm,k,igrd) )
cs   &   + tmp2(i) * ( gtpot(i,jp,k,igrd) - gtpot(i,jm,k,igrd) ) )
 
c  -- derivative in z direction --
 
        
      gadd = rd1zd * (
     &   + two * press(i,j,k,igrd) * (vxyz(i,j,kp,igrd,mxyz)
     &        -sg * vxyz(i,j,km,igrd,mxyz))
c     &   - two * vxyz(i,j,k,igrd,mxyz)
c     &          * (press(i,j,kp,igrd) - press(i,j,km,igrd))
     &   + tmp1(i) * (  gpot(i,j,kp,igrd) -  gpot(i,j,km,igrd) )
     &   + tmp2(i) * ( gtpot(i,j,kp,igrd) - gtpot(i,j,km,igrd) ) )
      g33 = g33  +  gadd

      enddo
      enddo
      enddo

C$OMP END PARALLEL DO
 
      gdij3dt(1,1,igrd) = g11
      gdij3dt(1,2,igrd) = g12
      gdij3dt(1,3,igrd) = g13
      gdij3dt(2,1,igrd) = g21
      gdij3dt(2,2,igrd) = g22
      gdij3dt(2,3,igrd) = g23
      gdij3dt(3,1,igrd) = g31
      gdij3dt(3,2,igrd) = g32
      gdij3dt(3,3,igrd) = g33

      do j = 1, 3
      do i = 1, 3
         dij3dt(i,j) = zero
         do lgrd = 1, ngrd
            dij3dt(i,j) = dij3dt(i,j) + two * gdij3dt(i,j,lgrd)
         enddo 
      enddo 
      enddo 
c     factor two from equation
 
cw    write(*,*) 'dij3dt:'
cw    write(*,'(1P,(3E13.5))') dij3dt
 
c     make dij3dt symmetric and trace free
 
      dijtrace = ( dij3dt(1,1) + dij3dt(2,2) + dij3dt(3,3) ) /3.D0
      dij3dt(1,1) = dij3dt(1,1) - dijtrace
      dij3dt(2,2) = dij3dt(2,2) - dijtrace
      dij3dt(3,3) = dij3dt(3,3) - dijtrace
      dij3dt(1,2) = ( dij3dt(1,2) + dij3dt(2,1) ) / two
      dij3dt(2,1) =   dij3dt(1,2)
      dij3dt(1,3) = ( dij3dt(1,3) + dij3dt(3,1) ) / two
      dij3dt(3,1) =   dij3dt(1,3)
      dij3dt(2,3) = ( dij3dt(2,3) + dij3dt(3,2) ) / two
      dij3dt(3,2) =   dij3dt(2,3)

      if ( mgonoff .eq. lzero )  then
         do j = 1, 3
         do i = 1, 3
            dij3dt(i,j) = zero       ! test
         enddo
         enddo
      endif

      if ( mgonoff .eq. ltwo ) then
         do j = 1, 3
         do i = 1, 3
            dij3dt(i,j) = dij3dt(i,j) * two
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------
      subroutine gpress(igrd)
c     find the pressure needed for the gravitational waves calculations
 
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
 
      real*8 rho(q),ek(q),e(q),ei(q),tmp(q),p(q),gamc(q),game(q),ye(q)
      common /tmp1/ rho,ek,e,ei,tmp,p,gamc,game,ye
C$OMP THREADPRIVATE (/tmp1/)

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+            PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,densty,velx,vely,velz,energy,
C$OMP+                   temper,press,chem,small)
      do k = 1, nz 
      do j = 1, ny
         do  i = 1, nx
            rho   (i) = densty(i,j,k,igrd)
            ek    (i) =  half *   ( velx(i,j,k,igrd)**2
     &                 + vely(i,j,k,igrd)**2 + velz(i,j,k,igrd)**2 )
 
            e     (i) = energy(i,j,k,igrd)
            ei    (i) = max( rho(i)*(e(i)-ek(i)), rho(i)*small*ek(i) )
            ye    (i) = chem(i,j,k,igrd,lone) / rho(i)
         enddo

         call eos( rho, tmp, ei, p, ye, gamc, game, nx, lzero )
 
         do  i = 1, nx
            temper(i,j,k,igrd)   = tmp(i)
            press (i,j,k,igrd)   = p(i)
            chem  (i,j,k,igrd,lone) = ye(i) * rho(i)
            energy(i,j,k,igrd) = max(ei(i)/rho(i),small*ek(i)) + ek(i)
         enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c --- reset pressure in the BH to a small value
c --- pressure in BH is not physical
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,bhdens,press)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         if ( bhdens(i,j,k,igrd) .ne. zero ) 
     &             press(i,j,k,igrd) = 1.D-5 * press(i,j,k,igrd)
        
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c-----------------------------------------------------------------------
      subroutine der2ro(igrd)
 
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
 
      real*8 rxyz(3), di, dj, dk
 
      d1x    = delx(igrd)

      cmgaz = dble( 0/2) + half    !  cmgax,cmgay given in compct.cmn
 
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+            PRIVATE(j,jp,jm,dj,k,kp,km,dk,i,ip,im,di,rxyz),
C$OMP+      SHARED(nx,ny,nz,d1x,cmgax,cmgay,cmgaz,ro,densty,dij3dt,igrd)
      do j = 1, ny
         jp = min ( j+lone, ny)
         jm = max ( j-lone, lone )
         dj = dble(jp - jm) * d1x
         rxyz(2) = (dble(j)-cmgay) * d1x
      do k = 1, nz
         kp = min ( k+lone, nz)
         km = max ( k-lone, lone )
         dk = dble(kp - km) * d1x
         if ( k .eq. lone ) dk = two * d1x !  reflection at k = 0.5
         rxyz(3) = (dble(k)-cmgaz) * d1x
      do i = 1, nx
         ip = min ( i+lone, nx)
         im = max ( i-lone, lone )
         di = dble(ip - im) * d1x
         rxyz(1) = (dble(i)-cmgax) * d1x
 
         ro(i,j,k) =
     &       (    dij3dt(1,1) * rxyz(1)
     &         +  dij3dt(2,1) * rxyz(2)
     &         +  dij3dt(3,1) * rxyz(3)  ) *
     &           ((densty(ip,j, k, igrd) - densty(im,j, k, igrd)) / di)
     &    +  (    dij3dt(1,2) * rxyz(1)
     &         +  dij3dt(2,2) * rxyz(2)
     &         +  dij3dt(3,2) * rxyz(3)  ) *
     &           ((densty(i, jp,k, igrd) - densty(i, jm,k, igrd)) / dj)
     &    +  (    dij3dt(1,3) * rxyz(1)
     &         +  dij3dt(2,3) * rxyz(2)
     &         +  dij3dt(3,3) * rxyz(3)  ) *
     &           ((densty(i, j, kp,igrd) - densty(i, j, km,igrd)) / dk)
 
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO 
 
      return
      end
c-----------------------------------------------------------------------
      subroutine adbakr(igrd,ifgr)
 
c     add background potential of coarse grid to fine grid
 
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
 
      noz = norgin(ifgr,lthree)
      noy = norgin(ifgr,ltwo)
      nox = norgin(ifgr,lone)
 
c do one less and more for tripo
c     do 10 k = noz, noz+1+nz/nfine
c exceptional: tripo contains refecting z=0 boundary

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i), 
C$OMP+            SHARED(igrd,nox,noy,noz,nx,ny,nz,po,delx,bakrg)
      do j = noy,   noy+1+ny/nfine
      do k = noz+1, noz+1+nz/nfine
      do i = nox,   nox+1+nx/nfine
         po(i,j,k) = po(i,j,k) * delx(igrd)**2  + bakrg(i,j,k,igrd)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
 
      return
      end
c-----------------------------------------------------------------------
      subroutine adbakcr(igrd)
 
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
 
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i), 
C$OMP+            SHARED(nx,ny,nz,po,rgarr,delx,bakrg,igrd)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         rgarr(i,j,k,igrd) = po(i,j,k) * delx(igrd)**2  
     &                                          + bakrg(i,j,k,igrd)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
 
      return
      end
c-----------------------------------------------------------------------
      subroutine reactphi(igrd)
 
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'
 
      real*8 rxyz(3), di, dj, dk
 
      d1x    =  delx(igrd)

      cmgaz = dble( 0/2) + half    !  cmgax,cmgay given in compct.cmn
 
c     frgx = zero
c     frgy = zero
 
c     do 11 k = 1, nz
c       kp = min ( k+lone, nz)
c       km = max ( k-lone, lone )
c     do 11 j = 1, ny
c       jp = min ( j+lone, ny)
c       jm = max ( j-lone, lone )
c     do 11 i = 1, nx/ltwo
c       ip = min ( i+lone, nx)
c       im = max ( i-lone, lone )
 
c       rpxtmp = rgarr(ip,j, k ) - rgarr(im,j, k )
c       rpytmp = rgarr(i, jp,k ) - rgarr(i, jm,k )
 
c       frgx = frgx + densty(i,j,k,igrd) * rpxtmp
c       frgy = frgy + densty(i,j,k,igrd) * rpytmp
 
c 11  continue
 
c     frgx = frgx*two/d1xd * 0.4D0 * g / cc**5 * d1x**3
c     frgy = frgy*two/d1xd * 0.4D0 * g / cc**5 * d1x**3
 
c     write(*,*)  'reactphi   '
c     write(*,*)  'frgx: ',frgx, '   frgy: ',frgy
 
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+            PRIVATE(j,jp,jm,dj,k,kp,km,dk,i,ip,im,di,rxyz),
C$OMP+            SHARED(nx,ny,nz,d1x,cmgax,cmgay,cmgaz,phireac,g,cc,
C$OMP+                   rgarr,dij3dt,gpot,igrd)
      do j = 1, ny
         jp = min ( j+lone, ny)
         jm = max ( j-lone, lone )
         dj = dble(jp - jm) * d1x
         rxyz(2) = (dble(j)-cmgay) * d1x
      do k = 1, nz    !k=1 symmetry
         kp = min ( k+lone, nz)
         km = max ( k-lone, lone )
         dk = dble(kp - km) * d1x
         if ( k .eq. lone ) dk = two * d1x
         rxyz(3) = (dble(k)-cmgaz) * d1x
      do i = 1, nx
         ip = min ( i+lone, nx)
         im = max ( i-lone, lone )
         di = dble(ip - im) * d1x
         rxyz(1) = (dble(i)-cmgax) * d1x
 
         phireac(i,j,k) = 0.4D0*g/cc**5 * ( rgarr(i,j,k,igrd) - 
     &    ( +  (    dij3dt(1,1) * rxyz(1)
     &           +  dij3dt(2,1) * rxyz(2)
     &           +  dij3dt(3,1) * rxyz(3)  ) * 
     &             ((gpot(ip,j, k, igrd) - gpot(im,j, k, igrd)) / di)
     &      +  (    dij3dt(1,2) * rxyz(1)
     &           +  dij3dt(2,2) * rxyz(2)
     &           +  dij3dt(3,2) * rxyz(3)  ) *
     &             ((gpot(i, jp,k, igrd) - gpot(i, jm,k, igrd)) / dj)
     &      +  (    dij3dt(1,3) * rxyz(1)
     &           +  dij3dt(2,3) * rxyz(2)
     &           +  dij3dt(3,3) * rxyz(3)  ) * 
     &             ((gpot(i, j, kp,igrd) - gpot(i, j, km,igrd)) / dk) ))
 
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO 
 
c --- the following is only to check
 
c     fpotx = zero
c     fpoty = zero
 
c     frgx = zero
c     frgy = zero
 
c     do 12 k = 1, nz
c       kp = min ( k+lone, nz)
c       km = max ( k-lone, lone )
c     do 12 j = 1, ny
c       jp = min ( j+lone, ny)
c       jm = max ( j-lone, lone )
c     do 12 i = 1, nx/ltwo
c       ip = min ( i+lone, nx)
c       im = max ( i-lone, lone )
 
c       gpxtmp =  gpot(ip,j, k, igrd) - gpot(im,j, k, igrd)
c       gpytmp =  gpot(i, jp,k, igrd) - gpot(i, jm,k, igrd)
 
c       fpotx = fpotx + densty(i,j,k,igrd) * gpxtmp
c       fpoty = fpoty + densty(i,j,k,igrd) * gpytmp
 
c       rpxtmp = rgarr(ip,j, k ) - rgarr(im,j, k )
c       rpytmp = rgarr(i, jp,k ) - rgarr(i, jm,k )
 
c       frgx = frgx + densty(i,j,k,igrd) * rpxtmp
c       frgy = frgy + densty(i,j,k,igrd) * rpytmp
 
c 12  continue
 
c     fpotx = fpotx*two/d1xd * d1x**3
c     fpoty = fpoty*two/d1xd * d1x**3
 
c     frgx = frgx*two/d1xd * d1x**3
c     frgy = frgy*two/d1xd * d1x**3
 
c     write(*,*)  'reactphi   '
c     write(*,*)  'fpotx:',fpotx,'   fpoty:',fpoty
c     write(*,*)  'frgx: ',frgx, '   frgy: ',frgy
c     stop
 
      return
      end
c-----------------------------------------------------------------------
      subroutine addpress(igrd)
 
c     calculate additional pressure term for energy update
 
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'
 
      real*8 di, dj, dk
 
      real*8 vxyz(qx,qy,qz,ngrd,3)
      equivalence ( vxyz(1,1,1,1,1), velx(1,1,1,1) ) ,
     &            ( vxyz(1,1,1,1,2), vely(1,1,1,1) ) ,
     &            ( vxyz(1,1,1,1,3), velz(1,1,1,1) )


      d1x    =  delx(igrd)

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+            PRIVATE(j,jp,jm,dj,k,kp,km,dk,i,ip,im,di,
C$OMP+                    pxtmp,pytmp,pztmp,rpxtm,rpytm,rpztm),
C$OMP+            SHARED(ny,nz,nx,d1x,igrd,cc,g,
C$OMP+                   press,densty,gpot,dij3dt,vxyz,adderg)
      do j = 1, ny
        jp = min ( j+lone, ny)
        jm = max ( j-lone, lone )
        dj = dble(jp-jm) * d1x
      do k = 1, nz
        kp = min ( k+lone, nz)
        km = max ( k-lone, lone )
        dk = dble(kp-km) * d1x
        if ( k .eq. lone ) dk = two * d1x
      do i = 1, nx
        ip = min ( i+lone, nx)
        im = max ( i-lone, lone )
        di = dble(ip-im) * d1x
          
        pxtmp = (press(ip,j, k ,igrd) - press(im,j, k ,igrd)) / di
        pytmp = (press(i, jp,k ,igrd) - press(i, jm,k ,igrd)) / dj
        pztmp = (press(i, j, kp,igrd) - press(i, j, km,igrd)) / dk
        
        rpxtm = densty(i,j,k,igrd) 
     &        * (gpot   (ip,j,k,igrd) - gpot   (im,j,k,igrd)) / di
        rpytm = densty(i,j,k,igrd) 
     &        * (gpot   (i,jp,k,igrd) - gpot   (i,jm,k,igrd)) / dj
        rpztm = densty(i,j,k,igrd) 
     &        * (gpot   (i,j,kp,igrd) - gpot   (i,j,km,igrd)) / dk

        adderg(i,j,k) = 0.8D0 * g / cc**5   * (
     &      (   dij3dt(1,1) * vxyz(i,j,k,igrd,1)
     &        + dij3dt(2,1) * vxyz(i,j,k,igrd,2)
     &        + dij3dt(3,1) * vxyz(i,j,k,igrd,3) ) *(pxtmp+rpxtm)
     &    + (   dij3dt(1,2) * vxyz(i,j,k,igrd,1)
     &        + dij3dt(2,2) * vxyz(i,j,k,igrd,2)
     &        + dij3dt(3,2) * vxyz(i,j,k,igrd,3) ) *(pytmp+rpytm)
     &    + (   dij3dt(1,3) * vxyz(i,j,k,igrd,1)
     &        + dij3dt(2,3) * vxyz(i,j,k,igrd,2)
     &        + dij3dt(3,3) * vxyz(i,j,k,igrd,3) ) *(pztmp+rpztm) )
 
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c --- addition 12/09/97
      ifgr = idxfgr(igrd)

      if ( ifgr .ne. -1 ) then

         noz = norgin(ifgr,lthree)
         noy = norgin(ifgr,ltwo)
         nox = norgin(ifgr,lone)
 
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nox,noy,noz,nx,ny,nz,adderg)
         do j = noy+1, noy+ny/nfine
         do k = noz+1, noz+nz/nfine
         do i = nox+1, nox+nx/nfine
            adderg(i,j,k) = zero
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine addchi(igrd)

c     calculate additional pressure terms for energy update
 
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'
 
      real*8 di, dj, dk, dgi, dgj, dgk
 
      real*8 vxyz(qx,qy,qz,ngrd,3)
      equivalence ( vxyz(1,1,1,1,1), velx(1,1,1,1) ) ,
     &            ( vxyz(1,1,1,1,2), vely(1,1,1,1) ) ,
     &            ( vxyz(1,1,1,1,3), velz(1,1,1,1) )
 
      d1x    =  delx(igrd)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,jp,jm,dj,k,kp,km,dk,
C$OMP+                   i,ip,im,di,pxtmp,pytmp,pztmp,dgi,dgj,dgk),
C$OMP+            SHARED(nx,ny,nz,d1x,press,igrd,gpot,adderg,dij3dt,
C$OMP+                   vxyz,corrmax,g,cc)
      do j = 1, ny
        jp = min ( j+lone, ny)
        jm = max ( j-lone, lone )
        dj = dble(jp-jm) * d1x
      do k = 1, nz
        kp = min ( k+lone, nz)
        km = max ( k-lone, lone )
        dk = dble(kp-km) * d1x
        if ( k .eq. lone )  dk = two * d1x      ! reflection 
      do i = 1, nx
        ip = min ( i+lone, nx)
        im = max ( i-lone, lone )
        di = dble(ip-im) * d1x
 
        pxtmp = (press(ip,j, k ,igrd) - press(im,j, k ,igrd)) / di
        pytmp = (press(i, jp,k ,igrd) - press(i, jm,k ,igrd)) / dj
        pztmp = (press(i, j, kp,igrd) - press(i, j, km,igrd)) / dk

        dgi = two / cc**2 * (gpot(ip,j,k,igrd) - gpot(im,j,k,igrd))
        dgj = two / cc**2 * (gpot(i,jp,k,igrd) - gpot(i,jm,k,igrd))
        dgk = two / cc**2 * (gpot(i,j,kp,igrd) - gpot(i,j,km,igrd))

        adderg(i,j,k) = 0.8D0 * g / cc**5   * (
     &         (   dij3dt(1,1) * vxyz(i,j,k,igrd,1)
     &           + dij3dt(2,1) * vxyz(i,j,k,igrd,2)
     &           + dij3dt(3,1) * vxyz(i,j,k,igrd,3) ) * pxtmp
     &       + (   dij3dt(1,2) * vxyz(i,j,k,igrd,1)
     &           + dij3dt(2,2) * vxyz(i,j,k,igrd,2)
     &           + dij3dt(3,2) * vxyz(i,j,k,igrd,3) ) * pytmp
     &       + (   dij3dt(1,3) * vxyz(i,j,k,igrd,1)
     &           + dij3dt(2,3) * vxyz(i,j,k,igrd,2)
     &           + dij3dt(3,3) * vxyz(i,j,k,igrd,3) ) * pztmp     )
     &  -  press(i,j,k,igrd)  *  (
     &     vxyz(i,j,k,igrd,1) * sign( min(abs(dgi),corrmax), dgi) / di
     &   + vxyz(i,j,k,igrd,2) * sign( min(abs(dgj),corrmax), dgj) / dj
     &   + vxyz(i,j,k,igrd,3) * sign( min(abs(dgk),corrmax), dgk) / dk )
 
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c --- addition 12/09/97
      ifgr = idxfgr(igrd)

      if ( ifgr .ne. -1 ) then

         noz = norgin(ifgr,lthree)
         noy = norgin(ifgr,ltwo)
         nox = norgin(ifgr,lone)
 
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,nox,noy,noz,adderg)
         do j = noy+1, noy+ny/nfine
         do k = noz+1, noz+nz/nfine
         do i = nox+1, nox+nx/nfine
            adderg(i,j,k) = zero
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO
 
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine trav2w(sg,igrd)
 
c     transform kinematic to dynamic velocity
 
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'
 
      real*8 sg, ronoff, wmat(3,3,q), winv(3,3,q)
      real*8 w1tmp(q), w2tmp(q), w3tmp(q), factmp(q)
 

      grt = 0.8D0 * g / cc**5
      ronoff = dble(mconoff)

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+          PRIVATE(j,i3,j3,i,k,factmp,wmat,winv,w1tmp,w2tmp,w3tmp),
C$OMP+          SHARED(nx,ny,nz,grt,dij3dt,h7ij,gpot,cc,corrmax,
C$OMP+                 sg,velx,vely,velz,energy,igrd,ronoff)
      do j = 1, ny
      do k = 1, nz

         do i3 = 1, 3
         do j3 = 1, 3
            do i = 1, nx 
               wmat(i3,j3,i) = grt * dij3dt(i3,j3)
            enddo
         enddo
         enddo

         do i = 1, nx
            factmp(i) =  one  +  ronoff * sign( min( abs(
     &        two*gpot(i,j,k,igrd)/cc**2), corrmax), gpot(i,j,k,igrd) )
            wmat(1,1,i) = wmat(1,1,i) + factmp(i)
            wmat(2,2,i) = wmat(2,2,i) + factmp(i)
            wmat(3,3,i) = wmat(3,3,i) + factmp(i)
         enddo

         if ( sg .lt. zero ) then
            call mat3inv(wmat,winv,lone,nx,q)
            do i3 = 1, 3
            do j3 = 1, 3
               do i = 1, nx 
                  wmat(i3,j3,i) = winv(i3,j3,i)
               enddo
            enddo
            enddo
         endif

      do i = 1, nx
         w1tmp(i) =     wmat(1,1,i) * velx(i,j,k,igrd)
     &                + wmat(2,1,i) * vely(i,j,k,igrd)
     &                + wmat(3,1,i) * velz(i,j,k,igrd)
 
         w2tmp(i) =     wmat(1,2,i) * velx(i,j,k,igrd)
     &                + wmat(2,2,i) * vely(i,j,k,igrd)
     &                + wmat(3,2,i) * velz(i,j,k,igrd)
 
         w3tmp(i) =     wmat(1,3,i) * velx(i,j,k,igrd)
     &                + wmat(2,3,i) * vely(i,j,k,igrd)
     &                + wmat(3,3,i) * velz(i,j,k,igrd)
 
         energy(i,j,k,igrd) = energy(i,j,k,igrd) - half * (
     &   velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2 +velz(i,j,k,igrd)**2)
     &       + half * ( w1tmp(i)**2 + w2tmp(i)**2 + w3tmp(i)**2 )
 
         velx(i,j,k,igrd) = w1tmp(i)
         vely(i,j,k,igrd) = w2tmp(i)
         velz(i,j,k,igrd) = w3tmp(i)
 
      enddo

      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c---------------------------------------------------------------
      subroutine strain(igrd)
c     calculate quantities for gravitational waves
 
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'
 
      integer*4 ierr
      real*8 rxyz(3)
 
      real*8 vxyz(qx,qy,qz,ngrd,3)
      equivalence ( vxyz(1,1,1,1,1), velx(1,1,1,1) ) ,
     &            ( vxyz(1,1,1,1,2), vely(1,1,1,1) ) ,
     &            ( vxyz(1,1,1,1,3), velz(1,1,1,1) )
 
      real*8 g11, g12, g13, g21, g22, g23, g31, g32, g33

      cmgaz = dble( 0/2) + half    !  cmgax,cmgay given in compct.cmn  

      do j = 1, 3
      do i = 1, 3
         gdij2dt(i,j,igrd) = zero
      enddo
      enddo

      g11=zero
      g12=zero
      g13=zero
      g21=zero
      g22=zero
      g23=zero
      g31=zero
      g32=zero
      g33=zero

      d1x   = delx(igrd)
      d3x   = d1x**3

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i), 
C$OMP+            SHARED(nx,ny,nz,ro,densty,igrd,d3x)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         ro(i,j,k) = densty(i,j,k,igrd) * two * d3x
c     one factor two because of symmetry about z-plane -> double volume
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c     zero area covered by finer grids
      ifgr = idxfgr(igrd)
      if( ifgr .ne. -lone ) call rozero(ifgr)
 
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+       PRIVATE(j,jp,jm,k,kp,km,i,ip,im,g13,g23,g31,g32,
C$OMP+               rxyz,rd1yd,rd1zd,rd1xd,mxyz), 
C$OMP+       SHARED(igrd,nx,ny,nz,cmgax,cmgay,cmgaz,d1x,
C$OMP+              ro,velx,vely,velz,gpot,vxyz),
C$OMP+       REDUCTION(+:g11,g12,g21,g22,g33) 
      do j = 1, ny
        jp = min ( j+lone, ny)
        jm = max ( j-lone, lone )
        rxyz(2) = (dble(j)-cmgay) * d1x
        rd1yd = one / ( d1x * dble(jp-jm) )
      do k = 1, nz
c        sg = +one
c        if ( k .eq. lone )  sg = -one    ! reflection of velz at k=0.5
        kp = min ( k+lone, nz)
        km = max ( k-lone, lone )
        rxyz(3) = (dble(k)-cmgaz) * d1x
        rd1zd = one / ( d1x * dble(kp-km) )
        if ( k .eq. lone) rd1zd = one / ( two*d1x ) ! reflection of velz at k=0.5
      do i = 1, nx
        ip = min ( i+lone, nx)
        im = max ( i-lone, lone )
        rxyz(1) = (dble(i)-cmgax) * d1x
        rd1xd = one / ( d1x * dble(ip-im) )
 
c ***
        mxyz = lone    ! unrolled loop
 
c  -- derivative in x direction --
 
        g11 = g11 + ro(i,j,k) * ( 
     &     vxyz(i,j,k,igrd,mxyz) * vxyz(i,j,k,igrd,1)
     &   - rxyz(mxyz) * (gpot(ip,j,k,igrd)-gpot(im,j,k,igrd))*rd1xd )
 
c  -- derivative in y direction --
 
        g12 = g12 + ro(i,j,k) * ( 
     &     vxyz(i,j,k,igrd,mxyz) * vxyz(i,j,k,igrd,2)
     &   - rxyz(mxyz) * (gpot(i,jp,k,igrd)-gpot(i,jm,k,igrd))*rd1yd )
 
c  -- derivative in z direction --
 
        g13 = zero
cs      gdij2dt(mxyz,3,igrd) = gdij2dt(mxyz,3,igrd) + ro(i,j,k) * ( 
csymmetry->0     &     vxyz(i,j,k,igrd,mxyz) * sg * vxyz(i,j,k,igrd,3)
cs   &   - rxyz(mxyz) * (gpot(i,j,kp,igrd)-gpot(i,j,km,igrd))*rd1zd )
 
c ***
        mxyz = ltwo    ! unrolled loop
  
c  -- derivative in x direction --
 
        g21 = g21 + ro(i,j,k) * ( 
     &     vxyz(i,j,k,igrd,mxyz) * vxyz(i,j,k,igrd,1)
     &   - rxyz(mxyz) * (gpot(ip,j,k,igrd)-gpot(im,j,k,igrd))*rd1xd )
  
c  -- derivative in y direction --
 
        g22 = g22 + ro(i,j,k) * ( 
     &     vxyz(i,j,k,igrd,mxyz) * vxyz(i,j,k,igrd,2)
     &   - rxyz(mxyz) * (gpot(i,jp,k,igrd)-gpot(i,jm,k,igrd))*rd1yd )
 
c  -- derivative in z direction --
 
        g23 = zero
cs      gdij2dt(mxyz,3,igrd) = gdij2dt(mxyz,3,igrd) + ro(i,j,k) * ( 
csymmetry->0     &     vxyz(i,j,k,igrd,mxyz) * sg * vxyz(i,j,k,igrd,3)
cs   &   - rxyz(mxyz) * (gpot(i,j,kp,igrd)-gpot(i,j,km,igrd))*rd1zd )
  
c ***
        mxyz = lthree    ! unrolled loop
c       because of symmetry about z=0.5, two parts have to be added:
c       sgp*vxyz wich cancel, except for the vz*vz part.
c       rxyz*gpot -> sgp*rayz*sgp*dgpot  => sum is double as large
c       however, this factor of two is taken care of further below

c  -- derivative in x direction --
 
        g31 = zero
cs      gdij2dt(mxyz,1,igrd) = gdij2dt(mxyz,1,igrd) + ro(i,j,k) * (
csymmetry->0     &   sg * vxyz(i,j,k,igrd,mxyz) * vxyz(i,j,k,igrd,1)
cs   &   - rxyz(mxyz) * (gpot(ip,j,k,igrd)-gpot(im,j,k,igrd))*rd1xd ) 
  
c  -- derivative in y direction --
 
        g32 = zero
cs      gdij2dt(mxyz,2,igrd) = gdij2dt(mxyz,2,igrd) + ro(i,j,k) * (
csymmetry->0     &   sg * vxyz(i,j,k,igrd,mxyz) * vxyz(i,j,k,igrd,2)
cs   &   - rxyz(mxyz) * (gpot(i,jp,k,igrd)-gpot(i,jm,k,igrd))*rd1yd )
 
c  -- derivative in z direction --
 
        g33 = g33 + ro(i,j,k) * (
     &     vxyz(i,j,k,igrd,mxyz) * vxyz(i,j,k,igrd,3)
     &   - rxyz(mxyz) * (gpot(i,j,kp,igrd)-gpot(i,j,km,igrd))*rd1zd )
  
c ***
 
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
 
c      dij2dt (3,:) = zero   ! z-symmetry

      gdij2dt(1,1,igrd) = g11
      gdij2dt(1,2,igrd) = g12
      gdij2dt(1,3,igrd) = g13
      gdij2dt(2,1,igrd) = g21
      gdij2dt(2,2,igrd) = g22
      gdij2dt(2,3,igrd) = g23
      gdij2dt(3,1,igrd) = g31
      gdij2dt(3,2,igrd) = g32
      gdij2dt(3,3,igrd) = g33

      do j = 1, 3
      do i = 1, 3
         dij2dt(i,j) = zero
         do lgrd = 1, ngrd
            dij2dt(i,j) = dij2dt(i,j) + two * gdij2dt(i,j,lgrd)
         enddo 
      enddo 
      enddo 

c     factor two from equation 2STF
c     make dij2dt symmetric and trace free
 
      dijtrace = ( dij2dt(1,1) + dij2dt(2,2) + dij2dt(3,3) ) /3.D0
      dij2dt(1,1) = dij2dt(1,1) - dijtrace
      dij2dt(2,2) = dij2dt(2,2) - dijtrace
      dij2dt(3,3) = dij2dt(3,3) - dijtrace
      dij2dt(1,2) = ( dij2dt(1,2) + dij2dt(2,1) ) / two
      dij2dt(2,1) =   dij2dt(1,2)
      dij2dt(1,3) = ( dij2dt(1,3) + dij2dt(3,1) ) / two
      dij2dt(3,1) =   dij2dt(1,3)
      dij2dt(2,3) = ( dij2dt(2,3) + dij2dt(3,2) ) / two
      dij2dt(3,2) =   dij2dt(2,3)
      if ( mgonoff .eq. lzero )  then 
         do j = 1, 3
         do i = 1, 3
            dij2dt(i,j) = zero       ! test
         enddo
         enddo
      endif

      if ( mgonoff .eq. ltwo ) then
         do j = 1, 3
         do i = 1, 3
            dij2dt(i,j) = dij2dt(i,j) * two
         enddo
         enddo
      endif

      return
      end
c---------------------------------------------------------------
      subroutine mat3inv(ma,iv,be,en,qq)
c     calculate inverse matrix inv of matrix mat
      include 'qparam.cmn'

      integer*4 be,en, n,i,j,qq
      real*8  ma(3,3,qq), iv(3,3,qq)
      real*8  det(q), div(q)
      
      if ( qq .gt. q ) stop ' mat3inv'

      do n = be, en

      det(n)=ma(1,1,n)*ma(2,2,n)*ma(3,3,n)+ma(1,2,n)*ma(2,3,n)*ma(3,1,n)
     &      +ma(1,3,n)*ma(3,2,n)*ma(2,1,n)-ma(1,3,n)*ma(2,2,n)*ma(3,1,n)
     &      -ma(2,3,n)*ma(3,2,n)*ma(1,1,n)-ma(2,1,n)*ma(1,2,n)*ma(3,3,n)

      if ( det(n) .eq. zero ) then
         write(*,*) 'matinv:  det .eq. 0 ',n,be,en
         write(*,'(1P,3E20.10)') ((ma(i,j,n),i=1,3),j=1,3)
         stop
      endif

      div(n) = one / det(n)

      iv(1,1,n) = ( ma(2,2,n)*ma(3,3,n) - ma(3,2,n)*ma(2,3,n) ) * div(n) 
      iv(2,1,n) = ( ma(2,3,n)*ma(3,1,n) - ma(2,1,n)*ma(3,3,n) ) * div(n) 
      iv(3,1,n) = ( ma(2,1,n)*ma(3,2,n) - ma(2,2,n)*ma(3,1,n) ) * div(n) 
      iv(1,2,n) = ( ma(3,2,n)*ma(1,3,n) - ma(1,2,n)*ma(3,3,n) ) * div(n) 
      iv(2,2,n) = ( ma(1,1,n)*ma(3,3,n) - ma(3,1,n)*ma(1,3,n) ) * div(n) 
      iv(3,2,n) = ( ma(1,2,n)*ma(3,1,n) - ma(1,1,n)*ma(3,2,n) ) * div(n) 
      iv(1,3,n) = ( ma(1,2,n)*ma(2,3,n) - ma(2,2,n)*ma(1,3,n) ) * div(n) 
      iv(2,3,n) = ( ma(2,1,n)*ma(1,3,n) - ma(1,1,n)*ma(2,3,n) ) * div(n) 
      iv(3,3,n) = ( ma(1,1,n)*ma(2,2,n) - ma(2,1,n)*ma(1,2,n) ) * div(n) 

      enddo

      return
      end
