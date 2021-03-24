c-----------------------------------------------------------------------
      subroutine video(istereo)
*--------------------------------------------------------------*
*     Voxel Intensity Density Emitter Encoded Output.          *
*     Variable  Density Emitter Encoded Output.                *
*     Plot subroutine to visualize three-dim. voxel data.      *
*     Based on density dependent emission and absorption.      *
*     Color attributed according to temperature.               *
*     Vector, Film, Multigrid, Adaptive Undersampling Version  *
*     Maximilian Ruffert,       (c)          28.Oct.1991       *
*--------------------------------------------------------------*

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'
      include 'rays.cmn'

      integer istereo

      real*8    denvo(vx,vy,vz), tempr(vx,vy,vz),
     &          rpix(0:npix,0:npix), gpix(0:npix,0:npix),
     &          bpix(0:npix,0:npix), temvo(vx)
      integer*4 ntemp(vx,vy,vz), intpix(0:npix,0:npix)
      integer*1 matrgb(3,0:npix,0:npix)   ! this should be unsigned integer
      common /coclock/ matrgb

      common /raypot/ denvo, tempr, rpix, gpix, bpix, ntemp, intpix

      integer*4   itmpnt, itmt, nrgb

      integer*4 i,j,k,n, kx,ky,kz
      integer*4 iptx, ipty

      real*8 ttt,tt(2)
      integer*4  lg
      integer*4  ntime, ltfs, mtr, mtg, mtb, mor, mog, mob
      real*8 ddm, ttm, ttmm, ddmm, ntmm
      character*4 ifra, wpix
      character*4 vpix
      character*13 filnamrgb
      character*14 filnamjpeg
      character*1  VW
      character*2  crest

C --- INITIALIZE PARAMETERS

      parameter (ddm = one, ttm = 1.D0)
      nttm = int(ttm*dble(len))

      rvx  = dble(vx)
      fdz  = dble(vz)/rvx/two
      rvzx = dble(vz)/rvx
      rvzxm = rvzx - 1.D-10
      rnpix = dble(npix)

      edis = 0.10D0 ! allowable variation between pixels to interpolate
      intval = 4
      opac = 0.3D0 * 64.D0 / rvx

      iste = istereo
      iframe = iframe + lone
      if ( iframe .gt. nview ) stop'iframe gt nview'
      tframe = time
      VW = 'V'
      opx = half
      opz = -one
      fby = 0.3D0
      fbz = one
      if ( istereo .eq. -1 ) opx = opx - 0.08D0
      if ( istereo .eq. +1 ) then
         iframe = iframe - lone
         VW = 'W'
         opx = opx + 0.08D0
      endif
c     only for +1, not for -1 !
      if ( istereo .ne. 0 ) then
         opz = -2.5D0
         fby = half
         fbz = zero
      endif

c - minima maxima of ray intensity

      frmax = one
      frmin = zero
c      frd = log10( frmax / frmin )
      fgmax = one
      fgmin = zero
c      fgd = log10( fgmax / fgmin )
      fbmax = one
      fbmin = zero
c      fbd = log10( fbmax / fbmin )

C --- INTERPOLATE NESTED GRIDS AROUND CENTER OF MASS

      rnf = 3.0D0

      vicenx = cenx
      viceny = ceny
      vicenz = cenz

      call gridnest(vicenx,viceny,vicenz,rnf)

* --- NORMALIZE DATA AND SHARPEN/DAMPEN

      ddmm = zero

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,temvo),
C$OMP+            SHARED(tpmin,tpmax,tempr,ntemp,denvo,dymin),
C$OMP+            REDUCTION(max:ddmm)
      do j = 1, vy
      do k = 1, vz
      do i = 1, vx
         temvo(i) = (tempr (i,j,k)-tpmin) / (tpmax-tpmin)
         temvo(i) = min( max( temvo(i), zero ), ttm )
c         temvo(i) = temvo(i) ** 2
         ntemp(i,j,k) = nint( 1000.D0 * temvo(i) )
c         denvo(i,j,k) = (denvo(i,j,k)-dymin) / (dymax-dymin)                   !Was commented out by MRR
c         denvo(i,j,k) = min( max( denvo(i,j,k), zero ), ddm )                  !!Including them results in `no image`
c         denvo(i,j,k) = denvo(i,j,k) ** 2                                      !!Unless parts below are commented out
         if (denvo(i,j,k) .lt. dymin) denvo(i,j,k) = dymin+half
         denvo(i,j,k) = denvo(i,j,k) - dble(nint(denvo(i,j,k)))
         denvo(i,j,k) = exp( -one * min( (denvo(i,j,k)/0.1D0)**4,50.D0))
         ddmm = max( ddmm, denvo(i,j,k) )
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO


* --- ADD EDGES FOR VISIBILITY

c      if ( iframe .lt. iendf+100 ) then
         call edgevisr ( denvo, 14.0D0 )  ! now any integer
         call edgevisi ( ntemp, lzero+1000 )
         call edgevisr ( tempr, tpmin+(tpmax-tpmin)*0.2D0)
c      endif
c      call cornvisr ( denvo, 1D3*ddmm )
c      call cornvisi ( ntemp, lzero+1000  )
c      call cornvisr ( tempr, tpmin+(tpmax-tpmin)*0.2D0)

* --- INITIALIZE VARIABLES

c     ttt = second(tt)
      scal = one / scale(iframe)
      if ( istereo .ne. 0 ) scal = 1.4D0
c     same 5000 as in inivid
c      if (mod(iframe,10).eq.0) write(*,*) 'video  ',iframe,scale(iframe)

      th = theta(iframe)/180.D0*pi
      cost = cos(th)
      sint = sin(th)
      ph = phi(iframe)/180.D0*pi
      cosp = cos(ph)
      sinp = sin(ph)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(iy,ix),
C$OMP+            SHARED(intpix,rpix,gpix,bpix,matrgb)
      do iy = 0, npix
      do ix = 0, npix
         intpix(ix,iy) =  lzero
         rpix  (ix,iy) =  zero
         gpix  (ix,iy) =  zero
         bpix  (ix,iy) =  zero
         matrgb(1,ix,iy) = 0
         matrgb(2,ix,iy) = 0
         matrgb(3,ix,iy) = 0
      enddo
      enddo
C$OMP END PARALLEL DO

* --- CALCULATE PLOT EDGES

      ixmax = lzero
      iymax = lzero
      ixmin = npix
      iymin = npix

      do kz = 0, 1
        dmz = (dble(kz)*rvzx-fdz) / scal
      do ky = 0, 1
        dmy = (dble(ky)-fdy) / scal
      do kx = 0, 1
        dmx = (dble(kx)-fdx) / scal

* --- ... FIND PIC POSITION FROM DATA POSITION

c            shift  and    turn   data cube
c              :            /:\
        blx = fbx +  cosp*     dmx + sinp*     dmy
        bly = fby   -sinp*sint*dmx + cosp*sint*dmy + cost*dmz
        blz = fbz   -sinp*cost*dmx + cosp*cost*dmy - sint*dmz

* --- ... FIND BEAM POSITION FROM PIC POSITION, AND FIND MAX,MIN

        plp  = one/ (one + blz/(-opz))
        iptx = nint ( ( (blx-opx)*plp + opx ) * rnpix )
        ipty = nint ( ( (bly-opy)*plp + opy ) * rnpix )

        ixmax = max( ixmax, iptx+1)
        ixmin = min( ixmin, iptx-1)

        iymax = max( iymax, ipty+1)
        iymin = min( iymin, ipty-1)

      enddo
      enddo
      enddo

      ixmin = max( ixmin, lzero )
      ixmax = min( ixmax, npix  )
      iymin = max( iymin, lzero )
      iymax = min( iymax, npix  )

      ltrace = lzero

      ixmin = lzero
      ixmax = (npix/intval)*intval
      iymin = lzero
      iymax = ixmax

c --- FIRST PASS: trace all rays in intervals of intval pixels

      do iy = iymin, iymax, intval*2
         do ix = ixmin, ixmax, intval
            ltrace = ltrace + lone
            xpix(ltrace) = dble(ix)
            ypix(ltrace) = dble(iy)
            ltrace = ltrace + lone
            xpix(ltrace) = dble(ix)
            ypix(ltrace) = dble(iy+intval)
         enddo
         !if ( ltrace .gt. ntra/ltwo ) call trace(lzero,lzero)          !DRW - Commented out
         if ( ltrace .gt. 10   ) call trace(lzero,lzero)                !DRW - Added to test
      enddo
      if ( ltrace .gt. lzero ) call trace(lzero,lzero)

c --- SECOND PASS: interpolate or trace diagonally every intval pixels
 222  continue

      intvl2 = intval/2

      do iy = iymin+intvl2, iymax, intval*2    ! shift intvl2 for
       call pixpol(ixmin+intvl2,ixmax,intval,intvl2,iy,       edis,ltwo)
       call pixpol(ixmin+intvl2,ixmax,intval,intvl2,iy+intval,edis,ltwo)
       !if ( ltrace .gt. ntra/ltwo ) call trace(lzero,lzero)            !DRW - Commented out
       if ( ltrace .gt. 10   ) call trace(lzero,lzero)                  !DRW - Added to test
      enddo
      if ( ltrace .gt. lzero ) call trace(lzero,lzero)

c --- THIRD PASS: interpolate or trace horiz/vert every intval pixels

c     in = lzero
      do iy = iymin,    iymax, intvl2*2   ! go to every line
c       in = in + lone
c       id = intvl2 * mod(in,ltwo)
        call pixpol(ixmin+intvl2,ixmax,intval,intvl2,iy,edis,lthree)
        call pixpol(ixmin,ixmax,intval,intvl2,iy+intvl2,edis,lthree)
        !if ( ltrace .gt. ntra/ltwo ) call trace(lzero,lzero)           !DRW - Commented out
        if ( ltrace .gt. 10   ) call trace(lzero,lzero)                 !DRW - Added to test
      enddo
      if ( ltrace .gt. lzero ) call trace(lzero,lzero)

c --- FOURTH PASS: check interpolated pixels and retrace

      intval = intvl2

 120  retrace = .false.
c               !  repeat retraceing until no more new pixels are traced
      do iy = iymin, iymax, intval*2
        call pixpol(ixmin,ixmax,intval,intvl2,iy,       edis,ltwo+ltwo)
        call pixpol(ixmin,ixmax,intval,intvl2,iy+intval,edis,ltwo+ltwo)
        !if ( ltrace .gt. ntra/ltwo ) call trace(lzero,lzero)           !DRW - Commented out
        if ( ltrace .gt. 10   ) call trace(lzero,lzero)                 !DRW - Added to test
      enddo
      if ( ltrace .gt. lzero ) call trace(lzero,lzero) ! clear up

      if ( retrace ) goto 120

c --- repeat PASSES TWO-FOUR until
      if ( intval .ne. lone ) goto 222

c --- LAST PASS: oversample where necessary                             !Was commented out by MRR
                                                                        !Including it seems to do nothing
                                                                        !to the initial image
c      intval = lone
c      edis = 0.60D0               ! only at very sharp changes
c
c      do 130 iy = iymin, iymax-1, intval                               !I don't know this logic, DRW
c         call pixpol(ixmin,ixmax-1,intval,intvl2,iy, edis,lthree+ltwo)
c         if ( ltrace .gt. ntra/3 ) call trace(lone,lzero)
c  130 continue
c      call trace(lone,lzero)             ! clear up
c
c ----------------------                                                !All of this commented out by MRR

c make random stellar background
c from random the,phi calculate positions on screen

      goto 115
      ltrace = lzero

      aist = abs(dble(istereo))

      iphi = lzero
      phish = zero

 116  do i = 1, lenstar
c positions of stars is not correct yet, just approximate
         xx = rnpix * ( opx - abs(opz) *
     &         tan(phistar(i)+phish-ph) )
         if ( istereo .eq. +1 ) xx = xx - shistar(i)
         yy = rnpix * ( opy - abs(opz) *
     &         tan(thestar(i)-th) )
         if ( (xx.gt.zero) .and. (xx.lt.rnpix-one) .and.
     &        (yy.gt.zero) .and. (yy.lt.rnpix-one)       ) then
            if (ltrace .ge. ntra-3) goto 114
c            if ( istereo .eq. 0 ) then
c               ltrace = ltrace + lone
c               xpix(ltrace) = xx
c               ypix(ltrace) = yy
c               ltrace = ltrace + lone
c               xpix(ltrace) = xx
c               ypix(ltrace) = yy+one
c               ltrace = ltrace + lone
c               xpix(ltrace) = xx+one
c               ypix(ltrace) = yy
c               ltrace = ltrace + lone
c               xpix(ltrace) = xx+one
c               ypix(ltrace) = yy+one
c            else
               ltrace = ltrace + lone
               xpix(ltrace) = xx
               ypix(ltrace) = yy
               if ( shistar(i) .gt. 57.5D0 ) then
                  ltrace = ltrace + lone
                  xpix(ltrace) = xx
                  ypix(ltrace) = yy+one
               endif
               if ( shistar(i) .gt. 55D0 ) then
                  ltrace = ltrace + lone
                  xpix(ltrace) = xx+one
                  ypix(ltrace) = yy
               endif
               if ( shistar(i) .gt. 52.5D0 ) then
                  ltrace = ltrace + lone
                  xpix(ltrace) = xx+one
                  ypix(ltrace) = yy+one
               endif
c            endif
         endif
      enddo

      if ( ltrace .lt. ntra-5 ) then
         phish = phish + 0.3D0
         iphi = iphi + lone
         if ( iphi .lt. 10 ) goto 116
      endif

  114 if ( ltrace .gt. lzero ) call trace(lzero,lone)

  115 continue
c ----------------------

c      rmax = one
c      gmax = one
c      bmax = one

C$OMP PARALLEL DEFAULT(NONE), PRIVATE(iy,ix),
C$OMP+         SHARED(rpix,gpix,bpix, matrgb,
C$OMP+                frmax,frmin,frd,fgmax,fgmin,fgd,fbmax,fbmin,fbd)

C$OMP DO
      do iy = 0, npix
      do ix = 0, npix
         rpix(ix,iy) = max(min(rpix(ix,iy),frmax),frmin)
         gpix(ix,iy) = max(min(gpix(ix,iy),fgmax),fgmin)
         bpix(ix,iy) = max(min(bpix(ix,iy),fbmax),fbmin)
      enddo
      enddo
C$OMP END DO

C$OMP DO
      do iy = 0, npix
      do ix = 0, npix
 	 matrgb (1,ix,iy) = max(min( nint( 255D0 *rpix(ix,iy) ), 255), 0)
 	 matrgb (2,ix,iy) = max(min( nint( 255D0 *gpix(ix,iy) ), 255), 0)
 	 matrgb (3,ix,iy) = max(min( nint( 255D0 *bpix(ix,iy) ), 255), 0)
      enddo
      enddo
C$OMP END DO

C$OMP END PARALLEL

c      ttt = second(tt)-ttt
c      ttt = etime(tt) - ttt
c      write(*,*) 'Time to calculate plot: ', ttt

      timeunit = 1.D-3
      ntime = nint(time*1000.D0/timeunit)
c      write(3) ntime
      call putclock(time/timeunit)

cC$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(iy,ix), SHARED(matrgb)
c      do iy = 0, npix
c      do ix = 0, npix
c 	 matrgb (1,ix,iy) = matrgb(1,ix,iy) * 256 ! *256 only for GRAND/SGI?
c 	 matrgb (2,ix,iy) = matrgb(2,ix,iy) * 256 ! *256 only for GRAND/SGI?
c 	 matrgb (3,ix,iy) = matrgb(3,ix,iy) * 256 ! *256 only for GRAND/SGI?
c      enddo
c      enddo
cC$OMP END PARALLEL DO

      ifram = 1000+iframe
      write(ifra,'(1I4)') ifram
      write(vpix,'(1I4)') npix+lone

c --- ppm format
c        filnam = '/home/mruffert/data/' // basenm // VW // ifra // '.ppm'
c        open(13, file=filnam, form='unformatted', access='direct',
c     &         recl=16+3*(npix+1)*(npix+1) )
c        write(wpix,'(1I4)') npix+1
c        write(13,rec=1) "P6",char(14),wpix,wpix," 255",char(14),matrgb

c --- rgb format
        filnamrgb  = basenm//VW//ifra// '.rgb'
        filnamjpeg = basenm//VW//ifra// '.jpeg'
        open(13, file=filnamrgb, form='unformatted', access='direct',
c     &       recl= 6*(npix+1)*(npix+1) ) ! only GRAND/SGI
     &       recl= 3*(npix+1)*(npix+1) )   ! for linux PC (gfortran)
c     &       recl=(3*(npix+1)*(npix+1))/4 )    ! for linux PC, ifort
        write(13,rec=1) matrgb ! , ntime
        close(13)

c      call system("/usr/local/bin/gzip -f "// filnamrgb)
      call system('/usr/bin/convert -depth 8 -size ' //
     &             vpix // 'x' // vpix // ' ' // filnamrgb // ' '
     &       // filnamjpeg // ' && /bin/rm -f ' // filnamrgb )
c     &       // filnamjpeg // ' && /usr/bin/rm -f ' // filnamrgb )
c      call system('convert -size ' // vpix//'x'//vpix // ' ' //
c     &        filnamrgb // '"[0]" ' // filnamjpeg //
c     &        ' && /usr/bin/rm -f ' // filnamrgb )

      mifra = mod(ifram,100)
      if ( mifra .eq. 99 ) then
         write(crest,'(1I2)') ifram/100
         call system('tar -cf '// basenm//'V' // crest // '.jpeg.tar '
     &                // basenm//'V' // crest // '??.jpeg' )
      endif

c      stop 'rays'
      return
      end
c-----------------------------------------------------------------------
      subroutine pixpol(mi,ma,md,m2,iy,edis,npass)

      include 'qparam.cmn'
      include 'rays.cmn'

      real*8    denvo(vx,vy,vz), tempr(vx,vy,vz),
     &          rgbpix(0:npix,0:npix,3)
      integer*4 ntemp(vx,vy,vz), intpix(0:npix,0:npix)

      common /raypot/ denvo, tempr, rgbpix, ntemp, intpix

      real*8 rgbint(0:npix,3), rgbdis(0:npix,3)
      equivalence ( rint, rgbint )


c-----
      if     ( npass .eq. ltwo ) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(ix,im,ip,v1,v2,v3,v4,irgb),
C$OMP+            SHARED(jm,jp,mi,ma,md,m2,iy,
C$OMP+                   ixmin,ixmax,iymin,iymax,rgbpix,rgbint,rgbdis)
      do ix = mi, ma, md

         im = ix - m2
         ip = ix + m2
         jm = iy - m2
         jp = iy + m2

         if ( im .lt. ixmin ) im = ip
         if ( ip .gt. ixmax ) ip = im
         if ( jm .lt. iymin ) jm = jp
         if ( jp .gt. iymax ) jp = jm

         do irgb = 1, 3

            v1 = rgbpix(im,jm,irgb)
            v2 = rgbpix(ip,jm,irgb)
            v3 = rgbpix(im,jp,irgb)
            v4 = rgbpix(ip,jp,irgb)

c         interpolate values
            rgbint(ix,irgb) = 0.25D0 * ( ( v1 + v2 ) + ( v3 + v4 ) )

c         distance between interpolated and given values
            rgbdis(ix,irgb) =
     &        max( abs(v1-rgbint(ix,irgb)), abs(v2-rgbint(ix,irgb)),
     &             abs(v3-rgbint(ix,irgb)), abs(v4-rgbint(ix,irgb))  )

         enddo
      enddo
C$OMP END PARALLEL DO

      do ix = mi, ma, md
         if (      (rgbdis(ix,1) .gt. edis*rgbint(ix,1))
     &        .or. (rgbdis(ix,2) .gt. edis*rgbint(ix,2))
     &        .or. (rgbdis(ix,3) .gt. edis*rgbint(ix,3))  ) then
            ltrace = ltrace + lone
            xpix(ltrace) = dble(ix)
            ypix(ltrace) = dble(iy)
         else
            intpix(ix,iy) = lone !  mark pixels that are interpolated
            rgbpix(ix,iy,1) = rgbint(ix,1)
            rgbpix(ix,iy,2) = rgbint(ix,2)
            rgbpix(ix,iy,3) = rgbint(ix,3)
         endif
      enddo

c-----
      elseif ( npass .eq. lthree ) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(ix,im,ip,v1,v2,v3,v4,irgb),
C$OMP+            SHARED(jm,jp,mi,ma,md,m2,iy,
C$OMP+                   ixmax,ixmin,iymax,iymin,rgbpix,rgbint,rgbdis)
      do ix = mi, ma, md

         im = ix - m2
         ip = ix + m2
         jm = iy - m2
         jp = iy + m2

         if ( im .lt. ixmin ) im = ip
         if ( ip .gt. ixmax ) ip = im
         if ( jm .lt. iymin ) jm = jp
         if ( jp .gt. iymax ) jp = jm

         do irgb = 1, 3

            v1 = rgbpix(im,iy,irgb)
            v2 = rgbpix(ip,iy,irgb)
            v3 = rgbpix(ix,jm,irgb)
            v4 = rgbpix(ix,jp,irgb)

c         interpolate values
            rgbint(ix,irgb) = 0.25D0 * ( ( v1 + v2 ) + ( v3 + v4 ) )

c         distance between interpolated and given values
            rgbdis(ix,irgb) =
     &        max( abs(v1-rgbint(ix,irgb)), abs(v2-rgbint(ix,irgb)),
     &             abs(v3-rgbint(ix,irgb)), abs(v4-rgbint(ix,irgb))  )

         enddo
      enddo
C$OMP END PARALLEL DO

      do ix = mi, ma, md
         if (      (rgbdis(ix,1) .gt. edis*rgbint(ix,1))
     &        .or. (rgbdis(ix,2) .gt. edis*rgbint(ix,2))
     &        .or. (rgbdis(ix,3) .gt. edis*rgbint(ix,3))  ) then
            ltrace = ltrace + lone
            xpix(ltrace) = dble(ix)
            ypix(ltrace) = dble(iy)
         else
            intpix(ix,iy) = lone !  mark pixels that are interpolated
            rgbpix(ix,iy,1) = rgbint(ix,1)
            rgbpix(ix,iy,2) = rgbint(ix,2)
            rgbpix(ix,iy,3) = rgbint(ix,3)
         endif
      enddo

c-----
      elseif ( npass .eq. 4 ) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(ix,im,ip,v1,v2,v3,v4,irgb),
C$OMP+            SHARED(mi,ma,md,m2,iy,jm,jp,
C$OMP+                   ixmax,ixmin,iymax,iymin,rgbpix,rgbint,rgbdis)
      do ix = mi, ma, md

         im = ix - m2
         ip = ix + m2
         jm = iy - m2
         jp = iy + m2

         if ( im .lt. ixmin ) im = ip
         if ( ip .gt. ixmax ) ip = im
         if ( jm .lt. iymin ) jm = jp
         if ( jp .gt. iymax ) jp = jm

         do irgb = 1, 3

            v1 = rgbpix(im,iy,irgb)
            v2 = rgbpix(ip,iy,irgb)
            v3 = rgbpix(ix,jm,irgb)
            v4 = rgbpix(ix,jp,irgb)

c         interpolate values
            rgbint(ix,irgb) = rgbpix(ix,iy,irgb)

c         distance between interpolated and given values
            rgbdis(ix,irgb) =
     &        max( abs(v1-rgbint(ix,irgb)), abs(v2-rgbint(ix,irgb)),
     &             abs(v3-rgbint(ix,irgb)), abs(v4-rgbint(ix,irgb))  )

         enddo
      enddo
C$OMP END PARALLEL DO

      do ix = mi, ma, md
         if (   (     (rgbdis(ix,1) .gt. edis*rgbint(ix,1))
     &           .or. (rgbdis(ix,2) .gt. edis*rgbint(ix,2))
     &           .or. (rgbdis(ix,3) .gt. edis*rgbint(ix,3)))
     &        .and. intpix(ix,iy).eq.1 ) then
c            at fourth pass check only interpolated pixels
            intpix(ix,iy) = lzero !  mark pixels to be traced
            retrace = .true.
            ltrace = ltrace + lone
            xpix(ltrace) = dble(ix)
            ypix(ltrace) = dble(iy)
         endif
      enddo

c-----
      elseif ( npass .eq. 5 ) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(ix,ip,v1,v2,v3,v4,irgb),
C$OMP+            SHARED(mi,ma,iy,jp,rgbint,rgbpix,rgbdis)
      do ix = mi, ma

         ip = ix + lone
         jp = iy + lone

         do irgb = 1, 3
c         interpolate values
            rgbint(ix,irgb) = ( rgbpix(ix,jp,irgb)
     &         + rgbpix(ip,iy,irgb) + rgbpix(ip,jp,irgb) ) / 3.D0
c       distance between pixels
            rgbdis(ix,irgb) =
     &        max( abs(rgbpix(ix,jp,irgb)-rgbpix(ix,iy,irgb)),
     &             abs(rgbpix(ip,iy,irgb)-rgbpix(ix,iy,irgb)),
     &             abs(rgbpix(ip,jp,irgb)-rgbpix(ix,iy,irgb))  )
         enddo
      enddo
C$OMP END PARALLEL DO

      do ix = mi, ma
c         if values differ too much then oversample by trace
         if (      (rgbdis(ix,1) .gt. edis*rgbint(ix,1))
     &        .or. (rgbdis(ix,2) .gt. edis*rgbint(ix,2))
     &        .or. (rgbdis(ix,3) .gt. edis*rgbint(ix,3)) ) then

            xpix(ltrace+1) = dble(ix)          !
            ypix(ltrace+1) = dble(iy)+0.45D0   !  0.45 and not 0.5
            xpix(ltrace+2) = dble(ix)+0.45D0   !  so that 'nint'
            ypix(ltrace+2) = dble(iy)          !  in trace will
            xpix(ltrace+3) = dble(ix)+0.45D0   !  choose the lower
            ypix(ltrace+3) = dble(iy)+0.45D0   !  integer value
            ltrace = ltrace + lthree
         endif
      enddo
c-----
      endif
c-----

      return
      end
c-----------------------------------------------------------------------
      subroutine trace(nover,lstar)
c     traces rays through data cube
c     input: arrays xpix, ypix, ltrace; positions of rays to be traced
c      as pixel number in x and y direction, number of rays to be traced
c     output: rsum, gsum, bsum

      include 'qparam.cmn'
      include 'rays.cmn'

      real*8    denvo(vx,vy,vz), tempr(vx,vy,vz),
     &          rpix(0:npix,0:npix), gpix(0:npix,0:npix),
     &          bpix(0:npix,0:npix)
      integer*4 ntemp(vx,vy,vz), intpix(0:npix,0:npix)

      common /raypot/ denvo, tempr, rpix, gpix, bpix, ntemp, intpix

      real*8           orix(ntra),  oriy(ntra),  oriz(ntra),
     &                 dirx(ntra),  diry(ntra),  dirz(ntra),
     &                 t1xy(ntra),  t2xy(ntra),  t1xz(ntra),
     &                 t2xz(ntra),  t1yz(ntra),  t2yz(ntra),
     &                  ax1(ntra),   ay1(ntra),   az1(ntra),
     &                  ax2(ntra),   ay2(ntra),   az2(ntra),
     &                  dex(ntra),   dey(ntra),   dez(ntra),
     &                tdepx(ntra), tdepy(ntra), tdepz(ntra),
     &                 rkdx(ntra),  rkdy(ntra),  rkdz(ntra),
     &                  pax(ntra),   pay(ntra),   paz(ntra),
     &                  rho(ntra),  semi(ntra),  emit(ntra),
     &                ridno(ntra),  xidn(ntra),  yidn(ntra),
     &                  tmi(ntra),   tma(ntra),   fac(ntra),
     &                   t0(ntra),    t1(ntra),  tpsi(ntra),
     &                   xt(ntra),    yt(ntra),    zt(ntra),
     &                 temp(ntra),  ekap(ntra)

      integer*4        sgnx(ntra),  sgny(ntra),  sgnz(ntra),
     &                  ipx(ntra),   ipy(ntra),   ipz(ntra),
     &                 ntau(ntra)

      logical        inside(ntra)

      real*8 twth, roex
      parameter ( twth = two/3.D0, roex = twth )

      integer*4 itma , itmi ,  left

      if (ltrace.gt.ntra) write(*,*) 'trace stop, ltrace.gt.ntra'
      if (ltrace .eq. lzero) return

* --- FOR EVERY RAY DO ...

* --- FIND ORIGINS AND DIRECTIONS OF LIGHT RAYS

C$OMP PARALLEL DEFAULT(NONE), PRIVATE(it),
C$OMP+         SHARED(ltrace,rnpix,scal,sinp,cosp,sint,cost,fdz,
C$OMP+                ypix,yidn,sidny,xpix,xidn,ridno,fac,xt,yt,zt,
C$OMP+                dirx,diry,dirz,orix,oriy,oriz,tmi,tma,inside,
C$OMP+                rsum,gsum,bsum,rkdz,rkdy,rkdx, rvzx,rvzxm,rvx,
C$OMP+                t1xy,t1xz,t1yz,t2xy,t2xz,t2yz,pax,pay,paz,opz,
C$OMP+                ax1,ay1,az1,ax2,ay2,az2,tdepx,tdepy,tdepz,opx,
C$OMP+                ipx,ipy,ipz,dex,dey,dez,sgnx,sgny,sgnz,t0,
C$OMP+                fby,fbz,iste)

C$OMP DO
      do it = 1, ltrace
         yidn(it) = ypix(it) / rnpix
         sidny = opz**2 + (yidn(it)-opy)**2

         xidn(it) = xpix(it) / rnpix

         ridno(it) = one / sqrt ( sidny + (xidn(it)-opx)**2 )
         fac (it) = scal * ridno(it)

         xt(it) = fac(it) * (xidn(it) - opx)
         yt(it) = fac(it) * (yidn(it) - opy)
         zt(it) = fac(it) * (         - opz)

         dirx(it) = cosp*xt(it) -sinp*sint*yt(it) -sinp*cost*zt(it)
         diry(it) = sinp*xt(it) +cosp*sint*yt(it) +cosp*cost*zt(it)
         dirz(it) =                   cost*yt(it) -     sint*zt(it)

         xt(it) =  scal * (xidn(it) - fbx)
         yt(it) =  scal * (yidn(it) - fby)
         zt(it) =  scal * (         - fbz)

c       origins of rays on screen
         orix(it) = fdx+cosp*xt(it)-sinp*sint*yt(it)-sinp*cost*zt(it)
         oriy(it) = fdy+sinp*xt(it)+cosp*sint*yt(it)+cosp*cost*zt(it)
         oriz(it) = fdz+                 cost*yt(it)-     sint*zt(it)

         tmi(it) = 1000.0D0
         tma(it) =    zero
         inside(it) = .false.

         fac(it)  = one         ! this is for later !
         rsum(it) = zero
         gsum(it) = zero
         bsum(it) = zero
      enddo
C$OMP END DO


* --- FIND BOUNDS OF RAY DEPTH TMI

C$OMP DO
      do it = 1, ltrace
         rkdz(it) = one/ (dirz(it)+1.D-10)
         t1xy(it) = -oriz(it) * rkdz(it)
         ax1(it) = dirx(it) * t1xy(it) + orix(it)
         ay1(it) = diry(it) * t1xy(it) + oriy(it)
         t2xy(it) = rvzx*rkdz(it) + t1xy(it)
         ax2(it) = dirx(it) * t2xy(it) + orix(it)
         ay2(it) = diry(it) * t2xy(it) + oriy(it)
         if (  (ax1(it).gt.zero) .and. (ax1(it).lt.one) .and.
     &         (ay1(it).gt.zero) .and. (ay1(it).lt.one)      )  then
            inside(it) = .true.
            tma(it) = max( tma(it), t1xy(it) )
            if ( t1xy(it) .lt. tmi(it) ) then
               tmi(it) = t1xy(it)
               pax(it) = ax1(it)
               pay(it) = ay1(it)
               paz(it) = zero
            endif
         endif
         if (  (ax2(it).gt.zero) .and. (ax2(it).lt.one) .and.
     &         (ay2(it).gt.zero) .and. (ay2(it).lt.one)      )  then
            inside(it) = .true.
            tma(it) = max( tma(it), t2xy(it) )
            if ( t2xy(it) .lt. tmi(it) ) then
               tmi(it) = t2xy(it)
               pax(it) = ax2(it)
               pay(it) = ay2(it)
               paz(it) = rvzxm
            endif
         endif
      enddo
C$OMP END DO

C$OMP DO
      do it = 1, ltrace
         rkdy(it) = one/ (diry(it)+1.D-10)
         t1xz(it) = -oriy(it) * rkdy(it)
         ax1(it) = dirx(it) * t1xz(it) + orix(it)
         az1(it) = dirz(it) * t1xz(it) + oriz(it)
         t2xz(it) = rkdy(it) + t1xz(it)
         ax2(it) = dirx(it) * t2xz(it) + orix(it)
         az2(it) = dirz(it) * t2xz(it) + oriz(it)
         if (  (ax1(it).gt.zero) .and. (ax1(it).lt.one) .and.
     &         (az1(it).gt.zero) .and. (az1(it).lt.rvzxm)    )  then
            inside(it) = .true.
            tma(it) = max( tma(it), t1xz(it) )
            if ( t1xz(it) .lt. tmi(it) ) then
               tmi(it) = t1xz(it)
               pax(it) = ax1(it)
               pay(it) = zero
               paz(it) = az1(it)
            endif
         endif
         if (  (ax2(it).gt.zero) .and. (ax2(it).lt.one) .and.
     &         (az2(it).gt.zero) .and. (az2(it).lt.rvzxm)    )  then
            inside(it) = .true.
            tma(it) = max( tma(it), t2xz(it) )
            if ( t2xz(it) .lt. tmi(it) ) then
               tmi(it) = t2xz(it)
               pax(it) = ax2(it)
               pay(it) = 0.99999999D0
               paz(it) = az2(it)
            endif
         endif
      enddo
C$OMP END DO

C$OMP DO
      do it = 1, ltrace
         rkdx(it) = one/ (dirx(it)+1.D-10)
         t1yz(it) = -orix(it) * rkdx(it)
         ay1(it) = diry(it) * t1yz(it) + oriy(it)
         az1(it) = dirz(it) * t1yz(it) + oriz(it)
         t2yz(it) = rkdx(it) + t1yz(it)
         ay2(it) = diry(it) * t2yz(it) + oriy(it)
         az2(it) = dirz(it) * t2yz(it) + oriz(it)
         if (  (ay1(it).gt.zero) .and. (ay1(it).lt.one) .and.
     &         (az1(it).gt.zero) .and. (az1(it).lt.rvzxm)    )  then
            inside(it) = .true.
            tma(it) = max( tma(it), t1yz(it) )
            if ( t1yz(it) .lt. tmi(it) ) then
               tmi(it) = t1yz(it)
               pax(it) = zero
               pay(it) = ay1(it)
               paz(it) = az1(it)
            endif
         endif
         if (  (ay2(it).gt.zero) .and. (ay2(it).lt.one) .and.
     &         (az2(it).gt.zero) .and. (az2(it).lt.rvzxm)    )  then
            inside(it) = .true.
            tma(it) = max( tma(it), t2yz(it) )
            if ( t2yz(it) .lt. tmi(it) ) then
               tmi(it) = t2yz(it)
               pax(it) = 0.99999999D0
               pay(it) = ay2(it)
               paz(it) = az2(it)
            endif
         endif
      enddo
C$OMP END DO

      if ( iste .eq. 0 ) then
c     check if part of screen is in data cube:
C$OMP DO
         do it = 1, ltrace
            if ( tmi(it) .lt. zero ) then
            if ( tma(it) .gt. zero ) then
               tmi(it) = zero
               pax(it) = orix(it)
               pay(it) = oriy(it)
               paz(it) = oriz(it)
            else
               inside(it) = .false.
            endif
            endif
         enddo
C$OMP END DO
      endif

c --- bright random stars

C$OMP DO
      do it = 1, ltrace
         if ( inside(it) ) then

* --- INITIALIZE INDEX-ARRAYS, find allowable index

            ipx(it) = int( pax(it) * rvx ) + lone
            ipy(it) = int( pay(it) * rvx ) + lone
            ipz(it) = int( paz(it) * rvx ) + lone

            dex(it)= abs(rkdx(it))
            dey(it)= abs(rkdy(it))
            dez(it)= abs(rkdz(it))

            sgnx(it) = int( sign( one, dirx(it) ) )
            sgny(it) = int( sign( one, diry(it) ) )
            sgnz(it) = int( sign( one, dirz(it) ) )

            t0(it) = tmi(it)* rvx

c            tdepx(it) = (dble(ipx(it)-lone)
c             tdepx(it) = (dble(ipx(it)-CVMGT(1,0,dirx(it).gt.zero))
            tdepx(it) = ( dble(ipx(it)-max(sgnx(it),lzero))
     &            -orix(it)*rvx) * rkdx(it) + dex(it)
c            tdepy(it) = (dble(ipy(it)-lone)
c             tdepy(it) = (dble(ipy(it)-CVMGT(1,0,diry(it).gt.zero))
            tdepy(it) = ( dble(ipy(it)-max(sgny(it),lzero))
     &            -oriy(it)*rvx) * rkdy(it) + dey(it)
c            tdepz(it) = (dble(ipz(it)-lone)
c             tdepz(it) = (dble(ipz(it)-CVMGT(1,0,dirz(it).gt.zero))
            tdepz(it) = ( dble(ipz(it)-max(sgnz(it),lzero))
     &            -oriz(it)*rvx) * rkdz(it) + dez(it)

         endif
      enddo
C$OMP END DO
C$OMP END PARALLEL

      itmi = lone
      itma = ltrace
      left = itma - itmi + lone

* ------
* --- ... INNER LOOP !
* ------

100   continue

      do it = itmi, itma
         if (inside(it)) goto 151
      enddo
      it = itma + lone
 151  itmt = it
      do it = itma, itmi, -1
         if (inside(it)) goto 161
      enddo
      it = itmi - lone
 161  itma = it
      itmi = itmt

      left = itma - itmi + lone

      if ( left .gt. lzero ) then
* --- ... GATHER DATA VALUES

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(it),
C$OMP+            SHARED(itmi,itma,ipx,ipy,ipz,inside, vdelx,
C$OMP+                   tempr,temp,ntemp,ntau,denvo,rho,semi,emit,
C$OMP+                   tdepx,tdepy,tdepz,tpsi,t0,t1,ekap,
C$OMP+                   rsum,gsum,bsum,fac,rtempr,gtempr,btempr,
C$OMP+                   dex,dey,dez,sgnx,sgny,sgnz,opac)
      do it = itmi, itma
      if ( inside(it) ) then

         ntau(it) = ntemp(ipx(it),ipy(it),ipz(it))
         rho (it) = denvo(ipx(it),ipy(it),ipz(it))

* --- ... SUM UP ALONG RAY

         tpsi(it) = min( tdepx(it), tdepy(it), tdepz(it))
         t1(it) = tpsi(it)

         ekap(it) = exp( -rho(it) * opac * ( t1(it)-t0(it) ) )
c         emit(it) = ( one - ekap(it)**2 )
         emit(it) = ( one - ekap(it) )
         semi(it) = emit(it) * fac(it)

         rsum(it) = rsum(it) + semi(it) * rtempr(ntau(it))
         gsum(it) = gsum(it) + semi(it) * gtempr(ntau(it))
         bsum(it) = bsum(it) + semi(it) * btempr(ntau(it))

         fac (it) = fac(it) * ekap(it)

* --- ... ADVANCE BEAM POSITION

         if (fac(it) .lt. 1D-3) inside(it) = .false.

         t0(it) = t1(it)
         if ( tpsi(it) .eq. tdepx(it) ) then
            tdepx(it) = tdepx(it) + dex(it)
            ipx(it) = ipx(it) + sgnx(it)
            if ((ipx(it).gt.vx).or.(ipx(it).lt.1)) inside(it)=.false.
         elseif ( tpsi(it) .eq. tdepy(it) ) then
            tdepy(it) = tdepy(it) + dey(it)
            ipy(it) = ipy(it) + sgny(it)
            if ((ipy(it).gt.vy).or.(ipy(it).lt.1)) inside(it)=.false.
         else
            tdepz(it) = tdepz(it) + dez(it)
            ipz(it) = ipz(it) + sgnz(it)
            if ((ipz(it).gt.vz).or.(ipz(it).lt.1)) inside(it)=.false.
         endif

      endif
      enddo
C$OMP END PARALLEL DO

      endif

      if (left .gt. lzero ) goto 100

* ------
* --- END OF INNER LOOP
* ------

c --- bright random stars

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(it), SHARED(ltrace,inside)
      do it = 1, ltrace
         inside(it) = .true.    ! reuse array inside for new purpose
      enddo
C$OMP END PARALLEL DO

      if (lstar .eq. lone) then
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(it,nitx,nity),
C$OMP+            SHARED(ltrace,inside,fac,xpix,ypix,bpix,
C$OMP+                   rsum,gsum,bsum,rpix,gpix)
         do it = 1, ltrace
            if ( fac(it) .gt. 1D-3 ) then
               nitx = nint(xpix(it))
               nity = nint(ypix(it))
               rsum(it) = max( fac(it) +rsum(it), rpix(nitx,nity))
               gsum(it) = max( fac(it) +gsum(it), gpix(nitx,nity))
               bsum(it) = max( fac(it) +bsum(it), bpix(nitx,nity))
            else
               inside(it) = .false.    ! reuse array inside for new purpose
            endif
         enddo
C$OMP END PARALLEL DO
      endif

c ---

      if ( nover .eq. lzero ) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(it,nitx,nity),
C$OMP+   SHARED(ltrace,xpix,ypix,rsum,bsum,gsum,rpix,gpix,bpix,intpix)
         do it = 1, ltrace
            nitx = nint(xpix(it))
            nity = nint(ypix(it))
            rpix(nitx,nity) = rsum(it)
            gpix(nitx,nity) = gsum(it)
            bpix(nitx,nity) = bsum(it)
            intpix(nitx,nity) = lzero !  mark traced pixels
         enddo
C$OMP END PARALLEL DO

         if ( lstar .eq. lone ) then
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(it,nitx,nity),
C$OMP+            SHARED(ltrace,xpix,ypix,rpix,gpix,bpix,intpix,
C$OMP+                   rsum,gsum,bsum)
         do it = 1, ltrace
            nitx = nint(xpix(it))
            nity = max( nint(ypix(it))-1, 0)
            rpix(nitx,nity) = rsum(it)
            gpix(nitx,nity) = gsum(it)
            bpix(nitx,nity) = bsum(it)
            intpix(nitx,nity) = lzero !  mark traced pixels
         enddo
C$OMP END PARALLEL DO
         endif

      else

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(it,nitx,nity)
C$OMP+            SHARED(ltrace,xpix,ypix,rpix,gpix,bpix,rsum,gsum,bsum)
         do it = 1, ltrace, 3
            nitx = nint(xpix(it))
            nity = nint(ypix(it))
            rpix(nitx,nity) = 0.25D0 * (  rpix(nitx,nity) +
     &                        rsum(it) + rsum(it+1) + rsum(it+2)  )
            gpix(nitx,nity) = 0.25D0 * (  gpix(nitx,nity) +
     &                        gsum(it) + gsum(it+1) + gsum(it+2)  )
            bpix(nitx,nity) = 0.25D0 * (  bpix(nitx,nity) +
     &                        bsum(it) + bsum(it+1) + bsum(it+2)  )
         enddo
C$OMP END PARALLEL DO

      endif

      ltrace = lzero


      return
      end
c -------------------------------------------------------------
      subroutine edgevisr( arr, val )
c  write val into edges of array arr

      include 'qparam.cmn'
      include 'rays.cmn'
      real*8  arr(vx,vy,vz), val

      integer*4 ix(4), iy(4), iz(4)
      data  ix / 1, vx, 2, vx /
      data  iy / 1, vy, 2, vy /
      data  iz / 1, vz, 2, vz /

      ix(4) = vx - 1         ! this cannot be done
      iy(4) = vy - 1         ! in the data-statements
      iz(4) = vz - 1         ! above

      do 100 j = 1, 2
      do 100 i = 1, 2

         do 20 n = 1, vx
 20         arr( n, iy(i), iz(j) ) = val

         do 25 n = 1, vy
 25         arr( ix(i), n, iz(j) ) = val

         do 30 n = 1, vz
 30         arr( ix(i), iy(j), n ) = val

 100  continue

      return
      end
c -------------------------------------------------------------
      subroutine edgevisi( iarr, ival )
c  write val into edges of array arr

      include 'qparam.cmn'
      include 'rays.cmn'
      integer*4  iarr(vx,vy,vz), ival

      integer*4 ix(4), iy(4), iz(4)
      data  ix / 1, vx, 2, vx /
      data  iy / 1, vy, 2, vy /
      data  iz / 1, vz, 2, vz /

      ix(4) = vx - 1         ! this cannot be done
      iy(4) = vy - 1         ! in the data-statements
      iz(4) = vz - 1         ! above

      do 100 j = 1, 2
      do 100 i = 1, 2

         do 20 n = 1, vx
 20         iarr( n, iy(i), iz(j) ) = ival

         do 25 n = 1, vy
 25         iarr( ix(i), n, iz(j) ) = ival

         do 30 n = 1, vz
 30         iarr( ix(i), iy(j), n ) = ival

 100  continue

      return
      end
c -------------------------------------------------------------
      subroutine cornvisr( arr, val )
c  write val into corners of array arr

      include 'qparam.cmn'
      include 'rays.cmn'
      real*8 arr(vx,vy,vz), val

      integer*4 ix(4), iy(4), iz(4)
      data  ix / 1, vx, 2, vx /
      data  iy / 1, vy, 2, vy /
      data  iz / 1, vz, 2, vz /

      ix(4) = vx - 1         ! this cannot be done
      iy(4) = vy - 1         ! in the data-statements
      iz(4) = vz - 1         ! above

      do 100 j = 1, 2
      do 100 i = 1, 2

         do 20 n = 1, 6
   20       arr( n, iy(i), iz(j) ) = val

         do 21 n = vx-5, vx
   21       arr( n, iy(i), iz(j) ) = val

         do 25 n = 1, 6
   25       arr( ix(i), n, iz(j) ) = val

         do 26 n = vy-5, vy
   26       arr( ix(i), n, iz(j) ) = val

         do 30 n = 1, 6
   30       arr( ix(i), iy(j), n ) = val

         do 31 n = vz-5, vz
   31       arr( ix(i), iy(j), n ) = val

 100  continue

      return
      end
c -------------------------------------------------------------
      subroutine cornvisi( iarr, ival )
c  write val into corners of array arr

      include 'qparam.cmn'
      include 'rays.cmn'
      integer*4 iarr(vx,vy,vz), ival

      integer*4 ix(4), iy(4), iz(4)
      data  ix / 1, vx, 2, vx /
      data  iy / 1, vy, 2, vy /
      data  iz / 1, vz, 2, vz /

      ix(4) = vx - 1         ! this cannot be done
      iy(4) = vy - 1         ! in the data-statements
      iz(4) = vz - 1         ! above

      do 100 j = 1, 2
      do 100 i = 1, 2

         do 20 n = 1, 6
   20       iarr( n, iy(i), iz(j) ) = ival

         do 21 n = vx-5, vx
   21       iarr( n, iy(i), iz(j) ) = ival

         do 25 n = 1, 6
   25       iarr( ix(i), n, iz(j) ) = ival

         do 26 n = vy-5, vy
   26       iarr( ix(i), n, iz(j) ) = ival

         do 30 n = 1, 6
   30       iarr( ix(i), iy(j), n ) = ival

         do 31 n = vz-5, vz
   31       iarr( ix(i), iy(j), n ) = ival

 100  continue

      return
      end
c -------------------------------------------------------------
      subroutine gridnest(vicenx,viceny,vicenz,rnf)
c     nest many grids into one to be rendered         June 1993

      include 'qparam.cmn'
      include 'aquants.cmn'
      include 'squants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'

      include 'rays.cmn'


      real*8    denvo(vx,vy,vz), tempr(vx,vy,vz),
     &          rpix(0:npix,0:npix), gpix(0:npix,0:npix),
     &          bpix(0:npix,0:npix)

      real*8    densbh(qx,qy,qz,ngrd)
      common /raygrnest/ densbh

      integer*4 ntemp(vx,vy,vz), intpix(0:npix,0:npix)

      common /raypot/ denvo, tempr, rpix, gpix, bpix, ntemp, intpix

      real*8   fx(vx,ngrd), fy(ngrd), fz(ngrd)
      integer*4 ixt(vx)
      real*8  wx0(vx), wx1(vx), wy0(vx), wy1(vx),
     &            wz0(vx), wz1(vx)
      integer*4 ix0(vx), ix1(vx), iy0(vx), iy1(vx), iz0(vx), iz1(vx)
      integer*4 izm(vx), igv(vx)

      integer*4 ilo, iyt, izt, izmt, vzh

      real*8 vicenx, viceny, vicenz, rnf, rnfin, rvxh, rdel,
     &       rcex, rcey, rcez

      ilo = mlev
c      ilo = indxgr(ltwo,maxlev(ltwo))
c      if ( ilo .lt. 6 ) then
c         write(*,*)  'change value of loop 11 and 13 below'
c         stop
c      endif

      rnfin = dble(nfine)

c      rvxh = dble(vx)/two
c      rcex = rvxh - rvxh/3.0D0
      rcex = dble(vx)/two
      rcey = dble(vy)/two
      vzh  = vz / ltwo
      rcez = zero
      rdel = dble(qx)/dble(vx)

c      vdelx= gridlx / rnf / dble(vx)   ! for ray integration

      do igrd = 1, ngrd
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(igrd,bhdens,densty,temper,densbh)
         do j = 1, qy
         do k = 1, qz
         do i = 1, qx
            if ( bhdens(i,j,k,igrd) .ne. zero ) then
               densbh(i,j,k,igrd) = 1D30
            else
               densbh(i,j,k,igrd) = densty(i,j,k,igrd)
            endif
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO
      enddo

c --- run over graphic cells

C$OMP PARALLEL DO PRIVATE(k,j,fy,fz,levl,num,igrd,ifgr,i,fx,igv),
C$OMP+            PRIVATE(ixt,iyt,izt,ix0,ix1,iy0,iy1,iz0,iz1),
C$OMP+            PRIVATE(wx0,wx1,wy0,wy1,wz0,wz1,izm,izmt)
      do j = 1, vy
      do k = 1, vzh

        fy(  1) = (dble(j)-half-rcey)*rdel / rnf  + viceny
        fz(  1) = (dble(k)-half-rcez)*rdel / rnf  + vicenz

        do levl = 1, mlev-1
        do num  = 1, indxgr(lzero,levl)
          igrd = indxgr(num,levl)
          ifgr = idxfgr(igrd)
          fy(  ifgr) = (fy(  igrd)-dble(norgin(ifgr,2)) -0.25D0)*rnfin
          fz(  ifgr) = (fz(  igrd)-dble(norgin(ifgr,3)) -0.25D0)*rnfin
        enddo
        enddo

      do i = 1, vx

c --- positions relative to topmost grid
        fx(i,1) = (dble(i)-half-rcex)*rdel / rnf  + vicenx

c       find position on finer grids, 1 to ilo
c       the following is copied from fipos

         do levl = 1, mlev-1
         do num  = 1, indxgr(lzero,levl)
            igrd = indxgr(num,levl)
            ifgr = idxfgr(igrd)
            fx(i,ifgr) = (fx(i,igrd)-dble(norgin(ifgr,1)) -0.25D0)*rnfin
         enddo
         enddo

         igv(i) = 0
      enddo

      do levl = ilo, 1, -1
      do num  = 1, indxgr(lzero,levl)
         igrd = indxgr(num,levl)

      do i = 1, vx

         ixt(i) = int(fx(i,igrd))
         iyt    = int(fy(  igrd))
         izt    = int(fz(  igrd))

	 izmt = max( izt, lone )

c         write(*,*) 'rays ',i,igrd,ixt(i),iyt,izt
         if ( (igv(i).eq.lzero) .and.
     &        (ixt(i).ge.lone ) .and. (ixt(i).le.qx-1) .and.
     &        (iyt   .ge.lone ) .and. (iyt   .le.qy-1) .and.
     &        (izmt  .ge.lone ) .and. (izmt  .le.qz-1)        ) then

            ix0(i) = ixt(i)
            iy0(i) = iyt
            iz0(i) = izt
            izm(i) = izmt

            ix1(i) = ix0(i) + lone
            iy1(i) = iy0(i) + lone
            iz1(i) = iz0(i) + lone

            wx1(i) = fx(i,igrd) - dble(ix0(i))
            wx0(i) = one - wx1(i)
            wy1(i) = fy(  igrd) - dble(iy0(i))
            wy0(i) = one - wy1(i)
            wz1(i) = fz(  igrd) - dble(iz0(i))
            wz0(i) = one - wz1(i)

            igv(i) = igrd

         endif
      enddo
      enddo
      enddo

      do i = 1, vx
         if (igv(i).eq.lzero) then
            write(*,*) 'video gridnest, igrd .eq. 0 ',i,j,k
            write(*,*) igv
            stop
         endif
      enddo

      do i = 1, vx
        denvo(i,j,vzh+k) = log10(
     &    wx0(i)*(wy0(i)* (wz0(i)*densbh(ix0(i),iy0(i),izm(i),igv(i))
     &                    +wz1(i)*densbh(ix0(i),iy0(i),iz1(i),igv(i)) )
     &           +wy1(i)* (wz0(i)*densbh(ix0(i),iy1(i),izm(i),igv(i))
     &                    +wz1(i)*densbh(ix0(i),iy1(i),iz1(i),igv(i)) ))
     &   +wx1(i)*(wy0(i)* (wz0(i)*densbh(ix1(i),iy0(i),izm(i),igv(i))
     &                    +wz1(i)*densbh(ix1(i),iy0(i),iz1(i),igv(i)) )
     &           +wy1(i)* (wz0(i)*densbh(ix1(i),iy1(i),izm(i),igv(i))
     &                    +wz1(i)*densbh(ix1(i),iy1(i),iz1(i),igv(i)))))

        tempr(i,j,vzh+k) = log10(
     &    wx0(i)*(wy0(i)* (wz0(i)*densbh(ix0(i),iy0(i),izm(i),igv(i))
     &                    +wz1(i)*densbh(ix0(i),iy0(i),iz1(i),igv(i)) )
     &           +wy1(i)* (wz0(i)*densbh(ix0(i),iy1(i),izm(i),igv(i))
     &                    +wz1(i)*densbh(ix0(i),iy1(i),iz1(i),igv(i)) ))
     &   +wx1(i)*(wy0(i)* (wz0(i)*densbh(ix1(i),iy0(i),izm(i),igv(i))
     &                    +wz1(i)*densbh(ix1(i),iy0(i),iz1(i),igv(i)) )
     &           +wy1(i)* (wz0(i)*densbh(ix1(i),iy1(i),izm(i),igv(i))
     &                    +wz1(i)*densbh(ix1(i),iy1(i),iz1(i),igv(i)))))

        denvo(i,j,vzh+1-k) = denvo(i,j,vzh+k)
        tempr(i,j,vzh+1-k) = tempr(i,j,vzh+k)

      enddo

      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine putclock(timu)
      include 'qparam.cmn'
      include 'rays.cmn'

      real*8  timu, ra, xm, ym, an, pi
      integer   i, j, ip, im, jp, jm
      integer*1 matrgb(3,0:npix,0:npix)
      common /coclock/ matrgb

      pi = 4.D0 * atan(1.D0)

      ra = dble(npix+1)/20.D0
      xm = 1.5D0*ra
      ym = 2.0D0*ra
      if ( iste .eq. +1 ) xm = xm - 7D0*float(npix)/511.D0

c --- draw circular border of clock

      do ian = 0, int(2.D0*pi/(0.8D0/ra))
         an = dble(ian) * 0.8D0/ra
         i = max( min( nint( xm + ra * sin(an) ), npix-1 ), 1 )
         j = max( min( nint( ym + ra * cos(an) ), npix-1 ), 1 )
         ip = i + 1
         im = i - 1
         jp = j + 1
         jm = j - 1
         matrgb(1,ip,j ) = 150
         matrgb(1,im,j ) = 150
         matrgb(1,i ,jp) = 150
         matrgb(1,i ,jm) = 150
         matrgb(2,ip,j ) = 150
         matrgb(2,im,j ) = 150
         matrgb(2,i ,jp) = 150
         matrgb(2,i ,jm) = 150
         matrgb(3,ip,j ) = 150
         matrgb(3,im,j ) = 150
         matrgb(3,i ,jp) = 150
         matrgb(3,i ,jm) = 150
      enddo

      do ian = 0, int(2.D0*pi/(0.4D0/ra))
         an = dble(ian) * 0.4D0/ra
         i = min( nint( xm + ra * sin(an) ), npix-1 )
         j = min( nint( ym + ra * cos(an) ), npix-1 )
         matrgb(1,i ,j ) = 255
         matrgb(2,i ,j ) = 255
         matrgb(3,i ,j ) = 255
      enddo

c --- draw ticks

      do ian = 0, int(0.901/0.1)
         an = dble(ian) * 0.1
         call clockhand(xm, ym, 0.8*ra, ra, an * 2.D0 * pi, 1.0D0 )
      enddo

c --- draw hands

      call clockhand(xm, ym, zero,       ra, timu * 2.D-1 * pi, 1.0D0 )
      call clockhand(xm, ym, 0.6D0*ra,   ra, timu * 2.D0  * pi, 1.0D0 )
      call clockhand(xm, ym, zero, 0.4D0*ra, timu * 2.D-2 * pi, 1.0D0 )

      return
      end
c-----------------------------------------------------------------------
      subroutine clockhand(xm,ym,ri,re,an,col)
c --- draw hand

      include 'qparam.cmn'
      include 'rays.cmn'

      real*8    ri, re, an, san, can, xm, ym, col
      integer   i, j, ip, im, jp, jm, icol1, icol2
      integer*1 matrgb(3,0:npix,0:npix)
      common /coclock/ matrgb

      san = sin(an)
      can = cos(an)

      icol1 = nint(col*150)
      icol2 = nint(col*255)

      do m = 0, nint((re-ri)/0.8D0)
         r = ri + 0.8D0 * m
         i = max( min( nint( xm + r * san ), npix-1), 1 )
         j = max( min( nint( ym - r * can ), npix-1), 1 )
         ip = i + 1
         im = i - 1
         jp = j + 1
         jm = j - 1
         matrgb(1,ip,j ) = icol1
         matrgb(1,im,j ) = icol1
         matrgb(1,i ,jp) = icol1
         matrgb(1,i ,jm) = icol1
         matrgb(2,ip,j ) = 150
         matrgb(2,im,j ) = 150
         matrgb(2,i ,jp) = 150
         matrgb(2,i ,jm) = 150
         matrgb(3,ip,j ) = 150
         matrgb(3,im,j ) = 150
         matrgb(3,i ,jp) = 150
         matrgb(3,i ,jm) = 150
      enddo

      do m = 0, nint((re-ri)/0.4D0)
         r = ri + 0.4D0 * m
         i = max( min( nint( xm + r * san ), npix-1), 1 )
         j = max( min( nint( ym - r * can ), npix-1), 1 )
         matrgb(1,i ,j ) = icol2
         matrgb(2,i ,j ) = 255
         matrgb(3,i ,j ) = 255
      enddo

      return
      end
