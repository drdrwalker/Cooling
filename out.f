c-----------------------------------------------------------------------
      subroutine detect(al,a,ar,smalla)

c     search for contact discontinuities in variable a and
c     steepen the zone structure if necessary

      include 'qparam.cmn'
      include 'squants.cmn'
      include   'ppms.cmn'

      real*8 scrch1(q), scrch2(q), scrch3(q), scrch4(q)

      real*8  al(q), a(q), ar(q), smalla

      print *, 'out, detect'


c------- the following parameters are set as in Colella and
c        Woodward (JCP, 54 (1984), 174)

      eta1  = 20.0D0
      eta2  = 0.05D0
      epsln = 0.01D0
      ak0   = 0.1D0

      do i = 3, np7
         scrch1(i) = dx(i) + dx(i-1)
         scrch2(i) = scrch1(i) + dx(i+1)
         scrch3(i) = a(i) - a(i-1)
         scrch1(i) = scrch3(i) / scrch1(i)
      enddo

      do i = 3, np6
         scrch2(i) = (scrch1(i+1)-scrch1(i)) / scrch2(i)
      enddo

      do i = 4, np6
c        scrch1(i) = x(i) - x(i-1)
         scrch1(i) = dx(i)
         scrch1(i) = scrch1(i)*scrch1(i)*scrch1(i)
      enddo

      do i = 4, np5
         scrch4(i) = a(i+1)-a(i-1)
         scrch4(i) = sign(max(abs(scrch4(i)),small*smalla),scrch4(i))
         scrch3(i) = (scrch2(i-1)-scrch2(i+1))*(scrch1(i)+scrch1(i+1))
         scrch3(i) = scrch3(i) / (two*dx(i)*scrch4(i))
c        scrch3(i) = scrch3(i) / ((x(i+1)-x(i-1))*scrch4(i))
      enddo

c     scrch2 and scrch3 now contain finite difference approximations
c     to the second and third derivatives of a.

      do i = 4, np5
         if ( scrch2(i-1)*scrch2(i+1) .ge. zero )   scrch3(i) = zero
      enddo

      do i = 3, np6
         scrch1(i) = abs(a(i))
         scrch2(i) = abs(a(i+1)-a(i-1))
      enddo

      do i = 4, np5
         scrch4(i) = epsln*min(scrch1(i+1),scrch1(i-1)) - scrch2(i)
         if ( scrch4(i) .ge. zero )  scrch3(i) = zero
         scrch3(i) = max(zero,min(eta1*(scrch3(i)-eta2),one))
      enddo

      do i = 4, np5
         scrch1(i) = abs(  p(i+1)-  p(i-1)) / min(  p(i+1),  p(i-1))
         scrch2(i) = abs(rho(i+1)-rho(i-1)) / min(rho(i+1),rho(i-1))
         if ( game(i)*ak0*scrch2(i)-scrch1(i) .lt. zero ) scrch3(i)=zero
      enddo

c     scrch3 now contains the contact steepening coefficient

      do i = 4, np5
         scrch1(i) = a(i-1) + half*dela(i-1)
         scrch2(i) = a(i+1) - half*dela(i+1)
      enddo

      do i = 4, np5
         al(i) = al(i)*(one-scrch3(i)) + scrch1(i)*scrch3(i)
         ar(i) = ar(i)*(one-scrch3(i)) + scrch2(i)*scrch3(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine flaten

c     flaten zone structure in regions where shocks are too thin

      include 'qparam.cmn'
      include 'squants.cmn'
      include   'ppms.cmn'

      real*8 scrch1(q), scrch2(q), scrch3(q), scrch4(q)

      print *, 'out, flaten'


c------- This version of subroutine FLATEN only uses the simplest
c        form of dissipation as described in the appendix of
c        Colella and Woodward (JCP, 54 (1984), 174). Therefore
c        the only constants required are omg1, omg2 and epsiln,
c        which are read in.
c        The "standard" values of the constants are:

c           epsiln = 0.33

c           omg1   = 0.75
c           omg2   = 10.0
c           sig1   = 0.50
c           sig2   = 1.00
c           ak1    = 2.00
c           ak2    = 0.01

c           wig1   = 2.00
c           wig2   = 0.00           for 1-d
c                    0.10           for 2-d
c           wig3   = 0.3333 - wig2


      do i = 2, np7
         dp(i)     = p(i+1) - p(i-1)
         du(i)     = u(i+1) - u(i-1)
         scrch1(i) = epsiln * min(p(i+1), p(i-1)) - abs( dp(i) )
      enddo

      do i = 2, np7
         scrch1(i) = max( zero, -sign(one,scrch1(i)) )
         if ( du(i) .ge. zero )  scrch1(i) = zero
      enddo

      do i = 3, np6
         dp2 = p(i+2) - p(i-2)
         dp2 = sign( max( abs(dp2), smallp ), dp2 )
         scrch2(i) = dp(i) / dp2 - omg1
         scrch3(i) = scrch1(i) * max (zero, scrch2(i) * omg2)
      enddo

      do i = 4, np5
         if ( dp(i) .lt. zero ) then
            scrch2(i) = scrch3(i+1)
         else
            scrch2(i) = scrch3(i-1)
         endif
      enddo

      do i = 4, np5
         flatn(i) = max (scrch3(i), scrch2(i))
         flatn(i) = max (zero, min(one, flatn(i)))
      enddo

      do i = 4, np5
         flatn (i) = flatn(i) * (one - igodu) + igodu
         flatn1(i) = one - flatn(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine monot(al,a,ar,da,a6)

c     apply monotonicity constraint to interpolation parabola
c     and calculate parabola coefficients da,a6

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 scrch1(q), scrch2(q), scrch3(q), scrch4(q)

      real*8  al(q), a(q), ar(q), da(q), a6(q)

      print *, 'out, monot'


      do i = 4, np5

         da(i)     =   ar(i) - al(i)
         scrch1(i) = ( ar(i) - a (i) )  *  ( al(i) - a(i) )

         if ( scrch1(i) .ge. zero ) then
            al(i) = a(i)
            ar(i) = a(i)
         endif

cc       scrch2(i) = 3.D0 * a(i)  -  two * ar(i)  !  cc
cc       scrch3(i) = 3.D0 * a(i)  -  two * al(i)  !  cc

c------- statements to avoid generation of spurious noise in case of
c        a constant state in one coordinate direction

         if ( (ar(i)-a(i)) * (al(i)-a(i)) .ne. zero ) then
            scrch2(i) = 3.D0*a(i)-two*ar(i)
            scrch3(i) = 3.D0*a(i)-two*al(i)
         else
            scrch2(i) = al(i)
            scrch3(i) = ar(i)
         endif

         if ( da(i) * (al(i)-scrch2(i)) .lt. zero )  al(i) = scrch2(i)
         if ( da(i) * (scrch3(i)-ar(i)) .lt. zero )  ar(i) = scrch3(i)

         da(i) = ar(i) - al(i)
         a6(i) = 6.D0 * a(i)  -  3.D0 * ( al(i) + ar(i) )

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine states

c     compute left and right states for input to riemann problem
c     input :  hydro , interp  , output :  reman

      include 'qparam.cmn'
      include 'squants.cmn'
      include   'ppms.cmn'

      real*8 scrch1(q), scrch2(q), scrch3(q), scrch4(q)

      print *, 'out, states'


      do  5  i = 4, np5
         urel  (i) = u(i) - ugrid(i)
c        dphidx(i) = ( gravr(i) - gravl(i) ) / dx(i)
         dphidx(i) = ( gravr(i) - gravl(i) ) / dx(i)  +  dgrav(i)
   5     continue


c------- calculation of left states

      do 10 i = 5, np5
      scrch1(i)=dtdx(i-1)*(urel(i-1)+ce(i-1))
      cflno(i)=max(zero,scrch1(i))
      scrch1(i)=half*min(one,cflno(i))
      scrch2(i)=one-forthd*scrch1(i)
      ppl  (i) = pr  (i-1) - scrch1(i)*(dp  (i-1)-scrch2(i)*p6  (i-1))
      ppl  (i) = max(smallp,ppl(i))
      upl  (i) = ur  (i-1) - scrch1(i)*(du  (i-1)-scrch2(i)*u6  (i-1))
      rhopl(i) = rhor(i-1) - scrch1(i)*(drho(i-1)-scrch2(i)*rho6(i-1))
   10 continue

      do 20 i = 5, np5
      scrch3(i)=dtdx(i-1)*(urel(i-1)-ce(i-1))
      scrch1(i)=half*min(one,max(zero,scrch3(i)))
      scrch2(i)=one-forthd*scrch1(i)
      pml  (i) = pr  (i-1) - scrch1(i)*(dp  (i-1)-scrch2(i)*p6  (i-1))
      pml  (i) = max(smallp,pml(i))
      uml  (i) = ur  (i-1) - scrch1(i)*(du  (i-1)-scrch2(i)*u6  (i-1))
      rhoml(i) = rhor(i-1) - scrch1(i)*(drho(i-1)-scrch2(i)*rho6(i-1))
  20  continue

      do 30 i = 5, np5
      scrch4(i) = dtdx(i-1)*urel(i-1)
      scrch1(i) = half*min(one,max(zero,scrch4(i)))
      scrch2(i) = one-forthd*scrch1(i)
      p0l  (i) = pr  (i-1) - scrch1(i)*(dp  (i-1)-scrch2(i)*p6  (i-1))
      p0l  (i) = max( p0l(i), smallp )
c?    u0l  (i) = ur  (i-1) - scrch1(i)*(du  (i-1)-scrch2(i)*u6  (i-1))
      rho0l(i) = rhor(i-1) - scrch1(i)*(drho(i-1)-scrch2(i)*rho6(i-1))
      ut0l (i) = utr (i-1) - scrch1(i)*(dut (i-1)-scrch2(i)*ut6 (i-1))
      utt0l(i) = uttr(i-1) - scrch1(i)*(dutt(i-1)-scrch2(i)*utt6(i-1))
      game0l(i) = gamer(i-1)-scrch1(i)*(dgame(i-1)-scrch2(i)*game6(i-1))
      gamc0l(i) = gamcr(i-1)-scrch1(i)*(dgamc(i-1)-scrch2(i)*gamc6(i-1))
  30  continue

      do m = 1, qc
      do i = 5, np5
         chelft(i,m) =
     &     cher(i-1,m) - scrch1(i) * (dche(i-1,m)-scrch2(i)*che6(i-1,m))
      enddo
      enddo

      do 40 i = 5, np5
      clft(i)  = sqrt( gamc(i-1) * ppl(i) * rhopl(i) )
      scrch1(i)= half * ( upl(i) - uml(i) - (ppl(i)-pml(i))/clft(i) )
      scrch2(i)= (ppl(i)-p0l(i)) / (clft(i)*clft(i)) + one/rhopl(i)
      scrch2(i)= scrch2(i) - one / rho0l(i)
      if ( scrch3(i) .le. zero ) scrch1(i) = zero
      if ( scrch4(i) .le. zero ) scrch2(i) = zero
c      scrch1(i)= cvmgm( scrch1(i), zero, -scrch3(i) )
c      scrch2(i)= cvmgm( scrch2(i), zero, -scrch4(i) )

      plft  (i) = ppl(i) + clft(i) * scrch1(i)
      plft  (i) = max( plft(i), smallp )
      ulft  (i) = upl(i) - scrch1(i)
      vlft  (i) = one / rhopl(i) - scrch2(i) - scrch1(i) / clft(i)
      utlft (i) = ut0l  (i)
      uttlft(i) = utt0l (i)
      gmelft(i) = game0l(i)
      gmclft(i) = gamc0l(i)
  40  continue

c      do 462  i = 4, np5
c 462     dloga(i) = zero

      do 49 i = 5, np5
cdloga         scrch1(i) = half*rho(i-1)*u(i-1)*dtppm*dloga(i-1)
cdloga         vlft(i)   = one/vlft(i)-scrch1(i)
cdloga         vlft(i)   = one/vlft(i)
cdloga         plft(i)   = plft(i)-scrch1(i)*ce(i-1)*ce(i-1)
         plft(i)   = max(plft(i),smallp)
c        ulft(i)   = ulft(i)  +  half * dt * ( fict(i-1) - dphidx(i-1) )
c                          planar geometry : fict = 0.
         ulft(i)   = ulft(i)  -  half * dtppm * dphidx(i-1)
  49  continue

c------- calculation of right states

      do 50 i = 5, np5
         scrch3(i)=-dtdx(i)*(urel(i)+ce(i))
         scrch1(i)=half*min(one,max(zero,scrch3(i)))
         scrch2(i)=one-forthd*scrch1(i)
         ppl  (i) = pl  (i) + scrch1(i) * (dp  (i) + scrch2(i)*p6  (i))
         upl  (i) = ul  (i) + scrch1(i) * (du  (i) + scrch2(i)*u6  (i))
c?       rhopl(i) = rhol(i) + scrch1(i) * (drho(i) + scrch2(i)*rho6(i))
  50  continue

      do 60 i = 5, np5
         scrch1(i)=-dtdx(i)*(urel(i)-ce(i))
         scrch1(i)=max(zero,scrch1(i))
         cflno(i)=max(cflno(i),scrch1(i))
         scrch1(i)=half*min(one,scrch1(i))
         scrch2(i)=one-forthd*scrch1(i)
         pml  (i) = pl  (i) + scrch1(i) * (dp  (i) + scrch2(i)*p6  (i))
         pml  (i) = max(pml(i),smallp)
         uml  (i) = ul  (i) + scrch1(i) * (du  (i) + scrch2(i)*u6  (i))
         rhoml(i) = rhol(i) + scrch1(i) * (drho(i) + scrch2(i)*rho6(i))
  60  continue

      do 70 i = 5, np5
      scrch4(i) = -dtdx(i) * urel(i)
      scrch1(i) = half * min ( one, max (zero, scrch4(i)) )
      scrch2(i) = one - forthd * scrch1(i)
      p0l  (i) = pl  (i) + scrch1(i) * (dp  (i) + scrch2(i)*p6  (i))
      p0l  (i) = max(p0l(i), smallp)
      rho0l(i) = rhol(i) + scrch1(i) * (drho(i) + scrch2(i)*rho6(i))
c?    u0l  (i) = ul  (i) + scrch1(i) * (du  (i) + scrch2(i)*u6  (i))
      ut0l (i) = utl (i) + scrch1(i) * (dut (i) + scrch2(i)*ut6 (i))
      utt0l(i) = uttl(i) + scrch1(i) * (dutt(i) + scrch2(i)*utt6(i))
      game0l(i) = gamel(i) + scrch1(i) * (dgame(i) + scrch2(i)*game6(i))
      gamc0l(i) = gamcl(i) + scrch1(i) * (dgamc(i) + scrch2(i)*gamc6(i))
  70  continue


      do m = 1, qc
      do i = 5, np5
         chergt(i,m) =
     &          chel(i,m) + scrch1(i) * (dche(i,m)+scrch2(i)*che6(i,m))
      enddo
      enddo

      do 80 i = 5, np5
      crght (i) = sqrt( gamc(i) * pml(i) * rhoml(i) )
      scrch1(i) = -half*( uml(i)-upl(i) + (pml(i)-ppl(i)) / crght(i) )
      scrch2(i) = (pml(i)-p0l(i)) / (crght(i)*crght(i)) + one/rhoml(i)
      scrch2(i) = scrch2(i) - one / rho0l(i)
      if ( scrch3(i) .le. zero ) scrch1(i) = zero
      if ( scrch4(i) .le. zero ) scrch2(i) = zero
c      scrch1(i) = cvmgm( scrch1(i), zero, -scrch3(i) )
c      scrch2(i) = cvmgm( scrch2(i), zero, -scrch4(i) )

      prght (i) = pml(i) + crght(i) * scrch1(i)
      prght (i) = max( prght(i), smallp )
      urght (i) = uml(i) + scrch1(i)
      vrght (i) = one / rhoml(i) - scrch2(i) - scrch1(i) / crght(i)
      utrght(i) = ut0l  (i)
      uttrgt(i) = utt0l (i)
      gmergt(i) = game0l(i)
      gmcrgt(i) = gamc0l(i)
  80  continue

      do 90 i = 5, np5
cdloga      scrch1(i) = half*rho(i)*u(i)*dtppm*dloga(i)
cdloga      vrght(i)  = one/vrght(i)-scrch1(i)
cdloga      vrght(i)  = one/vrght(i)
cdloga      prght(i)  = prght(i)-scrch1(i)*ce(i)*ce(i)
      prght(i)  = max(prght(i),smallp)
c     urght(i)  = urght(i)  +  half * dt * ( fict(i) - dphidx(i) )
c                      planar geometry : fict = 0.
      urght(i)  = urght(i)  -  half * dtppm * dphidx(i)
  90  continue

c-------
      if(igeom.eq.0) return
      if(igeom.eq.3) return

c?    urght(3) = two *u (3) * cflno(3)
c?    ulft (3) = - urght(3)

      return
      end
