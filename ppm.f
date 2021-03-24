c-----------------------------------------------------------------------
      subroutine hydrow (jj,kk,igrd,nxyzsw)
c     perform hydrow step on one row of zones
c     igrd is needed for the artificial viscosity

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      real*8  scrch1(q), scrch2(q), scrch3(q), scrch4(q), fudge
      real*8 wav(q), wtav(q), wttav(q), factmp(q), rontmp(q), sgn
      real*8 wnu(q), wtnu(q), wttnu(q), ronma
      real*8 ww (q), wt  (q), wtt  (q), wmat(3,3,q), winv(3,3,q)
      integer*4 lfive

      lfive = 5

      j = jj
      k = kk

c-------   j (and k) because of wiggling in FLATEN

      call intrfc
      call states
      call riemann

c------- fluxes

      if     ( nxyzsw .eq. lone ) then
         lpi = lone
         lpj = ltwo
         lpk = lthree
      elseif ( nxyzsw .eq. ltwo ) then
         lpi = ltwo
         lpj = lone
         lpk = lthree
      elseif ( nxyzsw .eq. lthree ) then
         lpi = lthree
         lpj = lone
         lpk = ltwo
      endif

c     no speeds are faster than light
c      cfud = ten*cc

c      do i = 5, np5
c         uav (i) = sign( min( abs(uav (i)), cfud ), uav (i) )
c         urel(i) = sign( min( abs(urel(i)), cfud ), urel(i) )
c      enddo

c     relativistic change:   uav = kinematic velocity  is not the
c     same as  wav = dynamic velocity = specific momentum !

      grt = 0.8D0 * g / cc**5
      ronoff = dble(mconoff)

      do i3 = 1, 3
      do j3 = 1, 3
      do i = 5, np5
         wmat(i3,j3,i) = grt * dij3dt(i3,j3)
      enddo
      enddo
      enddo

      do i = 5, np5
         rontmp(i) = one + ronoff
     &    * sign( min( abs(two * grav(i) / cc**2), corrmax ), grav(i) )
         wmat(1,1,i) = wmat(1,1,i) + rontmp(i)
         wmat(2,2,i) = wmat(2,2,i) + rontmp(i)
         wmat(3,3,i) = wmat(3,3,i) + rontmp(i)
      enddo

      call mat3inv(wmat,winv,lfive,np5,q)

      do i = 5, np5
         wav  (i) =   winv(lpi,lpi,i) * uav  (i)
     &              + winv(lpj,lpi,i) * utav (i)
     &              + winv(lpk,lpi,i) * uttav(i)

         wtav (i) =   winv(lpi,lpj,i) * uav  (i)
     &              + winv(lpj,lpj,i) * utav (i)
     &              + winv(lpk,lpj,i) * uttav(i)

         wttav(i) =   winv(lpi,lpk,i) * uav  (i)
     &              + winv(lpj,lpk,i) * utav (i)
     &              + winv(lpk,lpk,i) * uttav(i)
      enddo

      do i = 5, np5
         rhoflx(i) = rhoav (i) * urel (i)
      enddo

      fudge = 4D0

      do i = 5, np4
         rhoflx(i)   = max( rhoflx(i)  , -fudge*rho(i)/dtdx(i) )
      enddo
      do i = 5, np4
         rhoflx(i+1) = min( rhoflx(i+1),  fudge*rho(i)/dtdx(i) )
      enddo

      do m = 1, qc
         do i = 5, np5
            cheflx(i,m) = cheav(i,m) * urel(i)
         enddo
         do i = 5, np4
            cheflx(i,m)   = max( cheflx(i,m)  ,-fudge*che(i,m)/dtdx(i) )
         enddo
         do i = 5, np4
            cheflx(i+1,m) = min( cheflx(i+1,m), fudge*che(i,m)/dtdx(i) )
         enddo
      enddo

      do i = 5, np5
         uflx  (i) = rhoflx(i) * wav  (i)
         utflx (i) = rhoflx(i) * wtav (i)
         uttflx(i) = rhoflx(i) * wttav(i)
         scrch1(i) = pav(i) / ( rhoav(i) * (gameav(i)-one) )
         scrch1(i) = scrch1(i)
     &        + half * ( wav(i)**2 + wtav(i)**2 + wttav(i)**2 )
c        eflx(i)   = rhoflx(i) * scrch1(i) + urel(i) * pav(i)
         eflx(i)   = rhoflx(i) * scrch1(i) + wav (i) * pav(i)
         scrch2(i) = zero
      enddo

c-------
c--- artificial viscosity
c-------

      if (cvisc .gt. 1.D-5)  then


      if (nxyzsw .eq. lone)   then

         j = jj
         k = kk

         sigyb = +lone
         if (bndmny(1,1) .eq. one)  sigyb = -lone
         sigyt = +lone
         if (bndmxy(1,1) .eq. one)  sigyt = -lone

         sigzb = +lone
         if (bndmnz(1,1) .eq. one)  sigzb = -lone
         sigzt = +lone
         if (bndmxz(1,1) .eq. one)  sigzt = -lone

         jb = max ( j-1, lone )
         kb = max ( k-1, lone )

         jt = min ( j+1, ny)
         kt = min ( k+1, nz)

         do i = 6, np4
         divy= ( sigyb * ( vyold(i-4,jb,k,igrd) + vyold(i-5,jb,k,igrd) )
     &          -sigyt * ( vyold(i-4,jt,k,igrd) + vyold(i-5,jt,k,igrd)))
     &         * half * dx(i)           / (two*dx(j+2)         )
         divz= ( sigzb * ( vzold(i-4,j,kb,igrd) + vzold(i-5,j,kb,igrd) )
     &          -sigzt * ( vzold(i-4,j,kt,igrd) + vzold(i-5,j,kt,igrd)))
     &         * half * dx(i)           / (two*dx(k+2)         )
         scrch2(i) = u(i-1) - u(i)  +  divy  +  divz
         enddo

      elseif (nxyzsw .eq. ltwo)   then

         i = jj
         k = kk

         sigxb = +1
         if (bndmnx(1,1) .eq. one)  sigxb = -1
         sigxt = +1
         if (bndmxx(1,1) .eq. one)  sigxt = -1

         sigzb = +1
         if (bndmnz(1,1) .eq. one)  sigzb = -1
         sigzt = +1
         if (bndmxz(1,1) .eq. one)  sigzt = -1

         ib = max ( i-1, lone )
         kb = max ( k-1, lone )

         it = min ( i+1, nx)
         kt = min ( k+1, nz)

         do j = 6, np4
         divx= ( sigxb * ( vxold(ib,j-4,k,igrd) + vxold(ib,j-5,k,igrd) )
     &          -sigxt * ( vxold(it,j-4,k,igrd) + vxold(it,j-5,k,igrd)))
     &         * half * dx(j)           / (two*dx(i+2)         )
         divz= ( sigzb * ( vzold(i,j-4,kb,igrd) + vzold(i,j-5,kb,igrd) )
     &          -sigzt * ( vzold(i,j-4,kt,igrd) + vzold(i,j-5,kt,igrd)))
     &         * half * dx(j)           / (two*dx(k+2)         )
            scrch2(j) = u(j-1) - u(j)  +  divx  +  divz
         enddo

      elseif (nxyzsw .eq. lthree)   then

         i = jj
         j = kk

         sigxb = +1
         if (bndmnx(1,1) .eq. one)  sigxb = -1
         sigxt = +1
         if (bndmxx(1,1) .eq. one)  sigxt = -1

         sigyb = +1
         if (bndmny(1,1) .eq. one)  sigyb = -1
         sigyt = +1
         if (bndmxy(1,1) .eq. one)  sigyt = -1

         ib = max ( i-1, lone )
         jb = max ( j-1, lone )

         it = min ( i+1, nx)
         jt = min ( j+1, ny)

         do k = 6, np4
         divx= ( sigxb * ( vxold(ib,j,k-4,igrd) + vxold(ib,j,k-5,igrd) )
     &          -sigxt * ( vxold(it,j,k-4,igrd) + vxold(it,j,k-5,igrd)))
     &         * half * dx(k)           / (two*dx(i+2)         )
         divy= ( sigyb * ( vyold(i,jb,k-4,igrd) + vyold(i,jb,k-5,igrd) )
     &          -sigyt * ( vyold(i,jt,k-4,igrd) + vyold(i,jt,k-5,igrd)))
     &         * half * dx(k)           / (two*dx(j+2)         )
            scrch2(k) = u(k-1) - u(k)  +  divx  +  divy
         enddo


      end if
      end if

c-------
c--- end of artificial viscosity
c-------


      do i = 5, np5
         scrch2(i) = cvisc * max (scrch2(i), zero)
         rhoflx(i) = rhoflx(i) + scrch2(i) *
     &                        (rho(i-1)            - rho(i)         )
         uflx  (i) = uflx  (i) + scrch2(i) *
     &                        (rho(i-1) * u  (i-1) - rho(i) * u  (i))
         utflx (i) = utflx (i) + scrch2(i) *
     &                        (rho(i-1) * ut (i-1) - rho(i) * ut (i))
         uttflx(i) = uttflx(i) + scrch2(i) *
     &                        (rho(i-1) * utt(i-1) - rho(i) * utt(i))
         eflx  (i) = eflx  (i) + scrch2(i) *
     &                        (rho(i-1) * e  (i-1) - rho(i) * e  (i))
      enddo

      do m = 1, qc
      do i = 5, np5
         cheflx(i,m) = cheflx(i,m) + scrch2(i) *
     &                        (che(i-1,m)            - che(i,m)     )
      enddo
      enddo

c-------

c     relativistic change:   unu = kinematic velocity  is not the
c     same as  wnu = dynamic velocity = specific momentum !

      do i = 5, np4

         ww (i) =     winv(lpi,lpi,i) * u  (i)
     &              + winv(lpj,lpi,i) * ut (i)
     &              + winv(lpk,lpi,i) * utt(i)

         wt (i) =     winv(lpi,lpj,i) * u  (i)
     &              + winv(lpj,lpj,i) * ut (i)
     &              + winv(lpk,lpj,i) * utt(i)

         wtt(i) =     winv(lpi,lpk,i) * u  (i)
     &              + winv(lpj,lpk,i) * ut (i)
     &              + winv(lpk,lpk,i) * utt(i)

         e  (i) =  e  (i)
     &       - half * ( u (i)**2 + ut(i)**2 + utt(i)**2 )
     &       + half * ( ww(i)**2 + wt(i)**2 + wtt(i)**2 )

      enddo


      do i = 5, np5
         uflx  (i) = uflx(i) + pav(i)
      enddo

      do i = 5, np4

         rhonu(i) = max(smlrho,rho(i)-dtdx(i)*(rhoflx(i+1)-rhoflx(i)))

         wnu  (i) = ( rho(i) * ww(i) -
     &                dtdx(i) * ( uflx  (i+1) - uflx  (i) ) ) / rhonu(i)
         if(wnu(i).gt.half*cc) then
           print *, 'wnu big boi', igrd
           print *, dtdx(i), uflx(i+1), uflx(i), rhonu(i), i
           stop
         endif

         wtnu (i) = ( rho(i) * wt (i) -
     &                dtdx(i) * ( utflx (i+1) - utflx (i) ) ) / rhonu(i)

         wttnu(i) = ( rho(i) * wtt(i) -
     &                dtdx(i) * ( uttflx(i+1) - uttflx(i) ) ) / rhonu(i)

         enu  (i) = ( rho(i) * e(i) - dtdx(i)* ( eflx(i+1)- eflx(i)
     &                                                    ) ) / rhonu(i)

      enddo


C$OMP ATOMIC
      ergboun = ergboun + (delx(igrd)**2 * ( eflx(5) -  eflx(np4+1)
     &             + grav(5)*rhoflx(5) - grav(np4)*rhoflx(np4+1) ) )
c                   left rhoflx                right rhoflx

      do m = 1, qc
      do i = 5, np4
         chenu(i,m) = max( smlche(m),
     &          che(i,m)  - dtdx(i) * ( cheflx(i+1,m) - cheflx(i,m))   )
      enddo
      enddo

c     relativistic change:   unu = kinematic velocity  is not the
c     same as  wnu = dynamic velocity = specific momentum !

      do i = 5, np4

         unu  (i) =   wmat(lpi,lpi,i) * wnu  (i)
     &              + wmat(lpj,lpi,i) * wtnu (i)
     &              + wmat(lpk,lpi,i) * wttnu(i)

         utnu (i) =   wmat(lpi,lpj,i) * wnu  (i)
     &              + wmat(lpj,lpj,i) * wtnu (i)
     &              + wmat(lpk,lpj,i) * wttnu(i)

         uttnu(i) =   wmat(lpi,lpk,i) * wnu  (i)
     &              + wmat(lpj,lpk,i) * wtnu (i)
     &              + wmat(lpk,lpk,i) * wttnu(i)

         enu(i) =  enu(i)
     &       - half * ( wnu(i)**2 + wtnu(i)**2 + wttnu(i)**2 )
     &       + half * ( unu(i)**2 + utnu(i)**2 + uttnu(i)**2 )

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine intrfc

c     calculate interface values of all variables

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'ppms.cmn'

      call coeff

      call interp (rhol, rho, rhor)
      call detect (rhol, rho, rhor, smlrho)

      call interp (ul,   u,   ur   )
      call interp (utl,  ut,  utr  )
      call interp (uttl, utt ,uttr )
      call interp (pl,   p,   pr   )
      call interp (gamel,game,gamer)
      call interp (gamcl,gamc,gamcr)
      call interp (gravl,grav,gravr)

      do m = 1, qc
         call interp (chel(1,m), che(1,m), cher(1,m) )
         call detect (chel(1,m), che(1,m), cher(1,m), smlche(m))
      enddo

      call flaten

c      goto 111

      do i = 4, np5
         rhol (i) = flatn(i) * rho (i) + flatn1(i) * rhol (i)
         rhor (i) = flatn(i) * rho (i) + flatn1(i) * rhor (i)
         ul   (i) = flatn(i) * u   (i) + flatn1(i) * ul   (i)
         ur   (i) = flatn(i) * u   (i) + flatn1(i) * ur   (i)
         utl  (i) = flatn(i) * ut  (i) + flatn1(i) * utl  (i)
         utr  (i) = flatn(i) * ut  (i) + flatn1(i) * utr  (i)
         uttl (i) = flatn(i) * utt (i) + flatn1(i) * uttl (i)
         uttr (i) = flatn(i) * utt (i) + flatn1(i) * uttr (i)
         pl   (i) = flatn(i) * p   (i) + flatn1(i) * pl   (i)
         pr   (i) = flatn(i) * p   (i) + flatn1(i) * pr   (i)
         gamel(i) = flatn(i) * game(i) + flatn1(i) * gamel(i)
         gamer(i) = flatn(i) * game(i) + flatn1(i) * gamer(i)
         gamcl(i) = flatn(i) * gamc(i) + flatn1(i) * gamcl(i)
         gamcr(i) = flatn(i) * gamc(i) + flatn1(i) * gamcr(i)
      enddo

      do m = 1, qc
      do i = 4, np5
         chel (i,m) = flatn(i) * che (i,m) + flatn1(i) * chel (i,m)
         cher (i,m) = flatn(i) * che (i,m) + flatn1(i) * cher (i,m)
      enddo
      enddo

  111 continue

      call monot (rhol, rho, rhor, drho, rho6 )
      call monot (ul,   u,   ur,   du,   u6   )
      call monot (utl,  ut,  utr,  dut,  ut6  )
      call monot (uttl, utt, uttr, dutt, utt6 )
      call monot (pl,   p,   pr,   dp,   p6   )
      call monot (gamel,game,gamer,dgame,game6)
      call monot (gamcl,gamc,gamcr,dgamc,gamc6)
cmo   call monot (gravl,grav,gravr,dgrav,grav6)

      do m = 1, qc
         call monot(chel(1,m),che(1,m),cher(1,m),dche(1,m),che6(1,m))
      enddo

c      if ( relaxon .eq. lone ) then
c         call interp (urlxl, urlx, urlxr)
c         do i = 4, np5
c            urlxl(i) = flatn(i) * urlx(i) + flatn1(i) * urlxl(i)
c            urlxr(i) = flatn(i) * urlx(i) + flatn1(i) * urlxr(i)
c         enddo
c         call monot (urlxl,urlx,urlxr,durlx,urlx6)
c         do  i = 4, np5
c            ugrid(i) = half * ( urlxl(i) + urlxr(i-1) )
c         enddo
c      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine coeff

c     calculate coefficients of cubic interpolation polynomial
c     produces common/coeffs/ in ppms.cmn

      include 'qparam.cmn'
      include 'squants.cmn'
      include   'ppms.cmn'

      real*8 scrch1(q), scrch2(q), scrch3(q), scrch4(q)


      do i = 3, np7
         scrch1(i) = dx(i)     + dx(i-1)
         scrch2(i) = scrch1(i) + dx(i)
         scrch3(i) = scrch1(i) + dx(i-1)
      enddo

      do i = 3, np6
         scrch4(i) = dx(i)  /  ( scrch1(i) + dx(i+1) )
         coeff1(i) = scrch4(i) * scrch3(i)   / scrch1(i+1)
         coeff2(i) = scrch4(i) * scrch2(i+1) / scrch1(i)
      enddo

      do i = 3, np5
      scrch4(i) = one  /  ( scrch1(i) + scrch1(i+2) )
      coeff3(i) = -scrch4(i) * dx(i)   * scrch1(i)   / scrch3(i+1)
      coeff4(i) =  scrch4(i) * dx(i+1) * scrch1(i+2) / scrch2(i+1)
      coeff5(i) = dx(i) - two * (dx(i+1)*coeff3(i) + dx(i)*coeff4(i))
      coeff5(i) = coeff5(i) / scrch1(i+1)
      enddo

c --  for equidistant grids with dx(i) = dx follows:
c --  coeff1 = coeff2 = coeff5 = 1/2,  coeff3 = -1/6 = -coeff4

      return
      end
c-----------------------------------------------------------------------
      subroutine interp(al,a,ar)

c     interpolate interface values and add boundary values
c     uses common/coeffs/ in ppms.cmn

      include 'qparam.cmn'
      include 'squants.cmn'
      include   'ppms.cmn'

      real*8 scrch1(q), scrch2(q), scrch3(q), scrch4(q)

      real*8  al(q), a(q), ar(q), delabug(q)

      do i = 3, np7
         scrch1(i) = a(i) - a(i-1)
         scrch2(i) = two * abs(scrch1(i))
         scrch4(i) = sign (one, scrch1(i))
      enddo

      do i = 3, np6
         delabug(i) = coeff1(i)*scrch1(i+1) + coeff2(i)*scrch1(i)
         scrch3(i)  = sign(one, delabug(i))
         dela(i)= min(abs(delabug(i)),scrch2(i),scrch2(i+1)) * scrch3(i)
      enddo

      do i = 3, np6
         if ( scrch4(i) * scrch4(i+1) .le. zero )   dela(i) = zero
      enddo

      do i = 3, np5
cc      ar(i)  = a (i) + coeff5(i) * scrch1(i+1) +
c     instead of this... the following is for equidistant grids only !!
c     from  coeff5 = 1/2 follows
         ar(i)  = half * ( a(i) + a(i+1) ) +
     &             coeff3(i) * dela(i+1)  + coeff4(i) * dela(i)
      enddo

      do i = 4, np5
         al(i) = ar(i-1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine detect(al,a,ar,smalla)

c     search for contact discontinuities in variable a and
c     steepen the zone structure if necessary

      include 'qparam.cmn'
      include 'squants.cmn'
      include   'ppms.cmn'

      real*8 scrch1(q), scrch2(q), scrch3(q), scrch4(q)

      real*8  al(q), a(q), ar(q), smalla


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
         !if(p(i+1).eq.zero) print *, 'p(i+1)', i+1
         !if(p(i-1).eq.zero) print *, 'p(i-1)', i-1
         !if(rho(i+1).eq.zero) print *, 'rho(i+1)', i+1
         !if(rho(i-1).eq.zero) print *, 'rho(i+1)', i-1
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
c-----------------------------------------------------------------------
      subroutine riemann

c     solve riemann shock tube problem

      include 'qparam.cmn'
      include 'squants.cmn'
      include   'ppms.cmn'

      real*8 scrch1(q), scrch2(q), scrch3(q), scrch4(q)

      real*8 pstar1(q), pstar2(q), gmstrl(q), gmstrr(q),
     &     wlft1 (q), wrght1(q), gmin  (q), gmax  (q), aux(q)

      real*8 difu(q)

      integer*4 qi

c------- starting values

      do i = 5, np5
         pstar1(i)=prght(i)-plft(i)-crght(i)*(urght(i)-ulft(i))
         pstar1(i)=plft(i)+pstar1(i)*(clft(i)/(clft(i)+crght(i)))
         pstar1(i)=max(smallp,pstar1(i))

         ge=half*(gmelft(i)+gmergt(i))
         gc=half*(gmclft(i)+gmcrgt(i))
         gamfac=(one-ge/gc)*(ge-one)
         gmstrl(i)=gamfac*(pstar1(i)-plft(i))
         gmstrl(i)=gmelft(i)+two*gmstrl(i)/(pstar1(i)+plft(i))
         gmstrr(i)=gamfac*(pstar1(i)-prght(i))
         gmstrr(i)=gmergt(i)+two*gmstrr(i)/(pstar1(i)+prght(i))

         gmin(i)=min(game(i-1),game(i),game(i+1))
         gmax(i)=max(game(i-1),game(i),game(i+1))
         gmstrl(i)=max(gmin(i),min(gmstrl(i),gmax(i)))
         gmstrr(i)=max(gmin(i),min(gmstrr(i),gmax(i)))

         scrch1(i)=pstar1(i)-(gmstrl(i)-one)*plft(i)/(gmelft(i)-one)
         if ( scrch1(i) .eq. zero ) scrch1(i) = smallp
         wlft1(i)=pstar1(i)+half*(gmstrl(i)-one)*(pstar1(i)+plft(i))
         wlft1(i)=sqrt( abs(
     &       (pstar1(i)-plft(i)) * wlft1(i) / (vlft(i)*scrch1(i)) ) )
         scrch2(i)=pstar1(i)-(gmstrr(i)-one)*prght(i)/(gmergt(i)-one)
         if ( scrch2(i) .eq. zero ) scrch2(i) = smallp
         wrght1(i)=pstar1(i)+half*(gmstrr(i)-one)*(pstar1(i)+prght(i))
         wrght1(i)=sqrt( abs(
     &       (pstar1(i)-prght(i)) * wrght1(i) / (vrght(i)*scrch2(i)) ) )
         if ( abs(pstar1(i)/plft(i)-one)-small .lt. zero )
     &      wlft1(i) = clft(i)
         if ( abs(pstar1(i)/prght(i)-one)-small .lt. zero )
     &      wrght1(i) = crght(i)

         aux(i)=sqrt(half-half/game(i))
         wlft1(i) = max(wlft1(i), aux(i)*clft(i))
         wrght1(i)= max(wrght1(i),aux(i)*crght(i))

         pstar2(i)=prght(i)-plft(i)-wrght1(i)*(urght(i)-ulft(i))
         pstar2(i)=plft(i)+pstar2(i)*(wlft1(i)/(wlft1(i)+wrght1(i)))
         pstar2(i)=max(smallp,pstar2(i))

      enddo

c------- iteration
      n=0
   39 continue
      n=n+1
c      do 40 n = 1,nriem

      do 20 i = 5, np5
         ge=half*(gmelft(i)+gmergt(i))
         gc=half*(gmclft(i)+gmcrgt(i))
         gamfac=(one-ge/gc)*(ge-one)
         gmstrl(i)=gamfac*(pstar2(i)-plft(i))
         gmstrl(i)=gmelft(i)+two*gmstrl(i)/(pstar2(i)+plft(i))
         gmstrr(i)=gamfac*(pstar2(i)-prght(i))
         gmstrr(i)=gmergt(i)+two*gmstrr(i)/(pstar2(i)+prght(i))
         gmstrl(i)=max(gmin(i),min(gmax(i),gmstrl(i)))
         gmstrr(i)=max(gmin(i),min(gmax(i),gmstrr(i)))
 20   continue

       do 25 i = 5, np5
         scrch1(i)=pstar2(i)-(gmstrl(i)-one)*plft(i)/(gmelft(i)-one)
         if ( scrch1(i) .eq. zero ) scrch1(i) = smallp
c         scrch1(i)=cvmgz(smallp,scrch1(i),scrch1(i))
         wlft(i)=pstar2(i)+half*(gmstrl(i)-one)*(pstar2(i)+plft(i))
         wlft(i)   = sqrt( abs(
     &        (pstar2(i)-plft(i)) * wlft(i) / (vlft(i)*scrch1(i)) ) )
         scrch2(i)=pstar2(i)-(gmstrr(i)-one)*prght(i)/(gmergt(i)-one)
         if ( scrch2(i) .eq. zero ) scrch2(i) = smallp
c         scrch2(i)=cvmgz(smallp,scrch2(i),scrch2(i))
         wrght(i)=pstar2(i)+half*(gmstrr(i)-one)*(pstar2(i)+prght(i))
         wrght(i) = sqrt( abs(
     &        (pstar2(i)-prght(i)) * wrght(i) / (vrght(i)*scrch2(i)) ) )
         if ( abs(pstar2(i)/plft(i)-one)-small .lt. zero )
     &       wlft(i) = clft(i)
c         wlft(i)=cvmgp(wlft(i),clft(i),abs(pstar2(i)/plft(i)-one)-small)
         aux(i)=sqrt(half*(game(i)-one)/game(i))
         wlft(i)=max(wlft(i),aux(i)*clft(i))
         if ( abs(pstar2(i)/prght(i)-one)-small .lt. zero )
     &       wrght(i) = crght(i)
c         wrght(i)=cvmgp(wrght(i),crght(i),
c     &                   abs(pstar2(i)/prght(i)-one)-small)
         wrght(i)=max(wrght(i),aux(i)*crght(i))
 25   continue

      do 30 i = 5, np5
         ustrl1=ulft(i)-(pstar1(i)-plft(i))/wlft1(i)
         ustrr1=urght(i)+(pstar1(i)-prght(i))/wrght1(i)
         ustrl2=ulft(i)-(pstar2(i)-plft(i))/wlft(i)
         ustrr2=urght(i)+(pstar2(i)-prght(i))/wrght(i)
         delu1=ustrl1-ustrr1
         delu2=ustrl2-ustrr2
         scrch1(i)=delu2-delu1
         if ( scrch1(i) .eq. zero ) delu2 = zero
         if ( scrch1(i) .eq. zero ) scrch1(i) = smallu
c         scrch1(i)=cvmgz(smallu,scrch1(i),scrch1(i))
         pstar(i)=pstar2(i)-delu2*(pstar2(i)-pstar1(i))/scrch1(i)
         pstar(i)=max(smallp,pstar(i))

         difu(i) = abs((pstar2(i)-pstar(i))/pstar(i))
         pstar1(i)=pstar2(i)
         pstar2(i)=pstar(i)
         wlft1(i)=wlft(i)
         wrght1(i)=wrght(i)
  30  continue

      wmax = zero
      do i = 5, np5
	wmax = max ( wmax, difu(i) )
      enddo
      if ( ( wmax .gt. 1.D-5 ) .and. ( n .lt. nriem ) ) goto 39
c  40  continue
cwr   write (6,*) 'riemann iterations:',n
c------- end iteration

      do i = 5, np5
         scrch3(i) = ulft (i) - ( pstar(i) - plft (i) ) / wlft (i)
         scrch4(i) = urght(i) + ( pstar(i) - prght(i) ) / wrght(i)
         ustar (i) = half * ( scrch3(i) + scrch4(i) )

         urel(i)   = ustar(i) - ugrid(i)
         scrch1(i) = sign( one, urel(i) )
      enddo

      do 70 i = 5, np5
	 if ( scrch1(i) .ge. zero ) then
           ps   (i) = plft  (i)
           us   (i) = ulft  (i)
           uts  (i) = utlft (i)
           utts (i) = uttlft(i)
           vs   (i) = vlft  (i)
           games(i) = gmelft(i)
           gamcs(i) = gmclft(i)
	 else
           ps   (i) = prght  (i)
           us   (i) = urght  (i)
           uts  (i) = utrght (i)
           utts (i) = uttrgt(i)
           vs   (i) = vrght  (i)
           games(i) = gmergt(i)
           gamcs(i) = gmcrgt(i)
         endif
c         ps   (i) = CVMGP( plft  (i), prght (i), scrch1(i) )
c         us   (i) = CVMGP( ulft  (i), urght (i), scrch1(i) )
c         uts  (i) = CVMGP( utlft (i), utrght(i), scrch1(i) )
c         utts (i) = CVMGP( uttlft(i), uttrgt(i), scrch1(i) )
c         vs   (i) = CVMGP( vlft  (i), vrght (i), scrch1(i) )
c         games(i) = CVMGP( gmelft(i), gmergt(i), scrch1(i) )
c         gamcs(i) = CVMGP( gmclft(i), gmcrgt(i), scrch1(i) )
         rhos (i) = max( smlrho, one/vs(i) )
         vs   (i) = one / rhos(i)
	 if ( scrch1(i) .ge. zero ) then
           ws(i) = wlft(i)
         else
           ws(i) = wrght(i)
         endif
c         ws   (i) = CVMGP( wlft  (i), wrght (i), scrch1(i) )
         ces  (i) = sqrt( gamcs(i) * ps(i) * vs(i) )

         scrch2(i) = pstar(i) - ps(i)
         vstar(i)  = vs(i)-scrch2(i)/ws(i)/ws(i)
         rhostr(i) = one / vstar(i)
         rhostr(i) = max (smlrho, rhostr(i))
         vstar(i)  = one / rhostr(i)
         cestar(i) = sqrt(gamcs(i)*pstar(i)*vstar(i))

         wes   (i) = ces   (i) - scrch1(i)*us(i)
         westar(i) = cestar(i) - scrch1(i)*ustar(i)
c wrong  scrch4(i) = ws(i)*vstar(i)-scrch1(i)*ustar(i)
         scrch4(i) = ws(i)*vs(i)-scrch1(i)*us(i)
         if ( scrch2(i) .ge. zero ) then
	    wes(i)    = scrch4(i)
            westar(i) = scrch4(i)
         endif
c         wes   (i) = cvmgm(wes   (i),scrch4(i),scrch2(i))
c         westar(i) = cvmgm(westar(i),scrch4(i),scrch2(i))
         wes   (i) = wes   (i) + scrch1(i)*ugrid(i)
         westar(i) = westar(i) + scrch1(i)*ugrid(i)
  70  continue

      do 72 i = 5, np5
         gamfac = (one-games(i)/gamcs(i))*(games(i)-one)
         gmstar(i) = gamfac*(pstar(i)-ps(i))
         gmstar(i) = games(i)+two*gmstar(i)/(pstar(i)+ps(i))
         gmstar(i)=max(gmin(i),min(gmax(i),gmstar(i)))
 72   continue

c------- compute correct state (lin.interpol.) for rarefaction fan

      do i = 5, np5
         scrch3(i) = max( wes(i)-westar(i), wes(i)+westar(i), smallu )
         scrch3(i) = ( wes(i) + westar(i) )  /  scrch3(i)
         scrch3(i) = half * ( one + scrch3(i) )
         scrch2(i) =          one - scrch3(i)

         rhoav(i) = scrch3(i) * rhostr(i) + scrch2(i) * rhos(i)
         uav  (i) = scrch3(i) * ustar (i) + scrch2(i) * us  (i)
         utav (i) = uts(i)
         uttav(i) = utts(i)
         pav  (i) = scrch3(i) * pstar (i) + scrch2(i) * ps  (i)
         gameav(i)= scrch3(i) * gmstar(i) + scrch2(i) * games(i)

c------- check for supersonic flow

	 if ( westar(i) .gt. zero ) then
            rhoav(i) = rhostr(i)
            uav(i)   = ustar (i)
            pav(i)   = pstar (i)
            gameav(i)= gmstar(i)
         endif

	 if ( wes(i) .lt. zero ) then
            rhoav(i) = rhos(i)
            uav(i)   = us  (i)
            pav(i)   = ps  (i)
            gameav(i)= games(i)
         endif

c         rhoav(i) = CVMGM(rhoav (i),rhostr(i),westar(i))
c         uav(i)   = CVMGM(uav   (i),ustar (i),westar(i))
c         pav(i)   = CVMGM(pav   (i),pstar (i),westar(i))
c         gameav(i)= CVMGM(gameav(i),gmstar(i),westar(i))

c         rhoav(i) = CVMGP(rhoav (i), rhos(i), wes(i))
c         uav(i)   = CVMGP(uav   (i), us  (i), wes(i))
c         pav(i)   = CVMGP(pav   (i), ps  (i), wes(i))
c         gameav(i)= CVMGP(gameav(i),games(i), wes(i))

         urel(i)  = uav(i) - ugrid(i)
         ce(i) = sqrt(gameav(i) * pav(i) / rhoav(i) )
      enddo

c --- chemical composition
c BUGFIX T.PLEWA
c     -------------------------
c     upwind state is known now

      do i = 5, np5
         u_side    = sign (one, uav(i) - ugrid(i) )
         scrch2(i) = half * ( one + u_side )
         scrch3(i) = half * ( one - u_side )
      end do

      do m = 1, qc
      do i = 5, np5
         cheav(i,m) = chelft(i,m) * scrch2(i) + chergt(i,m) * scrch3(i)
      end do
      end do

c ---

c     do m = 1, qc
c     do i = 5, np5
c	 if ( scrch1(i) .ge. zero ) then
c          ches(i,m) = chelft(i,m)
c        else
c          ches(i,m) = chergt(i,m)
c        endif
cc        ches  (i,m) = CVMGP( chelft(i,m), chergt(i,m), scrch1(i) )
c        chestr(i,m) = ches(i,m) * (rhostr(i)/rhos(i))     ! guess
c        cheav (i,m) = scrch3(i) * chestr(i,m) + scrch2(i) * ches(i,m)
c	 if ( westar(i) .gt. zero ) cheav(i,m) = chestr(i,m)
cc        cheav (i,m) = CVMGM(cheav(i,m), chestr(i,m), westar(i))
c	 if ( wes(i) .lt. zero )    cheav(i,m) = ches(i,m)
cc        cheav (i,m) = CVMGP(cheav(i,m), ches  (i,m), wes   (i))
c     enddo
c     enddo

c     ---------
c     CMA start
c     ---------
c     corrections for consistency of sum over all mass fractions = 1
c     mass fractions are only in che(m=3..qc)

      if ( qc .ge. 3 ) then
         qi = qc

         do i = 5, np5
            scrch1(i) = zero
         end do

         do m = 3, qi
         do i = 5, np5
            scrch1(i) = scrch1(i) + cheav(i,m)
         end do
         end do

         do m = 3, qi
         do i = 5, np5
            cheav(i,m) = cheav(i,m) / scrch1(i)
         end do
         end do

      end if


c     -------
c     CMA end
c     -------

      return
      end
c-----------------------------------------------------------------------
      subroutine accel(igrd)
c     modify velocity and energy to account for
c     gravitational acceleration

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include 'compct.cmn'

      real*8 scrch1(q), scrch2(q), scrch3(q), scrch4(q)
      real*8 addegq,addphix,addphiy,addphiz, cfudge, ekinmax

      if (nsdim.ne.3) then
         write (6,*) 'nsdim must be 3 in this version of accel.'
        stop
      endif

c This is for the tstep-calculations
c velmax is set to zero in tstep

      vlmx = velmax(igrd)
C$OMP PARALLEL DO DEFAULT(NONE), REDUCTION(max:vlmx),PRIVATE(j,k,i),
C$OMP+    SHARED(nx,ny,nz,igrd,velx,vely,velz,ugridx,ugridy,ugridz)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         vlmx = max( vlmx,
     &               abs( velx(i,j,k,igrd) - ugridx(i,igrd) ) +
     &               abs( vely(i,j,k,igrd) - ugridy(j,igrd) ) +
     &               abs( velz(i,j,k,igrd) - ugridz(k,igrd) )   )
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
      velmax(igrd) = vlmx

c --- spread black hole mass for potentials

      call sucks(igrd,lone)

c -----
c first Newtonian gravitational potential
c -----
      call gravty(igrd)

      if ( mconoff .eq. lone )   call trav2w(-one,igrd)

c      iwdgr = indxgr(ltwo,maxlev(ltwo))
      iwdgr = indxgr(lone,maxlev(lone))
      itopgr= indxgr(lone,lone)
      rezd = one / delx(iwdgr)
      rkf  = one / (dble(nfine)**(levlgr(igrd,lone)-lone))
      rkf3 = rkf**3
c      topx2 = topos(posx2,lone,iwdgr)
c      topy2 = topos(posy2,ltwo,iwdgr)
c      topz2 = half * rkf
c      sop2 = ( soparp / (dble(nfine)**(maxlev(2)-levlgr(itopgr,1))) )**2
c      xopos= topos( zero, lone, igrd)
c      yopos= topos( zero, ltwo, igrd)

      ronoff = dble(mconoff)

c      print *, 'accel vel sample before:'
c      print *, 'grid', igrd
c      print *, velx(25,25,1,igrd),vely(25,25,1,igrd),velz(25,25,1,igrd)
c      print *, ''

c------ x-sweep --------------------

      nnn = nx
      np1 = nnn+1
      np2 = nnn+2
      np3 = nnn+3
      np4 = nnn+4
      np5 = nnn+5
      np6 = nnn+6
      np7 = nnn+7
      np8 = nnn+8


      call boundry(lzero,igrd,lone,ltwo)

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+         PRIVATE(jj,kk,ii,i4,i,dx,ytmp,xtmp,rtmpsq),
C$OMP+         SHARED(ny,nz,np4,nnn,igrd,ronoff,cc,corrmax,bxgpot,
C$OMP+                cmgax,cmgay,dtppm,velx,vely,energy,gpot,delx)
c#ifdef SGI
cC$OMP+        ,PRIVATE(unu,utnu,enu,grav,gravl,gravr,dphidx)
c#endif
      do jj = 1, ny
      do kk = 1, nz

       do ii = 1, nnn
          i4 = ii + 4
c          rho  (i4) = dold  (ii,jj,kk,igrd)
c          rhonu(i4) = densty(ii,jj,kk,igrd)
c          u    (i4) = vxold (ii,jj,kk,igrd)
          unu  (i4) = velx  (ii,jj,kk,igrd)
          utnu (i4) = vely  (ii,jj,kk,igrd)
          enu  (i4) = energy(ii,jj,kk,igrd)
          grav (i4) = gpot  (ii,jj,kk,igrd)
          dx   (i4) = delx(igrd)
       enddo

       call bndrwx( grav, bxgpot, jj, kk )
       call interp( gravl, grav, gravr)

       do i = 5, np4
          dphidx(i) = ( gravr(i) - gravl(i) )  /  dx(i)
     &      * ( one + ronoff * min( two * enu(i) / cc**2, corrmax ) )
       enddo

c       if ( relaxon .eq. lone ) then          ! centrifugal forces
c          ytmp = dble(jj) - cmgay
c          do i = 5, np4
c             xtmp = dble(i-4) - cmgax
c             rtmpsq = xtmp**2 + ytmp**2
c             dphidx(i) =dphidx(i) - xtmp/rtmpsq*((unu(i)**2+utnu(i)**2)
c     &           -(unu(i)*xtmp+utnu(i)*ytmp)**2/rtmpsq ) / dx(i)
c          enddo
c       endif

       do i = 5, np4
          unu(i) = unu(i) - dtppm * dphidx(i)
          enu(i) = enu(i) - dtppm * unu(i) * dphidx(i)
       enddo

       do ii = 1, nnn
          i4 = ii + 4
          velx  (ii,jj,kk,igrd) = unu (i4)
          energy(ii,jj,kk,igrd) = enu (i4)
       enddo

      enddo
      enddo
C$OMP END PARALLEL DO

c      print *, 'accel vel sample after x-sweep:'
c      print *, 'grid', igrd
c      print *, velx(25,25,1,igrd),vely(25,25,1,igrd),velz(25,25,1,igrd)
c      print *, ''

c------- y-sweep -------------------

      if (nsdim .gt. lone)  then

      nnn = ny
      np1 = nnn+1
      np2 = nnn+2
      np3 = nnn+3
      np4 = nnn+4
      np5 = nnn+5
      np6 = nnn+6
      np7 = nnn+7
      np8 = nnn+8

      call boundry(lzero,igrd,ltwo,ltwo)

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+         PRIVATE(ii,kk,jj,j4,j,dx,xtmp,ytmp,rtmpsq),
C$OMP+         SHARED(nx,nz,np4,nnn,igrd,ronoff,cc,corrmax,bygpot,
C$OMP+                cmgax,cmgay,dtppm,vely,velx,energy,gpot,delx)
c#ifdef SGI
cC$OMP+        ,PRIVATE(unu,utnu,enu,grav,gravl,gravr,dphidx)
c#endif
      do ii = 1, nx
      do kk = 1, nz

       do jj = 1, nnn
          j4 = jj + 4
c          rho  (j4) = dold  (ii,jj,kk,igrd)
c          rhonu(j4) = densty(ii,jj,kk,igrd)
c          u    (j4) = vyold (ii,jj,kk,igrd)
          unu  (j4) = vely  (ii,jj,kk,igrd)
          utnu (j4) = velx  (ii,jj,kk,igrd)
          enu  (j4) = energy(ii,jj,kk,igrd)
          grav (j4) = gpot  (ii,jj,kk,igrd)
          dx   (j4) = delx(igrd)
       enddo

       call bndrwy( grav, bygpot, ii, kk )
       call interp( gravl, grav, gravr)

       do j = 5, np4
          dphidx(j) = ( gravr(j) - gravl(j) )  /  dx(j)
     &      * ( one + ronoff * min( two * enu(j) / cc**2, corrmax ) )
       enddo

c       if ( relaxon .eq. lone ) then          ! centrifugal forces
c          xtmp = dble(ii) - cmgax
c          do j = 5, np4
c             ytmp = dble(j-4) - cmgay
c             rtmpsq = xtmp**2 + ytmp**2
c             if ( ( densty (ii,j-4,kk,igrd) .gt. 3.D14 ) .and.
c     &         ( igrd .eq. 3) .and. (vely(ii,j-4,kk,igrd).gt.1.D9)) then
c            write(*,*) dphidx(j), ytmp/rtmpsq*((unu(j)**2+utnu(j)**2)
c     &           -(unu(j)*ytmp+utnu(j)*xtmp)**2/rtmpsq ) / dx(j),
c     &             unu(j),utnu(j),grav(j)
c             endif
c             dphidx(j) =dphidx(j) - ytmp/rtmpsq*((unu(j)**2+utnu(j)**2)
c     &           -(unu(j)*ytmp+utnu(j)*xtmp)**2/rtmpsq ) / dx(j)
c          enddo
c       endif

       do j = 5, np4
          unu(j) = unu(j) - dtppm * dphidx(j)
          enu(j) = enu(j) - dtppm * unu(j) * dphidx(j)
       enddo

       do jj = 1, nnn
          j4 = jj + 4
          vely  (ii,jj,kk,igrd) = unu (j4)
          energy(ii,jj,kk,igrd) = enu (j4)
       enddo

      enddo
      enddo
C$OMP END PARALLEL DO

      endif

c------- z-sweep -------------------

      if (nsdim .ge. lthree)  then

      nnn = nz
      np1 = nnn+1
      np2 = nnn+2
      np3 = nnn+3
      np4 = nnn+4
      np5 = nnn+5
      np6 = nnn+6
      np7 = nnn+7
      np8 = nnn+8

      call boundry(lzero,igrd,lthree,ltwo)

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+         PRIVATE(jj,ii,kk,k4,k,dx),
C$OMP+         SHARED(ny,nx,np4,nnn,igrd,ronoff,cc,corrmax,bzgpot,
C$OMP+                dtppm,velz,energy,gpot,delx)
c#ifdef SGI
cC$OMP+        ,PRIVATE(unu,utnu,enu,grav,gravl,gravr,dphidx)
c#endif
      do jj = 1, ny
      do ii = 1, nx

       do kk = 1, nnn
          k4 = kk + 4
c          rho  (k4) = dold  (ii,jj,kk,igrd)
c          rhonu(k4) = densty(ii,jj,kk,igrd)
c          u    (k4) = vzold (ii,jj,kk,igrd)
          unu  (k4) = velz  (ii,jj,kk,igrd)
          enu  (k4) = energy(ii,jj,kk,igrd)
          grav (k4) = gpot  (ii,jj,kk,igrd)
          dx   (k4) = delx(igrd)
       enddo

       call bndrwz( grav, bzgpot, ii, jj )
       call interp( gravl, grav, gravr)

       do k = 5, np4
          dphidx(k) = ( gravr(k) - gravl(k) )  /  dx(k)
     &      * ( one + ronoff * min( two * enu(k) / cc**2, corrmax ) )
       enddo

       do k = 5, np4
          unu(k) = unu(k) - dtppm * dphidx(k)
          enu(k) = enu(k) - dtppm * unu(k) * dphidx(k)
       enddo

       do kk = 1, nnn
          k4 = kk + 4
          velz  (ii,jj,kk,igrd) = unu (k4)
          if( velz(ii,jj,kk,igrd).gt.half*cc) then
            print *, 'velz mad crazy ppm'
            print *, ii,jj,kk,igrd
          endif
          energy(ii,jj,kk,igrd) = enu (k4)
       enddo

      enddo
      enddo
C$OMP END PARALLEL DO

      end if

c      print *, 'accel vel sample after z-sweep:'
c      print *, 'grid', igrd
c      print *, velx(25,25,1,igrd),vely(25,25,1,igrd),velz(25,25,1,igrd)
c      print *, ''

      if ( mconoff .eq. lone )   call trav2w(+one,igrd)
c -------- more general relativity: transform back
c          this is opposite of call trav2w further above

c ===============================================

c -----
c Take care of the gravitational wave corrections
c -----

      call gwaves  ( igrd )
      call rgpot   ( igrd )
      call reactphi( igrd )
      call strain  ( igrd )

      if     ( mconoff .eq. lzero )   then
         call addpress( igrd )
      elseif ( mconoff .eq. lone )   then
         call chigravty( igrd )      !   gchi is set to zero in rstzero
         call trav2w(-one,igrd)
         call addchi( igrd )
      else
         write(*,*)  'accel, something wrong with mconoff ', mconoff
      endif

c ===============================================

      iwdgr = indxgr(lone,maxlev(lone))
      itopgr= indxgr(lone,lone)
      d3x = delx(igrd)**3
      rezd = one / delx(iwdgr)
      rkf  = one / (dble(nfine)**(levlgr(igrd,lone)-lone))
      rkf3 = rkf**3

c -----

      ifgr = idxfgr(igrd)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i), SHARED(nx,ny,nz,d3x,po)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         po(i,j,k) = -two * d3x
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      if ( ifgr .ne. -lone )  then
         noz = norgin(ifgr,lthree)
         noy = norgin(ifgr,ltwo)
         nox = norgin(ifgr,lone)
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i),
C$OMP+            SHARED(nox,noy,noz,nx,ny,nz,po)
         do k = noz+1, noz+nz/nfine
         do j = noy+1, noy+ny/nfine
         do i = nox+1, nox+nx/nfine
            po(i,j,k) = zero
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO
      endif

c------ x-sweep --------------------

      nnn = nx
      np1 = nnn+1
      np2 = nnn+2
      np3 = nnn+3
      np4 = nnn+4
      np5 = nnn+5
      np6 = nnn+6
      np7 = nnn+7
      np8 = nnn+8

      call boundry(ltwo,igrd,lone,ltwo)

C$OMP PARALLEL DO DEFAULT(NONE),REDUCTION(-:egquel),
C$OMP+            PRIVATE(jj,kk,ii,i4,i,dx,dunu)
C$OMP+            SHARED(delx,ronoff,dtppm,igrd,np4,nnn,nz,ny,
C$OMP+                   densty,velx,vely,velz,energy,gpot,po,press,
C$OMP+                   phireac,gchi,bxvelx,corrmax,cc)
c#ifdef SGI
cC$OMP+        ,PRIVATE(dphidx,rhonu,enu,unu,grav,gravr,gravl,p,pl,pr)
c#endif
      do jj = 1, ny
      do kk = 1, nz

       do ii = 1, nnn
          i4 = ii + 4
c          rho  (i4) = dold  (ii,jj,kk,igrd)
          rhonu(i4) = densty(ii,jj,kk,igrd)
c          u    (i4) = vxold (ii,jj,kk,igrd)
          unu  (i4) = velx  (ii,jj,kk,igrd)
          enu  (i4) = energy(ii,jj,kk,igrd)
          grav (i4) = phireac(ii,jj,kk) + gchi(ii,jj,kk,igrd)
          dx   (i4) = delx(igrd)
          p    (i4) = press(ii,jj,kk,igrd)
     &      * sign( min( abs(two*gpot(ii,jj,kk,igrd)/cc**2), corrmax ),
     &                         gpot(ii,jj,kk,igrd)  )
       enddo

       call bndrwx( grav, bxvelx, jj, kk )
       call interp(gravl, grav, gravr)

       do i = 1, 4
          p(i)     = p(5)
          p(np4+i) = p(np4)
       enddo

       call interp(pl, p, pr)

       do i = 5, np4
          dphidx(i) = ( gravr(i) - gravl(i) )  /  dx(i)
       enddo

       do i = 5, np4
          egquel= egquel - rhonu(i)*dtppm*unu(i)*dphidx(i)*po(i-4,jj,kk)
       enddo

       do i = 5, np4
c          unu(i) = unu(i) - dtppm *
c     &       ( dphidx(i)  +  ronoff * (pr(i)-pl(i)) / dx(i) / rhonu(i) )
          dunu   = dtppm * ronoff * (pr(i)-pl(i)) / dx(i) / rhonu(i)
          unu(i) = unu(i) - dtppm * dphidx(i)
     &                   - sign( min( half*corrmax*cc, abs(dunu)), dunu)
          enu(i) = enu(i) - dtppm * unu(i) * dphidx(i)
c1     &          - unu(i) * sign( min( half*corrmax*cc, abs(dunu) ), dunu)
c1  line above instead of adderg
       enddo

       do ii = 1, nnn
          i4 = ii + 4
          velx  (ii,jj,kk,igrd) = unu (i4)
          energy(ii,jj,kk,igrd) = enu (i4)
       enddo

      enddo
      enddo
C$OMP END PARALLEL DO

c------- y-sweep -------------------

      if (nsdim .gt. lone)  then

      nnn = ny
      np1 = nnn+1
      np2 = nnn+2
      np3 = nnn+3
      np4 = nnn+4
      np5 = nnn+5
      np6 = nnn+6
      np7 = nnn+7
      np8 = nnn+8

      call boundry(ltwo,igrd,ltwo,ltwo)

C$OMP PARALLEL DO DEFAULT(NONE),REDUCTION(-:egquel),
C$OMP+            PRIVATE(ii,kk,jj,j4,j,i,dx,dunu),
C$OMP+            SHARED(delx,ronoff,dtppm,igrd,np4,nnn,nx,nz,
C$OMP+                   densty,velx,vely,velz,energy,gpot,po,press,
C$OMP+                   phireac,gchi,byvely,corrmax,cc)
c#ifdef SGI
cC$OMP+        ,PRIVATE(dphidx,rhonu,enu,unu,grav,gravr,gravl,p,pl,pr)
c#endif
      do ii = 1, nx
      do kk = 1, nz

       do jj = 1, nnn
          j4 = jj + 4
c          rho  (j4) = dold  (ii,jj,kk,igrd)
          rhonu(j4) = densty(ii,jj,kk,igrd)
c          u    (j4) = vyold (ii,jj,kk,igrd)
          unu  (j4) = vely  (ii,jj,kk,igrd)
          enu  (j4) = energy(ii,jj,kk,igrd)
          grav (j4) = phireac(ii,jj,kk) + gchi(ii,jj,kk,igrd)
          dx   (j4) = delx(igrd)
          p    (j4) = press(ii,jj,kk,igrd)
     &       * sign( min( abs(two*gpot(ii,jj,kk,igrd)/cc**2), corrmax ),
     &                         gpot(ii,jj,kk,igrd)  )
       enddo

       call bndrwy( grav, byvely, ii, kk )
       call interp( gravl, grav, gravr)

       do j = 1, 4
          p(j)     = p(5)
          p(np4+j) = p(np4)
       enddo

       call interp(pl, p, pr)

       do i = 5, np4
          dphidx(i) = ( gravr(i) - gravl(i) )  /  dx(i)
       enddo

       do i = 5, np4
          egquel= egquel - rhonu(i)*dtppm*unu(i)*dphidx(i)*po(ii,i-4,kk)
       enddo

       do i = 5, np4
c          unu(i) = unu(i) - dtppm *
c     &       ( dphidx(i)  +  ronoff * (pr(i)-pl(i)) / dx(i) / rhonu(i) )
          dunu   = dtppm * ronoff * (pr(i)-pl(i)) / dx(i) / rhonu(i)
          unu(i) = unu(i) - dtppm * dphidx(i)
     &                   - sign( min( half*corrmax*cc, abs(dunu)), dunu)
          enu(i) = enu(i) - dtppm * unu(i) * dphidx(i)
c1     &           - unu(i) * sign( min( half*corrmax*cc, abs(dunu) ), dunu)
c1  line above instead of adderg
       enddo

       do jj = 1, nnn
          j4 = jj + 4
          vely  (ii,jj,kk,igrd) = unu (j4)
          energy(ii,jj,kk,igrd) = enu (j4)
       enddo

      enddo
      enddo
C$OMP END PARALLEL DO

      end if

c------- z-sweep -------------------

      if (nsdim .ge. lthree)  then

      nnn = nz
      np1 = nnn+1
      np2 = nnn+2
      np3 = nnn+3
      np4 = nnn+4
      np5 = nnn+5
      np6 = nnn+6
      np7 = nnn+7
      np8 = nnn+8

      call boundry(ltwo,igrd,lthree,ltwo)

C$OMP PARALLEL DO DEFAULT(NONE),REDUCTION(-:egquel),
C$OMP+            PRIVATE(jj,ii,kk,k4,k,i,dx,dunu),
C$OMP+            SHARED(delx,ronoff,dtppm,igrd,np4,nnn,ny,nx,
C$OMP+                   densty,velx,vely,velz,energy,gpot,po,press,
C$OMP+                   phireac,gchi,bzvelz,corrmax,cc)
c#ifdef SGI
cC$OMP+        ,PRIVATE(dphidx,rhonu,enu,unu,grav,gravr,gravl,p,pl,pr)
c#endif
      do jj = 1, ny
      do ii = 1, nx

       do kk = 1, nnn
          k4 = kk + 4
c          rho  (k4) = dold  (ii,jj,kk,igrd)
          rhonu(k4) = densty(ii,jj,kk,igrd)
c          u    (k4) = vzold (ii,jj,kk,igrd)
          unu  (k4) = velz  (ii,jj,kk,igrd)
          enu  (k4) = energy(ii,jj,kk,igrd)
          grav (k4) = phireac(ii,jj,kk) + gchi(ii,jj,kk,igrd)
          dx   (k4) = delx(igrd)
          p    (k4) = press(ii,jj,kk,igrd)
     &       * sign( min( abs(two*gpot(ii,jj,kk,igrd)/cc**2), corrmax ),
     &                         gpot(ii,jj,kk,igrd)  )
       enddo

       call bndrwz( grav, bzvelz, ii, jj )
       call interp( gravl, grav, gravr)

       do k = 1, 4
          p(k)     = p(9-k) ! z reflection
          p(np4+k) = p(np4)
       enddo

       call interp(pl, p, pr)

       do i = 5, np4
          dphidx(i) = ( gravr(i) - gravl(i) )  /  dx(i)
       enddo

       do i = 5, np4
          egquel= egquel - rhonu(i)*dtppm*unu(i)*dphidx(i)*po(ii,jj,i-4)
       enddo

       do i = 5, np4
c          unu(i) = unu(i) - dtppm *
c     &       ( dphidx(i)  +  ronoff * (pr(i)-pl(i)) / dx(i) / rhonu(i) )
          dunu   = dtppm * ronoff * (pr(i)-pl(i)) / dx(i) / rhonu(i)
          unu(i) = unu(i) - dtppm * dphidx(i)
     &                   - sign( min( half*corrmax*cc, abs(dunu)), dunu)
          enu(i) = enu(i) - dtppm * unu(i) * dphidx(i)
c1     &          - unu(i) * sign( min( half*corrmax*cc, abs(dunu) ), dunu)
c1  line above instead of adderg
       enddo

       do kk = 1, nnn
          k4 = kk + 4
          velz  (ii,jj,kk,igrd) = unu (k4)
          if( velz(ii,jj,kk,igrd).gt.half*cc) then
            print *, 'velz mad crazy second sweep'
            print *, ii,jj,kk,igrd
          endif
          energy(ii,jj,kk,igrd) = enu (k4)
       enddo

      enddo
      enddo
C$OMP END PARALLEL DO

      end if

c --- add back reaction pressure source term

C$OMP PARALLEL DO DEFAULT(NONE), REDUCTION(+:egquel),
C$OMP+            PRIVATE(jj,kk,ii,derg),
C$OMP+            SHARED(ny,nz,nx,igrd,dtppm,cc,corrmax,
C$OMP+                   adderg,densty,energy,po)

      do jj = 1, ny
      do kk = 1, nz
      do ii = 1, nx
         derg = dtppm * adderg(ii,jj,kk) / densty(ii,jj,kk,igrd)
         derg = sign( min( abs(derg),
     &       max( (half*corrmax*cc)**2, energy(ii,jj,kk,igrd) ) ), derg)
         energy(ii,jj,kk,igrd) = energy(ii,jj,kk,igrd) + derg
         egquel = egquel + derg*densty(ii,jj,kk,igrd) * po(ii,jj,kk)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      if ( mconoff .eq. lone )   call trav2w(+one,igrd)
c -------- more general relativity: transform back
c     this is opposite of call trav2w further above

      call kickbh(igrd,lone)
      if ( npawi .eq. 1 ) call PaWi(igrd)

c ---

      cfudge = ten
      !ekinmax = half * (cfudge*cc)**2                                  !DRW
      ekinmax = half * cc**2                                            !DRW - Edited

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i,j4,k4,i4,ip,im),
C$OMP+         SHARED(nx,ny,nz,igrd,dtppm,cc,vfac,densty,chem,dtgr,
C$OMP+                velx,vely,velz,energy,press,temper,bhdens,small,
C$OMP+                ekinmax,trelax) ! ,relaxvx,relaxvy,relaxvz)
c#ifdef SGI
cC$OMP+        ,PRIVATE(rho,enu,tmp,ek,ei,p,gamc,game,ye),
c#endif
      do j = 1, ny
      do k = 1, nz

         do i = 1, nx
            ek (i) = half*(velx(i,j,k,igrd)**2
     &             + vely(i,j,k,igrd)**2 + velz(i,j,k,igrd)**2)
            enu(i) = energy(i,j,k,igrd)
         enddo

         do i = 1, nx

            if ( ( ek(i) .gt. half*ekinmax )  .and.
     &           ( bhdens(i,j,k,igrd) .eq. zero )   )  then
               write(*,'(x,a,4i3,a,3F11.6,A,F11.6,A,F11.6)') 'accel,  ',
     &              i,j,k,igrd, '  velx,y,z: ',  velx(i,j,k,igrd)/cc,
     &              vely(i,j,k,igrd)/cc, velz(i,j,k,igrd)/cc,
     &              '  rho:',  log10(densty(i,j,k,igrd)),'  Ye:',
     &              chem(i,j,k,igrd,1) / densty(i,j,k,igrd)
            endif

            if ( ( ek(i) .gt. ekinmax ) .and.
     &           ( bhdens(i,j,k,igrd) .eq. zero )   )  then
               write(*,'(x,a,4i3,a,3F11.6)') 'accel, panic!,  ',
     &              i,j,k,igrd, '  velx,y,z: ',  velx(i,j,k,igrd)/cc,
     &              vely(i,j,k,igrd)/cc, velz(i,j,k,igrd)/cc
               call exit(1)
            endif

            im = max( i-lone, ltwo-i )
            ip = min( i+lone, ltwo*nx-i)
            if ( enu(i) .lt. (one+small)*ek(i) ) then
               enu(i) = half*( enu(im)-ek(im) + enu(ip)-ek(ip) )
               energy(i,j,k,igrd) = max( enu(i), small*ek(i) ) + ek(i)
            endif

         enddo

c         if ( relaxon .eq. lone ) then
c            j4 = j+4
c            k4 = k+4
c            trelax = 0.2D-3     ! for neutron stars only !!!
c            vfac = one -  dtgr(igrd) / trelax
c            do i = 1, nx        !   trelax <= toscill (cf. Rasio et al.)
c               i4 = i + 4
c               energy(i,j,k,igrd) = energy(i,j,k,igrd) - ek(i)
c               velx  (i,j,k,igrd) =       relaxvx(i4,j4,k4,igrd) +
c     &             ( velx(i,j,k,igrd) - relaxvx(i4,j4,k4,igrd) ) * vfac
c               vely  (i,j,k,igrd) =       relaxvy(i4,j4,k4,igrd) +
c     &             ( vely(i,j,k,igrd) - relaxvy(i4,j4,k4,igrd) ) * vfac
c               velz  (i,j,k,igrd) =       relaxvz(i4,j4,k4,igrd) +
c     &             ( velz(i,j,k,igrd) - relaxvz(i4,j4,k4,igrd) ) * vfac
c               ek (i) = half * ( velx(i,j,k,igrd)**2
c     &             + vely(i,j,k,igrd)**2 + velz(i,j,k,igrd)**2 )
c               energy(i,j,k,igrd) = energy(i,j,k,igrd) + ek(i)
c            enddo
c         endif

         do i = 1, nx
            rho(i) = densty(i,j,k,igrd)
            ei (i) = (energy(i,j,k,igrd)-ek(i)) * rho(i)
            ye (i) = chem(i,j,k,igrd,1) / rho(i)
         enddo

         call eos( rho, tmp, ei, p, ye, gamc, game, nx, lzero )
c *** watch out: eos changes also ei and ye where necessary !!!


         do i = 1, nx
            energy(i,j,k,igrd)   = ei(i)/rho(i) + ek(i)
            chem  (i,j,k,igrd,1) = ye(i) * rho(i)
            temper(i,j,k,igrd)   = tmp(i)
            press (i,j,k,igrd)   = p(i)
         enddo

      enddo
      enddo
C$OMP END PARALLEL DO
c-----------------------------------

      call nospike3d(igrd)

      return
      end
c-----------------------------------------------------------------------
      subroutine viscos
! STILL TO BE CHECKED!!!
c     include physcial viscosity
c     input: cnvisc, u(), ut(), utt() rho()
c     output: uflx(), utflx(), uttflx()
c     PARALLEL OK

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      real*8 cnvi, rholim, rhounit, duflx, dutflx, duttflx, chamax

      itopgr=  indxgr(lone,lone)
      cnvi = cnvisc * delx(itopgr) * cc

      rhounit = (one/timeunit)**2 /pi/g
      rholim = rhounit * 3D-2

      chamax = 0.2D0   ! 20% maximum change in flux
      sufl = 1.D-20

      do i = 5, np5
         im = i-1

         if ( rhoav(i) .gt. rholim ) then

c         duflx   =  -  cnvi   / dx(i) * (
c     &         rhoav(i) * 4D0/3D0 * (u(i)-u(im))  - half*half * (
c     &    + (rhogt(i )+rhogb(i ))*2D0/3D0*(utgt (i )-utgb (i ))/two
c     &    + (rhogd(i )+rhogs(i ))*2D0/3D0*(uttgd(i )-uttgs(i ))/two
c     &    + (rhogt(im)+rhogb(im))*2D0/3D0*(utgt (im)-utgb (im))/two
c     &    + (rhogd(im)+rhogs(im))*2D0/3D0*(uttgd(im)-uttgs(im))/two))
c         dutflx  =  -  cnvi   / dx(i) * (
c     &         rhoav(i) * (ut(i)-ut(im)) +  half*half* (
c     &    + (rhogt(i )+rhogb(i )) *(ugt (i )-ugb (i ))/two
c     &    + (rhogt(im)+rhogb(im)) *(ugt (im)-ugb (im))/two   ))
c         duttflx =  -  cnvi   / dx(i) * (
c     &         rhoav(i) * (utt(i)-utt(im)) +  half*half* (
c     &    + (rhogd(i )+rhogs(i )) *(ugd (i )-ugs (i ))/two
c     &    + (rhogd(im)+rhogs(im)) *(ugd (im)-ugs (im))/two   ))
c
c         uflx  (i) = uflx  (i)  +
c     &         sign( (abs(uflx(i))  +sufl)*min( chamax,
c     &                  abs(duflx)  /(abs(uflx(i))  +sufl) ), duflx  )
c         utflx (i) = utflx (i)  +
c     &         sign( (abs(utflx(i)) +sufl)*min( chamax,
c     &                  abs(dutflx) /(abs(utflx(i)) +sufl) ), dutflx )
c         uttflx(i) = uttflx(i)  +
c     &         sign( (abs(uttflx(i))+sufl)*min( chamax,
c     &                  abs(duttflx)/(abs(uttflx(i))+sufl) ), duttflx)

         endif

      enddo

      return
      end

c=======================================================================
c                      B E T A   C O O L I N G                         =
c=======================================================================

      subroutine betacool(igrd)
c
c     B E T A   C O O L I N G
c     Danny Rhys Walker
c     - beta-cooling is used if bcool>0 and sbyn!=1. Crank-Nicolson is
c       used if cnyn==1, else forward Euler is used.
c
c     Input: igrd.
c
c     Output: none.
c
c-----------------------------------------------------------------------

      include 'qparam.cmn'  ! qx and qy
      include 'ppdisk.cmn'  ! bcool, r_in, r_out, cnyn, sbyn, ergirr, tirr4, xy_rs, omegak
      ! include 'squants.cmn' ! Don't need this?
      include 'aquants.cmn' ! energy, densty, all vels
      include 'grnest.cmn'  ! dtgr

      real*8 ergloss, ! Loss in internal energy density, Euler only
     &       ergnew, ! New internal energy
     &       cncoeff, ! Crank-Nicolson coefficient
     &       ek, ei, ! specific kinetic and internal energies
     &       tcool, ! Cooling time
     &       ergav, ! Mean specific internal energy on midplane
     &       cnirr, ! Crank-Nicolson irradiation coefficient
     &       tmptst, ! Test the temperature
     &       sigma(qx,qy,ngrd), ! Column/surface density
     &       ergrad(qx), ! Radiative contribution to energy change
     &       erghyd(qx), ! Hydrodynamic contirbution to energy change
     &       tau(qx), ! Optical depth
     &       tauden(qx), ! Denominator involving tau
     &       opac(qx), ! Opacity
     &       teq(qx), ergeq(qx), ! Equilibrium temp and energy
     &       tmp, ! Temperature
     &       rho(qx), ! Density
     &       ttherm ! Thermalisation timescale

      ergav = zero

      if(bcool .gt. zero) then
c   E X P O N E N T I A L    C O O L I N G
c   --------------------------------------

        if(cnyn .eq. zero) then
c   F O R W A R D   E U L E R   M E T H O D
c   ---------------------------------------
          if(igrd .eq. lone) print *, 'C O O L I N G; EULER'

          do i = 1, nx
            do j = 1, ny
              tcool = bcool / omegak(i,j,igrd)
              do k = 1, nz
                ek = half*(velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2 +
     &                        velz(i,j,k,igrd)**2)
                ei = energy(i,j,k,igrd) - ek
                tmptst = (gamma-one) * ei / specr
                if( ei.lt.ergirr ) then
                  energy(i,j,k,igrd) = ergirr + ek
                ! Perform cooling only within disk annulus
                elseif(xy_rs(i,j,igrd).gt.r_in .and.
     &                 xy_rs(i,j,igrd).lt.r_out     ) then

                  ! Make sure the temperature isn't cray
                  !if(tmptst.gt.1E3) stop'massive temp, cooling'

                  ergloss = - dtgr(igrd)*ei / tcool

                  ! Loss in internal energy should not exceed the internal
                  ! energy
                  if (abs(ergloss) .gt. ei) then
                    stop'ergloss wild, something gone wrong, cooling'
                  endif

                  energy(i,j,k,igrd) = ei + ergloss +
     &                               ( dtgr(igrd)/tcool ) * ergirr + ek

                  !For testing
                  if( k.eq.one .and. igrd.eq.lone ) then
                    ergav = ergav + energy(i,j,k,igrd)
                  endif
                endif
              enddo
            enddo
          enddo

          ergav = ergav / (nx * ny)

          if(igrd.eq.lone) then
            open(unit=777,file='ergseuler.txt')
            write(777,*) ergav
          endif

        else
c   C R A N K - N I C O L S O N   M E T H O D
c   -----------------------------------------
          if(igrd .eq. lone) print *, 'C O O L I N G; CRANK-NICOLSON'

          do i = 1, nx
            do j = 1, ny
              tcool = bcool / omegak(i,j,igrd)
              if(dtgr(igrd).gt.tcool) stop'dtgr too large, cooling'
              cncoeff = (2*tcool - dtgr(igrd))/(2*tcool + dtgr(igrd))
              cnirr = 2*dtgr(igrd) / (2*tcool + dtgr(igrd))
              do k = 1, nz
                ek = half*(velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2 +
     &                        velz(i,j,k,igrd)**2)
                ei = energy(i,j,k,igrd) - ek

                ! Put irradiation energy at all points
                if( ei.le.ergirr ) then
                  energy(i,j,k,igrd) = ergirr + ek
                  ! Update other state variables
c                  temper(i,j,k,igrd) = (gamma-one) * ergirr / specr
c                  press(i,j,k,igrd) = specr * densty(i,j,k,igrd) *
c     &                                temper(i,j,k,igrd)
                ! Perform cooling only within disk radial boundaries
                elseif(xy_rs(i,j,igrd).gt.r_in .and.
     &                 xy_rs(i,j,igrd).lt.r_out       ) then
                  ergnew = cnirr*ergirr + cncoeff*ei

                  ! Temperature should not be very high
c                  tmptst = (gamma-one) * ergnew / specr
c                  if(tmptst.gt.1E3) then
c                    print *, i,j,k,igrd, tmptst, cncoeff
c                    stop'temperature wild, cooling'
c                  endif

                  ! Should be cooling, not heating very much. Heating
                  ! Will happen often close to ergirr, but should only
                  ! be slight.
                  if (ergnew .gt. 1.01*ei) then
                    stop'ergnew growth, something gone wrong, cooling'
                  endif
                  ! Check we've not dipped below irradiation energy
                  if(ergnew .lt. ergirr) stop'ergnew below ergirr'

                  ! New energy with cooled internal
                  energy(i,j,k,igrd) = ergnew + ek
                  ! Update also the other state variables
c                  temper(i,j,k,igrd) = (gamma-one) * ergnew / specr
c                  press(i,j,k,igrd)  = specr * densty(i,j,k,igrd) *
c     &                                 temper(i,j,k,igrd)

                  ! Update other state variables
c                  press(i,j,k,igrd) = (gamma-one) * ergnew *
c     &                                              densty(i,j,k,igrd)
c                  temper(i,j,k,igrd) = (gamma-one) * ergnew / specr

                  !For purposes of testing this subroutine, we check the
                  !mean midplane internal energy.
                  if( k.eq.one .and. igrd.eq.lone) then
                    !if(i.eq.32) print *, (gamma-one) * ergnew / specr
                    ergav = ergav + ergnew
                  endif
                ! Make sure all points are at least irradiation energy
                endif
              enddo
            enddo
          enddo

          ergav = ergav / (nx * ny)

          if(igrd.eq.lone) then
            open(unit=777,file='ergs.txt')
            write(777,*) ergav, dtgr(igrd)
          endif

        endif


      elseif(bcool .lt. zero) then
c   E X P O N E N T I A L   H E A T I N G ! ? ! ?
c   ---------------------------------------------
! This shouldn't happen. Maybe in the future heating will be implemented
! but probably not in this subroutine
        stop'exponential heating! bcool negative.'

      else
c   C O O L I N G   O F F
c   ---------------------
        if(igrd .eq. lone) then
          print *, 'Cooling off'

          do i = 1,nx
            do j = 1,ny
              ek = half * (velx(i,j,1,igrd)**2 + vely(i,j,1,igrd)**2 +
     &                     velz(i,j,1,igrd)**2)
              ei = energy(i,j,1,igrd) - ek
              if(ei .lt. small*ek) then
                ei = small*ek
                energy(i,j,1,igrd) = ei + ek
              endif
              ergav = ergav + ei
            enddo
          enddo

          ergav = ergav / (nx * ny)

          open(unit=777,file='cooloff.txt')
          write(777,*) ergav

        endif

      endif

      return
      end

c=======================================================================
c               S T A M A T E L L O S   C O O L I N G                  =
c=======================================================================

      subroutine stamcool(igrd)
c
c     S T A M A T E L L O S   R A D I A T I V E   T R A N S F E R
c     Based on the method from Stamatellos et al. (2007).
c     Danny Rhys Walker
c
c     Input: igrd.
c
c     Output: none.
c
c-----------------------------------------------------------------------

      include 'qparam.cmn'  ! qx, qy, ngrd
      include 'ppdisk.cmn'  ! ergirr, tirr4
      ! include 'squants.cmn' ! Don't need this?
      include 'aquants.cmn' ! energy, densty
      include 'grnest.cmn'  ! dtgr

      real*8 ergnew,             ! New internal energy
     &       ek, ei,             ! specific kinetic and internal energies
     &       sigma(qx,qy,ngrd),  ! Column/surface density
     &       ergrad(qx),         ! Radiative contribution to energy change
     &       erghyd(qx),         ! Hydrodynamic contirbution to energy change
     &       tau(qx,qy,ngrd),    ! Optical depth
     &       tauden(qx), ! Denominator involving tau
     &       opac(qx),           ! Opacity
     &       teq(qx), ergeq(qx), ! Equilibrium temp and energy
     &       tmp,                ! Temperature
     &       rho(qx),            ! Density
     &       ttherm              ! Thermalisation timescale

      if(igrd .eq. lone) print *, 'C O O L I N G; STAMATELLOS'

      ! First we need optical depth
      ! call eosopac(igrd,opac) ! This is now in optdepz()
      ! call colden(sigma) ! Previous version used this
      call optdepz(igrd,tau)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ! tau(i) = opac(i,j,k,igrd) * sigma(i,j,igrd) ! Old version

            ! Denominator that appears in cooling law
            tauden(i) = tau(i,j,igrd) + 1/tau(i,j,igrd)

            ! Needed for eos
            rho(i) = densty(i,j,k,igrd)

            ! Radiative contribution to energy change
            tmp = temper(i,j,k,igrd)
            ergrad(i) = stefan * (tirr4 - tmp**4) / tauden(i)

            ! Specific kinetic energy
            ek = half*( velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2 +
     &                    velz(i,j,k,igrd)**2 )
            ! Specific internal energy
            ei = energy(i,j,k,igrd) - ek

            ! Hydrodynamic contribution to energy change
            erghyd(i) = (ei - eiprior(i,j,k,igrd)) / dtgr(igrd)

            ! Equilibrium temperature
            teq(i) = ( (tauden(i)/stefan)*erghyd(i) + tirr4 )**(1/4)

          enddo

          ! Get the equilibrium energy
          call eos(rho, teq, ergeq, peq, yeeq, gceq, geeq, nx, lone)

          do i = 1,nx
            ! Thermalisation timescale
            ttherm = (ergeq(i) - ei) / (ergrad(i) + erghyd(i))

            ! Finally evolve the energy
            ergnew = ergeq(i)+(ei-ergeq(i))*exp(-dtgr(igrd)/ttherm)

            if( ergnew .le. ergirr ) then
              ! At least irradiation energy
              energy(i,j,k,igrd) = ergirr + ek
            else
              energy(i,j,k,igrd) = ergnew + ek
            endif

          enddo
        enddo
      enddo

      return
      end

c=======================================================================
c                    C O L U M N   D E N S I T Y                       =
c=======================================================================

      subroutine colden(igrd,sig)
c
c     C O L U M N   (S U R F A C E)   D E N S I T Y
c     Danny Rhys Walker
c
c     Calculate surface density by integrating volume density in z.
c     Computationally this is approximated as a Riemann sum using the
c     grid length.
c
c     Input: igrd.
c
c     Output: sig.
c
c-----------------------------------------------------------------------

      include 'qparam.cmn'  ! qx, qy, ngrd
      include 'aquants.cmn' ! densty
      include 'grnest.cmn'  ! delx

      real*8 sig(qx,qy,ngrd)

      ! Include igrd for posterity. Could remove.
      do i = 1,nx
        do j = 1,ny
          sig(i,j,igrd) = 0.0
          do k = 1,nz
            sig(i,j,igrd) = sig(i,j,igrd)+delx(igrd)*densty(i,j,k,igrd)
          enddo
        enddo
      enddo

      return
      end

c=======================================================================
c                     O P T I C A L   D E P T H                        =
c=======================================================================

      subroutine optdepz(igrd,tauz)
c
c     O P T I C A L   D E P T H   (z-direction)
c     Danny Rhys Walker
c
c     For each (x,y) coordinate on the current grid compute the optical
c     depth in the z-direction and store the result in tauz(i,j,igrd).
c     The optical depth is the integral of kappa*rho wrt z. We simply
c     compute this as a Riemman sum using the grid length.
c
c     Input: igrd.
c
c     Output: tauz.
c
c-----------------------------------------------------------------------

      include 'qparam.cmn'  ! qx, qy, ngrd
      include 'aquants.cmn' ! densty

      real*8 taudif ! Differential
      real*8 tauz(qx,qy,ngrd)
      real*8 opac(qx,qy,qz,ngrd)

      ! Get the opacity
      call eosopac(igrd,opac)

      ! Do the Riemann sum
      do i = 1,nx
        do j = 1,ny
          tauz(i,j,igrd) = 0.0
          do k = 1,nz
            taudif = delx(igrd)*opac(i,j,k,igrd)*densty(i,j,k,igrd)
            tauz(i,j,igrd) = tauz(i,j,igrd) + taudif
          enddo
        enddo
      enddo

      return
      end
