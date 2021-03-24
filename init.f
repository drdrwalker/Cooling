c-------------------------------------------------------------------
      subroutine inippm

c     Define initial model: Protoplanetary disk

      include 'qparam.cmn'

      integer*4    qx1, qz1
      parameter ( qx1 = qx+1, qz1 = qz+1 )

      include 'ppdisk.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'
      include 'eos.cmn'
      include 'rays.cmn'  ! for iframe

      real*8 grdidi, gradrho, renv, dromax,
     &     solmas, totmas,       gasmas, dumas,
     &     dum1, dum2, rmua(300),  vlx(q),vly(q),
     &     rmum(q), rho(q), tau(q), drpos(q), pr(q), p(q), tmp(q),
     &     erg(q), gamc(q), game(q), vel(q), ek(q), e(q), ei(q), ye(q),
     &     tent(q), ent(q), dis1, dis2

      parameter( nn=10000, mm = 11000 )

      real*8    rhons(mm), radns(mm), ergns(mm), yens(mm), drns,
     &        tmpns(mm), prsns(mm),  mmns(mm), tmpns0(mm),
     &         dumm(mm), nergn(mm)

      real*8    rhotab(nn), ergtab(nn), yetab(nn), radtab(nn),
     &        mmtab(nn), prstab(nn),  amunu(nn), temtab(nn),
     &        cce(nn), ccn(nn), ccp(nn)

      real*8    rho2(2), tem2(2), yee2(2), cpn2(2), cpp2(2), ete2(2)

      integer*4 idum, iptr, nmax, num, menv, icggg, qi
      integer*4 ndi, ndip1, ndip2, ndim1,ndim2, nd1, nd2
ct    real*8 grav(q), d2phdx(q), dtga(q)

      integer*4  numin
      real*4     rhotin(nn),  ergtin(nn), yetain(nn), radtin(nn),
     &           mmtabin(nn), prstin(nn), amunin(nn)

c DRW - New quantities

      real*8 diskm

      ! Variables for testing
      real*8 eiini, ekini, tmptst

      ! Variables for averaging state quantities
      integer*4 xp, xm, yp, ym, zp, zm
      real*8 avdens(qx,qy,qz,ngrd)

      real*8 vkep(qx,qy),        !Keplerian orbital velocity
     &       vertsh,             !Vertical scale height
c     &       xy_rs(qx,qy,ngrd),  !xy-distance from centre of star
     &       xy_rd(qx,qy,ngrd),  !xy-distance from centre of disk
     &       z_dist(qz),         !z-distance from the midplane
     &       sigma0,             !Principal surface density
     &       sigma(qx,qy,ngrd),  !Surface density (qx == qy usually)
     &       vtheta,             !Initial azimuthal velocity
     &       tmp0,               !Principal temperature
     &       midrho              !Density at midplane

      sigma0 = pmass2 / ( 2*pi*(r_out - r_in) ) !Principal sigma based on double integral to disk mass
      !tmp0   = 50 * sqrt(au) !based on observed temperatures at 1AU.
      ! Principle temperature based on outer Q = 2
      tmp0 = ( 4*(pi**2)*g*(sigma0**2)*(r_out**1.5) ) /
     &       ( gamma*specr*pmass1 )

c Back to MRR

      if ( relaxon .eq. lone ) then
         write(*,*)  '*********************'
         write(*,*)  '**** RELAX PHASE ****'
         write(*,*)  '*********************'
      endif

c DRW - Following kept from Max's initialisation =======================
      dt  = dtini
      itopgr=  indxgr(lone,lone)
      dto = dt / dble( nfine ** (maxlev(lone)-lone) )
      d1x = delx(itopgr)

      cenx = dble(nxh) + half  !  for subroutines  'grdvel'
      ceny = dble(nyh) + half  !  and 'hydro'
      cenz =             half

      cmgax = cenx
      cmgay = ceny

c------- boundary conditions
c------- types of boundary, cf. 'getrwx'
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,i),
C$OMP+            SHARED(bndmnx,bndmxx,bndmny,bndmxy,bndmnz,bndmxz)
      do j = 1, qx
      do i = 1, qx
         bndmnx(i,j) = two
         bndmxx(i,j) = two
         bndmny(i,j) = two
         bndmxy(i,j) = two
         bndmnz(i,j) = one
         bndmxz(i,j) = two
      enddo
      enddo

      xcm = dble(nxh) + half
      ycm = dble(nxh) + half
      open(21, file=homedata//'ns_shen.data',
     &     form='unformatted', err=1002)
      read(21) numin, rhotin, yetain, radtin, mmtabin, prstin
      close(21)
      num    = numin
      do i = 1, num-1
         !radns(i) = radtab(i+1)
         radns(i) = dble(radtin(i+1))
      enddo
      num = num - 1
      drns = (radns(num)-radns(1)) / dble(num-lone)
      drns = drns       / d1x
      num = mm

c=======================================================================

!      print * , ' --------------------------------------------'
!      write(6,102)  rhons(1), tmpns(1)
!  102 format(1x,' rho0  = ', 1pe11.4,'      temp0 = ', e11.4)
!      write(6,105)  ergns(1), zero
!  105 format(1x,' erg0  = ', 1pe11.4,'       phi0 = ', e11.4)
!      print * , ' --------------------------------------------'
!      print * , ' '

c -- center of star and disk (was BH and NS)
      !pmass1 = solmas !Mass of star
      !pmass2 = 0.1*solmas !Mass of disk.
      dis1 = bidist * pmass2 / (pmass1+pmass2) !Barycentre radii
      dis2 = bidist * pmass1 / (pmass1+pmass2) !!
      ttx1 = xcm
      tty1 = ycm                  ! BH position
      ttx2 = xcm                  ! NS1 position
      tty2 = ycm                  ! NS1 position

      !bidist = sqrt ( (ttx1-ttx2)**2 + (tty1-tty2)**2 ) * d1x
      bidist = au                                                       !DRW -Temporary fix so I can centre 'BH'
      nd1 = nint(0.1D0/(drns*d1x/delx(indxgr(lone,mlev))))
      nd2 = two * nd1

      posx1 = ttx1
      posy1 = tty1
      polx1 = posx1
      poly1 = posy1

      posx2 = ttx2
      posy2 = tty2
      polx2 = posx2
      poly2 = posy2

      angm1 = zero
      angm2 = zero
      angcc1= zero
      lzgrav= zero
      eint1 = zero

c------- fill arrays with array values

      do 100 igrd = 1, ngrd

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(igrd,nx,ny,nz,tty1,ttx1,tty2,
C$OMP+       drns,nume,nd1,nd2,rhons,rhoran,ode,yens,tmpns,temper,ttx2,
C$OMP+       velx,vely,velz,densty,energy,tmpent,chem,gamma,smallu,
C$OMP+       d1x,nspin,ergns),
C$OMP+    PRIVATE(j,tmy1,tmy2,ty1,ty2,tmpy1,tmpy2,k,qi,
C$OMP+            i,tmx1,diss1,ndst1,tmx2,diss2,ndst2,ndi,erg,pr,gamc,
C$OMP+            game,ndip1,ndip2,ndim1,ndim2,rho,ye,tmp,tent,m,
C$OMP+            tmpy,tmpx,tmpz,ent)

      ! Radial distances
      do i=1,nx
        tmx1  = (ttx1 - topos(dble(i),lone,igrd))
        tmx2  = (ttx2 - topos(dble(i),lone,igrd))
      do j=1,ny
        tmy1 = (tty1 - topos(dble(j),ltwo,igrd))
        tmy2 = (tty2 - topos(dble(j),ltwo,igrd))

        ! xy distance from disk centre
        xy_rd(i,j,igrd)  = delx(igrd) * sqrt(tmy2**2 + tmx2**2) ! delx factor to convert to 'real' units
        ! xy distance from star centre
        xy_rs(i,j,igrd)  = delx(igrd) * sqrt(tmy1**2 + tmx1**2)
      enddo
      enddo

      do k=1,nz
        tmpz  = half - topos(dble(k),lthree,igrd)

        ! Distance from midplane
        z_dist(k) = delx(igrd) * tmpz
      enddo

      open(unit=909,file='midrho.txt')

      do j=1,ny
         tmy1 = (tty1 - topos(dble(j),ltwo,igrd))
         tmy2 = (tty2 - topos(dble(j),ltwo,igrd))
         ty1  = tmy1**2
         ty2  = tmy2**2
      do k=1,nz
         !tmz2  = abs(half - topos(dble(k),lthree,igrd))                 !DRW - I added this but not using it currently
         tmpz  = ((half - topos(dble(k),lthree,igrd)))**2
         tmpy1 = tmpz +  ty1
         tmpy2 = tmpz +  ty2
      do i=1,nx
         tmx1  = (ttx1 - topos(dble(i),lone,igrd))
         diss1 = tmpy1 + tmx1**2
         ndst1 = nint ( sqrt(diss1) / drns )
         tmx2  = (ttx2 - topos(dble(i),lone,igrd))
         diss2 = tmpy2 + tmx2**2
         ndst2 = nint ( sqrt(diss2) / drns )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are now calculated in previous loops
         !These quantities are calculated in grid units, so we require
         !a factor of the grid space, delx(igrd)
         ! xy distance from disk centre
         !xy_rd(i,j,igrd)  = delx(igrd) * sqrt(ty2 + tmx2**2)
         ! xy distance from star centre
         !xy_rs(i,j,igrd)  = delx(igrd) * sqrt(ty1 + tmx1**2)
         ! Distance from midplane
         !z_dist = delx(igrd) * sqrt(tmpz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ndim1 = max( ndi-nd1, lone )
         ndim2 = max( ndi-nd2, lone )

         vkep(i,j)   = sqrt( g * pmass1 / xy_rs(i,j,igrd) )
         omegak(i,j,igrd) = vkep(i,j) / xy_rs(i,j,igrd)

         ! If we are inside the radial disk boundaries
         if( xy_rs(i,j,igrd).gt.r_in  .and.
     &       xy_rs(i,j,igrd).lt.r_out       ) then

           ! Temperature profile
           tmp(i) = tmp0 / sqrt( xy_rs(i,j,igrd) )
           if(tmp(i).lt.smltem) tmp(i) = smltem
           ! Vertical scale height
           vertsh = sqrt(specr * tmp(i)) / omegak(i,j,igrd)

!!!!!!!!!!!!!!!!!!!!!!!! I have removed this check, it seems unnecessary
           ! More than 99% of matter is within 2rt(2)h of the midplane
!           if( z_dist.lt.(2*sqrt(2.0)*vertsh) ) then
!!!!!!!!!!!!!!!!!!!!!!!!

           ! Surface density profile
           sigma(i,j,igrd) = sigma0 / xy_rs(i,j,igrd)

           ! Volume density
           midrho = sigma(i,j,igrd) / (vertsh*sqrt(2*pi))

           ! Record midplane density
           if(i.eq.32 .and. k.eq.1 .and. igrd.eq.1) then
             write(909,*) midrho
           endif

           rho(i) = midrho * exp( -(z_dist(k)**2) / (2*(vertsh**2)) )
           if(rho(i).lt.smlrho) then
!             if(igrd.eq.one) print *, 'small rho:', i, j, k
             rho(i) = smlrho
           endif

           ye (i) = one !DRW. This electron fraction is a relic.

           ! Keplerian flow, i.e. orbital velocities
           velx(i,j,k,igrd)= -vkep(i,j)*(delx(igrd)*tmy1 /
     &                                    xy_rs(i,j,igrd))
           if(abs(velx(i,j,k,igrd)).lt.smallu) then
             velx(i,j,k,igrd) = smallu
           endif

           vely(i,j,k,igrd) = vkep(i,j)*(delx(igrd)*tmx1/
     &                                    xy_rs(i,j,igrd))
           if(abs(vely(i,j,k,igrd)).lt.smallu) then
             vely(i,j,k,igrd) = smallu
           endif

           velz(i,j,k,igrd) =  smallu

!           else
           ! We are above/below the main disk
!           tmp(i) = smltem
!           velx(i,j,k,igrd) = smallu
!           vely(i,j,k,igrd) = smallu
!           velz(i,j,k,igrd) = smallu
!           rho(i) = smlrho
!           ye(i) = 1D-20

!           endif

         else
           ! We are outside the radial disk boundaries
           tmp(i) = smltem
           velx(i,j,k,igrd) = smallu
           vely(i,j,k,igrd) = smallu
           velz(i,j,k,igrd) = smallu
           rho(i) = smlrho
           ye(i) = 1D-20

         endif

         if(velx(i,j,k,igrd).gt.cc .or.
     &      vely(i,j,k,igrd).gt.cc .or.
     &      velz(i,j,k,igrd).gt.cc) stop'exceeding light speed, init'

         ! Specific kinetic energy
         ek(i) = half * ( velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2 +
     &                    velz(i,j,k,igrd)**2 )
         if(ek(i).eq.zero) stop'ek(i) zero!'

c Back to MRR.

           tent(i) = 0.05D0     !  small temperature
c           tent(i) = tmp(i)

      enddo

      ! Get energy from eos
      call eos(rho,tmp,erg,dumm,dumm,dumm,dumm,nx,lone)
      ! The eos gives energy density so we need to convert
      do i = 1,nx
        erg(i) = erg(i) / rho(i)
        if(erg(i).lt.smalle) erg(i) = smalle
        ! Check these energies are reasonable
        tmptst = (gamma-one) * erg(i) / specr
        !if(tmptst.gt.1E3) stop'massive temp, inippm'
      enddo

c------- fill arrays
      do i=1,nx
        densty(i, j, k, igrd) = rho(i)
        ! Sort kinetic and internal energies and make sure internal
        ! is not too much smaller than kinetic
        energy(i, j, k, igrd) = ek(i) + max(erg(i), small*ek(i))
        temper(i, j, k, igrd) = tmp(i)
        ! This is a relic that's yet to be dealt with
        tmpent(i, j, k, igrd) = tent(i)
      enddo

      call eosent( rho, tent, ye, ent, nx, lone )     !  tempr -> entropy

      do i = 1, nx
         chem  (i, j, k, igrd, lone) = ye (i) * rho(i)
         chem  (i, j, k, igrd, ltwo) = ent(i) * rho(i)
      enddo

      qi = qc ! just to trick compiler not to complain
      do m = 3, qi
      do i = 1, nx
         chem  (i, j, k, igrd, m) = zero
      enddo
      enddo

      enddo
      enddo ! End of big loop on this grid


C$OMP END PARALLEL DO

 100  continue ! This is where the grid loop continues/ends

      print *, 'max en init', maxval(energy)/1E10
      print *, 'max tmp init', maxval(temper)

      open(unit=341,file='rho1.txt')
      open(unit=342,file='rho2.txt')
      open(unit=343,file='rho3.txt')
      open(unit=344,file='rho4.txt')
      open(unit=441,file='tmp1.txt')
      open(unit=442,file='tmp2.txt')
      open(unit=443,file='tmp3.txt')
      open(unit=444,file='tmp4.txt')
      open(unit=541,file='xyd1.txt')
      open(unit=542,file='xyd2.txt')
      open(unit=543,file='xyd3.txt')
      open(unit=544,file='xyd4.txt')
      open(unit=641,file='xys1.txt')
      open(unit=642,file='xys2.txt')
      open(unit=643,file='xys3.txt')
      open(unit=644,file='xys4.txt')
      do i = 1,nx
        write(341,*) densty(i,32,1,1)
        write(342,*) densty(i,32,1,2)
        write(343,*) densty(i,32,1,3)
        write(344,*) densty(i,32,1,4)
        write(441,*) temper(i,32,1,1)
        write(442,*) temper(i,32,1,2)
        write(443,*) temper(i,32,1,3)
        write(444,*) temper(i,32,1,4)
        write(541,*) xy_rd(i,32,1)/au
        write(542,*) xy_rd(i,32,2)/au
        write(543,*) xy_rd(i,32,3)/au
        write(544,*) xy_rd(i,32,4)/au
        write(641,*) xy_rs(i,32,1)/au
        write(642,*) xy_rs(i,32,2)/au
        write(643,*) xy_rs(i,32,3)/au
        write(644,*) xy_rs(i,32,4)/au
      enddo

      print *, 'densty:', densty(32,32,1,1)

c------- inflow boundary values

      !This will do for now
      rhoe = smlrho
      erge = smalle

      rhoin  = rhoe
      uin    = smallu
      utin   = smallu
      uttin  = smallu
      gin    = -1.D15
      gamein = gamma
      gamcin = gamma
      ein    = erge
c      write(*,*)'inippm:pin ',pin
      pin = smallp
      tin = 1.0
      do m = 1, qc
         chein(m) = 1.D-20
      enddo
      chein(1)= ye(1) * rhoe                                            !DRW - edited from yens(nume) to ye(1)
      chein(2)= ten**entl(nroent,1,1) * rhoe
c  smallest entropies are at largest densities
      smlrho = rhoin
!      do m = 1, qc
!         smlche(m) = chein(m)
!      enddo

!      write(*,'(A,1I4,A,1P,1E13.6,A,0P,1F8.5,A,1F8.5)')
!     &      ' init edge of ns: ',nume, '   rhoe: ',rhoe,
!     &      '   tin: ', tin

c------- unify grids, without point masses

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,igrd),
C$OMP+            SHARED( nx,ny,nz,gpot,gtpot,rgarr)
      do j = 1, ny
      do igrd = 1, ngrd
      do k = 1, nz
      do i = 1, nx
         gpot(i,j,k,igrd) = zero
         gtpot(i,j,k,igrd) = zero
         rgarr(i,j,k,igrd) = zero
      enddo
      enddo
      enddo
      enddo

C$OMP END PARALLEL DO

      do lev = mlev, 2, -1
         call fineup(lev,lzero)
      enddo

c      do i = 1, nx
c        do j = 1, ny
c          do k = 1, nz
c            ekini = half * (velx(i,j,k,1)**2 + vely(i,j,k,1)**2
c     &                      + velz(i,j,k,1)**2)
c            eiini = energy(i,j,k,1) - ekini
c            tmptst = (gamma-one)*eiini / specr
c            if(tmptst.gt.1E3) then
c              print *, eiini, ekini, xy_rs(i,j,1)/au, i, j, k
c              stop'high temp inippm, fineup'
c            endif
c          enddo
c        enddo
c      enddo

c------- initial mass of point masses

      do igrd = 1, ngrd

         gasmas = zero
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,densty), REDUCTION(+:gasmas)
         do j = 1, ny-1
         do k = 1, nz-1
         do i = 1, nx-1
            gasmas = gasmas + densty(i,j,k,igrd)
         enddo
         enddo
         enddo

C$OMP END PARALLEL DO
         gasmas = two*gasmas * delx(igrd)**3

         gasmas = zero
         do i = 1,nx-1
         do j = 1,ny-1
         do k = 1,nz-1
           gasmas = gasmas + densty(i,j,k,igrd) + densty(i+1,j,k,igrd)
     &                 + densty(i,j+1,k,igrd) + densty(i+1,j+1,k,igrd)
     &                 + densty(i,j,k+1,igrd) + densty(i+1,j,k+1,igrd)
     &             + densty(i,j+1,k+1,igrd) + densty(i+1,j+1,k+1,igrd)
         enddo
         enddo
         enddo
         gasmas = 0.25*(gasmas * delx(igrd)**3)
         !print *, 'new gasmas:', igrd, gasmas/solmas

c         gasmas = zero
c         do i = 1, nx
c           do j = 1, ny
c             gasmas = gasmas + sigma(i,j,igrd)
c           enddo
c         enddo
c         gasmas = two * gasmas * delx(igrd)**2

c         print *, 'gasmas 2:', gasmas/solmas
c         print *, ''
c         stop


         write(6,'(a,i2,a,f9.6)')
     &           ' gasmas on grid', igrd, ' :', gasmas/solmas
      enddo

      gasmas = zero
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,itopgr,densty), REDUCTION(+:gasmas)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         gasmas = gasmas + densty(i,j,k,itopgr)
      enddo
      enddo
      enddo

C$OMP END PARALLEL DO
      gasmas = two*gasmas * d1x**3
      totmas = gasmas + pmass1 + pmass2


c------- unify grids, now with point masses

      do 161 lev = mlev, 2, -1
  161    call fineup(lev,lzero)

      !
      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            ekini = half * (velx(i,j,k,1)**2 + vely(i,j,k,1)**2
     &                      + velz(i,j,k,1)**2)
            eiini = energy(i,j,k,1) - ekini
            tmptst = (gamma-one)*eiini / specr
c            if(tmptst.gt.1E3) stop'high temp inippm, fineup 2'
          enddo
        enddo
      enddo

      sumass = zero

      do 500 levl = 1, mlev
      numgr = indxgr(lzero,levl)

      do 500 num = 1, numgr
         lgr = indxgr(num,levl)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,ro,densty,lgr)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            ro(i,j,k) = densty(i,j,k,lgr)
         enddo
         enddo
         enddo

C$OMP END PARALLEL DO

         do 120 ifgr = 1, ngrd
 120        if ( idxcgr(ifgr) .eq. lgr ) call rozero(ifgr)
         !
         do i = 1, nx
           do j = 1, ny
             do k = 1, nz
               ekini = half * (velx(i,j,k,1)**2 + vely(i,j,k,1)**2
     &                      + velz(i,j,k,1)**2)
               eiini = energy(i,j,k,1) - ekini
               tmptst = (gamma-one)*eiini / specr
c               if(tmptst.gt.1E3) stop'high temp inippm, rozero'
             enddo
           enddo
         enddo

         sum = zero

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,ro), REDUCTION(+:sum)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            sum = sum + ro(i,j,k)
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      sumass = sumass + two * sum * delx(lgr)**3

 500  continue

c      write(*,*) 'gasmas+pmass1: ', gasmas/solmas,
c     &           '    sumass+pmass1: ',(sumass+pmass1)/solmas
      print *, ''
      write(*,'(a,F10.6)') ' sumass: ',sumass/solmas
      print *, ''

      absakerr = abs(akerr)
      if ( (absakerr .ge. one) .and. (absakerr.lt.two) ) then
         angm1 = (absakerr-one) / cc * (g * pmass1**2)
      elseif ( (absakerr .gt. zero) .and. (absakerr.lt.one) ) then
         angm1 = absakerr / cc * (g * pmass1**2)
      elseif ( absakerr .ge. two ) then
         write(*,*) 'inippm something wrong,  akerr',akerr
         stop
      endif
      angm1 = sign(angm1,akerr)

      if ( akerr .eq. zero ) then
         rota = zero
      elseif ( absakerr .ge. one ) then
         rota = abs(angm1) * cc / (g * pmass1**2)
      else
         rota = absakerr
      endif

      thrd = one/3D0
      z1 = one  +  (one-rota**2)**thrd
     &                    *( (one+rota)**thrd + (one-rota)**thrd )
      z2 = sqrt( 3D0*rota**2 + z1**2 )
      rin = 3D0 + z2 - sign( sqrt( (3D0-z1)*(3D0+z1+two*z2) ), akerr )
c       rin = 3D0 + z2 - sqrt( (3D0-z1)*(3D0+z1+two*z2) )   ! prograde
c       rin = 3D0 + z2 + sqrt( (3D0-z1)*(3D0+z1+two*z2) )   ! retrograde

      r1 = one + sqrt(one-rota**2)
      rs = two * g * pmass1 / cc**2
      rbhkerr = ( rin + r1 ) / two   *rs/two
c                        average      rin,r1 in units of rs/2
      beta = rin / r1  - one
      beta2 = two - beta
      r1 = r1/two * rs       ! change to  Schwarzschild Radius at rs

c------- velocities

      vcmx = smallu
      vcmy = smallu

c DRW - I will comment all of this out, except bits that seem like they
c could be needed. Max said this omega averaging is unneeded.

c      if     ( npawi .eq. 0 ) then
c         omega= sqrt( g * (pmass1+pmass2) / bidist**3 ) !Newt pot
c         write(*,*) 'Newt init omega ',omega
c      elseif ( npawi .eq. 1 ) then
c     omega= sqrt( g * (pmass1+pmass2) /bidist /(bidist-rs)**2 ) ! PaWi pot
c         omega= sqrt( g * (pmass1+pmass2) /bidist       ! KerrPaWi pot
c     &                      / ( bidist**beta2 * (bidist-r1)**beta ) )
c         write(*,*) 'PaWi init omega ',omega
c      else
c         write(*,*) 'something wrong with npawi variable'
c         stop 'init'
c      endif

c--- try to find more appropriate omega (extended NS instead of point mass)
c    by mass averaging omegas

c      avromega = zero
c      wgt = zero

c      do 502 levl = 1, mlev                ! loop over all levels
c      do 502 num = 1, indxgr(lzero,levl)   ! loop over all grids on levl
c         igrd = indxgr(num,levl)
c         dgx3 = delx(igrd)**3

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,ro,densty)
         do 112 j = 1, ny
         do 112 k = 1, nz
         do 112 i = 1, nx
            ro(i,j,k) = densty(i,j,k,igrd)
  112    continue

C$OMP END PARALLEL DO

c        zero area covered by finer grids
c         ifgr = idxfgr(igrd)
c         if( ifgr .ne. -lone ) call rozero(ifgr)

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(igrd,nx,ny,nz,tty1,ttx1,
C$OMP+    ro,d1x,npawi,beta,beta2,r1,g,pmass1,pmass2,dgx3),
C$OMP+    PRIVATE(j,tmy1,tmpy1,k,tmpz,i,tmx1,rrr),
C$OMP+    REDUCTION(+:avromega,wgt)

c      do j=1,ny
c         tmy1  = ( tty1 - topos(dble(j),ltwo,igrd) )**2
c      do k=1,nz
c         tmpz  = ( half - topos(dble(k),lthree,igrd) )**2
c         tmpy1 =   tmpz +  tmy1
c      do i=1,nx
c           tmx1  = ( ttx1 - topos(dble(i),lone,igrd) )**2
c           rrr = sqrt(tmpy1 + tmx1) * d1x

c           if ( ro(i,j,k) .gt. 2.D10 ) then                            !DRW - Commented out since ro() won't ever be that high
!           if ( ro(i,j,k) .gt. 1D-20 ) then                             !DRW - New ro() condition
c              wgt = wgt + ro(i,j,k)*dgx3
c              if     ( npawi .eq. 0 ) then   ! Newt pot
c                 avromega = avromega + ro(i,j,k)*dgx3 *
c     &                  sqrt( g * (pmass1+pmass2) / rrr**3 )
c              elseif ( npawi .eq. 1 ) then   ! KerrPaWi pot
c                 avromega = avromega + ro(i,j,k)*dgx3 *
c     &                  sqrt( g * (pmass1+pmass2) /rrr
c     &                   / ( rrr**beta2 * (rrr-r1)**beta ) )
c              endif
c           endif      ! of density check

c      enddo
c      enddo
c      enddo

c  502 continue

c      if (.not.(avromega.eq.avromega)) then   ! avromega = NaN
c         write(*,*) 'too big rhoe, >2E10?', avromega, rhoe
c         stop
c      endif

c      if (wgt .ne. zero) then
c         avromega = avromega / wgt
c      else
c         write(*,*) 'something wrong with wgt ', wgt
c      endif

c      write(*,*) 'init, omega average ', avromega

c      omega = avromega

c DRW - Back to business after this.

c --- end omega average

      if ( nspin .eq. 4 ) then   !colli

        print *, 'Should not be here, nspin==4, init'

         velx1 = - vrns * dis1 / bidist
     &        - sqrt( two * g * sumass**2 / bidist / (pmass1+sumass) )
         vely1 = zero
         velx2 = + vrns * dis2 / bidist
     &        + sqrt( two * g * pmass1**2 / bidist / (pmass1+sumass) )
         vely2 = zero
      else
         !velx1 = +vrns  * dis1 / bidist                                !DRW - Commented out
         !vely1 = +omega * dis1                                         !!
         !velx2 = -vrns  * dis2 / bidist                                !!
         !vely2 = -omega * dis2                                         !!
         velx1 = zero
         vely1 = zero
         velx2 = zero
         vely2 = zero
      endif

      vnsx = velx2
      vnsy = vely2

      omegaspin = -two * pi / 2D-3    !  spin period = 2 ms

c      write(*,'(1P,3(A,1E12.4))') 'init vrns: ',vrns,'   omega: ',omega,!DRW - Commented out.
c     &        '   omegaspin: ',omegaspin                                !!
c      write(*,'(1P,2(a,1E12.4))') 'init velx1 ',velx1,'  vely1 ',vely1
c      write(*,'(1P,2(a,1E12.4))') 'init velx2 ',velx2,'  vely2 ',vely2

      do 301 igrd = 1, ngrd

      if ( nspin .eq. 4 ) then
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,erg),
C$OMP+    SHARED(nx,ny,nz,igrd,velx,vely,velz,vnsx,vnsy,energy,small)

      print *, 'Should not be here, nspin==4, init'

      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         velx(i,j,k,igrd) = velx(i,j,k,igrd) * vnsx
         vely(i,j,k,igrd) = vely(i,j,k,igrd) * vnsy
         erg(i) = half * (  velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2
     &                   + velz(i,j,k,igrd)**2 )
         energy(i,j,k,igrd) = erg(i) +
     &             max( energy(i,j,k,igrd), small*erg(i) )
         stop'we should not be here, inippm 1'
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
      else
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,sgt,erg),
C$OMP+            SHARED(nx,ny,nz,igrd,velx,vely,velz,temper,energy,
C$OMP+                   vnsx,vnsy,omegaspin,small)
      do j = 1, ny                                                      !DRW -  EDIT HERE?
      do k = 1, nz
      do i = 1, nx
         if ( velx(i,j,k,igrd) .ne. zero ) then  ! inside NS
            sgt = sign( one, temper(i,j,k,igrd) )
            !velx(i,j,k,igrd) = sgt*( vnsx + velx(i,j,k,igrd)*omegaspin )!DRW - Commented out
            !vely(i,j,k,igrd) = sgt*( vnsy + vely(i,j,k,igrd)*omegaspin )!!
c            erg(i) = half * (  velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2
c     &                   + velz(i,j,k,igrd)**2 )
c            energy(i,j,k,igrd) = erg(i) +
c     &             max( energy(i,j,k,igrd), small*erg(i) )
         endif
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
      endif

c      if ( relaxon .eq. lone ) then
c         do kk = 1, nz+qb
c            k = min( max( kk-qb/2, lone), nz )
c         do jj = 1, ny+qb
c            j = min( max( jj-qb/2, lone), ny )
c         do ii = 1, nx+qb
c            i = min( max( ii-qb/2, lone), nx )
c
c            relaxvx(ii,jj,kk,igrd) = relaxvx(ii,jj,kk,igrd) * vnsy
c     &       + sign(vnsx,temper(i,j,k,igrd))*abs(relaxvy(ii,jj,kk,igrd))
c            relaxvy(ii,jj,kk,igrd) = relaxvy(ii,jj,kk,igrd) * vnsy
c
c         enddo
c         enddo
c         enddo
c      endif


C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,temper)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         temper(i, j, k, igrd) = abs(temper(i, j, k, igrd))
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c *** now only for relax mode
      if ( relaxon .eq. lone ) then

        stop'relaxon, init.f'

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,tmy1,tmx1,omega,erg),
C$OMP+            SHARED(nx,ny,nz,igrd,xcm,ycm,rhoin,d1x,
C$OMP+                   velx,vely,energy,densty)
         do j=1,ny                                                      !DRW - Notice that this loop does nothing
            tmy1  = topos(dble(j),ltwo,igrd) - ycm
         do k=1,nz
         do i=1,nx
c            if ( velx(i,j,k,igrd) .ne. zero ) then ! inside NS
            if ( densty(i,j,k,igrd) .ge. 4.D0*rhoin ) then
c           matter inside NS as taken in subroutine calm

               tmx1  = topos(dble(i),lone,igrd) - xcm

c               erg(i)= half*( velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2)

c               velx(i,j,k,igrd) = velx(i,j,k,igrd)! + omega*tmy1*d1x    !DRW - commented out omega part.
c               vely(i,j,k,igrd) = vely(i,j,k,igrd) !- omega*tmx1*d1x    !!

c               energy(i,j,k,igrd) = energy(i,j,k,igrd)  - erg(i)
c     &               + half*( velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2)

            endif
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

         velx1 = zero
         vely1 = zero
         velx2 = zero
         vely2 = zero

      endif

  301 continue

      vnsx = velx2
      vnsy = vely2

      vlox1 = velx1
      vloy1 = vely1
      vlox2 = velx2
      vloy2 = vely2

      dtgs = 1.D10*dtmax

      do igrd = 1, ngrd
         dtgg(igrd) = 1.D10*dtmax
      enddo

      tempo = one
      awdprx = zero
      awdpry = zero

      enuloss = zero
      eneloss = zero
      enaloss = zero
      enxloss = zero
      anelo   = one
      analo   = one
      anxlo   = one

      egquel  = zero

      iframe = lzero

c --- first without fineup

      !
      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            ekini = half * (velx(i,j,k,1)**2 + vely(i,j,k,1)**2
     &                      + velz(i,j,k,1)**2)
            eiini = energy(i,j,k,1) - ekini
            tmptst = (gamma-one)*eiini / specr
c            if(tmptst.gt.1E3) stop'high temp inippm, presucks'
          enddo
        enddo
      enddo

      do igrd = ngrd, 1, -1
         call sucks(igrd,lzero)
      enddo

      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            ekini = half * (velx(i,j,k,1)**2 + vely(i,j,k,1)**2
     &                      + velz(i,j,k,1)**2)
            eiini = energy(i,j,k,1) - ekini
            tmptst = (gamma-one)*eiini / specr
c            if(tmptst.gt.1E3) stop'high temp inippm, sucks 1'
          enddo
        enddo
      enddo

c     now with fineup

      do igrd = ngrd, 1, -1
         call sucks(igrd,lone)
      enddo

      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            ekini = half * (velx(i,j,k,1)**2 + vely(i,j,k,1)**2
     &                      + velz(i,j,k,1)**2)
            eiini = energy(i,j,k,1) - ekini
            tmptst = (gamma-one)*eiini / specr
c            if(tmptst.gt.1E3) stop'high temp inippm, sucks 2'
          enddo
        enddo
      enddo

      pmass2 = zero   !  not needed any more

      do igrd = ngrd, 1, -1
         call kickbh(igrd,lzero)    ! reset to zero
      enddo

      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            ekini = half * (velx(i,j,k,1)**2 + vely(i,j,k,1)**2
     &                      + velz(i,j,k,1)**2)
            eiini = energy(i,j,k,1) - ekini
            tmptst = (gamma-one)*eiini / specr
c            if(tmptst.gt.1E3) stop'high temp inippm, kickbh'
          enddo
        enddo
      enddo

c---  find machine precision
      te = two
      ib = -lone
 20   te = te/two
      ib = ib+lone
      if ( one+te .ne. one ) goto 20
      signif = two**dble(lthree-ib)
c---                 3 : leave 3 bits for rounding
c      write(*,'(a,i3,a)') 'machine precision: ',ib,' bits'

      return

c File I/O errors *****************************************************

 1002 write(*,*) "Couldn't find ns_shen.data"
      stop

      end
c-----------------------------------------------------------------------
      subroutine inispec
c Produces the array r,g,b,tempr(0:1000)    output: from 0.0 to 1.0
c The array is an approximation for the spectrum.

      include 'qparam.cmn'
      include 'rays.cmn'

        integer*4  rma1, rmi1, rma2, rmi2
        integer*4  gma1, gmi1, gma2, gmi2, bmi1, bma1

        parameter(rma1=250, rmi1=400, rmi2=950 ,rma2=1050)
        parameter(gmi1=150, gma1=250, gma2=650 ,gmi2=1000)
        parameter(bmi1=600, bma1=750)

        integer*4    i, ia, ib

      !print *, 'init, inispec'

        do 5 i = 0, len
           rtempr(i)=zero
           gtempr(i)=zero
           btempr(i)=zero
  5     continue

* red
        do 10 i=0,rma1
 10        rtempr(i) = one

        ia = rma1+1
        ib = rmi1
        do 20 i=ia,ib
 20        rtempr(i) = one -  dble(i-ia) / dble(ib-ia)

        ia = rmi2+1
        ib = len
        do 30 i=ia,ib
 30        rtempr(i) = half + dble(i-ia) / dble(rma2-ia) / two


* green
        ia = gmi1+1
        ib = gma1
        do 40 i=ia,ib
  40       gtempr(i) =    dble(i-ia) / dble(ib-ia)

        do 50 i=gma1+1,gma2
  50      gtempr(i) = one

        ia = gma2+1
        ib = gmi2
        do 60 i=ia,ib
 60       gtempr(i) = one  - dble(i-ia) / dble(ib-ia)

* blue
        ia = bmi1+1
        ib = bma1
        do 70 i=ia,ib
 70       btempr(i) =     dble(i-ia) / dble(ib-ia)

        do 80 i=bma1+1,len
 80       btempr(i) = one

* gray
        rtempr(len) = 0.3D0
        gtempr(len) = 0.3D0
        btempr(len) = 0.3D0

        return
        end
c-----------------------------------------------------------------------
      subroutine inivid


      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'rays.cmn'
      include 'grnest.cmn'

      !print *, 'init, inivid'

* --- INITIALIZE VARIABLES

      iendf = 250*nx/32 * 2/ifilm

      do 5 i = 0, nview
!         theta(i) = -44.D0                                             !DRW - Commented out
         theta(i) = 90.D0
c        phi  (i) = 160.0D0 + 1.D0*i                                    !DRW - The appearance of i, a coordinate?
!         phi  (i) = 95.0D0 + 90.D0                                     !DRW - Commented out
         phi(i) = 0.D0                                                  !DRW - "Birds-eye view"
c         if ( i .le. iendf ) then
c            ri = dble(i)/(dble(iendf)/pi)
c            scale(i) = 2.0D0 + 0.5 * cos(ri)
c         else
c            scale(i) = 2.0D0 + 0.5 * cos(pi)
c         endif
         scale(i) = 2.0D0
  5   continue

      call inispec

* --- FIND MIN and MAX OF DENSTY AND TEMPR

c     nxy  = nx  * ny
c     nxyz = nxy * nz

c     irhmax = ismax(nxyz, denvec, lone)
c     dymax   = denvec(irhmax)
c     irhmin = ismin(nxyz, denvec, lone)
c     dymin   = denvec(irhmin)
c     itemax = ismax(nxyz, temvec, lone)
c     tpmax   = temvec(itemax)
c     itemin = ismin(nxyz, temvec, lone)
c     tpmin   = temvec(itemin)

      dymin = 10D0
c     dymax = half
c from input
c      tpmin = zero
      tpmin =  10D0
c      tpmax = 30D0
c from input

      if (lenstar .gt. 4*ntra) then
         write(*,*) 'lenstar should not be too big because of trace'
         stop
      endif

c      rand = 0.6D0
c      do i = 1, lenstar
c         rand = (rand+5D0)**ten
c         rand = rand-dble(int(rand))
c         phistar(i) = rand* two*pi
c         rand = (rand+5D0)**ten
c         rand = rand-dble(int(rand))
c         thestar(i) = acos(two*rand-one) - pi/two
c         rand = (rand+5D0)**ten
c         rand = rand-dble(int(rand))
c         shistar(i) = rand * ten + 50D0*float(npix)/511.D0
c      enddo

      itopgr=  indxgr(lone,lone)
      tframe = time - two*dtgr(itopgr)

      return
      end
c-----------------------------------------------------------------------
      subroutine inigrid

c     one time /nest/ initialisation
c     only equidistant cartesian grids are used.
c     the topmost grid has to be grid number 1 !! This is used elsewhere
c     the position of WD is implied in origins of grids 2,4,6

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'

      integer*4 levcnt(mlev)

      real*8 orgin(ngrd,3)

      !print *, 'init, inigrid'
c      equivalence( norgin, orgin )


c --- planar geometry, cf. 'geom'
      igeom = lzero

      maxlev(1) = mlev
      maxlev(2) = mlev

c --- allocate grids to levels: is done in blokda !!

c --- check conistancy

      nsg = lzero
      do 20 lev = 1, mlev
         nsg = nsg + indxgr(lzero,lev)
 20   continue

      if ( (nsg.ne.ngrd) .or. (nsg.lt.mlev) ) then
         write(*,*) 'error in consitancy of indxgr!'
         stop'=> grids'
      endif

      do 30 igr = 1, ngrd
         levlgr(igr,1) = lzero
 30   continue

c --- given grid, return (1) level and (2) column in indxgr
      do 40 levl = 1, mlev
       do 40 igr = 1, indxgr(lzero,levl)
          levlgr(indxgr(igr,levl),lone) = levl
          levlgr(indxgr(igr,levl),ltwo) = igr
 40   continue

c     in which fine grid lies
      do 50 igr = 2, ngrd
          idxcgr(igr) = indxgr( levlgr(igr,ltwo), levlgr(igr,lone)-lone)
 50   continue

c --- topmost grid does not have any coarser grid

      idxcgr(indxgr(lone,lone)) = -lone

      do 55 igr = 1, ngrd-1
          idxfgr(igr) = indxgr( levlgr(igr,ltwo), levlgr(igr,lone)+lone)
 55   continue

c --- finest grid does not have any finer grid

      idxfgr(indxgr(lone,maxlev(lone))) = -lone
c     idxfgr(indxgr(ltwo,maxlev(ltwo))) = -lone

c --- check consistancy

      do 60 igr = 1, ngrd
         if ( levlgr(igr,lone) .eq. lzero ) then
            write(*,*) 'error, a grid is missing in levlgr!'
            stop'=> grids'
         endif
 60   continue


      nx4 = nx/4

      orgin(1,1) =  zero
      orgin(1,2) =  zero
      orgin(1,3) =  zero

      do 70 ig = 2, indxgr(lone,maxlev(lone))

      norgin(ig,1) = nx4
      norgin(ig,2) = nx4
      norgin(ig,3) = lzero

 70   continue


c --- topmost grid (index 1) is absolute and needs no origin

c --- all cell lengths are equal

      deltax = gridlx / nx

      do 10 igr = 1, ngrd
         delx(igr) = deltax / dble(nfine**(levlgr(igr,lone)-lone))
 10   continue

c --- initialize time

      do 80 igr = 1, ngrd
         dtgr(igr) = dtini  / dble(nfine**(levlgr(igr,lone)-lone))
 80   continue

c --- initialize scheduling of levels

      do 110 levl = 1, mlev
         levcnt(levl) = lone
110   continue

      do 111 node = 1, 2000
         ifinup(node) = .false.
111   continue

      levl = lone
      node = lzero

120   continue

         if ( levcnt(levl) .le. nfine ) then
            node = node + 1
            ishedul(node) = levl
            if ( levcnt(levl) .eq. nfine ) ifinup(node) = .true.
            levcnt(levl) = levcnt(levl) + 1
            levl = levl + 1
         else
            levcnt(levl) = 1
            levl = levl - 1
         endif

         if (levl .gt. mlev) then
            levl = mlev
         endif

         nodenr = node

      if (levl .ne. 1 ) goto 120

      if (nodenr .gt. 2000) then
         write(*,*) 'inigrid,  2000 .lt. nodenr:',nodenr
         stop
      endif

      return
      end
c-----------------------------------------------------------------------
      block data blokda

c     to initialize indxgr, instead of in inigrid.

      include 'qparam.cmn'
      include 'grnest.cmn'

      data indxgr / 1,  ngrd * 1 ,
     &              1,  ngrd * 2 ,
     &              1,  ngrd * 3 ,
     &              1,  ngrd * 4  /
c     &              1,  ngrd * 5    /

      end
c-----------------------------------------------------------------------
      subroutine inigrav
c------- initialize background potentials

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'
      include 'aquants.cmn'

      !print *, 'init, inigrav'

c------- background of other levels

      do igrd = 1, ngrd
         call fillold(igrd)
         do j = 1, 3
         do i = 1, 3
            gdij3dt(i,j,igrd) = zero
            gdij2dt(i,j,igrd) = zero
         enddo
         enddo
      enddo

      do 172 lev = 1, mlev
      do 172 num = 1, indxgr(lzero,lev)
         igrd = indxgr(num,lev)
         call gpress      ( igrd )
         call bhdenup     ( igrd )
         call gravty      ( igrd )
         call gwaves      ( igrd )
         call rgpot       ( igrd )
         call strain      ( igrd )
         if ( mconoff .eq. lone )  call chigravty( igrd )
c         call bhdendown   ( igrd )
  172 continue

c------- set back densty to small values for hydro
c      do k = 1, nz
c      do j = 1, ny
c      do i = 1, nx
c         if ( bhdens(i,j,k,igrd) .ne. zero )
c     &             densty(i,j,k,igrd) = 5.D0*rhoin  ! 5 because of calm
c      enddo
c      enddo
c      enddo

      return
      end
