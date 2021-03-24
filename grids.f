c-----------------------------------------------------------------------
      subroutine plotieee

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'

      real*4  dty1(qx,qy), dty2(qx,qy), dty3(qx,qy), dty4(qx,qy),
     &        dty5(qx,qy), dty6(qx,qy)
      real*8  e(q),ek(q),ei(q),tmp(q),rho(q),gamc(q),game(q),p(q),ye(q)

      integer*4 qx2, qy2
      parameter( qx2=qx/2, qy2=qy/2 )

      !print *, 'grids, plotieee'

      write(44) time,posx1,posy1,pmass1,gridlx

      do nngr = 1, ngrd

c         call bhdenup(nngr)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,i,ye),
C$OMP+            SHARED(densty,velx,vely,tmpent,gpot,bhdens,
C$OMP+                   chem,dty1,dty2,dty3,dty4,dty5,dty6,nngr)
         do j = 1, qy
         do i = 1, qx
            dty1(i,j) = densty(i,j,lone,nngr)
            if ( bhdens(i,j,lone,nngr) .ne. zero ) dty1(i,j)=-one
            dty2(i,j) = sign(
     &        max(1.D-30,abs(velx(i,j,lone,nngr))),velx(i,j,lone,nngr))
            dty3(i,j) = sign(
     &        max(1.D-30,abs(vely(i,j,lone,nngr))),vely(i,j,lone,nngr))
            dty4(i,j) =  max(1.D-30, abs(tmpent(i,j,lone,nngr)))
            dty5(i,j) = sign(
     &        max(1.D-30,abs(gpot(i,j,lone,nngr))),gpot(i,j,lone,nngr))
         if (bhdens(i,j,lone,nngr).ne.zero)chem(i,j,lone,nngr,lone)=zero
            ye(i) = chem(i,j,lone,nngr,lone) / densty(i,j,lone,nngr)
            dty6(i,j) = sign( max(1.D-30, abs(ye(i))), ye(i) )
         enddo
         enddo
C$OMP END PARALLEL DO

         write(44) dty1
         write(44) dty2
         write(44) dty3
         write(44) dty4
         write(44) dty5
         write(44) dty6

c         call bhdendown(nngr)

      enddo

c ---

      do nngr = 1, ngrd

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,k,ye),
C$OMP+            SHARED(densty,velx,vely,velz,tmpent,gpot,bhdens,
C$OMP+                   chem,dty1,dty2,dty3,dty4,dty5,dty6,nngr)
       do i = 1, qx
       do k = 1, qz
          dty1(i,k)    = densty(i,qy2,k,nngr)
          dty1(i,k+qz) = densty(qx2,i,k,nngr)
          if ( bhdens(i,qy2,k,nngr) .ne. zero ) dty1(i,k)=-one
          if ( bhdens(qx2,i,k,nngr) .ne. zero ) dty1(i,k+qz)=-one
          dty2(i,k) = sign(
     &      max(1.D-30,abs(velx(i,qy2,k,nngr))),velx(i,qy2,k,nngr) )
          dty2(i,k+qz) = sign(
     &      max(1.D-30,abs(vely(qx2,i,k,nngr))),vely(qx2,i,k,nngr) )
          dty3(i,k) = sign(
     &      max(1.D-30,abs(velz(i,qy2,k,nngr))),velz(i,qy2,k,nngr) )
          dty3(i,k+qz) = sign(
     &      max(1.D-30,abs(velz(qx2,i,k,nngr))),velz(qx2,i,k,nngr) )
          dty4(i,k)    =  max(1.D-30, abs(tmpent(i,qy2,k,nngr)))
          dty4(i,k+qz) =  max(1.D-30, abs(tmpent(qx2,i,k,nngr)))
          dty5(i,k) = sign(
     &      max(1.D-30,abs(gpot(i,qy2,k,nngr))),gpot(i,qy2,k,nngr) )
          dty5(i,k+qz) = sign(
     &      max(1.D-30,abs(gpot(qx2,i,k,nngr))),gpot(qx2,i,k,nngr) )
          if (bhdens(i,qy2,k,nngr).ne.zero) chem(i,qy2,k,nngr,lone)=zero
          ye(i) = chem(i,qy2,k,nngr,lone) / densty(i,qy2,k,nngr)
          dty6(i,k) = sign( max(1.D-30, abs(ye(i))), ye(i) )
          if (bhdens(qx2,i,k,nngr).ne.zero) chem(qx2,i,k,nngr,lone)=zero
          ye(i) = chem(qx2,i,k,nngr,lone) / densty(qx2,i,k,nngr)
          dty6(i,k+qz) = sign( max(1.D-30, abs(ye(i))), ye(i) )
       enddo
       enddo
C$OMP END PARALLEL DO

       write(44) dty1
       write(44) dty2
       write(44) dty3
       write(44) dty4
       write(44) dty5
       write(44) dty6

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine fineup(levl,nfix)

c --- bring up to date all grids one level coarser
c     which contain finer grids on level levl

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include 'ppdisk.cmn'

c      call gasumas('infine')

co    logical overlap

      !print *, 'grids, fineup'

c --- take all grids on level levl
      do 300 inum = 1, indxgr(lzero,levl)
         ifgr = indxgr(inum,levl)

c --- take all grids one level coarser
      do 300 inuf = 1, indxgr(lzero,levl-lone)
         icgr = indxgr(inuf,levl-lone)

c --- copy overlap of grids ifgr and icgr
co    if (overlap(ifgr,icgr)) then

c         call avrage(ifgr,icgr)
c         do i = 1, nx
c           do j = 1, ny
c             do k = 1, nz
c               ekgrb = half * (velx(i,j,k,1)**2 + vely(i,j,k,1)**2
c     &                         + velz(i,j,k,1)**2)
c               eigrb = energy(i,j,k,1) - ekgrb
c               tmpgrb = (gamma-one)*eigrb / specr
c               if(tmpgrb.gt.1E3) then
c                 print *, i,j,k, tmpgrb, ekgrb, eigrb
c                 stop'high temp avrage'
c               endif
c             enddo
c           enddo
c         enddo

         if ( nfix .eq. lone ) then
            call flupfix(ifgr,icgr)
            call fluzero(ifgr)
         endif
co    endif

 300  continue

c      call gasumas('exfine')

      return
      end
c-----------------------------------------------------------------------
      subroutine fillold(igrd)

c --- fill arrays in common block vold with values from vnew

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,velx,vely,velz,
C$OMP+                   densty,energy,vxold,vyold,vzold,dold,eold)

      !print *, 'grids, fillold'

      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         vxold(i,j,k,igrd) = velx  (i,j,k,igrd)
         vyold(i,j,k,igrd) = vely  (i,j,k,igrd)
         vzold(i,j,k,igrd) = velz  (i,j,k,igrd)
         dold (i,j,k,igrd) = densty(i,j,k,igrd)
         eold (i,j,k,igrd) = energy(i,j,k,igrd)
      enddo
      enddo
      enddo

C$OMP END PARALLEL DO

      ergboun = zero
      ergwave = zero

      return
      end
c-----------------------------------------------------------------------
      subroutine tstep(new)

c     compute new timestep value, search all levels
c     use this routine only for (nsdim .eq. 3) !!!

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'

      real*8  scrch1(q), scrch2(q), scrch3(q), scrch4(q)
      real*8  iscrch(q), dphi(q)
      real*8  dtfac(ngrd), ceul(ngrd), ceuligr, dtggigr
      real*8  rho(q),ek(q),e(q),ei(q),tmp(q),p(q),gamc(q),game(q),ye(q)
      integer*4  new

      !!print *, 'grids, tstep'

      if (nsdim.ne.lthree) then
         write (6,*) 'nsdim must be 3 in this version of tstep.'
         stop'=> tstep'
      endif

      olddt  = dtgr(lone)

c-------  check CFL-condition

      tempo = one
      dtemp = zero
      do 20 igr = 1, ngrd

         ceul(igr) = zero
         ceuligr   = zero
         dtgg(igr) = 1D10*dtmax
         dtggigr   = 1D10*dtmax
         dex = delx(igr)
         dtfac(igr) = dble( nfine ** ( levlgr(igr,lone) -lone ) )
         dphi(nx) = 1.D-99

C$OMP PARALLEL DO DEFAULT(NONE), REDUCTION(+:dtemp,tempo),
C$OMP+         REDUCTION(max:ceuligr),REDUCTION(min:dtggigr),
C$OMP+         SHARED(densty,velx,vely,velz,energy,chem,
C$OMP+                temper,dtfac,ceul,velmax,gpot,press,
C$OMP+                nx,ny,nz,igr,small,dex),
C$OMP+         PRIVATE(j,k,i,ye,jj,kk,dphi),
C$OMP+         PRIVATE(rho,ek,e,ei,tmp,p,gamc,game)
         do j = 1, ny
         do k = 1, nz

            do i = 1, nx
               rho   (i) = densty(i,j,k,igr)
               ek    (i) =  half *
     &        (velx(i,j,k,igr)**2+vely(i,j,k,igr)**2+velz(i,j,k,igr)**2)
               e     (i) = energy(i,j,k,igr)
               ei    (i) = rho(i)  *  ( e(i) - ek(i) )
               ei    (i) = max( ei(i), rho(i)*small*ek(i) )
               ye    (i) = chem(i,j,k,igr,lone) / rho(i)
               tmp   (i) = temper(i,j,k,igr)
            enddo

            do i = 1, nx
               if( e(i) .lt. 0.999D0*ek(i) )
     &         write(*,'(A,4I3,A,1P,1E14.5,A,1E14.5)')' tstep  energy:',
     &            i,j,k,igr, '  ek:',ek(i), '  e:',e(i)
            enddo

c      call eos(rho(1),tmp(1),ei(1),p(1),ye(1),gamc(1),game(1),nx,lzero)
            call eos( rho, tmp, ei, p, ye, gamc, game, nx, lone)

            do i = 1, nx
               press (i,j,k,igr) = p(i)
               ceuligr = max ( ceuligr, sqrt(gamc(i)*p(i)/rho(i)) )
            enddo

            jj = min( j, ny-1 )
            kk = min( k, nz-1 )
            do i = 1, nx-1
               dphi(i) = abs( gpot(i,jj,kk,igr)-gpot(i+1,jj,kk,igr) ) +
     &                   abs( gpot(i,jj,kk,igr)-gpot(i,jj+1,kk,igr) ) +
     &                   abs( gpot(i,jj,kk,igr)-gpot(i,jj,kk+1,igr) )
               dtggigr = min( dtggigr, dtfac(igr) * dex *
     &                 (ceul(igr)+velmax(igr)) / max(dphi(i), 1.D-99)  )
            enddo

            if ( igr .eq. ngrd ) then
               do i = 1, nx
                  dtemp = dtemp + rho(i)
                  tempo = tempo + tmp(i)*rho(i)
               enddo
            endif

         enddo
         enddo
C$OMP END PARALLEL DO

         ceul(igr) = ceuligr
         dtgg(igr) = dtggigr
         dtcc(igr) = dtfac(igr)*dex/(ceul(igr)+velmax(igr)) / sqrt(3.D0)

  20  continue

      tempo = tempo / dtemp

      do igr = 1, ngrd
         dtc = min( dtc, cfl*dtcc(igr) )
      enddo

c-------  check motion of two point-masses:
c         less than 0.4 cells at the finest level per timestep

c     dtp = dtmax
      dtp = 0.4D0 / (! max (
     &  sqrt(1.D-3 + velx1**2 + vely1**2) / (
     &delx(indxgr(lone,maxlev(lone)))*dble(nfine**(maxlev(lone)-lone))))
c    &, sqrt(1.D-3 + velx2**2 + vely2**2) / (
c    &delx(indxgr(ltwo,maxlev(ltwo)))*dble(nfine**(maxlev(ltwo)-lone))))


c-------  check gravitational oscillation periods:
c         omega = sqrt(4 pi g rho) = sqrt(laplace phi)  !not any more!
c         This ( dtg = dtfac*min(1 / omega) ) is calculated in accel !

c     do 120 igr = 1, ngrd
c        dtg = min( dtg, cgr*dtgg(igr) )
         dtg = min( dtg, cgr*dtgg(indxgr(lone,maxlev(lone))) )
c        dtgg(igr) = 1.D10*dtmax
 120  continue


c-------  take smallest of all timesteps

!      dt = min( dtc, dtp, dtg, half*dtn, half*dte )                            !Commented out - DRW
!      if ( dtg .lt. dtmin )  dt = min( dtc, dtp, half*dtn, half*dte )          !!!
!      if ( half*dtn .lt. dtmin )  dt = min( dtc, dtp, dtg, half*dte )          !!!

!      dt = min( dt, 1.2D0*olddt )                                              !!!
!      if (dt .lt. dtmin)  then                                                 !!!
!         write(*,*) 'dt:',dt,'  dtc:',dtc,'   dtp:',dtp,'   dtg:',dtg,         !!!
!     &                       '  dtn:',dtn,'   dte:',dte                        !!!

c DRW  - Remove neutrino timesteps, dtn and dte

      dt = min( dtc, dtp, dtg )
      if ( dtg .lt. dtmin )  dt = min( dtc, dtp )

      dt = min( dt, 1.2D0*olddt )
      if (dt .lt. dtmin)  then
         write(*,*) 'dt:',dt,'  dtc:',dtc,'   dtp:',dtp,'   dtg:',dtg

!      dt = min( dtc, dtg )                                                      !DRW - Removing dtp since NaN
!      if ( dtg .lt. dtmin )  dt = dtc

!      dt = min( dt, 1.2D0*olddt )
!      if (dt .lt. dtmin)  then
!         write(*,*) 'dt:',dt,'  dtc:',dtc,'   dtp:',dtp,'   dtg:',dtg

c Back to MRR

         call cloplt
c         call clograv
         call restrt(lzero)

         stop'=> tstep'
      endif

      dt = min(dt,dtmax)

!      if ( itstp .ne. lzero ) then                                             !Commented out - DRW
!      if ( mod(nstep,itstp) .eq. lzero ) then                                  !!!
!      write(*,1001) 'nstep: ',nstep, '    dt:',dt,                             !!!
!      &   '  dtc:', dtc, '  dtp:', dtp, '  dtg:', dtg,                         !!!
!      &   '  dte:', dte, '  dtn:', dtn                                         !!!
!      write(*,'(x,a,1P,10(10G11.3/))') 'tstep   cfl:  ', dtcc                  !!!
!      write(*,'(x,a,1P,10(10G11.3/))') 'tstep   pot:  ', dtgg                  !!!
!      endif                                                                    !!!
!      endif                                                                    !!!
!      do igr = 1, ngrd                                                         !!!
!         dtgg(igr) = 1.D10*dtmax                                               !!!
!      enddo                                                                    !!!

c DRW - Remove neutrino timestep, dtn and dte, numero dos

      if ( itstp .ne. lzero ) then
      if ( mod(nstep,itstp) .eq. lzero ) then
      write(*,1001) 'nstep: ',nstep, '    dt:',dt,
     &   '  dtc:', dtc, '  dtp:', dtp, '  dtg:', dtg
      write(*,'(x,a,1P,10(10G11.3/))') 'tstep   cfl:  ', dtcc
      write(*,'(x,a,1P,10(10G11.3/))') 'tstep   pot:  ', dtgg
      endif
      endif
      do igr = 1, ngrd
         dtgg(igr) = 1.D10*dtmax
      enddo

c Back to MRR

 1001 format(1x,a,1i5,10(a,1pe8.1))

      if ( new .eq. lone ) then
         do igr = 1, ngrd
            dtgr(igr) = dt / dtfac(igr)
            velmax(igr) = one
         enddo
         dtc = 1.D10*dtmax
         dtg = 1.D10*dtmax
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine getrwx (j,k,igrd)

c     get row j,k from /vnew/ of igrd,
c     add boundary values and put in /hydro/
c     PARALLEL OK

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include 'compct.cmn'
      include    'eos.cmn'

      real*8  tme(q)

      !print *, 'grids, getrwx'

      nnn = nx
      np1 = nnn + 1
      np2 = nnn + 2
      np3 = nnn + 3
      np4 = nnn + 4
      np5 = nnn + 5
      np6 = nnn + 6
      np7 = nnn + 7
      np8 = nnn + 8

      do i = 1, q
         dgrav(i) = zero
         tme  (i) = smltme
      enddo

c --- take one row from data arrays

      do i = 1, nnn
         i4 = i + 4
         rho (i4) = densty(i,j,k,igrd)
         u   (i4) = velx  (i,j,k,igrd)
         if(u(i4).ne.u(i4)) then
           print *, 'velx:', velx(i,j,k,igrd)
           print *, i,j,k,igrd
           stop'u NaN in grids, velx'
         endif
         ut  (i4) = vely  (i,j,k,igrd)
         utt (i4) = velz  (i,j,k,igrd)
         e   (i4) = energy(i,j,k,igrd)
         grav(i4) = gpot  (i,j,k,igrd)
         ugridx(i,igrd) = zero
         tme (i4) = tmpent(i,j,k,igrd)
      enddo

      !print *, 'getrwx vels'
      !print *, maxval(velx)
      !print *, maxval(vely)
      !print *, maxval(velz)

      do m = 1, qc
      do i = 1, nnn
         i4 = i + 4
         che(i4,m) = chem(i,j,k,igrd,m)
      enddo
      enddo

      do i4 = 1, np8                ! non conforming to
         ugrid (i4) = zero
      enddo
c9    ugrid(i) = ugridx(lone,igrd)        ! anything else !!!

c      if ( relaxon .eq. lone ) then
c         do i4 = 1, np8
c            urlx (i4) = relaxvx(i4,j+4,k+4,igrd)
c         enddo
c      endif

c --- put boundary values from coarser grids

      call bndrwx( rho (1), bxdens(1,1,1), j, k )
      call bndrwx( u   (1), bxvelx(1,1,1), j, k )
      call bndrwx( ut  (1), bxvely(1,1,1), j, k )
      call bndrwx( utt (1), bxvelz(1,1,1), j, k )
      call bndrwx( e   (1), bxener(1,1,1), j, k )
      call bndrwx( grav(1), bxgpot(1,1,1), j, k )

      do m = 1, qc
         call bndrwx( che(1,m), bxchem(1,1,1,m), j, k )
      enddo

      do ii = 1, np8
         ek  (ii) =  half * ( u(ii)**2 + ut(ii)**2 + utt(ii)**2 )
         ei  (ii) = ( e(ii) - ek(ii) ) * rho(ii)
         ye  (ii) = che(ii,1) / rho(ii)
         dx  (ii) = delx(igrd)
         dtdx(ii) = dtgr(igrd) / delx(igrd)
      enddo

      do ii = 1, np8
         if(ei(ii).lt.zero) write(*,*) 'getrwx  ',ii,j,k,igrd
      enddo

      call eos( rho, tmp, ei, p, ye, gamc, game, np8, lzero )
c *** watch out: eos changes also ei and ye where necessary !!!

!      call nospike( rho, tmp, tme, ei, p, ye, gamc, game, np8 )

      do i = 1, np8
         che(i,1) = ye(i) * rho(i)
         e(i) = ei(i)/rho(i) + ek(i)
c        v(i)  = one / rho(i)
         c(i)  = sqrt( gamc(i) * p(i) * rho(i) )
         ce(i) = c(i) / rho(i)
      enddo

      do i = 1, nnn
         i4 = i + 4
         temper(i,j,k,igrd) = tmp(i4)
         press (i,j,k,igrd) = p  (i4)
      enddo

c     needed in states and in accel
      dtppm = dtgr(igrd)

      !print *, 'getrwx'
      !print *, velx(1,5,8,igrd), vely(6,3,2,igrd), vely(10,54,6,igrd)

      return
      end
c-----------------------------------------------------------------------
      subroutine bndrwx( row, bxplane, j, k )
c     add boundary values and put in /hydro/
c     PARALLEL OK: PRIVATE(row,j,k), SHARED(bxplane)

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 row(q), bxplane(qb, qy, qz)

      !print *, 'grids, bndrwx'

      row   (1) = bxplane (1,j,k)
      row   (2) = bxplane (2,j,k)
      row   (3) = bxplane (3,j,k)
      row   (4) = bxplane (4,j,k)
      row (np5) = bxplane (5,j,k)
      row (np6) = bxplane (6,j,k)
      row (np7) = bxplane (7,j,k)
      row (np8) = bxplane (8,j,k)

      return
      end
c-----------------------------------------------------------------------
      subroutine putrwx (j,k,igrd)
c     get /hydro/..nu values and put them in /vnew/
c     PARALLEL OK

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include 'compct.cmn'

      integer*4 j,k,igrd,i,i4,m,i4m,i4p
c     these must all be private

      !print *, 'grids, putrwx'

      do i = 1, nnn
         i4 = i + 4
         ek (i4) = half*(unu(i4)**2 + utnu(i4)**2 + uttnu(i4)**2)
      enddo

!      call lightning(unu(5),utnu(5),uttnu(5),enu(5),ek(5), nnn)

      do i = 1, nx
         i4 = i + 4
         densty(i,j,k,igrd) = rhonu(i4)
         velx  (i,j,k,igrd) = unu  (i4)
         vely  (i,j,k,igrd) = utnu (i4)
         velz  (i,j,k,igrd) = uttnu(i4)
         energy(i,j,k,igrd) = enu  (i4)
      enddo

      do m = 1, qc
      do i = 1, nx
         i4 = i + 4
         chem(i,j,k,igrd,m) = max( chenu(i4,m), smlche(m) )
      enddo
      enddo

      enu(4)   = max(enu(5)  -ek(5),   small*ek(5)  )
      enu(np5) = max(enu(np4)-ek(np4), small*ek(np4))
      ek (4)   = ek(5)
      ek (np5) = ek(np4)

      do i = 1, nx
         i4 = i + 4
         enu(i4) = enu(i4) - ek(i4)
      enddo

      do i = 1, nx
         i4 = i + 4
         i4m = max( i4-1, 11-i4 )
         i4p = min( i4+1, 2*nx+7-i4)
         if ( enu(i4) .le. small*ek(i4) ) then
            enu(i4) = half * ( enu(i4m) + max(enu(i4p),zero) )
            energy(i,j,k,igrd) = max( enu(i4), small*ek(i4) ) + ek(i4)
cc        if( igrd.ne.indxgr(ltwo,maxlev(ltwo)) )
c            if ( bhdens(i,j,k,igrd) .eq. zero )
c     &         write(*,'(a,6i4,1P,3E19.10)')
c     &             'putrwx, energy:',i,j,k,igrd,i4p,i4m,
c     &                energy(i,j,k,igrd), enu(i4), ek(i4)
         endif
      enddo

      !print *, 'putrwx'
      !print *, velx(1,5,8,igrd), vely(6,3,2,igrd), vely(10,54,6,igrd)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine getrwy (i,k,igrd)
c     get row i,k from /vnew/ of igrd,
c     add boundary values and put in /hydro/
c     PARALLEL OK

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include 'compct.cmn'
      include    'eos.cmn'

      real*8  tme(q)

      !print *, 'grids, getrwy'

      nnn = ny
      np1 = nnn + 1
      np2 = nnn + 2
      np3 = nnn + 3
      np4 = nnn + 4
      np5 = nnn + 5
      np6 = nnn + 6
      np7 = nnn + 7
      np8 = nnn + 8

      do j = 1, q
         dgrav(j) = zero
         tme  (j) = smltme
      enddo

c --- take one row from data arrays

      do j = 1, nnn
         j4 = j + 4
         rho (j4) = densty(i,j,k,igrd)
         u   (j4) = vely  (i,j,k,igrd)
         if(u(j4).ne.u(j4)) stop'u, vely'
         ut  (j4) = velx  (i,j,k,igrd)
         utt (j4) = velz  (i,j,k,igrd)
         e   (j4) = energy(i,j,k,igrd)
         grav(j4) = gpot  (i,j,k,igrd)
         ugridy(j,igrd) = zero
         tme (j4) = tmpent(i,j,k,igrd)
      enddo

      do m = 1, qc
      do j = 1, nnn
         j4 = j + 4
         che(j4,m) = chem(i,j,k,igrd,m)
      enddo
      enddo

      do j4 = 1, np8                ! non conforming to
         ugrid (j4) = zero
c     ugrid(j) = ugridy(lone,igrd)        ! anything else !!!
      enddo

c      if ( relaxon .eq. lone ) then
c         do j4 = 1, np8
c            urlx (j4) = relaxvy(i+4,j4,k+4,igrd)
c         enddo
c      endif

      call bndrwy( rho (1), bydens(1,1,1), i, k )
      call bndrwy( u   (1), byvely(1,1,1), i, k )
      call bndrwy( ut  (1), byvelx(1,1,1), i, k )
      call bndrwy( utt (1), byvelz(1,1,1), i, k )
      call bndrwy( e   (1), byener(1,1,1), i, k )
      call bndrwy( grav(1), bygpot(1,1,1), i, k )

      do m = 1, qc
         call bndrwy( che(1,m), bychem(1,1,1,m), i, k )
      enddo

      do jj = 1, np8
         ek  (jj) =  half * ( u(jj)**2 + ut(jj)**2 + utt(jj)**2 )
         ei  (jj) = ( e(jj) - ek(jj) ) * rho(jj)
         ye  (jj) = che(jj,1) / rho(jj)
         dx  (jj) = delx(igrd)
         dtdx(jj) = dtgr(igrd) / delx(igrd)
      enddo

      do jj = 1, np8
         if(ei(jj).lt.zero) write(*,*) 'getrwy  ',i,jj,k,igrd
      enddo

      call eos( rho, tmp, ei, p, ye, gamc, game, np8, lzero )

      ! print *, 'max en getrwy', maxval(energy)/1E10
c *** watch out: eos changes also ei and ye where necessary !!!

!      call nospike( rho, tmp, tme, ei, p, ye, gamc, game, np8 )

      do j = 1, np8
         che(j,1) = ye(j) * rho(j)
         e(j) = ei(j)/rho(j) + ek(j)
c        v (j) = one / rho(j)
         c (j) = sqrt( gamc(j) * p(j) * rho(j) )
         ce(j) = c(j) / rho(j)
      enddo

      do j = 1, nnn
         j4 = j + 4
         temper(i,j,k,igrd) = tmp(j4)
         press (i,j,k,igrd) = p  (j4)
      enddo

      dtppm = dtgr(igrd)

      return
      end
c-----------------------------------------------------------------------
      subroutine bndrwy( row, byplane, i, k )
c     add boundary values and put in /hydro/
c     PARALLEL OK: PRIVATE(row,i,k), SHARED(byplane)

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 row(q), byplane(qx, qb, qz)

      !print *, 'grids, bndrwy'

      row   (1) = byplane (i,1,k)
      row   (2) = byplane (i,2,k)
      row   (3) = byplane (i,3,k)
      row   (4) = byplane (i,4,k)
      row (np5) = byplane (i,5,k)
      row (np6) = byplane (i,6,k)
      row (np7) = byplane (i,7,k)
      row (np8) = byplane (i,8,k)

      return
      end
c-----------------------------------------------------------------------
      subroutine putrwy (i,k,igrd)
c     get /hydro/..nu values and put them in /vnew/
c     PARALLEL OK

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include 'compct.cmn'

      integer*4 i,k,igrd,j,j4,m,j4m,j4p
c     these must all be private

      !print *, 'grids, putrwy'

      do j = 1, nnn
         j4 = j + 4
         ek (j4) = half*(unu(j4)**2 + utnu(j4)**2 + uttnu(j4)**2)
      enddo

!      call lightning(unu(5),utnu(5),uttnu(5),enu(5),ek(5), nnn)

      do j = 1, ny
         j4 = j + 4
         densty(i,j,k,igrd) = rhonu(j4)
         vely  (i,j,k,igrd) = unu  (j4)
         velx  (i,j,k,igrd) = utnu (j4)
         if(utnu(j4).ne.utnu(j4)) stop'utnu NaN putrwy'
         velz  (i,j,k,igrd) = uttnu(j4)
         energy(i,j,k,igrd) = enu  (j4)
      enddo

      do m = 1, qc
      do j = 1, ny
         j4 = j + 4
         chem(i,j,k,igrd,m) = max( chenu(j4,m), smlche(m) )
      enddo
      enddo

      enu(4)   = max(enu(5)  -ek(5),   small*ek(5)  )
      enu(np5) = max(enu(np4)-ek(np4), small*ek(np4))
      ek (4)   = ek(5)
      ek (np5) = ek(np4)

      do j = 1, ny
         j4 = j + 4
         enu(j4) = enu(j4) - ek(j4)
      enddo

      do j = 1, ny
         j4 = j + 4
         j4m = max( j4-1, 11-j4 )
         j4p = min( j4+1, 2*ny+7-j4)
         if ( enu(j4) .le. small*ek(j4) ) then
            enu(j4) = half * ( enu(j4m) + max(enu(j4p),zero) )
            energy(i,j,k,igrd) = max( enu(j4), small*ek(j4) ) + ek(j4)
cc        if( igrd.ne.indxgr(ltwo,maxlev(ltwo)) )
c            if ( bhdens(i,j,k,igrd) .eq. zero )
c     &         write(*,'(a,6i4,1P,3E19.10)')
c     &             'putrwy, energy:',i,j,k,igrd,j4p,j4m,
c     &                energy(i,j,k,igrd), enu(j4), ek(j4)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine getrwz (i,j,igrd)
c     get row i,j from /vnew/ of igrd,
c     add boundary values and put in /hydro/
c     PARALLEL OK

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include 'compct.cmn'
      include    'eos.cmn'

      real*8  tme(q)

      !print *, 'grids, getrwz'

      nnn = nz
      np1 = nnn + 1
      np2 = nnn + 2
      np3 = nnn + 3
      np4 = nnn + 4
      np5 = nnn + 5
      np6 = nnn + 6
      np7 = nnn + 7
      np8 = nnn + 8

      do k = 1, q
         dgrav(k) = zero
         tme  (k) = smltme
      enddo

c --- take one row from data arrays

      do k = 1, nnn
         k4 = k + 4
         rho (k4) = densty(i,j,k,igrd)
         u   (k4) = velz  (i,j,k,igrd)
         ut  (k4) = velx  (i,j,k,igrd)
         utt (k4) = vely  (i,j,k,igrd)
         e   (k4) = energy(i,j,k,igrd)
         grav(k4) = gpot  (i,j,k,igrd)
         ugridz(k,igrd) = zero
         tme (k4) = tmpent(i,j,k,igrd)
      enddo

      do m = 1, qc
      do k = 1, nnn
         k4 = k + 4
         che(k4,m) = chem(i,j,k,igrd,m)
      enddo
      enddo

      do k4 = 1, np8                ! non conforming to
         ugrid (k4) = zero
c     ugrid(k) = ugridz(lone,igrd)        ! anything else !!!
      enddo

c      if ( relaxon .eq. lone ) then
c         do k4 = 1, np8
c            urlx (k4) = relaxvz(i+4,j+4,k4,igrd)
c         enddo
c      endif

      call bndrwz( rho (1), bzdens(1,1,1), i, j )
      call bndrwz( u   (1), bzvelz(1,1,1), i, j )
      call bndrwz( ut  (1), bzvelx(1,1,1), i, j )
      call bndrwz( utt (1), bzvely(1,1,1), i, j )
      call bndrwz( e   (1), bzener(1,1,1), i, j )
      call bndrwz( grav(1), bzgpot(1,1,1), i, j )

      do m = 1, qc
         call bndrwz( che(1,m), bzchem(1,1,1,m), i, j )
      enddo

      do kk = 1, np8
         ek  (kk) =  half * ( u(kk)**2 + ut(kk)**2 + utt(kk)**2 )
         ei  (kk) = ( e(kk) - ek(kk) ) * rho(kk)
         ye  (kk) = che(kk,1) / rho(kk)
         dx  (kk) = delx(igrd)
         dtdx(kk) = dtgr(igrd) / delx(igrd)
      enddo

      do kk = 1, np8
         if(ei(kk).lt.zero) write(*,*) 'getrwz  ',i,j,kk,igrd
      enddo

      call eos( rho, tmp, ei, p, ye, gamc, game, np8, lzero )
c *** watch out: eos changes also ei and ye where necessary !!!

!      call nospike( rho, tmp, tme, ei, p, ye, gamc, game, np8 )

      do k = 1 , np8
         che(k,1) = ye(k) * rho(k)
         e(k) = ei(k)/rho(k) + ek(k)
c        v(k)  = one / rho(k)
         c(k)  = sqrt( gamc(k) * p(k) * rho(k) )
         ce(k) = c(k) / rho(k)
      enddo

      do k = 1, nnn
         k4 = k + 4
         temper(i,j,k,igrd) = tmp(k4)
         press (i,j,k,igrd) = p  (k4)
      enddo

      dtppm = dtgr(igrd)

      return
      end
c-----------------------------------------------------------------------
      subroutine bndrwz( row, bzplane, i, j )
c     add boundary values and put in /hydro/
c     PARALLEL OK: PRIVATE(row,i,j), SHARED(bzplane)

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 row(q), bzplane(qx, qy, qb)

      !print *, 'grids, bndrwz'

      row   (1) = bzplane (i,j,1)
      row   (2) = bzplane (i,j,2)
      row   (3) = bzplane (i,j,3)
      row   (4) = bzplane (i,j,4)
      row (np5) = bzplane (i,j,5)
      row (np6) = bzplane (i,j,6)
      row (np7) = bzplane (i,j,7)
      row (np8) = bzplane (i,j,8)

      return
      end
c-----------------------------------------------------------------------
      subroutine putrwz (i,j,igrd)
c     get /hydro/..nu values and put them in /vnew/
c     PARALLEL OK

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include 'compct.cmn'

      integer*4 i,j,igrd,k,k4,m,k4m,k4p
c     these must all be private

      do k = 1, nnn
         k4 = k + 4
         ek (k4) = half*(unu(k4)**2 + utnu(k4)**2 + uttnu(k4)**2)
      enddo

!      call lightning(unu(5),utnu(5),uttnu(5),enu(5),ek(5), nnn)

      do k = 1, nz
         k4 = k + 4
         densty(i,j,k,igrd) = rhonu(k4)
         velz  (i,j,k,igrd) = unu  (k4)
         velx  (i,j,k,igrd) = utnu (k4)
         vely  (i,j,k,igrd) = uttnu(k4)
         energy(i,j,k,igrd) = enu  (k4)
      enddo

      do m = 1, qc
      do k = 1, nz
         k4 = k + 4
         chem(i,j,k,igrd,m) = max( chenu(k4,m), smlche(m) )
      enddo
      enddo

      enu(4)   = max(enu(5)  -ek(5),   small*ek(5)  )
      enu(np5) = max(enu(np4)-ek(np4), small*ek(np4))
      ek (4)   = ek(5)
      ek (np5) = ek(np4)

      do k = 1, nz
         k4 = k + 4
         enu(k4) = enu(k4) - ek(k4)
      enddo

      do k = 1, nz
         k4 = k + 4
         k4m = max( k4-1, 11-k4 )
         k4p = min( k4+1, 2*nz+7-k4)
         if ( enu(k4) .le. small*ek(k4) ) then                          !DRW - ask Max about this check
            tmpp = enu(k4)                                              !!
            enu(k4) = half * ( enu(k4m) + max(enu(k4p),zero) )          !!
            energy(i,j,k,igrd) = max( enu(k4), small*ek(k4) ) + ek(k4)  !!
c!c        if( igrd.ne.indxgr(ltwo,maxlev(ltwo)) )                      !!
            if ( bhdens(i,j,k,igrd) .eq. zero )                         !!
     &         write(*,'(a,6i4,1P,4E19.10)')                            !!
     &             'putrwz, energy:',i,j,k,igrd,k4p,k4m,                !!
     &                energy(i,j,k,igrd), enu(k4), ek(k4),tmpp          !!
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine boundry (ifgrv,igrd,nxyzsw,icycl)
c     put boundary values around /vnew/

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include 'compct.cmn'

      logical overlap

      real*8 ek(q)

      !print *, 'grids, boundry'

      itgr = indxgr(lone,lone)
      if ( igrd .eq. itgr ) then
c --- is igrd the topmost level? Then use the bndm.. values.

c     bndm.. = 1  ====>  reflecting boundary
c     bndm.. = 2  ====>  flow out   boundary
c     bndm.. = 3  ====>  flow in    boundary
c     these are initialized in 'inippm',
c     but those that are not 'reflecting' are changed here again

c --------
c x boundary
c --------
      if (nxyzsw .eq. lone) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,j), SHARED(bndmnx,bndmxx,
C$OMP+                   velx,refl,refr,floutl,floutr,flinl,flinr,
C$OMP+                   igrd,nx,ny,nz)
      do j = 1, ny
      do k = 1, nz
         if ( bndmnx(j,k) .ne. one ) then
            if ( velx(1,j,k,igrd) .lt. zero ) then
               bndmnx(j,k) = two
            else
               bndmnx(j,k) = one
            endif
         endif
         if ( bndmxx(j,k) .ne. one ) then
            if ( velx(nx,j,k,igrd) .gt. zero ) then
               bndmxx(j,k) = two
            else
               bndmxx(j,k) = one
            endif
         endif
         refl  (j,k) = (half*bndmnx(j,k)-2.5D0)*bndmnx(j,k)+3.D0
         refr  (j,k) = (half*bndmxx(j,k)-2.5D0)*bndmxx(j,k)+3.D0
         floutl(j,k) = (-bndmnx(j,k)+4.D0)*bndmnx(j,k)-3.D0
         floutr(j,k) = (-bndmxx(j,k)+4.D0)*bndmxx(j,k)-3.D0
         flinl (j,k) = one-refl(j,k)-floutl(j,k)
         flinr (j,k) = one-refr(j,k)-floutr(j,k)
      enddo
      enddo
C$OMP END PARALLEL DO

      nnn = nx

      if (ifgrv .eq. lone) then
         call bndtox ( densty(1,1,1,itgr), rhoin, +one, bxdens(1,1,1) )
         call bndtox (   gpot(1,1,1,itgr),   gin, +one, bxgpot(1,1,1) )
         call bndtox (   velx(1,1,1,itgr),   uin, -one, bxvelx(1,1,1) )
         call bndtox (   vely(1,1,1,itgr),  utin, +one, bxvely(1,1,1) )
         call bndtox (   velz(1,1,1,itgr), uttin, +one, bxvelz(1,1,1) )
         call bndtox ( energy(1,1,1,itgr),   ein, +one, bxener(1,1,1) )
         do m = 1, qc
         call bndtox (chem(1,1,1,itgr,m),chein(m),+one, bxchem(1,1,1,m))
         enddo
      elseif ( ifgrv .eq. lzero ) then
c1       call bndtox (  gpold(1,1,1,itgr),   gin, +one, bxgold(1,1,1) )
         call bndtox (   gpot(1,1,1,itgr),   gin, +one, bxgpot(1,1,1) )
      elseif ( ifgrv .eq. ltwo ) then
         call bndtox (phireac(1,1,1     ),   gin, +one, bxvelx(1,1,1) )
      endif

c --------
c y boundary
c --------
      elseif (nxyzsw .eq. ltwo) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,i), SHARED(bndmny,bndmxy,
C$OMP+                   vely,refl,refr,floutl,floutr,flinl,flinr,
C$OMP+                   igrd,nx,ny,nz)
      do i = 1, nx
      do k = 1, nz
         if ( bndmny(i,k) .ne. one ) then
            if ( vely(i,1,k,igrd) .lt. zero ) then
               bndmny(i,k) = two
            else
               bndmny(i,k) = one
            endif
         endif
         if ( bndmxy(i,k) .ne. one ) then
            if ( vely(i,ny,k,igrd) .gt. zero ) then
               bndmxy(i,k) = two
            else
               bndmxy(i,k) = one
            endif
         endif
         refl  (i,k) = (half*bndmny(i,k)-2.5D0)*bndmny(i,k)+3.D0
         refr  (i,k) = (half*bndmxy(i,k)-2.5D0)*bndmxy(i,k)+3.D0
         floutl(i,k) = (-bndmny(i,k)+4.D0)*bndmny(i,k)-3.D0
         floutr(i,k) = (-bndmxy(i,k)+4.D0)*bndmxy(i,k)-3.D0
         flinl (i,k) = one-refl(i,k)-floutl(i,k)
         flinr (i,k) = one-refr(i,k)-floutr(i,k)
      enddo
      enddo
C$OMP END PARALLEL DO

      nnn = ny

      if (ifgrv .eq. lone) then
         call bndtoy ( densty(1,1,1,itgr), rhoin, +one, bydens(1,1,1) )
         call bndtoy (   gpot(1,1,1,itgr),   gin, +one, bygpot(1,1,1) )
         call bndtoy (   velx(1,1,1,itgr),   uin, +one, byvelx(1,1,1) )
         call bndtoy (   vely(1,1,1,itgr),  utin, -one, byvely(1,1,1) )
         call bndtoy (   velz(1,1,1,itgr), uttin, +one, byvelz(1,1,1) )
         call bndtoy ( energy(1,1,1,itgr),   ein, +one, byener(1,1,1) )
         do m = 1, qc
         call bndtoy (chem(1,1,1,itgr,m),chein(m),+one, bychem(1,1,1,m))
         enddo
      elseif ( ifgrv .eq. lzero ) then
c1       call bndtoy (  gpold(1,1,1,itgr),   gin, +one, bygold(1,1,1) )
         call bndtoy (   gpot(1,1,1,itgr),   gin, +one, bygpot(1,1,1) )
      elseif ( ifgrv .eq. ltwo ) then
         call bndtoy (phireac(1,1,1     ),   gin, +one, byvelx(1,1,1) )
      endif

c --------
c z boundary
c --------
      elseif (nxyzsw .eq. lthree) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j), SHARED(bndmnz,bndmxz,
C$OMP+                   velz,refl,refr,floutl,floutr,flinl,flinr,
C$OMP+                   igrd,nx,ny,nz)
      do j = 1, ny
      do i = 1, nx
         if ( bndmnz(i,j) .ne. one ) then
            if ( velz(i,j,1,igrd) .lt. zero ) then
               bndmnz(i,j) = two
            else
               bndmnz(i,j) = one
            endif
         endif
         if ( bndmxz(i,j) .ne. one ) then
            if ( velz(i,j,nz,igrd) .gt. zero ) then
               bndmxz(i,j) = two
            else
               bndmxz(i,j) = one
            endif
         endif
         refl  (i,j) = (half*bndmnz(i,j)-2.5D0)*bndmnz(i,j)+3.D0
         refr  (i,j) = (half*bndmxz(i,j)-2.5D0)*bndmxz(i,j)+3.D0
         floutl(i,j) = (-bndmnz(i,j)+4.D0)*bndmnz(i,j)-3.D0
         floutr(i,j) = (-bndmxz(i,j)+4.D0)*bndmxz(i,j)-3.D0
         flinl (i,j) = one-refl(i,j)-floutl(i,j)
         flinr (i,j) = one-refr(i,j)-floutr(i,j)
      enddo
      enddo
C$OMP END PARALLEL DO

      nnn = nz

      if (ifgrv .eq. lone) then
         call bndtoz ( densty(1,1,1,itgr), rhoin, +one, bzdens(1,1,1) )
         call bndtoz (   gpot(1,1,1,itgr),   gin, +one, bzgpot(1,1,1) )
         call bndtoz (   velx(1,1,1,itgr),   uin, +one, bzvelx(1,1,1) )
         call bndtoz (   vely(1,1,1,itgr),  utin, +one, bzvely(1,1,1) )
         call bndtoz (   velz(1,1,1,itgr), uttin, -one, bzvelz(1,1,1) )
         call bndtoz ( energy(1,1,1,itgr),   ein, +one, bzener(1,1,1) )
         do m = 1, qc
         call bndtoz (chem(1,1,1,itgr,m),chein(m),+one, bzchem(1,1,1,m))
         enddo
      elseif (ifgrv .eq. lzero ) then
c1       call bndtoz (  gpold(1,1,1,itgr),   gin, +one, bzgold(1,1,1) )
         call bndtoz (   gpot(1,1,1,itgr),   gin, +one, bzgpot(1,1,1) )
      elseif (ifgrv .eq. ltwo ) then
         call bndtoz (phireac(1,1,1     ),   gin, +one, bzvelx(1,1,1) )
      endif

      endif
c --- of nxyzsw

      else
c --------
c --- of special case for topmost level
c --------

c --- igrd lies within which coarser grid icog?

      icog = idxcgr(igrd)


c --- get values from coarser grid (level-1) and interpolate
c --- use /vold/ values because coarser grid is already updated.


      if (nxyzsw .eq. lone) then
c ------- x boundary

      nnn = nx

      if ( ifgrv .ne. ltwo ) then
      if ( icycl .eq. lone ) then
         call btrialm ( lone, igrd, icog,-lthree, lzero, ifgrv )
         call btrialm ( lone, igrd, icog,  nnn+1, nnn+4, ifgrv )
      else
         call btrialn ( lone, igrd, icog,-lthree, lzero, ifgrv )
         call btrialn ( lone, igrd, icog,  nnn+1, nnn+4, ifgrv )
      endif
      endif

      if (ifgrv .eq. lone) then
      do i = 1, qb
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,j,ek),
C$OMP+            SHARED(i,bxvelx,bxvely,bxvelz,bxener,small,
C$OMP+            igrd,ny,nz)
      do j = 1, ny
      do k = 1, nz
         ek(j) = (one+small) * half *
     &      ( bxvelx(i,j,k)**2 + bxvely(i,j,k)**2 + bxvelz(i,j,k)**2 )
         if ( bxener(i,j,k) .lt. ek(j) ) then
             write(*,*) 'boundry,x,  energy:',i,j,k,igrd,
     &                   bxener(i,j,k),ek(j)
              stop 'boundry'
         endif
         bxener(i,j,k) = max( bxener(i,j,k), ek(j) )
      enddo
      enddo
C$OMP END PARALLEL DO
      enddo
      endif

      elseif (nxyzsw .eq. ltwo) then
c ------- y boundary

      nnn = ny

      if ( ifgrv .ne. ltwo ) then
      if ( icycl .eq. lone ) then
         call btrialm ( ltwo, igrd, icog,-lthree, lzero, ifgrv )
         call btrialm ( ltwo, igrd, icog,  nnn+1, nnn+4, ifgrv )
      else
         call btrialn ( ltwo, igrd, icog,-lthree, lzero, ifgrv )
         call btrialn ( ltwo, igrd, icog,  nnn+1, nnn+4, ifgrv )
      endif
      endif

      if (ifgrv .eq. lone) then
      do j = 1, qb
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,i,ek),
C$OMP+            SHARED(j,byvelx,byvely,byvelz,byener,small,
C$OMP+            igrd,nx,nz)
      do i = 1, nx
      do k = 1, nz
         ek(i) = (one+small) * half *
     &      ( byvelx(i,j,k)**2 + byvely(i,j,k)**2 + byvelz(i,j,k)**2 )
         if ( byener(i,j,k) .lt. ek(i) )
     &        write(*,*) 'boundry,y,  energy:',i,j,k,igrd,
     &                   byener(i,j,k),ek(i)
         byener(i,j,k) = max( byener(i,j,k), ek(i) )
      enddo
      enddo
C$OMP END PARALLEL DO
      enddo
      endif

      elseif (nxyzsw .eq. lthree) then
c ------- z boundary

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j),
C$OMP+      SHARED(refl,refr,floutl,floutr,flinl,flinr,bndmnz,bndmxz,
C$OMP+             igrd,nx,ny)
      do j = 1, ny
      do i = 1, nx
         refl  (i,j) = (0.5D0*bndmnz(i,j)-2.5D0)*bndmnz(i,j)+3.D0
         refr  (i,j) = (0.5D0*bndmxz(i,j)-2.5D0)*bndmxz(i,j)+3.D0
         floutl(i,j) = (-bndmnz(i,j)+4.D0)*bndmnz(i,j)-3.D0
         floutr(i,j) = (-bndmxz(i,j)+4.D0)*bndmxz(i,j)-3.D0
         flinl (i,j) = one-refl(i,j)-floutl(i,j)
         flinr (i,j) = one-refr(i,j)-floutr(i,j)
      enddo
      enddo
C$OMP END PARALLEL DO

      nnn = nz

c     reflecting z-boundry at z=0.5D0
      if (ifgrv .eq. lone) then
         call bndtoz ( densty(1,1,1,igrd), rhoin, +one, bzdens(1,1,1) )
         call bndtoz (   gpot(1,1,1,igrd),   gin, +one, bzgpot(1,1,1) )
         call bndtoz (   velx(1,1,1,igrd),   uin, +one, bzvelx(1,1,1) )
         call bndtoz (   vely(1,1,1,igrd),  utin, +one, bzvely(1,1,1) )
         call bndtoz (   velz(1,1,1,igrd), uttin, -one, bzvelz(1,1,1) )
         call bndtoz ( energy(1,1,1,igrd),   ein, +one, bzener(1,1,1) )
         do m = 1, qc
         call bndtoz (chem(1,1,1,igrd,m),chein(m),+one, bzchem(1,1,1,m))
         enddo
      elseif (ifgrv.ne.ltwo) then
c1       call bndtoz (  gpold(1,1,1,igrd),   gin, +one, bzgold(1,1,1) )
         call bndtoz (   gpot(1,1,1,igrd),   gin, +one, bzgpot(1,1,1) )
      endif

      if ( ifgrv .ne. ltwo ) then
      if ( icycl .eq. lone ) then
c        call btrialm ( lthree, igrd, icog,-lthree, lzero, ifgrv )
         call btrialm ( lthree, igrd, icog,  nnn+1, nnn+4, ifgrv )
      else
c        call btrialn ( lthree, igrd, icog,-lthree, lzero, ifgrv )
         call btrialn ( lthree, igrd, icog,  nnn+1, nnn+4, ifgrv )
      endif
      endif

      if (ifgrv .eq. lone) then
      do k = 5, qb
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,i,ek),
C$OMP+            SHARED(k,bzvelx,bzvely,bzvelz,bzener,small,
C$OMP+                   igrd,nx,ny)
      do j = 1, ny
      do i = 1, nx
         ek(i) = (one+small) * half *
     &      ( bzvelx(i,j,k)**2 + bzvely(i,j,k)**2 + bzvelz(i,j,k)**2 )
         if ( bzener(i,j,k) .lt. ek(i) )
     &        write(*,*) 'boundry,z,  energy:',i,j,k,igrd,
     &                   bzener(i,j,k),ek(i)
         bzener(i,j,k) = max( bzener(i,j,k), ek(i) )
      enddo
      enddo
C$OMP END PARALLEL DO
      enddo
      endif

      endif
c --- of nxyzsw

      if (nxyzsw .ne. lthree) then
      if (ifgrv .eq. lone) then
c --- take values from grid of equal level, if present

      levl  = levlgr(igrd,lone)
      numgr = indxgr(lzero,levl)
      ncol  = levlgr(igrd,ltwo)

      if ( numgr .gt. lone ) then

         do inum = 1, ncol-1
         igr2 = indxgr(inum,levl)
         if (overlap(igrd,igr2)) then
            call takval ( bxdens, bydens,   dold(1,1,1,igr2) )
c1          call takval ( bxgpot, bygpot,  gpold(1,1,1,igr2) )
            call takval ( bxvelx, byvelx,  vxold(1,1,1,igr2) )
            call takval ( bxvely, byvely,  vyold(1,1,1,igr2) )
            call takval ( bxvelz, byvelz,  vzold(1,1,1,igr2) )
            call takval ( bxener, byener,   eold(1,1,1,igr2) )
         endif
         enddo

         do inum = ncol+1, numgr
         igr2 = indxgr(inum,levl)
         if (overlap(igrd,igr2)) then
            call takval ( bxdens, bydens, densty(1,1,1,igr2) )
            call takval ( bxgpot, bygpot,   gpot(1,1,1,igr2) )
            call takval ( bxvelx, byvelx,   velx(1,1,1,igr2) )
            call takval ( bxvely, byvely,   vely(1,1,1,igr2) )
            call takval ( bxvelz, byvelz,   velz(1,1,1,igr2) )
            call takval ( bxener, byener, energy(1,1,1,igr2) )
         endif
         enddo

      endif
c --- of values of equal level
      endif

      endif
c --- of nxyzsw

c --------
      endif
c --------
c --- of other than topmost level


      return
      end
c-----------------------------------------------------------------------
      subroutine bndtox (qty,qin,sg,bqty)

c     put boundary values around /vnew/
c     x-boundry only for topmost level (igrd = 1)

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 qty( qx, qy, qz ),  bqty(qb, qy, qz)

      !print *, 'grids, bndtox'

      nm1 = nnn - lone
      nm2 = nnn - ltwo
      nm3 = nnn - lthree

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,j), SHARED(qty,bqty,
C$OMP+                       refl,refr,floutl,floutr,flinr,flinl,
C$OMP+                       qin,sg,ny,nz,nnn,nm1,nm2,nm3)
      do j = 1, ny
      do k = 1, nz

      bqty(1,j,k)= sg * refl(j,k) * qty(  4,j,k)
     &            + floutl(j,k)*qty(  1,j,k) + flinl(j,k)*qin
      bqty(2,j,k)= sg * refl(j,k) * qty(  3,j,k)
     &            + floutl(j,k)*qty(  1,j,k) + flinl(j,k)*qin
      bqty(3,j,k)= sg * refl(j,k) * qty(  2,j,k)
     &            + floutl(j,k)*qty(  1,j,k) + flinl(j,k)*qin
      bqty(4,j,k)= sg * refl(j,k) * qty(  1,j,k)
     &            + floutl(j,k)*qty(  1,j,k) + flinl(j,k)*qin
      bqty(5,j,k)= sg * refr(j,k) * qty(nnn,j,k)
     &            + floutr(j,k)*qty(nnn,j,k) + flinr(j,k)*qin
      bqty(6,j,k)= sg * refr(j,k) * qty(nm1,j,k)
     &            + floutr(j,k)*qty(nnn,j,k) + flinr(j,k)*qin
      bqty(7,j,k)= sg * refr(j,k) * qty(nm2,j,k)
     &            + floutr(j,k)*qty(nnn,j,k) + flinr(j,k)*qin
      bqty(8,j,k)= sg * refr(j,k) * qty(nm3,j,k)
     &            + floutr(j,k)*qty(nnn,j,k) + flinr(j,k)*qin

      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine bndtoy (qty,qin,sg,bqty)

c     put boundary values around /vnew/
c     y-boundry only for topmost level (igrd = 1)

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 qty( qx, qy, qz ),  bqty(qx, qb, qz)

      !print *, 'grids, bndtoy'

      nm1 = nnn - lone
      nm2 = nnn - ltwo
      nm3 = nnn - lthree

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,i), SHARED(qty,bqty,
C$OMP+                      refl,refr,floutl,floutr,flinr,flinl,
C$OMP+                      qin,sg,nx,nz,nnn,nm1,nm2,nm3)
      do i = 1, nx
      do k = 1, nz

      bqty(i,1,k)= sg * refl(i,k) * qty(i,  4,k)
     &            + floutl(i,k)*qty(i,  1,k) + flinl(i,k)*qin
      bqty(i,2,k)= sg * refl(i,k) * qty(i,  3,k)
     &            + floutl(i,k)*qty(i,  1,k) + flinl(i,k)*qin
      bqty(i,3,k)= sg * refl(i,k) * qty(i,  2,k)
     &            + floutl(i,k)*qty(i,  1,k) + flinl(i,k)*qin
      bqty(i,4,k)= sg * refl(i,k) * qty(i,  1,k)
     &            + floutl(i,k)*qty(i,  1,k) + flinl(i,k)*qin
      bqty(i,5,k)= sg * refr(i,k) * qty(i,nnn,k)
     &            + floutr(i,k)*qty(i,nnn,k) + flinr(i,k)*qin
      bqty(i,6,k)= sg * refr(i,k) * qty(i,nm1,k)
     &            + floutr(i,k)*qty(i,nnn,k) + flinr(i,k)*qin
      bqty(i,7,k)= sg * refr(i,k) * qty(i,nm2,k)
     &            + floutr(i,k)*qty(i,nnn,k) + flinr(i,k)*qin
      bqty(i,8,k)= sg * refr(i,k) * qty(i,nm3,k)
     &            + floutr(i,k)*qty(i,nnn,k) + flinr(i,k)*qin

      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine bndtoz (qty,qin,sg,bqty)

c     put boundary values around /vnew/
c     z-boundry for all levels

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 qty( qx, qy, qz ),  bqty(qx, qy, qb)

      !print *, 'grids, bndtoz'

      nm1 = nnn - lone
      nm2 = nnn - ltwo
      nm3 = nnn - lthree

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,i), SHARED(qty,bqty,
C$OMP+                      refl,refr,floutl,floutr,flinr,flinl,
C$OMP+                      qin,sg,nx,ny,nnn,nm1,nm2,nm3)
      do j = 1, ny
      do i = 1, nx

      bqty(i,j,1)= sg * refl(i,j) * qty(i,j,  4)
     &            + floutl(i,j)*qty(i,j,  1) + flinl(i,j)*qin
      bqty(i,j,2)= sg * refl(i,j) * qty(i,j,  3)
     &            + floutl(i,j)*qty(i,j,  1) + flinl(i,j)*qin
      bqty(i,j,3)= sg * refl(i,j) * qty(i,j,  2)
     &            + floutl(i,j)*qty(i,j,  1) + flinl(i,j)*qin
      bqty(i,j,4)= sg * refl(i,j) * qty(i,j,  1)
     &            + floutl(i,j)*qty(i,j,  1) + flinl(i,j)*qin
      bqty(i,j,5)= sg * refr(i,j) * qty(i,j,nnn)
     &            + floutr(i,j)*qty(i,j,nnn) + flinr(i,j)*qin
      bqty(i,j,6)= sg * refr(i,j) * qty(i,j,nm1)
     &            + floutr(i,j)*qty(i,j,nnn) + flinr(i,j)*qin
      bqty(i,j,7)= sg * refr(i,j) * qty(i,j,nm2)
     &            + floutr(i,j)*qty(i,j,nnn) + flinr(i,j)*qin
      bqty(i,j,8)= sg * refr(i,j) * qty(i,j,nm3)
     &            + floutr(i,j)*qty(i,j,nnn) + flinr(i,j)*qin

      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine bxther(be, bvx, bvy, bvz, nb, ne, nn)

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 be(nn,qy,qz), bvx(nn,qy,qz), bvy(nn,qy,qz), bvz(nn,qy,qz)

      !print *, 'grids, bxther'

      np = lzero
      if (ne.eq. lzero)   np = qb/ltwo
      if (nb.eq.nx+lone)  np = qb/ltwo - nx
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,i,k),
C$OMP+            SHARED(be,bvx,bvy,bvz,nb,ne,np)
      do j = 1, qy
      do i = nb+np, ne+np
      do k = 1, qz
         be(i,j,k) = be(i,j,k) + half *
     &       ( bvx(i,j,k)**2 + bvy(i,j,k)**2 + bvz(i,j,k)**2 )
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine byther(be, bvx, bvy, bvz, nb, ne, nn)

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 be(qx,nn,qz), bvx(qx,nn,qz), bvy(qx,nn,qz), bvz(qx,nn,qz)

      !print *, 'grids, byther'

      np = lzero
      if (ne.eq. lzero)   np = qb/ltwo
      if (nb.eq.ny+lone)  np = qb/ltwo - ny
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(be,bvx,bvy,bvz,nb,ne,np)
      do i = 1, qx
      do k = 1, qz
      do j = nb+np, ne+np
         be(i,j,k) = be(i,j,k) + half *
     &       ( bvx(i,j,k)**2 + bvy(i,j,k)**2 + bvz(i,j,k)**2 )
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine bzther(be, bvx, bvy, bvz, nb, ne, nn)

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 be(qx,qy,nn), bvx(qx,qy,nn), bvy(qx,qy,nn), bvz(qx,qy,nn)

      !print *, 'grids, bzther'

      np = lzero
      if (ne.eq. lzero)   np = qb/ltwo
      if (nb.eq.nz+lone)  np = qb/ltwo - nz
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i),
C$OMP+            SHARED(be,bvx,bvy,bvz,nb,ne,np)
      do j = 1, qy
      do k = nb+np, ne+np
      do i = 1, qx
         be(i,j,k) = be(i,j,k) + half *
     &       ( bvx(i,j,k)**2 + bvy(i,j,k)**2 + bvz(i,j,k)**2 )
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine etherm(erg, vx, vy, vz)

c     fine internal energy and put into ro

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8 erg(qx,qy,qz), vx(qx,qy,qz), vy(qx,qy,qz), vz(qx,qy,qz)

      !print *, 'grids, etherm'

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i),
C$OMP+            SHARED(ro,erg,vx,vy,vz,nx,ny,nz)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         ro(i,j,k) = erg(i,j,k) - half *
     &               ( vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2 )
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine btrialm(ixyz,igrd,icog,ifpl,ifpr,ifgrv)

c     call tri x,y,z, for all physical quantities

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8  dem(qx,qy,qz), gpm(qx,qy,qz), enm(qx,qy,qz),
     &        vxm(qx,qy,qz), vym(qx,qy,qz), vzm(qx,qy,qz)
      common /raypot/ dem, gpm, enm, vxm, vym, vzm

c calculate time average of old and new values
c this is done for all cells; waste!

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i),
C$OMP+            SHARED(dem,gpm,enm,vxm,vym,vzm,
C$OMP+                   densty,gpot,energy,velx,vely,velz,
C$OMP+                   dold,eold,vxold,vyold,vzold,
C$OMP+                   nx,ny,nz,icog)

      !print *, 'grids, btrialm'

      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
        dem(i,j,k) = half * ( densty(i,j,k,icog) +  dold(i,j,k,icog) )
        gpm(i,j,k) = half * (   gpot(i,j,k,icog) +  gpot(i,j,k,icog) )
c       gpm(i,j,k) = half * (   gpot(i,j,k,icog) + gpold(i,j,k,icog) )
        enm(i,j,k) = half * ( energy(i,j,k,icog) +  eold(i,j,k,icog) )
        vxm(i,j,k) = half * (   velx(i,j,k,icog) + vxold(i,j,k,icog) )
        vym(i,j,k) = half * (   vely(i,j,k,icog) + vyold(i,j,k,icog) )
        vzm(i,j,k) = half * (   velz(i,j,k,icog) + vzold(i,j,k,icog) )
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c -----
      if (ixyz .eq. lone) then
c -----
      if (ifgrv .eq. lone) then

         call trix (bxdens(1,1,1), dem(1,1,1),
     &              bxdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call trix (bxgpot(1,1,1), gpm(1,1,1),
     &              bxdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call etherm( enm(1,1,1), vxm(1,1,1), vym(1,1,1), vzm(1,1,1) )

         call trix (bxener(1,1,1),  ro(1,1,1),
     &              bxdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr,  lone , qb, +one )

         call trix (bxvelx(1,1,1), vxm(1,1,1),
     &              bxdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call trix (bxvely(1,1,1), vym(1,1,1),
     &              bxdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call trix (bxvelz(1,1,1), vzm(1,1,1),
     &              bxdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, -one )

         call bxther( bxener(1,1,1), bxvelx(1,1,1),
     &                bxvely(1,1,1), bxvelz(1,1,1),
     &                ifpl, ifpr, qb )

         do m = 1, qc
         call trix (bxchem(1,1,1,m),  chem(1,1,1,icog,m),
     &              bxdens(1,1,1),  densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )
         enddo

      else

         call trix (bxgpot(1,1,1), gpm(1,1,1),
     &              bxdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

c1       call trix (bxgold(1,1,1), gpm(1,1,1),
c1   &              bxdens(1,1,1), dem(1,1,1),
c1   &               igrd, ifpl, ifpr, lzero , qb, +one )

         write(*,*) 'something wrong in btrialm'
         stop'=>btrialm'
      endif

c -----
      else if (ixyz .eq. ltwo) then
c -----
      if (ifgrv .eq. lone) then

         call triy (bydens(1,1,1), dem(1,1,1),
     &              bydens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triy (bygpot(1,1,1), gpm(1,1,1),
     &              bydens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call etherm( enm(1,1,1), vxm(1,1,1), vym(1,1,1), vzm(1,1,1) )

         call triy (byener(1,1,1),  ro(1,1,1),
     &              bydens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr,  lone , qb, +one )

         call triy (byvelx(1,1,1), vxm(1,1,1),
     &              bydens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triy (byvely(1,1,1), vym(1,1,1),
     &              bydens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triy (byvelz(1,1,1), vzm(1,1,1),
     &              bydens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, -one )

         call byther( byener(1,1,1), byvelx(1,1,1),
     &                byvely(1,1,1), byvelz(1,1,1),
     &                ifpl, ifpr, qb )

         do m = 1, qc
         call triy (bychem(1,1,1,m),  chem(1,1,1,icog,m),
     &              bydens(1,1,1),  densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )
         enddo

      else

         call triy (bygpot(1,1,1), gpm(1,1,1),
     &              bydens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

c1       call triy (bygold(1,1,1), gpm(1,1,1),
c1   &              bydens(1,1,1), dem(1,1,1),
c1   &               igrd, ifpl, ifpr, lzero , qb, +one )

         write(*,*) 'something wrong in btrialm'
         stop'=>btrialm'
      endif

c -----
      else if (ixyz .eq. lthree) then
c -----
      if (ifgrv .eq. lone) then

         call triz (bzdens(1,1,1), dem(1,1,1),
     &              bzdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triz (bzgpot(1,1,1), gpm(1,1,1),
     &              bzdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call etherm( enm(1,1,1), vxm(1,1,1), vym(1,1,1), vzm(1,1,1) )

         call triz (bzener(1,1,1),  ro(1,1,1),
     &              bzdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr,  lone , qb, +one )

         call triz (bzvelx(1,1,1), vxm(1,1,1),
     &              bzdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triz (bzvely(1,1,1), vym(1,1,1),
     &              bzdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triz (bzvelz(1,1,1), vzm(1,1,1),
     &              bzdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, -one )

         call bzther( bzener(1,1,1), bzvelx(1,1,1),
     &                bzvely(1,1,1), bzvelz(1,1,1),
     &                ifpl, ifpr, qb )

         do m = 1, qc
         call triz (bzchem(1,1,1,m),  chem(1,1,1,icog,m),
     &              bzdens(1,1,1),  densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )
         enddo

      else

         call triz (bzgpot(1,1,1), gpm(1,1,1),
     &              bzdens(1,1,1), dem(1,1,1),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

c1       call triz (bzgold(1,1,1), gpm(1,1,1),
c1   &              bzdens(1,1,1), dem(1,1,1),
c1   &               igrd, ifpl, ifpr, lzero , qb, +one )

         write(*,*) 'something wrong in btrialm'
         stop'=>btrialm'
      endif

      else

         write(*,*) 'something wrong in btrialm'
         stop'=>btrialm'

      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine btrialn(ixyz,igrd,icog,ifpl,ifpr,ifgrv)

c     call tri x,y,z, for all physical quantities

      include 'qparam.cmn'
      include 'grnest.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'bounds.cmn'

      !print *, 'grids, btrialn'

c -----
      if (ixyz .eq. lone) then
c -----
      if (ifgrv .eq. lone) then

         call trix (bxdens(1,1,1),densty(1,1,1,icog),
     &              bxdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call trix (bxgpot(1,1,1),  gpot(1,1,1,icog),
     &              bxdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call etherm(energy(1,1,1,icog),  velx(1,1,1,icog),
     &                 vely(1,1,1,icog),  velz(1,1,1,icog) )

         call trix (bxener(1,1,1),  ro  (1,1,1     ),
     &              bxdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr,  lone , qb, +one )

         call trix (bxvelx(1,1,1),  velx(1,1,1,icog),
     &              bxdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call trix (bxvely(1,1,1),  vely(1,1,1,icog),
     &              bxdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call trix (bxvelz(1,1,1),  velz(1,1,1,icog),
     &              bxdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, -one )

         call bxther( bxener(1,1,1), bxvelx(1,1,1),
     &                bxvely(1,1,1), bxvelz(1,1,1),
     &                ifpl, ifpr, qb )

         do m = 1, qc
         call trix (bxchem(1,1,1,m),  chem(1,1,1,icog,m),
     &              bxdens(1,1,1),  densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )
         enddo

      else

         call trix (bxgpot(1,1,1),  gpot(1,1,1,icog),
     &              bxdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

c1       call trix (bxgold(1,1,1), gpold(1,1,1,icog),
c1   &              bxdens(1,1,1),densty(1,1,1,icog),
c1   &               igrd, ifpl, ifpr, lzero , qb, +one )

      endif

c -----
      else if (ixyz .eq. ltwo) then
c -----
      if (ifgrv .eq. lone) then

         call triy (bydens(1,1,1),densty(1,1,1,icog),
     &              bydens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triy (bygpot(1,1,1),  gpot(1,1,1,icog),
     &              bydens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call etherm(energy(1,1,1,icog),  velx(1,1,1,icog),
     &                 vely(1,1,1,icog),  velz(1,1,1,icog) )

         call triy (byener(1,1,1),  ro  (1,1,1     ),
     &              bydens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr,  lone , qb, +one )

         call triy (byvelx(1,1,1),  velx(1,1,1,icog),
     &              bydens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triy (byvely(1,1,1),  vely(1,1,1,icog),
     &              bydens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triy (byvelz(1,1,1),  velz(1,1,1,icog),
     &              bydens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, -one )

         call byther( byener(1,1,1), byvelx(1,1,1),
     &                byvely(1,1,1), byvelz(1,1,1),
     &                ifpl, ifpr, qb )

         do m = 1, qc
         call triy (bychem(1,1,1,m),  chem(1,1,1,icog,m),
     &              bydens(1,1,1),  densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )
         enddo

      else

         call triy (bygpot(1,1,1),  gpot(1,1,1,icog),
     &              bydens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

c1       call triy (bygold(1,1,1), gpold(1,1,1,icog),
c1   &              bydens(1,1,1),densty(1,1,1,icog),
c1   &               igrd, ifpl, ifpr, lzero , qb, +one )

      endif

c -----
      else if (ixyz .eq. lthree) then
c -----
      if (ifgrv .eq. lone) then

         call triz (bzdens(1,1,1),densty(1,1,1,icog),
     &              bzdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triz (bzgpot(1,1,1),  gpot(1,1,1,icog),
     &              bzdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call etherm(energy(1,1,1,icog),  velx(1,1,1,icog),
     &                 vely(1,1,1,icog),  velz(1,1,1,icog) )

         call triz (bzener(1,1,1),  ro  (1,1,1     ),
     &              bzdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr,  lone , qb, +one )

         call triz (bzvelx(1,1,1),  velx(1,1,1,icog),
     &              bzdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triz (bzvely(1,1,1),  vely(1,1,1,icog),
     &              bzdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

         call triz (bzvelz(1,1,1),  velz(1,1,1,icog),
     &              bzdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, -one )

         call bzther( bzener(1,1,1), bzvelx(1,1,1),
     &                bzvely(1,1,1), bzvelz(1,1,1),
     &                ifpl, ifpr, qb )

         do m = 1, qc
         call triz (bzchem(1,1,1,m),  chem(1,1,1,icog,m),
     &              bzdens(1,1,1),  densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )
         enddo

      else

         call triz (bzgpot(1,1,1),  gpot(1,1,1,icog),
     &              bzdens(1,1,1),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qb, +one )

c1       call triz (bzgold(1,1,1), gpold(1,1,1,icog),
c1   &              bzdens(1,1,1),densty(1,1,1,icog),
c1   &               igrd, ifpl, ifpr, lzero , qb, +one )

      endif

      else

         write(*,*) 'something wrong in btrialn'
         stop'=>btrian'

      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine triall(ixyz,igrd,icog,ifpl,ifpr)

c     call tri x,y,z, for all physical quantities

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'

      !print *, 'grids, triall'

c -----
      if (ixyz .eq. lone) then
c -----

         call trix (densty(1,1,1,igrd),densty(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qx, +one )

         call trix (  gpot(1,1,1,igrd),  gpot(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qx, +one )

         call etherm(energy(1,1,1,icog),  velx(1,1,1,icog),
     &                 vely(1,1,1,icog),  velz(1,1,1,icog) )

         call trix (energy(1,1,1,igrd),  ro  (1,1,1     ),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr,  lone , qx, +one )

         call trix (  velx(1,1,1,igrd),  velx(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qx, +one )

         call trix (  vely(1,1,1,igrd),  vely(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qx, +one )

         call trix (  velz(1,1,1,igrd),  velz(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qx, -one )

         call bxther( energy(1,1,1,igrd), velx(1,1,1,igrd),
     &                vely  (1,1,1,igrd), velz(1,1,1,igrd),
     &                ifpl, ifpr, qx )

         do m = 1, qc
         call trix (  chem(1,1,1,igrd,m),  chem(1,1,1,icog,m),
     &              densty(1,1,1,igrd),  densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qx, +one )
         enddo

c -----
      else if (ixyz .eq. ltwo) then
c -----

         call triy (densty(1,1,1,igrd),densty(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qy, +one )

         call triy (  gpot(1,1,1,igrd),  gpot(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qy, +one )


         call etherm(energy(1,1,1,icog),  velx(1,1,1,icog),
     &                 vely(1,1,1,icog),  velz(1,1,1,icog) )

         call triy (energy(1,1,1,igrd),  ro  (1,1,1     ),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr,  lone , qy, +one )

         call triy (  velx(1,1,1,igrd),  velx(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qy, +one )

         call triy (  vely(1,1,1,igrd),  vely(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qy, +one )

         call triy (  velz(1,1,1,igrd),  velz(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qy, -one )

         call byther( energy(1,1,1,igrd), velx(1,1,1,igrd),
     &                vely  (1,1,1,igrd), velz(1,1,1,igrd),
     &                ifpl, ifpr, qy )

         do m = 1, qc
         call triy (  chem(1,1,1,igrd,m),  chem(1,1,1,icog,m),
     &              densty(1,1,1,igrd),  densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qy, +one )
         enddo

c -----
      else if (ixyz .eq. lthree) then
c -----

         call triz (densty(1,1,1,igrd),densty(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qz, +one )

         call triz (  gpot(1,1,1,igrd),  gpot(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qz, +one )

         call etherm(energy(1,1,1,icog),  velx(1,1,1,icog),
     &                 vely(1,1,1,icog),  velz(1,1,1,icog) )

         call triz (energy(1,1,1,igrd),  ro  (1,1,1     ),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr,  lone , qz, +one )

         call triz (  velx(1,1,1,igrd),  velx(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qz, +one )

         call triz (  vely(1,1,1,igrd),  vely(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qz, +one )

         call triz (  velz(1,1,1,igrd),  velz(1,1,1,icog),
     &              densty(1,1,1,igrd),densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qz, -one )

         call bzther( energy(1,1,1,igrd), velx(1,1,1,igrd),
     &                vely  (1,1,1,igrd), velz(1,1,1,igrd),
     &                ifpl, ifpr, qz )

         do m = 1, qc
         call triz (  chem(1,1,1,igrd,m),  chem(1,1,1,icog,m),
     &              densty(1,1,1,igrd),  densty(1,1,1,icog),
     &               igrd, ifpl, ifpr, lzero , qz, +one )
         enddo

      else

         write(*,*) 'something wrong in triall'
         stop'=> trial'

      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine trix(qtyfin,qtycor,dfin,dcor,igrd,ifpl,ifpr,ifd,nbx,sg)

c     triquadratic interpolation of qtycor onto qtyfin
c     y-z-slices from icpl to icpr on coarse grid
c     this version works only for nfine = 2 !

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8   qtycor(  qx, qy, qz ),  dcor(  qx, qy, qz ),
     &         qtyfin( nbx, qy, qz ),  dfin( nbx, qy, qz )
      integer*4 ifd

      !print *, 'grids, trix'

      nox = norgin(igrd,lone)
      noy = norgin(igrd,ltwo)
      noz = norgin(igrd,lthree)
      noq = noy
      nq = ny

      ione = lone
      if (ifpr .eq. lzero) ione = lzero

      do if = ifpl, ifpr
         ic = nox + (if+ione)/nfine
         icp = ic - ltwo*abs(mod(if,ltwo)) + lone
c        icp = ic - ltwo*iand(if,lone) + lone
         icm = ic - (icp-ic)
         ifp = if
         if (ifpr.eq.lzero) ifp = ifp + qb/2
         if (ifpl.eq.nx+1)  ifp = ifp - nx + qb/2

C$OMP PARALLEL DO DEFAULT(NONE), ! COPYIN(ic,icp,icm),
C$OMP+            PRIVATE(kf,jf),
c#ifdef SGI
C$OMP+            PRIVATE(kc,kcm,kcp,qty,qtymin,qtymax,sgnp,sgnm),
c#endif
C$OMP+            SHARED(nx,ny,nz,nox,noy,noz,noq,nq,sg,
C$OMP+                   qtycor,qtyfin,dcor,dfin,ifd,ifp)
      do kf = 1, nz
         kc = noz + (kf+1)/nfine
c        kcp = kc - ltwo*abs(mod(kf,ltwo)) + lone
         kcp = kc - ltwo*iand(kf,lone) + lone
         kcm = kc - (kcp-kc)
         sgnp = +one
         sgnm = +one
         if (kcp .lt. lone) then      ! reflecting boundary !
             kcp = lone
             sgnp = sg
         endif
         if (kcm .lt. lone) then      ! reflecting boundary !
             kcm = lone
             sgnm = sg
         endif

         call qintrj( qtycor, dcor, ifd )

         if (ifd .eq. lone) then
            do 60 jf = 1, ny
               qtyfin(ifp,jf,kf) =
     &            max( qtymin(jf), min( qtymax(jf),
     &                 qty(jf,1) / dfin(ifp,jf,kf)    ))
 60      continue
         else
            do 61 jf = 1, ny
 61            qtyfin(ifp,jf,kf) = qty(jf,1)
         endif

      enddo
C$OMP END PARALLEL DO
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine triy(qtyfin,qtycor,dfin,dcor,igrd,jfpl,jfpr,ifd,nby,sg)

c     triquadratic interpolation of qtycor onto qtyfin
c     x-z-slices from icpl to icbr on coarse grid
c     this version works only for nfine = 2 !

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8   qtycor( qx, qy, qz ),  dcor( qx, qy, qz ),
     &         qtyfin( qx,nby, qz ),  dfin( qx,nby, qz )
      integer*4 ifd

      !print *, 'grids, triy'

      nox = norgin(igrd,lone)
      noy = norgin(igrd,ltwo)
      noz = norgin(igrd,lthree)
      noq = nox
      nq = nx

      jone = lone
      if (jfpr .eq. lzero) jone = lzero

      do jf = jfpl, jfpr
         jc = noy + (jf+jone)/nfine
         jcp = jc - ltwo*iand(jf,lone) + lone
         jcm = jc - (jcp-jc)
         jfp = jf
         if (jfpr.eq.lzero) jfp = jfp + qb/2
         if (jfpl.eq.ny+1)  jfp = jfp - ny + qb/2

C$OMP PARALLEL DO DEFAULT(NONE), ! COPYIN(jc,jcp,jcm),
C$OMP+            PRIVATE(kf,if),
c#ifdef SGI
C$OMP+            PRIVATE(kc,kcm,kcp,qty,qtymin,qtymax,sgnp,sgnm),
c#endif
C$OMP+            SHARED(nx,ny,nz,nox,noy,noz,noq,nq,sg,
C$OMP+                   qtycor,qtyfin,dcor,dfin,ifd,jfp)
      do kf = 1, nz
         kc = noz + (kf+1)/nfine
         kcp = kc - ltwo*iand(kf,lone) + lone
         kcm = kc - (kcp-kc)
         sgnp = +one
         sgnm = +one
         if (kcp .lt. lone) then      ! reflecting boundary !
             kcp = lone
             sgnp = sg
         endif
         if (kcm .lt. lone) then      ! reflecting boundary !
             kcm = lone
             sgnm = sg
         endif

         call qintri( qtycor, dcor, ifd )

         if (ifd .eq. lone) then
            do 60 if = 1, nx
               qtyfin(if,jfp,kf) =
     &            max( qtymin(if), min( qtymax(if),
     &                 qty(if,1) / dfin(if,jfp,kf)    ))
 60         continue
         else
            do 61 if = 1, nx
 61            qtyfin(if,jfp,kf) = qty(if,1)
         endif

      enddo
C$OMP END PARALLEL DO
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine triz(qtyfin,qtycor,dfin,dcor,igrd,kfpl,kfpr,ifd,nbz,sg)

c     triquadratic interpolation of qtycor onto qtyfin
c     x-y-slices from ifpl to ifpr on fine grid
c     this version works only for nfine = 2 !

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8   qtycor( qx, qy, qz ),  dcor( qx, qy, qz ),
     &         qtyfin( qx, qy,nbz ),  dfin( qx, qy,nbz )
      integer*4 ifd

      !!print *, 'grids, triz'

      nox = norgin(igrd,lone)
      noy = norgin(igrd,ltwo)
      noz = norgin(igrd,lthree)
      noq = nox
      nq = nx

      kone = lone
      if (kfpr .eq. lzero) kone = lzero

      do kf = kfpl, kfpr
         kc = noz + (kf+kone)/nfine
         kcp = kc - ltwo*iand(kf,lone) + lone
         kcm = kc - (kcp-kc)
         sgnp = +one
         sgnm = +one
         if (kcp .lt. lone) then      ! reflecting boundary !
             kcp = lone
             sgnp = sg
         endif
         if (kcm .lt. lone) then      ! reflecting boundary !
             kcm = lone
             sgnm = sg
         endif
         kfp = kf
         if (kfpr.eq.lzero) kfp = kfp + qb/2
         if (kfpl.eq.nz+1)  kfp = kfp - nz + qb/2

C$OMP PARALLEL DO DEFAULT(NONE), ! COPYIN(kc,kcp,kcm,sgnm,sgnp),
C$OMP+            PRIVATE(jf,if),
c#ifdef SGI
C$OMP+            PRIVATE(jc,jcm,jcp,qty,qtymin,qtymax),
c#endif
C$OMP+            SHARED(nx,ny,nz,nox,noy,noz,noq,nq, kfpl,kfpr,
C$OMP+                          qtycor,qtyfin,dcor,dfin,ifd,kfp)
      do jf = 1, ny
         jc = noy + (jf+1)/nfine
         jcp = jc - ltwo*iand(jf,lone) + lone
         jcm = jc - (jcp-jc)

         call qintri( qtycor, dcor, ifd )

         if (ifd .eq. lone) then
            do 60 if = 1, nx
               qtyfin(if,jf,kfp) =
     &            max( qtymin(if), min( qtymax(if),
     &                 qty(if,1) / dfin(if,jf,kfp)    ))
 60         continue
         else
            do 61 if = 1, nx
 61            qtyfin(if,jf,kfp) = qty(if,1)
         endif

      enddo
C$OMP END PARALLEL DO
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine qintri (qtycor,dcor,ifd)

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8     qtycor(qx,qy,qz), dcor(qx,qy,qz), ftemp
      integer*4  ifd, if, nf, ls, m3

      !print *, 'grids, qintri'

      do if = 1, nx
         ic = noq + (if+1)/nfine
         icp = ic - ltwo*iand(if,lone) + lone
         icm = ic - (icp-ic)

         qty(if, 1) = qtycor(icm,jcm,kcm) * sgnm
         qty(if, 2) = qtycor(ic ,jcm,kcm)
         qty(if, 3) = qtycor(icp,jcm,kcm) * sgnp
         qty(if, 4) = qtycor(icm,jc ,kcm) * sgnm
         qty(if, 5) = qtycor(ic ,jc ,kcm)
         qty(if, 6) = qtycor(icp,jc ,kcm) * sgnp
         qty(if, 7) = qtycor(icm,jcp,kcm) * sgnm
         qty(if, 8) = qtycor(ic ,jcp,kcm)
         qty(if, 9) = qtycor(icp,jcp,kcm) * sgnp
         qty(if,10) = qtycor(icm,jcm,kc ) * sgnm
         qty(if,11) = qtycor(ic ,jcm,kc )
         qty(if,12) = qtycor(icp,jcm,kc ) * sgnp
         qty(if,13) = qtycor(icm,jc ,kc ) * sgnm
         qty(if,14) = qtycor(ic ,jc ,kc )
         qty(if,15) = qtycor(icp,jc ,kc ) * sgnp
         qty(if,16) = qtycor(icm,jcp,kc ) * sgnm
         qty(if,17) = qtycor(ic ,jcp,kc )
         qty(if,18) = qtycor(icp,jcp,kc ) * sgnp
         qty(if,19) = qtycor(icm,jcm,kcp) * sgnm
         qty(if,20) = qtycor(ic ,jcm,kcp)
         qty(if,21) = qtycor(icp,jcm,kcp) * sgnp
         qty(if,22) = qtycor(icm,jc ,kcp) * sgnm
         qty(if,23) = qtycor(ic ,jc ,kcp)
         qty(if,24) = qtycor(icp,jc ,kcp) * sgnp
         qty(if,25) = qtycor(icm,jcp,kcp) * sgnm
         qty(if,26) = qtycor(ic ,jcp,kcp)
         qty(if,27) = qtycor(icp,jcp,kcp) * sgnp
      enddo

c -----
      if (ifd .eq. lone) then
c -----

c         call mamiq
         do nf = 1, nq
            qtymin(nf) = min (
     &              qty(nf, 1), qty(nf, 2), qty(nf, 3),
     &              qty(nf, 4), qty(nf, 5), qty(nf, 6),
     &              qty(nf, 7), qty(nf, 8), qty(nf, 9),
     &              qty(nf,10), qty(nf,11), qty(nf,12),
     &              qty(nf,13), qty(nf,14), qty(nf,15),
     &              qty(nf,16), qty(nf,17), qty(nf,18),
     &              qty(nf,19), qty(nf,20), qty(nf,21),
     &              qty(nf,22), qty(nf,23), qty(nf,24),
     &              qty(nf,25), qty(nf,26), qty(nf,27)   )

            qtymax(nf) = max (
     &              qty(nf, 1), qty(nf, 2), qty(nf, 3),
     &              qty(nf, 4), qty(nf, 5), qty(nf, 6),
     &              qty(nf, 7), qty(nf, 8), qty(nf, 9),
     &              qty(nf,10), qty(nf,11), qty(nf,12),
     &              qty(nf,13), qty(nf,14), qty(nf,15),
     &              qty(nf,16), qty(nf,17), qty(nf,18),
     &              qty(nf,19), qty(nf,20), qty(nf,21),
     &              qty(nf,22), qty(nf,23), qty(nf,24),
     &              qty(nf,25), qty(nf,26), qty(nf,27)   )
         enddo

         do if = 1, nx
            ic = noq + (if+1)/nfine
            icp = ic - ltwo*iand(if,lone) + lone
            icm = ic - (icp-ic)

            qty(if, 1) = qty(if, 1) * dcor(icm,jcm,kcm)
            qty(if, 2) = qty(if, 2) * dcor(ic ,jcm,kcm)
            qty(if, 3) = qty(if, 3) * dcor(icp,jcm,kcm)
            qty(if, 4) = qty(if, 4) * dcor(icm,jc ,kcm)
            qty(if, 5) = qty(if, 5) * dcor(ic ,jc ,kcm)
            qty(if, 6) = qty(if, 6) * dcor(icp,jc ,kcm)
            qty(if, 7) = qty(if, 7) * dcor(icm,jcp,kcm)
            qty(if, 8) = qty(if, 8) * dcor(ic ,jcp,kcm)
            qty(if, 9) = qty(if, 9) * dcor(icp,jcp,kcm)
            qty(if,10) = qty(if,10) * dcor(icm,jcm,kc )
            qty(if,11) = qty(if,11) * dcor(ic ,jcm,kc )
            qty(if,12) = qty(if,12) * dcor(icp,jcm,kc )
            qty(if,13) = qty(if,13) * dcor(icm,jc ,kc )
            qty(if,14) = qty(if,14) * dcor(ic ,jc ,kc )
            qty(if,15) = qty(if,15) * dcor(icp,jc ,kc )
            qty(if,16) = qty(if,16) * dcor(icm,jcp,kc )
            qty(if,17) = qty(if,17) * dcor(ic ,jcp,kc )
            qty(if,18) = qty(if,18) * dcor(icp,jcp,kc )
            qty(if,19) = qty(if,19) * dcor(icm,jcm,kcp)
            qty(if,20) = qty(if,20) * dcor(ic ,jcm,kcp)
            qty(if,21) = qty(if,21) * dcor(icp,jcm,kcp)
            qty(if,22) = qty(if,22) * dcor(icm,jc ,kcp)
            qty(if,23) = qty(if,23) * dcor(ic ,jc ,kcp)
            qty(if,24) = qty(if,24) * dcor(icp,jc ,kcp)
            qty(if,25) = qty(if,25) * dcor(icm,jcp,kcp)
            qty(if,26) = qty(if,26) * dcor(ic ,jcp,kcp)
            qty(if,27) = qty(if,27) * dcor(icp,jcp,kcp)
         enddo

c -----
      endif
c -----

c      call interq
      do 30 ls = 3, 1, -1
      do 30 m3 = 2, 3**ls, 3
      do 30 nf = 1, nq
        ftemp = qty(nf,m3+1) - qty(nf,m3-1)
        qty(nf,(m3-2)/3+1) = qty(nf,m3)  +  half * max(
     &     (qty(nf,m3)-qty(nf,m3-1)) * (qty(nf,m3+1)-qty(nf,m3)), zero )
     &           / (ftemp + 1.D-50)
 30   continue

      return
      end
c-----------------------------------------------------------------------
      subroutine qintrj (qtycor,dcor,ifd)

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8     qtycor(qx,qy,qz), dcor(qx,qy,qz), ftemp
      integer*4  ifd, if, nf, ls, m3

      !print *, 'grids, qintrj'

      do jf = 1, ny
         jc = noq + (jf+1)/nfine
         jcp = jc - ltwo*iand(jf,lone) + lone
         jcm = jc - (jcp-jc)

         qty(jf, 1) = qtycor(icm,jcm,kcm) * sgnm
         qty(jf, 2) = qtycor(ic ,jcm,kcm)
         qty(jf, 3) = qtycor(icp,jcm,kcm) * sgnp
         qty(jf, 4) = qtycor(icm,jc ,kcm) * sgnm
         qty(jf, 5) = qtycor(ic ,jc ,kcm)
         qty(jf, 6) = qtycor(icp,jc ,kcm) * sgnp
         qty(jf, 7) = qtycor(icm,jcp,kcm) * sgnm
         qty(jf, 8) = qtycor(ic ,jcp,kcm)
         qty(jf, 9) = qtycor(icp,jcp,kcm) * sgnp
         qty(jf,10) = qtycor(icm,jcm,kc ) * sgnm
         qty(jf,11) = qtycor(ic ,jcm,kc )
         qty(jf,12) = qtycor(icp,jcm,kc ) * sgnp
         qty(jf,13) = qtycor(icm,jc ,kc ) * sgnm
         qty(jf,14) = qtycor(ic ,jc ,kc )
         qty(jf,15) = qtycor(icp,jc ,kc ) * sgnp
         qty(jf,16) = qtycor(icm,jcp,kc ) * sgnm
         qty(jf,17) = qtycor(ic ,jcp,kc )
         qty(jf,18) = qtycor(icp,jcp,kc ) * sgnp
         qty(jf,19) = qtycor(icm,jcm,kcp) * sgnm
         qty(jf,20) = qtycor(ic ,jcm,kcp)
         qty(jf,21) = qtycor(icp,jcm,kcp) * sgnp
         qty(jf,22) = qtycor(icm,jc ,kcp) * sgnm
         qty(jf,23) = qtycor(ic ,jc ,kcp)
         qty(jf,24) = qtycor(icp,jc ,kcp) * sgnp
         qty(jf,25) = qtycor(icm,jcp,kcp) * sgnm
         qty(jf,26) = qtycor(ic ,jcp,kcp)
         qty(jf,27) = qtycor(icp,jcp,kcp) * sgnp
      enddo

c -----
      if (ifd .eq. lone) then
c -----

c         call mamiq
         do nf = 1, nq
            qtymin(nf) = min (
     &              qty(nf, 1), qty(nf, 2), qty(nf, 3),
     &              qty(nf, 4), qty(nf, 5), qty(nf, 6),
     &              qty(nf, 7), qty(nf, 8), qty(nf, 9),
     &              qty(nf,10), qty(nf,11), qty(nf,12),
     &              qty(nf,13), qty(nf,14), qty(nf,15),
     &              qty(nf,16), qty(nf,17), qty(nf,18),
     &              qty(nf,19), qty(nf,20), qty(nf,21),
     &              qty(nf,22), qty(nf,23), qty(nf,24),
     &              qty(nf,25), qty(nf,26), qty(nf,27)   )

            qtymax(nf) = max (
     &              qty(nf, 1), qty(nf, 2), qty(nf, 3),
     &              qty(nf, 4), qty(nf, 5), qty(nf, 6),
     &              qty(nf, 7), qty(nf, 8), qty(nf, 9),
     &              qty(nf,10), qty(nf,11), qty(nf,12),
     &              qty(nf,13), qty(nf,14), qty(nf,15),
     &              qty(nf,16), qty(nf,17), qty(nf,18),
     &              qty(nf,19), qty(nf,20), qty(nf,21),
     &              qty(nf,22), qty(nf,23), qty(nf,24),
     &              qty(nf,25), qty(nf,26), qty(nf,27)   )
         enddo

         do jf = 1, ny
            jc = noq + (jf+1)/nfine
c           jcp = jc - ltwo*abs(mod(jf,ltwo)) + lone
            jcp = jc - ltwo*iand(jf,lone) + lone
            jcm = jc - (jcp-jc)

            qty(jf, 1) = qty(jf, 1) * dcor(icm,jcm,kcm)
            qty(jf, 2) = qty(jf, 2) * dcor(ic ,jcm,kcm)
            qty(jf, 3) = qty(jf, 3) * dcor(icp,jcm,kcm)
            qty(jf, 4) = qty(jf, 4) * dcor(icm,jc ,kcm)
            qty(jf, 5) = qty(jf, 5) * dcor(ic ,jc ,kcm)
            qty(jf, 6) = qty(jf, 6) * dcor(icp,jc ,kcm)
            qty(jf, 7) = qty(jf, 7) * dcor(icm,jcp,kcm)
            qty(jf, 8) = qty(jf, 8) * dcor(ic ,jcp,kcm)
            qty(jf, 9) = qty(jf, 9) * dcor(icp,jcp,kcm)
            qty(jf,10) = qty(jf,10) * dcor(icm,jcm,kc )
            qty(jf,11) = qty(jf,11) * dcor(ic ,jcm,kc )
            qty(jf,12) = qty(jf,12) * dcor(icp,jcm,kc )
            qty(jf,13) = qty(jf,13) * dcor(icm,jc ,kc )
            qty(jf,14) = qty(jf,14) * dcor(ic ,jc ,kc )
            qty(jf,15) = qty(jf,15) * dcor(icp,jc ,kc )
            qty(jf,16) = qty(jf,16) * dcor(icm,jcp,kc )
            qty(jf,17) = qty(jf,17) * dcor(ic ,jcp,kc )
            qty(jf,18) = qty(jf,18) * dcor(icp,jcp,kc )
            qty(jf,19) = qty(jf,19) * dcor(icm,jcm,kcp)
            qty(jf,20) = qty(jf,20) * dcor(ic ,jcm,kcp)
            qty(jf,21) = qty(jf,21) * dcor(icp,jcm,kcp)
            qty(jf,22) = qty(jf,22) * dcor(icm,jc ,kcp)
            qty(jf,23) = qty(jf,23) * dcor(ic ,jc ,kcp)
            qty(jf,24) = qty(jf,24) * dcor(icp,jc ,kcp)
            qty(jf,25) = qty(jf,25) * dcor(icm,jcp,kcp)
            qty(jf,26) = qty(jf,26) * dcor(ic ,jcp,kcp)
            qty(jf,27) = qty(jf,27) * dcor(icp,jcp,kcp)
         enddo

c -----
      endif
c -----

c      call interq
      do 30 ls = 3, 1, -1
      do 30 m3 = 2, 3**ls, 3
      do 30 nf = 1, nq
         ftemp = qty(nf,m3+1) - qty(nf,m3-1)
         qty(nf,(m3-2)/3+1) = qty(nf,m3)  +  half * max(
     &     (qty(nf,m3)-qty(nf,m3-1)) * (qty(nf,m3+1)-qty(nf,m3)), zero )
     &           / (ftemp + 1.D-50)
 30   continue

      return
      end
c-----------------------------------------------------------------------
      subroutine interq

      include 'qparam.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

cq    facq = one / 32768.D0
cq
cq    do 20 nf = 1, nq
cq
c --  triquadratic interpolation
cq       qtyout(nf) =  facq * (
cq   &   -   27.* qty(nf, 1)  +  270.* qty(nf, 2)  +   45.* qty(nf, 3)
cq   &   +  270.* qty(nf, 4)  - 2700.* qty(nf, 5)  -  450.* qty(nf, 6)
cq   &   +   45.* qty(nf, 7)  -  450.* qty(nf, 8)  -   75.* qty(nf, 9)
cq   &   +  270.* qty(nf,10)  - 2700.* qty(nf,11)  -  450.* qty(nf,12)
cq   &   - 2700.* qty(nf,13)  +27000.* qty(nf,14)  + 4500.* qty(nf,15)
cq   &   -  450.* qty(nf,16)  + 4500.* qty(nf,17)  +  750.* qty(nf,18)
cq   &   +   45.* qty(nf,19)  -  450.* qty(nf,20)  -   75.* qty(nf,21)
cq   &   -  450.* qty(nf,22)  + 4500.* qty(nf,23)  +  750.* qty(nf,24)
cq   &   -   75.* qty(nf,25)  +  750.* qty(nf,26)  +  125.* qty(nf,27))
cq
c --  trilinear interpolation
cq       qtylin(nf) =  0.015625 * (
cq   &         27. * qty(nf,14) +  9. * qty(nf,15)
cq   &       +  9. * qty(nf,23) +  3. * qty(nf,18)
cq   &       +  9. * qty(nf,17) +  3. * qty(nf,26)
cq   &       +  3. * qty(nf,24) +       qty(nf,27) )
cq
c20   continue

c --  Bram van Leer -interpolation  (JCP 23, 267)

      !print *, 'grids, interq'

      do 30 ls = 3, 1, -1
      do 30 m3 = 2, 3**ls, 3
      do 30 nf = 1, nq
        qty(nf,(m3-2)/3+1) = qty(nf,m3)  +  half * max(
     &     (qty(nf,m3)-qty(nf,m3-1)) * (qty(nf,m3+1)-qty(nf,m3)), zero )
     &           / (qty(nf,m3+1) - qty(nf,m3-1) + 1.D-90)
 30   continue
cq
cq    do 40 nf = 1, nq
cq        if ( (qtyout(nf).lt.qtymin(nf) ) .or.
cq   &         (qtyout(nf).gt.qtymax(nf) ) .or.
cq   &         (abs(qtymax(nf)).gt.ten*abs(qtymin(nf))) )
cq   &          qtyout(nf) = qtylin(nf)
c40   continue

      return
      end
c-----------------------------------------------------------------------
      subroutine mamiq

      include 'qparam.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      !print *, 'grids, mamiq'


      do nf = 1, nq

            qtymin(nf) = min (
     &              qty(nf, 1), qty(nf, 2), qty(nf, 3),
     &              qty(nf, 4), qty(nf, 5), qty(nf, 6),
     &              qty(nf, 7), qty(nf, 8), qty(nf, 9),
     &              qty(nf,10), qty(nf,11), qty(nf,12),
     &              qty(nf,13), qty(nf,14), qty(nf,15),
     &              qty(nf,16), qty(nf,17), qty(nf,18),
     &              qty(nf,19), qty(nf,20), qty(nf,21),
     &              qty(nf,22), qty(nf,23), qty(nf,24),
     &              qty(nf,25), qty(nf,26), qty(nf,27)   )

            qtymax(nf) = max (
     &              qty(nf, 1), qty(nf, 2), qty(nf, 3),
     &              qty(nf, 4), qty(nf, 5), qty(nf, 6),
     &              qty(nf, 7), qty(nf, 8), qty(nf, 9),
     &              qty(nf,10), qty(nf,11), qty(nf,12),
     &              qty(nf,13), qty(nf,14), qty(nf,15),
     &              qty(nf,16), qty(nf,17), qty(nf,18),
     &              qty(nf,19), qty(nf,20), qty(nf,21),
     &              qty(nf,22), qty(nf,23), qty(nf,24),
     &              qty(nf,25), qty(nf,26), qty(nf,27)   )

         enddo

      return
      end
c-----------------------------------------------------------------------
      logical function overlap (igr1,igr2)

c   if grid1 and grid2 overlap, then .true. .
c   comparison is done on the topmost grid.  here: checks only x and y !
c   nox1,noy1,nox2,noy2,overx1, overx2, overy1, overy2 are also output !

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'


      write(*,*) 'overlap   igr1:',igr1,'   igr2:',igr2
      nox1 = norgin(igr1,lone)
      noy1 = norgin(igr1,ltwo)
      noz1 = norgin(igr1,lthree)

      nox2 = norgin(igr2,lone)
      noy2 = norgin(igr2,ltwo)
      noz2 = norgin(igr2,lthree)

      toxl1 = topos (  one   , lone  , igr1)
      toyl1 = topos (  one   , ltwo  , igr1)
c     tozl1 = topos (  one   , lthree, igr1)
      toxr1 = topos (dble(nx), lone  , igr1)
      toyr1 = topos (dble(ny), ltwo  , igr1)
c     tozr1 = topos (dble(nz), lthree, igr1)

      toxl2 = topos (  one   , lone  , igr2)
      toyl2 = topos (  one   , ltwo  , igr2)
c     tozl2 = topos (  one   , lthree, igr2)
      toxr2 = topos (dble(nx), lone  , igr2)
      toyr2 = topos (dble(ny), ltwo  , igr2)
c     tozr2 = topos (dble(nz), lthree, igr2)
      write(*,*) 'toxl1:',toxl1,'  toxr1:',toxr1,'  toyl1:',toyl1,
     &           '  toyr1:',toyr1
      write(*,*) 'toxl2:',toxl2,'  toxr2:',toxr2,'  toyl2:',toyl2,
     &           '  toyr2:',toyr2

      overx1 =       ( (toxl1-toxl2) .gt. -1.D-9*abs(toxr1-toxl1) )
     &         .and. ( (toxr2-toxl1) .gt. -1.D-9*abs(toxr1-toxl1) )
      overx2 =       ( (toxl2-toxl1) .gt. -1.D-9*abs(toxr2-toxl2) )
     &         .and. ( (toxr1-toxl2) .gt. -1.D-9*abs(toxr2-toxl2) )

      overy1 =       ( (toyl1-toyl2) .gt. -1.D-9*abs(toyr1-toyl1) )
     &         .and. ( (toyr2-toyl1) .gt. -1.D-9*abs(toyr1-toyl1) )
      overy2 =       ( (toyl2-toyl1) .gt. -1.D-9*abs(toyr2-toyl2) )
     &         .and. ( (toyr1-toyl2) .gt. -1.D-9*abs(toyr2-toyl2) )

      write(*,*) 'overx1:',overx1,'  overy1:',overy1,
     &           '  overx2:',overx2,'  overy2:',overy2
      overlap = (overx1 .or. overx2) .and. (overy1 .or. overy2)
c    &          .and. (idxcgr(igr1) .eq. idxcgr(igr2))

      lenx1 = min(nint(dble(2**(levlgr(igr1,1)-1))*(toxr1-toxl1)) ,
     &            nint(dble(2**(levlgr(igr1,1)-1))*(toxr2-toxl1)-.25D0))

      lenx2 = min(nint(dble(2**(levlgr(igr2,1)-1))*(toxr2-toxl2)) ,
     &            nint(dble(2**(levlgr(igr2,1)-1))*(toxr1-toxl2))  )
      write(*,*) 'lenx1:',lenx1,'  lenx2:',lenx2

      if ( overx1 ) then
        il1 = lone
        ir1 = il1 + lenx1
        il2 = 1 + nint( dble(2**(levlgr(igr2,1)-1))*(toxl1-toxl2) )
        ir2 = il2 + lenx1/(2**(levlgr(igr1,1)-levlgr(igr2,1)))
      else
        il2 = lone
        ir2 = il2 + lenx2
        il1 =lone+nint(dble(2**(levlgr(igr1,1)-1))*(toxl2-toxl1)-.25D0)
        ir1 = il1 + lenx2*(2**(levlgr(igr1,1)-levlgr(igr2,1)))
      endif

      leny1 = min(nint(dble(2**(levlgr(igr1,1)-1))*(toyr1-toyl1)) ,
     &            nint(dble(2**(levlgr(igr1,1)-1))*(toyr2-toyl1)-.25D0))

      leny2 = min(nint(dble(2**(levlgr(igr2,1)-1))*(toyr2-toyl2)) ,
     &            nint(dble(2**(levlgr(igr2,1)-1))*(toyr1-toyl2))  )

      write(*,*) 'leny1:',leny1,'  leny2:',leny2
      if ( overy1 ) then
        jl1 = lone
        jr1 = jl1 + leny1
        jl2 = 1 + nint( dble(2**(levlgr(igr2,1)-1))*(toyl1-toyl2) )
        jr2 = jl2 + leny1/(2**(levlgr(igr1,1)-levlgr(igr2,1)))
      else
        jl2 = lone
        jr2 = jl2 + leny2
        jl1 =lone+nint(dble(2**(levlgr(igr1,1)-1))*(toyl2-toyl1)-.25D0)
        jr1 = jl1 + leny2*(2**(levlgr(igr1,1)-levlgr(igr2,1)))
      endif

      write(*,*)  'overlap   ', il1, ir1, jl1, jr1, il2, ir2, jl2, jr2
      return
      end
c-----------------------------------------------------------------------
      subroutine takval (qtakx,qtaky,qtgiv)

c     take boundary values for qtytak from qtgiv
c     this version works only for nfine = 2 !

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8 qtakx( qb, qy, qz ), qtaky( qx, qb, qz )
      real*8 qtgiv( qx, qy, qz )

      !print *, 'grids, takval'


      if ( overx1 ) then
        ixmi = nox1
        ixma = nox2 + nx/nfine
        nb1x = lone
        ng1x = nfine * (nox1 + (  -1+nfine-1)/nfine -1 - nox2) + 1
      else
        ixmi = nox2
        ixma = nox1 + nx/nfine
        nb1x = 5
        ng1x = nfine * (nox1 + (nx+1+nfine-1)/nfine -1 - nox2) + 1
      endif

      nb2x = nb1x + 1
      nb3x = nb1x + 2
      nb4x = nb1x + 3
      ng2x = ng1x + 1
      ng3x = ng1x + 2
      ng4x = ng1x + 3

      if ( overy1 ) then
        iymi = noy1
        iyma = noy2 + ny/nfine
        nb1y = 1
        ng1y = nfine * (noy1 + ( -1+nfine-1)/nfine -1 - noy2) + 1
      else
        iymi = noy2
        iyma = noy1 + ny/nfine
        nb1y = 5
        ng1y = nfine * (noy1 + (nx+1+nfine-1)/nfine -1 - noy2) + 1
      endif

      nb2y = nb1y + 1
      nb3y = nb1y + 2
      nb4y = nb1y + 3
      ng2y = ng1y + 1
      ng3y = ng1y + 2
      ng4y = ng1y + 3

      if ((ng4y.gt.qy).or.(ng4x.gt.qx))
     &      write(*,*) 'takval   ng4y:',ng4y,'  ng4x:',ng4x

      is1 = (ixmi-nox1)*nfine
      js1 = (iymi-noy1)*nfine
      is2 = (ixmi-nox2)*nfine
      js2 = (iymi-noy2)*nfine
      write(*,*) 'takval     1, il1:',il1,ir1,jl1,jr1,il2,ir2,jl2,jr2

C$OMP PARALLEL DO PRIVATE(i,k)
      do i = 1, (ixma-ixmi)*nfine
      do k = 1, nz
         qtaky( i+is1, nb1y, k) = qtgiv( i+is2, ng1y, k)
         qtaky( i+is1, nb2y, k) = qtgiv( i+is2, ng2y, k)
         qtaky( i+is1, nb3y, k) = qtgiv( i+is2, ng3y, k)
         qtaky( i+is1, nb4y, k) = qtgiv( i+is2, ng4y, k)
      enddo
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO PRIVATE(i,k)
      do j = 1, (iyma-iymi)*nfine
      do k = 1, nz
         qtakx( nb1x, j+js1, k) = qtgiv( ng1x, j+js2, k)
         qtakx( nb2x, j+js1, k) = qtgiv( ng2x, j+js2, k)
         qtakx( nb3x, j+js1, k) = qtgiv( ng3x, j+js2, k)
         qtakx( nb4x, j+js1, k) = qtgiv( ng4x, j+js2, k)
      enddo
      enddo
C$OMP END PARALLEL DO


      return
      end
c-----------------------------------------------------------------------
      integer*4 function ifico (ipos,ixyz,ifgr)

c     given position ipos on coarse grid ,
c     return position on fine grid ifgr (smallest value)

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'

      ifico = nfine * (ipos - norgin(ifgr,ixyz) - lone ) + lone

      return
      end
c-----------------------------------------------------------------------
      integer*4 function icofi (ipos,ixyz,ifgr)

c     given position ipos on fine grid ifgr,
c     return position on grid one step coarser

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'


      icofi = norgin(ifgr,ixyz) + (ipos+nfine-lone)/nfine


      return
      end
c-----------------------------------------------------------------------
      real*8 function topos (pos,ixyz,igrd)

c     given position pos on grid igrd,
c     return position on topmost grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'

      apos = pos
      icgr = igrd
      rkfin = one / dble(nfine)

 10   if (levlgr(icgr,lone) .eq. lone ) goto 20
         apos = dble(norgin(icgr,ixyz)) + 0.25D0   +   apos * rkfin
         icgr = idxcgr(icgr)
      goto 10

 20   topos = apos

      return
      end
c-----------------------------------------------------------------------
      real*8 function copos (pos,ixyz,igrd,itgr)

c     given position pos on grid igrd,
c     return position on coarser grid itgr

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'

      apos = pos
      icgr = igrd
      rkfin = one / dble(nfine)

 10   if (levlgr(icgr,lone) .eq. levlgr(itgr,lone) ) goto 20
         apos = dble(norgin(icgr,ixyz)) + 0.25D0   +   apos * rkfin
         icgr = idxcgr(icgr)
      goto 10

 20   copos = apos

      return
      end
c-----------------------------------------------------------------------
      real*8 function fipos (pos,ixyz,igrd)

c     given position pos on topmost grid,
c     return position on finer grid igrd.

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'

      apos = pos
      ifgr = indxgr(lone,lone)
      rnfin = dble(nfine)

 10   if ( ifgr .eq. igrd ) goto 20
         ifgr = idxfgr(ifgr)
         apos = ( apos - dble(norgin(ifgr,ixyz)) - 0.25D0 ) * rnfin
      goto 10

 20   fipos = apos

      return
      end
c-----------------------------------------------------------------------
      subroutine regrid

c     compute grid origin coordinates of moved grid,
c     only equidistant cartesian grids are used.

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'
      include 'bounds.cmn'

      !print *, 'grids, regrid'

cf    logical flag4

cw    write(*,*) 'begin of regrid'
c --- z coordinates of origins are constant !

c --- did the WD move too much from the centre
c     of the finest grid with number nwdgrd?

      nwdgrd = indxgr(lone,maxlev(lone))
c     nwdgrd = indxgr(ltwo,maxlev(ltwo))

      id2 = nint( posx2 - (dble(nx/2)+half) )
      jd2 = nint( posy2 - (dble(ny/2)+half) )
      mwdsx= nx/2-nx/3
      mwdsy= ny/2-ny/3

      if ( ( abs(id2) .gt. mwdsx ) .or.
     &     ( abs(jd2) .gt. mwdsy )   )  then
         call movegrid(nwdgrd,id2/nfine,jd2/nfine)
         posx2 = posx2 - dble(nfine*(id2/nfine))
         posy2 = posy2 - dble(nfine*(jd2/nfine))
      endif


c --- topmost grid (index 1) is absolute and needs no origin.
c --- all odd level grids are fixed.

cf--- flag4 = .true. if grid 4 moves
cf    flag4 = .false.

c --- did some intermediate grid igr approach the border
c     of the next coarser grid icog ?

      do 10 igr = ngrd, 2, -1

         icog = idxcgr(igr)

         id2 = norgin(igr,1) - nx/4
         jd2 = norgin(igr,2) - ny/4

         if ( ( abs(id2) .gt. nx/4-4-mwdsx/2 ) .or.
     &        ( abs(jd2) .gt. ny/4-4-mwdsy/2 ) )   then
            if (icog .eq. lone) then
                write(*,*) 'grid 2 moved too far: ',id2,jd2
                stop
            endif
            call movegrid(icog,id2/nfine,jd2/nfine)
         endif

  10  continue


cw    write(*,*) 'end of regrid'

      return
      end
c-----------------------------------------------------------------------
      subroutine movegrid(igrd,id,jd)

c     move grid igrd by id in x and ij in y direction
c     id, jd are the number of cells on the coarse grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      write(*,*) '** movegrid **    igrd:',igrd,'  id:',id,'  jd:',jd

      !print *, 'grids, movegrid'

      do 160 lev = mlev, 2, -1
 160      call fineup(lev,lzero)

      idf = id * nfine
      jdf = jd * nfine

c --- change origin of igrd relative to coarser grid

      norgin(igrd,lone) = norgin(igrd,lone) + id
      norgin(igrd,ltwo) = norgin(igrd,ltwo) + jd

      norx = norgin(igrd,lone)
      nory = norgin(igrd,ltwo)
      if((norx.lt.1).or.(norx.gt.nx).or.(nory.lt.1).or.(nory.gt.ny))then
         write(*,*) 'regrid  norgin x,y : ', norx, nory
         stop
      endif

c --- change origin of all grids one level finer
c     and contained in igrd

      do 10 ifgr = 1, ngrd
      if ( idxcgr(ifgr) .eq. igrd ) then
         norgin(ifgr,lone) = norgin(ifgr,lone) - idf
         norgin(ifgr,ltwo) = norgin(ifgr,ltwo) - jdf
      endif
 10   continue


c --- save as much data as possible by copying it

      if ( id .lt. lzero ) then
         ifl =     lone
         ifr =  -idf
         ibeg=  nx - ifr
         iend=  lone
      else
         ifl =  nx -idf +1
         ifr =  nx
         ibeg=  idf + 1
         iend=  nx
      endif

      if ( jd .lt. lzero ) then
         jfl =     lone
         jfr =  -jdf
         jbeg=  ny - jfr
         jend=  lone
      else
         jfl =  ny - jdf +1
         jfr =  ny
         jbeg=  jdf + 1
         jend=  ny
      endif

      call ficopy(igrd,ibeg,iend,idf,jbeg,jend,jdf)


c --- igrd lies within which coarser grid icog?

      icog = idxcgr(igrd)

c --- interpolate missing data on grid from coarser grid

c --- y-z-slice

c     icpl = icofi(ifl,igrd,lone)
c     icpr = icofi(ifr,igrd,lone)

      write(*,*) 'movegrid    y-z-slice  ifl:',ifl,'  ifr:',ifr
      call triall ( lone, igrd, icog, ifl, ifr )


c --- x-z-slice

c     icpl = icofi(jfl,igrd,lone)
c     icpr = icofi(jfr,igrd,lone)

      write(*,*) 'movegrid    x-z-slice  jfl:',jfl,'  jfr:',jfr
      call triall ( ltwo, igrd, icog, jfl, jfr )

c --- x-y-slice
c     not needed, since the grids only move in x and y direction.

      return
      end
c-----------------------------------------------------------------------
      subroutine ficopy(igrd,ibeg,iend,idf,jbeg,jend,jdf)

c     shift contents of igrd by idf and jdf;
c     copy all quantities

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'

      !print *, 'grids, ficopy'


      jsd = isign(lone,jdf)
      isd = isign(lone,idf)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,j,js,i,is),
C$OMP+     SHARED(velx,vely,velz,densty,energy,gpot,gtpot,rgarr,
C$OMP+            jbeg,jend,jsd,jdf,ibeg,iend,isd,idf,igrd)
      do k = 1, qz
      do j = jbeg, jend, jsd
      js = j - jdf
      do i = ibeg, iend, isd
      is = i - idf
           velx(is,js,k,igrd) =   velx(i,j,k,igrd)
           vely(is,js,k,igrd) =   vely(i,j,k,igrd)
           velz(is,js,k,igrd) =   velz(i,j,k,igrd)
         densty(is,js,k,igrd) = densty(i,j,k,igrd)
         energy(is,js,k,igrd) = energy(i,j,k,igrd)
           gpot(is,js,k,igrd) =   gpot(i,j,k,igrd)
          gtpot(is,js,k,igrd) =  gtpot(i,j,k,igrd)
          rgarr(is,js,k,igrd) =  rgarr(i,j,k,igrd)
       enddo
       enddo
       enddo
C$OMP END PARALLEL DO


      return
      end
c-----------------------------------------------------------------------
      subroutine detoro(igrd,mm)

c     copy densty to ro

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8 di, dj ,dk
      integer*4 mm

      !print *, 'grids, detoro'

      d1x = delx(igrd)

      if ( mm .eq. lone ) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,j,i),
C$OMP+            SHARED(ro,densty,nx,ny,nz,igrd)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            ro(i,j,k) = densty(i,j,k,igrd)
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      elseif ( mm .eq. lthree ) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,i,k),
C$OMP+            SHARED(ro,densty,energy,nx,ny,nz,igrd,cc,corrmax)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            ro(i,j,k) = densty(i,j,k,igrd)
     &            * min( two * energy(i,j,k,igrd) / cc**2, corrmax )
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      elseif ( mm .eq. ltwo ) then

c     copy divergence ( densty * velocity ) to ro
c     this ro is also used in gwaves !! cf. dtgravty, inigrav !

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(ro,densty,velx,vely,velz),
C$OMP+            PRIVATE(j,jp,jm,dj,i,ip,im,di,k,kp,km,dk,sg),
C$OMP+            SHARED(nx,ny,nz,igrd,d1x)
         do j = 1, ny
            jp = min( j+1, ny )
            jm = max( j-1, lone  )
            dj = dble(jp - jm) * d1x
         do k = 1, nz
            sg = +one
            if ( k .eq. lone )  sg = -one ! reflection of velz at k=0.5
            kp = min( k+1, nz )
            km = max( k-1, lone  )
            dk = two  * d1x
         do i = 1, nx
            ip = min( i+1, nx )
            im = max( i-1, lone  )
            di = dble(ip - im) * d1x
        ro(i,j,k) =
     &    (     densty(ip,j, k, igrd) * velx(ip,j, k, igrd)
     &       -  densty(im,j, k, igrd) * velx(im,j, k, igrd) ) / di
     &  + (     densty(i, jp,k, igrd) * vely(i, jp,k, igrd)
     &       -  densty(i, jm,k, igrd) * vely(i, jm,k, igrd) ) / dj
     &  + (     densty(i, j, kp,igrd) * velz(i, j, kp,igrd)
     &      -sg*densty(i, j, km,igrd) * velz(i, j, km,igrd) ) / dk
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      else
         write(*,*) 'something wrong in detoro ',mm
         stop
      endif

      return

      end
c-----------------------------------------------------------------------
      subroutine rozero(ifgr)

c     zero the part of the coarse grid igrd in which ifgr lies

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      !print *, 'grids, rozero'


      noz = norgin(ifgr,lthree)
      noy = norgin(ifgr,ltwo)
      nox = norgin(ifgr,lone)

C$OMP PARALLEL DO PRIVATE(k,j,i)
      do j = noy+1, noy+ny/nfine
      do k = noz+1, noz+nz/nfine
      do i = nox+1, nox+nx/nfine
         ro(i,j,k) = zero
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine adbakf(igrd,ifgr,mm)

c     add background potential of coarse grid to fine grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      !print *, 'grids, adbakf'

      noz = norgin(ifgr,lthree)
      noy = norgin(ifgr,ltwo)
      nox = norgin(ifgr,lone)

c do one less and more for tripo
c     do 10 k = noz, noz+1+nz/nfine
c exceptional: tripo contains refecting z=0 boundary

      if ( mm .eq. lone ) then

C$OMP PARALLEL DO PRIVATE(k,j,i)
         do j = noy,   noy+1+ny/nfine
         do k = noz+1, noz+1+nz/nfine
         do i = nox,   nox+1+nx/nfine
            po(i,j,k) = po(i,j,k) * delx(igrd)**2   + bakpot(i,j,k,igrd)
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      elseif ( mm .eq. ltwo ) then

C$OMP PARALLEL DO PRIVATE(k,j,i)
         do j = noy,   noy+1+ny/nfine
         do k = noz+1, noz+1+nz/nfine
         do i = nox,   nox+1+nx/nfine
            po(i,j,k) = -po(i,j,k) * delx(igrd)**2  + baktpo(i,j,k,igrd)
                      ! minus from equation
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      elseif ( mm .eq. lthree ) then

C$OMP PARALLEL DO PRIVATE(k,j,i)
         do j = noy,   noy+1+ny/nfine
         do k = noz+1, noz+1+nz/nfine
         do i = nox,   nox+1+nx/nfine
            po(i,j,k) = po(i,j,k) * delx(igrd)**2   + bakchi(i,j,k,igrd)
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      else
         write (*,*) 'something wrong in adbakf ', mm
         stop
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine tripo (igrd,ifgr,mm)

c     triquadratic interpolation of po on igrd to bakpot on ifgr
c     in x-y-slices along the whole z direction
c     this version works only for nfine = 2 !
c !!! only for reflecting z=0 boundary !!!

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      !print *, 'grids, tripo'


      nox = norgin(ifgr,lone)
      noy = norgin(ifgr,ltwo)
      noz = norgin(ifgr,lthree)

      ifd = lzero
      noq = nox
      nq = nx

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+       PRIVATE(if,jf,kf),
c#ifdef SGI
C$OMP+       PRIVATE(jc,jcp,jcm,qty,kc,kcm,kcp,sgnp,sgnm),  ! threadprivate
c#endif
C$OMP+       SHARED(mm,ifgr,nx,ny,nz,nox,noy,noz,noq,nq,po,ifd,bakpot,
C$OMP+              baktpo,bakchi,bakrg,bakttp,bakvpx,bakvpy,bakvpz)

      do jf = 1, ny
         jc = noy + (jf+1)/nfine
         jcp = jc - ltwo*iand(jf,lone) + lone
         jcm = jc - (jcp-jc)
         sgnp = +one
         sgnm = +one

      do kf = 1, nz
         kc = noz + (kf+1)/nfine
         kcp = kc - ltwo*iand(kf,lone) + lone
         kcm = kc - (kcp-kc)
         kcp = max( kcp, lone )       ! reflecting boundary !
         kcm = max( kcm, lone )       ! reflecting boundary !

         call qintri( po, po, ifd )

         if ( mm .eq. lone ) then
            do 60 if = 1, nx
   60          bakpot(if,jf,kf,ifgr) = qty(if,1)
         elseif (mm .eq. ltwo) then
            do 61 if = 1, nx
   61          baktpo(if,jf,kf,ifgr) = qty(if,1)
         elseif (mm .eq. lthree) then
            do 63 if = 1, nx
   63          bakchi(if,jf,kf,ifgr) = qty(if,1)
         elseif (mm .eq. 4) then
            do 64 if = 1, nx
   64          bakrg (if,jf,kf,ifgr) = qty(if,1)
         else
            write(*,*) 'something wrong in tripo'
         endif

      enddo
      enddo
C$OMP END PARALLEL DO


      return
      end
c-----------------------------------------------------------------------
      subroutine adbakcf(igrd,mm)

c     add background potential bakpot to po and put into gpot

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'

      !print *, 'grids, adbakcf'

      if ( mm .eq. lone ) then

C$OMP PARALLEL DO PRIVATE(i,j,k)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            gpot(i,j,k,igrd) = po(i,j,k) * delx(igrd)**2
     &                                  + bakpot(i,j,k,igrd)
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      elseif ( mm .eq. ltwo ) then

C$OMP PARALLEL DO PRIVATE(i,j,k)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            gtpot(i,j,k,igrd) = -po(i,j,k) * delx(igrd)**2
     &                                  + baktpo(i,j,k,igrd)
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      elseif ( mm .eq. lthree ) then

C$OMP PARALLEL DO PRIVATE(i,j,k)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            gchi(i,j,k,igrd) = po(i,j,k) * delx(igrd)**2
     &                                  + bakchi(i,j,k,igrd)
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine avrqty(qtyfin,qtycor,ifgr,dfin,dcor,ll)

c     average values of fine grid onto coarse grid
c     this version works only for nfine = 2 !

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8 qtyfin( qx, qy, qz ), qtycor( qx, qy, qz ),
     &         dfin( qx, qy, qz ),   dcor( qx, qy, qz )

      noz = norgin(ifgr,lthree)
      noy = norgin(ifgr,ltwo)
      nox = norgin(ifgr,lone)

      refin = one / dble(nfine**3)

      if (ll .eq. lzero) then

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,j1,k,k1,i,i1),
C$OMP+            SHARED(qtycor,qtyfin,nx,ny,nz,nox,noy,noz,refin)
c#ifdef SGI
C$OMP+            ,PRIVATE(jc,kc,ic)     ! already in threadprivate
c#endif
         do k = 1, nz-1, nfine
            kc = noz + (k+nfine-1)/nfine
            k1 = k + 1
         do j = 1, ny-1, nfine
            jc = noy + (j+nfine-1)/nfine
            j1 = j + 1
         do i = 1, nx-1, nfine
            ic = nox + (i+nfine-1)/nfine
            i1 = i + 1
            qtycor(ic,jc,kc) = refin * (
     &                 ((qtyfin(i ,j ,k ) + qtyfin(i1,j ,k ))
     &               +  (qtyfin(i ,j1,k ) + qtyfin(i1,j1,k )))
     &               + ((qtyfin(i ,j ,k1) + qtyfin(i1,j ,k1))
     &               +  (qtyfin(i ,j1,k1) + qtyfin(i1,j1,k1))))
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      else

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,j1,k,k1,i,i1),
C$OMP+            SHARED(dcor,qtycor,qtyfin,dfin,
C$OMP+                   nx,ny,nz,nox,noy,noz,refin)
c#ifdef SGI
C$OMP+            ,PRIVATE(jc,kc,ic)   ! already in threadprivate
c#endif
         do k = 1, nz-1, nfine
            kc = noz + (k+nfine-1)/nfine
            k1 = k + 1
         do j = 1, ny-1, nfine
            jc = noy + (j+nfine-1)/nfine
            j1 = j + 1
         do i = 1, nx-1, nfine
            ic = nox + (i+nfine-1)/nfine
            i1 = i + 1
            qtycor(ic,jc,kc) = refin / dcor(ic,jc,kc) * (
     &             ( ( qtyfin(i ,j ,k ) * dfin(i ,j, k )
     &                +qtyfin(i1,j ,k ) * dfin(i1,j ,k ))
     &              +( qtyfin(i ,j1,k ) * dfin(i ,j1,k )
     &                +qtyfin(i1,j1,k ) * dfin(i1,j1,k )))
     &            +( ( qtyfin(i ,j ,k1) * dfin(i ,j ,k1)
     &                +qtyfin(i1,j ,k1) * dfin(i1,j ,k1))
     &              +( qtyfin(i ,j1,k1) * dfin(i ,j1,k1)
     &                +qtyfin(i1,j1,k1) * dfin(i1,j1,k1)))  )
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine avrage(ifgr,icgr)

c --- average values of fine grid onto coarse grid
c --- add point mass if ifgr is finest level

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include 'ppdisk.cmn'

      integer*4 ic
      real*8 ek(q)

      if ( icgr .eq. idxcgr(ifgr) ) then

      call avrqty ( densty(1,1,1,ifgr), densty(1,1,1,icgr), ifgr,
     &              densty(1,1,1,ifgr), densty(1,1,1,icgr), lzero)
      call avrqty (   gpot(1,1,1,ifgr),   gpot(1,1,1,icgr), ifgr,
     &              densty(1,1,1,ifgr), densty(1,1,1,icgr), lzero)

      call avrqty (   velx(1,1,1,ifgr),   velx(1,1,1,icgr), ifgr,
     &              densty(1,1,1,ifgr), densty(1,1,1,icgr), lone )
      call avrqty (   vely(1,1,1,ifgr),   vely(1,1,1,icgr), ifgr,
     &              densty(1,1,1,ifgr), densty(1,1,1,icgr), lone )
      call avrqty (   velz(1,1,1,ifgr),   velz(1,1,1,icgr), ifgr,
     &              densty(1,1,1,ifgr), densty(1,1,1,icgr), lone )
      call avrqty ( energy(1,1,1,ifgr), energy(1,1,1,icgr), ifgr,
     &              densty(1,1,1,ifgr), densty(1,1,1,icgr), lone )

      do m = 1, qc
      call avrqty (   chem(1,1,1,ifgr,m), chem(1,1,1,icgr,m),ifgr,
     &              densty(1,1,1,ifgr), densty(1,1,1,icgr), lzero)
      enddo

      else

      write(*,*)  'avrage  something wrong '

      endif

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,i,k,ek),
C$OMP+            SHARED(velx,vely,velz,energy,nx,ny,nz,small,icgr)

      do j = 1, ny
      do k = 1, nz
         do i = 1, nx
            ek(i) = (one+small) * half *  ( velx(i,j,k,icgr)**2
     &                + vely(i,j,k,icgr)**2 + velz(i,j,k,icgr)**2 )
            energy(i,j,k,icgr) = max( energy(i,j,k,icgr), ek(i) )
            if ( energy(i,j,k,icgr) .eq. ek(i) )
     &        write(*,*) 'avrage,     energy:',i,j,k,icgr
         enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c --- add point mass if ifgr is finest level

cp1   if     ( ifgr .eq. indxgr(lone,maxlev(lone)) ) then
cp1 for collisions no pmass1 is used

c     write(*,*) 'avrage    add point mass 1'
cp1      icox = icofi( int(posx1), lone, ifgr )
cp1      icoy = icofi( int(posy1), ltwo, ifgr )
cp1      icox1= icox + 1
cp1      icoy1= icoy + 1

cp1      densty( icox , icoy , 1, icgr) =
cp1  &   densty( icox , icoy , 1, icgr) + pmass1 / delx(icgr)**3 / 8.D0
cp1      densty( icox1, icoy , 1, icgr) =
cp1  &   densty( icox1, icoy , 1, icgr) + pmass1 / delx(icgr)**3 / 8.D0
cp1      densty( icox , icoy1, 1, icgr) =
cp1  &   densty( icox , icoy1, 1, icgr) + pmass1 / delx(icgr)**3 / 8.D0
cp1      densty( icox1, icoy1, 1, icgr) =
cp1  &   densty( icox1, icoy1, 1, icgr) + pmass1 / delx(icgr)**3 / 8.D0

c        because of symmetry around z=0: /2.D0
c     elseif ( ifgr .eq. indxgr(ltwo,maxlev(ltwo)) ) then
c      ** this part is now done in subroutine detoro ??

c        icox = icofi( nint(posx2), lone, ifgr )
c        icoy = icofi( nint(posy2), ltwo, ifgr )

c        densty( icox, icoy, lone, icgr) =
c    &   densty( icox, icoy, lone, icgr) + pmass2 / delx(icgr)**3 *half


c        velx  ( icox, icoy, lone, icgr) = velx2
c        vely  ( icox, icoy, lone, icgr) = vely2
c        energy( icox, icoy, lone, icgr) =
c    &   energy( icox, icoy, lone, icgr) + half*(velx2 ** 2 + vely2 ** 2)

cp1   endif

      return
      end
c-----------------------------------------------------------------------
      subroutine finzero(qtty,icgr)

c     zero all parts of the coarse grid qty with index icgr
c     which are overlayed by fine grids ifgr

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8 qtty(qx, qy, qz)

      !print *, 'grids, finzero'

      do 20 ifgr = 1, ngrd
         if ( idxcgr(ifgr) .eq. icgr ) then

            noz = norgin(ifgr,lthree)
            noy = norgin(ifgr,ltwo)
            nox = norgin(ifgr,lone)

            do 10 k = noz+1, noz+nz/nfine
            do 10 j = noy+1, noy+ny/nfine
            do 10 i = nox+1, nox+nx/nfine
 10            qtty(i,j,k) = zero

         endif
 20   continue

      return
      end
c-----------------------------------------------------------------------
      subroutine fluzero(igr)

      include 'qparam.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      !print *, 'grids, fluzero'

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,i),
C$OMP+            SHARED(flxrho,flxe,flxux,flxuy,flxuz,igr)
      do j = 1, qx
      do i = 1, qx
         flxrho(i,j,1,igr,1) =  zero
         flxrho(i,j,2,igr,1) =  zero
         flxe  (i,j,1,igr,1) =  zero
         flxe  (i,j,2,igr,1) =  zero
         flxux (i,j,1,igr,1) =  zero
         flxux (i,j,2,igr,1) =  zero
         flxuy (i,j,1,igr,1) =  zero
         flxuy (i,j,2,igr,1) =  zero
         flxuz (i,j,1,igr,1) =  zero
         flxuz (i,j,2,igr,1) =  zero
         flxrho(i,j,1,igr,2) =  zero
         flxrho(i,j,2,igr,2) =  zero
         flxe  (i,j,1,igr,2) =  zero
         flxe  (i,j,2,igr,2) =  zero
         flxux (i,j,1,igr,2) =  zero
         flxux (i,j,2,igr,2) =  zero
         flxuy (i,j,1,igr,2) =  zero
         flxuy (i,j,2,igr,2) =  zero
         flxuz (i,j,1,igr,2) =  zero
         flxuz (i,j,2,igr,2) =  zero
         flxrho(i,j,1,igr,3) =  zero
         flxrho(i,j,2,igr,3) =  zero
         flxe  (i,j,1,igr,3) =  zero
         flxe  (i,j,2,igr,3) =  zero
         flxux (i,j,1,igr,3) =  zero
         flxux (i,j,2,igr,3) =  zero
         flxuy (i,j,1,igr,3) =  zero
         flxuy (i,j,2,igr,3) =  zero
         flxuz (i,j,1,igr,3) =  zero
         flxuz (i,j,2,igr,3) =  zero
      enddo
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,i,m), SHARED(flxche,igr)
      do j = 1, qx
      do m = 1, qc
      do i = 1, qx
         flxche(i,j,1,igr,1,m) =  zero
         flxche(i,j,2,igr,1,m) =  zero
         flxche(i,j,1,igr,2,m) =  zero
         flxche(i,j,2,igr,2,m) =  zero
         flxche(i,j,1,igr,3,m) =  zero
         flxche(i,j,2,igr,3,m) =  zero
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine flupmem(ii,jj,igr,nsw,nswi,nswj)
c     PARALLEL OK

      include 'qparam.cmn'
      include 'squants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      integer*4 ii, jj, igr, nsw, nswi, nswj,
     &          nn(3), ifgr, ic1, ic2, jc1, jc2, kc1, kc2, i, j

      !print *, 'grids, flupmem'

c all these should be private

c --- save coarse fluxes at position of finer grid boundry onto
c     finer grid flxbnd

      ifgr = idxfgr(igr)
      nn(1) = nx
      nn(2) = ny
      nn(3) = nz

      if ( ifgr .ne. -lone ) then    ! so igr is not the finest

      ic1 = icofi(     lone, nswi, ifgr)
      ic2 = icofi( nn(nswi), nswi, ifgr)
      jc1 = icofi(     lone, nswj, ifgr)
      jc2 = icofi( nn(nswj), nswj, ifgr)

      if ( (ii.ge.ic1).and.(ii.le.ic2).and.
     &     (jj.ge.jc1).and.(jj.le.jc2) ) then

c         i = ii - ic1 + 1
c         j = jj - jc1 + 1
         i = (ii - ic1)*nfine + 1     ! specifically for parallelisation !
         j = (jj - jc1)*nfine + 1

         kc1 = icofi(    lone, nsw, ifgr)       + 4  ! the rows of PPM
         kc2 = icofi( nn(nsw), nsw, ifgr) + 1   + 4  ! are shifted by 4

         flxrho(i,j,1,ifgr,nsw) = flxrho(i,j,1,ifgr,nsw) + rhoflx(kc1)
         flxe  (i,j,1,ifgr,nsw) = flxe  (i,j,1,ifgr,nsw) +   eflx(kc1)

         flxrho(i,j,2,ifgr,nsw) = flxrho(i,j,2,ifgr,nsw) + rhoflx(kc2)
         flxe  (i,j,2,ifgr,nsw) = flxe  (i,j,2,ifgr,nsw) +   eflx(kc2)

         do m = 1, qc
         flxche(i,j,1,ifgr,nsw,m)=flxche(i,j,1,ifgr,nsw,m)+cheflx(kc1,m)
         flxche(i,j,2,ifgr,nsw,m)=flxche(i,j,2,ifgr,nsw,m)+cheflx(kc2,m)
         enddo

         if    ( nsw .eq. lone ) then

          flxux (i,j,1,ifgr,nsw) = flxux (i,j,1,ifgr,nsw) +   uflx(kc1)
          flxuy (i,j,1,ifgr,nsw) = flxuy (i,j,1,ifgr,nsw) +  utflx(kc1)
          flxuz (i,j,1,ifgr,nsw) = flxuz (i,j,1,ifgr,nsw) + uttflx(kc1)

          flxux (i,j,2,ifgr,nsw) = flxux (i,j,2,ifgr,nsw) +   uflx(kc2)
          flxuy (i,j,2,ifgr,nsw) = flxuy (i,j,2,ifgr,nsw) +  utflx(kc2)
          flxuz (i,j,2,ifgr,nsw) = flxuz (i,j,2,ifgr,nsw) + uttflx(kc2)

         elseif( nsw .eq. ltwo ) then

          flxuy (i,j,1,ifgr,nsw) = flxuy (i,j,1,ifgr,nsw) +   uflx(kc1)
          flxux (i,j,1,ifgr,nsw) = flxux (i,j,1,ifgr,nsw) +  utflx(kc1)
          flxuz (i,j,1,ifgr,nsw) = flxuz (i,j,1,ifgr,nsw) + uttflx(kc1)

          flxuy (i,j,2,ifgr,nsw) = flxuy (i,j,2,ifgr,nsw) +   uflx(kc2)
          flxux (i,j,2,ifgr,nsw) = flxux (i,j,2,ifgr,nsw) +  utflx(kc2)
          flxuz (i,j,2,ifgr,nsw) = flxuz (i,j,2,ifgr,nsw) + uttflx(kc2)

         elseif( nsw .eq. lthree ) then

          flxuz (i,j,1,ifgr,nsw) = flxuz (i,j,1,ifgr,nsw) +   uflx(kc1)
          flxux (i,j,1,ifgr,nsw) = flxux (i,j,1,ifgr,nsw) +  utflx(kc1)
          flxuy (i,j,1,ifgr,nsw) = flxuy (i,j,1,ifgr,nsw) + uttflx(kc1)

          flxuz (i,j,2,ifgr,nsw) = flxuz (i,j,2,ifgr,nsw) +   uflx(kc2)
          flxux (i,j,2,ifgr,nsw) = flxux (i,j,2,ifgr,nsw) +  utflx(kc2)
          flxuy (i,j,2,ifgr,nsw) = flxuy (i,j,2,ifgr,nsw) + uttflx(kc2)

         endif

      endif
      endif

c --- add difference coarse boundry fluxes onto coarse grid flxbnd

c     nfine timesteps of fine grid per coarse grid
c     nfine**2 fine fluxes for each coarse flux
      rfin3 = dble(nfine**3)

c     if ( idxcgr(igr) .ne. -lone ) then   ! so igr is not coarsest

c      i = (ii+lone)/nfine
c      j = (jj+lone)/nfine
      i = ii     ! specifically for parallelisation !
      j = jj     ! is then corrected in subroutine flupar !


      flxrho(i,j,1,igr,nsw) = flxrho(i,j,1,igr,nsw) - rhoflx(5) / rfin3
      flxe  (i,j,1,igr,nsw) = flxe  (i,j,1,igr,nsw) -   eflx(5) / rfin3

      flxrho(i,j,2,igr,nsw) = flxrho(i,j,2,igr,nsw) - rhoflx(np5)/rfin3
      flxe  (i,j,2,igr,nsw) = flxe  (i,j,2,igr,nsw) -   eflx(np5)/rfin3

      do m = 1, qc
      flxche(i,j,1,igr,nsw,m)=flxche(i,j,1,igr,nsw,m)-cheflx(5,m)/rfin3
      flxche(i,j,2,igr,nsw,m)=flxche(i,j,2,igr,nsw,m)
     &                                             -cheflx(np5,m)/rfin3
      enddo

      if    ( nsw .eq. lone ) then

      flxux (i,j,1,igr,nsw) = flxux (i,j,1,igr,nsw) -   uflx(5) / rfin3
      flxuy (i,j,1,igr,nsw) = flxuy (i,j,1,igr,nsw) -  utflx(5) / rfin3
      flxuz (i,j,1,igr,nsw) = flxuz (i,j,1,igr,nsw) - uttflx(5) / rfin3

      flxux (i,j,2,igr,nsw) = flxux (i,j,2,igr,nsw) -   uflx(np5)/rfin3
      flxuy (i,j,2,igr,nsw) = flxuy (i,j,2,igr,nsw) -  utflx(np5)/rfin3
      flxuz (i,j,2,igr,nsw) = flxuz (i,j,2,igr,nsw) - uttflx(np5)/rfin3

      elseif( nsw .eq. ltwo ) then

      flxuy (i,j,1,igr,nsw) = flxuy (i,j,1,igr,nsw) -   uflx(5) / rfin3
      flxux (i,j,1,igr,nsw) = flxux (i,j,1,igr,nsw) -  utflx(5) / rfin3
      flxuz (i,j,1,igr,nsw) = flxuz (i,j,1,igr,nsw) - uttflx(5) / rfin3

      flxuy (i,j,2,igr,nsw) = flxuy (i,j,2,igr,nsw) -   uflx(np5)/rfin3
      flxux (i,j,2,igr,nsw) = flxux (i,j,2,igr,nsw) -  utflx(np5)/rfin3
      flxuz (i,j,2,igr,nsw) = flxuz (i,j,2,igr,nsw) - uttflx(np5)/rfin3

      elseif( nsw .eq. lthree ) then

      flxuz (i,j,1,igr,nsw) = flxuz (i,j,1,igr,nsw) -   uflx(5) / rfin3
      flxux (i,j,1,igr,nsw) = flxux (i,j,1,igr,nsw) -  utflx(5) / rfin3
      flxuy (i,j,1,igr,nsw) = flxuy (i,j,1,igr,nsw) - uttflx(5) / rfin3

      flxuz (i,j,2,igr,nsw) = flxuz (i,j,2,igr,nsw) -   uflx(np5)/rfin3
      flxux (i,j,2,igr,nsw) = flxux (i,j,2,igr,nsw) -  utflx(np5)/rfin3
      flxuy (i,j,2,igr,nsw) = flxuy (i,j,2,igr,nsw) - uttflx(np5)/rfin3

      endif

c     endif  ! if igr is not coarsest grid

      return
      end
c-----------------------------------------------------------------------
      subroutine flupar(igr,nsw,nii,njj)
c     sum up 4 fluxes into one; this is not done in flupmem,
c     because awkward to parallelise

      include 'qparam.cmn'
      include 'squants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      integer*4 igr, nsw, nii, njj, ii, jj, i, j, im, jm, m

      real*8    trho(qx,qx,2),  te(qx,qx,2), tche(qx,qx,2,ngrd),
     &           tux(qx,qx,2), tuy(qx,qx,2),  tuz(qx,qx,2)

      !print *, 'grids, flupar'

c --- save coarse fluxes at position of finer grid boundry onto
c     finer grid flxbnd
c     nfine**2 fine fluxes for each coarse flux

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(jj,ii,i,j,im,jm,m),
C$OMP+            SHARED(igr,nsw,nii,njj,flxrho,flxe,flxche,
C$OMP+                   flxux,flxuy,flxuz,trho,te,tche,tux,tuy,tuz)
      do ii = 1, nii
         i = ii*nfine
         im = i - lone
      do jj = 1, njj
         j = jj*nfine
         jm = j - lone

      trho(ii,jj,1)=   flxrho(i ,j ,1,igr,nsw) + flxrho(i ,jm,1,igr,nsw)
     &               + flxrho(im,j ,1,igr,nsw) + flxrho(im,jm,1,igr,nsw)
      te  (ii,jj,1)=   flxe  (i ,j ,1,igr,nsw) + flxe  (i ,jm,1,igr,nsw)
     &               + flxe  (im,j ,1,igr,nsw) + flxe  (im,jm,1,igr,nsw)

      trho(ii,jj,2)=   flxrho(i ,j ,2,igr,nsw) + flxrho(i ,jm,2,igr,nsw)
     &               + flxrho(im,j ,2,igr,nsw) + flxrho(im,jm,2,igr,nsw)
      te  (ii,jj,2)=   flxe  (i ,j ,2,igr,nsw) + flxe  (i ,jm,2,igr,nsw)
     &               + flxe  (im,j ,2,igr,nsw) + flxe  (im,jm,2,igr,nsw)

      do m = 1, qc
      tche(ii,jj,1,m)=flxche(i ,j,1,igr,nsw,m)+flxche(i ,jm,1,igr,nsw,m)
     &               +flxche(im,j,1,igr,nsw,m)+flxche(im,jm,1,igr,nsw,m)
      tche(ii,jj,2,m)=flxche(i ,j,2,igr,nsw,m)+flxche(i ,jm,2,igr,nsw,m)
     &               +flxche(im,j,2,igr,nsw,m)+flxche(im,jm,2,igr,nsw,m)
      enddo

      tux(ii,jj,1)=   flxux(i ,j ,1,igr,nsw) + flxux(i ,jm,1,igr,nsw)
     &              + flxux(im,j ,1,igr,nsw) + flxux(im,jm,1,igr,nsw)
      tuy(ii,jj,1)=   flxuy(i ,j ,1,igr,nsw) + flxuy(i ,jm,1,igr,nsw)
     &              + flxuy(im,j ,1,igr,nsw) + flxuy(im,jm,1,igr,nsw)
      tuz(ii,jj,1)=   flxuz(i ,j ,1,igr,nsw) + flxuz(i ,jm,1,igr,nsw)
     &              + flxuz(im,j ,1,igr,nsw) + flxuz(im,jm,1,igr,nsw)

      tux(ii,jj,2)=   flxux(i ,j ,2,igr,nsw) + flxux(i ,jm,2,igr,nsw)
     &              + flxux(im,j ,2,igr,nsw) + flxux(im,jm,2,igr,nsw)
      tuy(ii,jj,2)=   flxuy(i ,j ,2,igr,nsw) + flxuy(i ,jm,2,igr,nsw)
     &              + flxuy(im,j ,2,igr,nsw) + flxuy(im,jm,2,igr,nsw)
      tuz(ii,jj,2)=   flxuz(i ,j ,2,igr,nsw) + flxuz(i ,jm,2,igr,nsw)
     &              + flxuz(im,j ,2,igr,nsw) + flxuz(im,jm,2,igr,nsw)

      enddo
      enddo
C$OMP END PARALLEL DO


C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(jj,ii,m),
C$OMP+            SHARED(igr,nsw,nii,njj,flxrho,flxe,flxche,
C$OMP+                   flxux,flxuy,flxuz,trho,te,tche,tux,tuy,tuz)
      do ii = 1, nii
      do jj = 1, njj

         flxrho(ii,jj,1,igr,nsw) = trho(ii,jj,1)
         flxrho(ii,jj,2,igr,nsw) = trho(ii,jj,2)

         flxe  (ii,jj,1,igr,nsw) = te  (ii,jj,1)
         flxe  (ii,jj,2,igr,nsw) = te  (ii,jj,2)

         do m = 1, qc
            flxche(ii,jj,1,igr,nsw,m) = tche(ii,jj,1,m)
            flxche(ii,jj,2,igr,nsw,m) = tche(ii,jj,2,m)
         enddo

         flxux(ii,jj,1,igr,nsw) = tux(ii,jj,1)
         flxuy(ii,jj,1,igr,nsw) = tuy(ii,jj,1)
         flxuz(ii,jj,1,igr,nsw) = tuz(ii,jj,1)

         flxux(ii,jj,2,igr,nsw) = tux(ii,jj,2)
         flxuy(ii,jj,2,igr,nsw) = tuy(ii,jj,2)
         flxuz(ii,jj,2,igr,nsw) = tuz(ii,jj,2)

      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine flupfix(ifgr,icgr)

c --- fix coarse grid cells surrunding fine grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      !print *, 'grids, flupfix'

      call flupar (ifgr, lone,   nyh, nzh)
      call flupar (ifgr, ltwo,   nxh, nzh)
      call flupar (ifgr, lthree, nxh, nyh)

      do i = 1, q
         dtdx(i) = dtgr(icgr) / delx(icgr)
      enddo

      dtdxs = dtdx(1)

      icf(1) = icofi( lone, lone,   ifgr) - 1
      icf(2) = icofi(   qx, lone,   ifgr) + 1
      jcf(1) = icofi( lone, ltwo,   ifgr) - 1
      jcf(2) = icofi(   qy, ltwo,   ifgr) + 1
      kcf(1) = icofi( lone, lthree, ifgr) - 1
      kcf(2) = icofi(   qz, lthree, ifgr) + 1

      call dexfix( icgr, ifgr, dtdxs )
      call deyfix( icgr, ifgr, dtdxs )
      call dezfix( icgr, ifgr, dtdxs )

c     flxrho contains old density values

      call qtxfix( energy(1,1,1,icgr), flxe  (1,1,1,ifgr,1), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,1)        )
      call qtyfix( energy(1,1,1,icgr), flxe  (1,1,1,ifgr,2), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,2)        )
      call qtzfix( energy(1,1,1,icgr), flxe  (1,1,1,ifgr,3), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,3)        )

      call qtxfix( velx  (1,1,1,icgr), flxux (1,1,1,ifgr,1), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,1)        )
      call qtyfix( velx  (1,1,1,icgr), flxux (1,1,1,ifgr,2), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,2)        )
      call qtzfix( velx  (1,1,1,icgr), flxux (1,1,1,ifgr,3), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,3)        )

      call qtxfix( vely  (1,1,1,icgr), flxuy (1,1,1,ifgr,1), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,1)        )
      call qtyfix( vely  (1,1,1,icgr), flxuy (1,1,1,ifgr,2), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,2)        )
      call qtzfix( vely  (1,1,1,icgr), flxuy (1,1,1,ifgr,3), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,3)        )

      call qtxfix( velz  (1,1,1,icgr), flxuz (1,1,1,ifgr,1), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,1)        )
      call qtyfix( velz  (1,1,1,icgr), flxuz (1,1,1,ifgr,2), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,2)        )
      call qtzfix( velz  (1,1,1,icgr), flxuz (1,1,1,ifgr,3), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,3)        )

      call ergfix(icgr)

      return
      end
c-----------------------------------------------------------------------
      subroutine qtxfix( qtty, flx, dtx, dnw, dol )

c --- fix coarse grid cells surrunding fine grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8 rest
      parameter ( rest = 0.9D0 )

      real*8 qtty( qx, qy, qz ), flx ( qx, qx, 2 )
      real*8 dnw ( qx, qy, qz ), dol ( qx, qx, 2 )

      real*8 qtd(q,q), cor(q,q)

      !print *, 'grids, qtxfix'

      do mp = 1, 2
         msgn = lthree-ltwo*mp

C$OMP PARALLEL DEFAULT(NONE), PRIVATE(j,k)
C$OMP+         SHARED(icf,jcf,kcf,qtty,flx,dnw,dol,qtd,cor,
C$OMP+                nyh,nzh,mp,msgn,dtx)
C$OMP DO
         do j = 1, nyh
         do k = 1, nzh
            qtd(j,k) =  qtty(icf(mp),jcf(1)+j,kcf(1)+k) * dol(j,k,mp)
            cor(j,k) =  dtx * msgn * flx(j,k,mp)
            if ( abs(qtd(j,k)) .lt. abs(cor(j,k)) )
     &           cor(j,k) = rest * sign( abs(qtd(j,k)), cor(j,k) )
         enddo
         enddo
C$OMP END DO
C$OMP DO
         do j = 1, nyh
         do k = 1, nzh
            qtty(icf(mp),jcf(1)+j,kcf(1)+k)
     &          = (qtd(j,k) + cor(j,k)) / dnw(icf(mp),jcf(1)+j,kcf(1)+k)
         enddo
         enddo
C$OMP END DO
C$OMP END PARALLEL

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine qtyfix( qtty, flx, dtx, dnw, dol )

c --- fix coarse grid cells surrunding fine grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8 rest
      parameter ( rest = 0.9D0 )

      real*8 qtty( qx, qy, qz ), flx ( qx, qx, 2 )
      real*8 dnw ( qx, qy, qz ), dol ( qx, qx, 2 )

      real*8 qtd(q,q), cor(q,q)

      !print *, 'grids, qtyfix'

      do mp = 1, 2
         msgn = lthree-ltwo*mp

C$OMP PARALLEL DEFAULT(NONE), PRIVATE(i,k)
C$OMP+         SHARED(icf,jcf,kcf,qtty,flx,dnw,dol,qtd,cor,
C$OMP+                nxh,nzh,mp,msgn,dtx)
C$OMP DO
         do i = 1, nxh
         do k = 1, nzh
            qtd(i,k) =  qtty(icf(1)+i,jcf(mp),kcf(1)+k) * dol(i,k,mp)
            cor(i,k) =  dtx * msgn * flx(i,k,mp)
            if ( abs(qtd(i,k)) .lt. abs(cor(i,k)) )
     &           cor(i,k) = rest * sign( abs(qtd(i,k)), cor(i,k) )
         enddo
         enddo
C$OMP END DO
C$OMP DO
         do i = 1, nxh
         do k = 1, nzh
            qtty(icf(1)+i,jcf(mp),kcf(1)+k)
     &          = (qtd(i,k) + cor(i,k)) / dnw(icf(1)+i,jcf(mp),kcf(1)+k)
         enddo
         enddo
C$OMP END DO
C$OMP END PARALLEL

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine qtzfix( qtty, flx, dtx, dnw, dol )

c --- fix coarse grid cells surrunding fine grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8 rest
      parameter ( rest = 0.9D0 )

      real*8 qtty( qx, qy, qz ), flx ( qx, qx, 2 )
      real*8 dnw ( qx, qy, qz ), dol ( qx, qx, 2 )

      real*8 qtd(q,q), cor(q,q)

      !print *, 'grids, qtzfix'

      do mp = 2, 2      !  1 is not needed: reflecting boundry
         msgn = lthree-ltwo*mp

C$OMP PARALLEL DEFAULT(NONE), PRIVATE(i,j)
C$OMP+         SHARED(icf,jcf,kcf,qtty,flx,dnw,dol,qtd,cor,
C$OMP+                nxh,nyh,mp,msgn,dtx)
C$OMP DO
         do j = 1, nyh
         do i = 1, nxh
            qtd(i,j) =  qtty(icf(1)+i,jcf(1)+j,kcf(mp)) * dol(i,j,mp)
            cor(i,j) =  dtx * msgn * flx(i,j,mp)
            if ( abs(qtd(i,j)) .lt. abs(cor(i,j)) )
     &           cor(i,j) = rest * sign( abs(qtd(i,j)), cor(i,j) )
         enddo
         enddo
C$OMP END DO
C$OMP DO
         do j = 1, nyh
         do i = 1, nxh
            qtty(icf(1)+i,jcf(1)+j,kcf(mp))
     &          = (qtd(i,j) + cor(i,j)) / dnw(icf(1)+i,jcf(1)+j,kcf(mp))
         enddo
         enddo
C$OMP END DO
C$OMP END PARALLEL

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dexfix( icgr, ifgr, dtx )

c --- fix coarse grid cells surrunding fine grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8 den(q,q), cor(q,q), fac(q,q), rest
      parameter ( rest = 0.9D0 )

      !print *, 'grids, dexfix'


      do 10 mp = 1, 2
         msgn = lthree-ltwo*mp

C$OMP PARALLEL DEFAULT(NONE), PRIVATE(j,k),
C$OMP+         SHARED(nyh,nzh,mp,msgn,dtx,icgr,ifgr,icf,jcf,kcf,
C$OMP+               densty,den,cor,fac,flxrho,flxe,flxux,flxuy,flxuz)
C$OMP DO
         do j = 1, nyh
         do k = 1, nzh
            den(j,k) = densty(icf(mp),jcf(1)+j,kcf(1)+k,icgr)
            cor(j,k) = dtx * msgn * flxrho(j,k,mp,ifgr,1)
            fac(j,k) = one
         enddo
         enddo
C$OMP END DO
C$OMP DO
         do j = 1, nyh
         do k = 1, nzh
            if( rest*den(j,k) .lt. -cor(j,k) ) then
               fac(j,k) = rest*den(j,k) / (-cor(j,k))
               flxe  (j,k,mp,ifgr,1) = fac(j,k) * flxe  (j,k,mp,ifgr,1)
               flxux (j,k,mp,ifgr,1) = fac(j,k) * flxux (j,k,mp,ifgr,1)
               flxuy (j,k,mp,ifgr,1) = fac(j,k) * flxuy (j,k,mp,ifgr,1)
               flxuz (j,k,mp,ifgr,1) = fac(j,k) * flxuz (j,k,mp,ifgr,1)
               cor(j,k) = fac(j,k) * cor(j,k)
            endif
            if(fac(j,k).ne.one)write(*,*) 'dexfix ',j,k,mp,icgr,fac(j,k)
            flxrho(j,k,mp,ifgr,1) = den(j,k)
            densty(icf(mp),jcf(1)+j,kcf(1)+k,icgr) = den(j,k) + cor(j,k)
         enddo
         enddo
C$OMP END DO
C$OMP END PARALLEL

c --- now chem

      do m = 1, qc
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k),
C$OMP+            SHARED(icf,jcf,kcf,chem,den,cor,fac,flxche,
C$OMP+                   nyh,nzh,mp,msgn,dtx,icgr,ifgr,m)
         do j = 1, nyh
         do k = 1, nzh
            den(j,k) = chem(icf(mp),jcf(1)+j,kcf(1)+k,icgr,m)
            cor(j,k) = fac(j,k) * dtx * msgn * flxche(j,k,mp,ifgr,1,m)
            if(rest*den(j,k) .lt. -cor(j,k)) cor(j,k) = -rest * den(j,k)
            chem(icf(mp),jcf(1)+j,kcf(1)+k,icgr,m) = den(j,k) + cor(j,k)
         enddo
         enddo
C$OMP END PARALLEL DO
      enddo

c---

 10   continue


      return
      end
c-----------------------------------------------------------------------
      subroutine deyfix( icgr, ifgr, dtx )

c --- fix coarse grid cells surrunding fine grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8 den(q,q), cor(q,q), fac(q,q), rest
      parameter ( rest = 0.9D0 )

      !print *, 'grids, deyfix'


      do 10 mp = 1, 2
         msgn = lthree-ltwo*mp

C$OMP PARALLEL DEFAULT(NONE), PRIVATE(i,k),
C$OMP+         SHARED(icf,jcf,kcf,
C$OMP+               densty,den,cor,fac,flxrho,flxe,flxux,flxuy,flxuz,
C$OMP+               nxh,nzh,mp,msgn,dtx,icgr,ifgr)
C$OMP DO
         do i = 1, nxh
         do k = 1, nzh
            den(i,k) = densty(icf(1)+i,jcf(mp),kcf(1)+k,icgr)
            cor(i,k) = dtx * msgn * flxrho(i,k,mp,ifgr,2)
            fac(i,k) = one
         enddo
         enddo
C$OMP END DO
C$OMP DO
         do i = 1, nxh
         do k = 1, nzh
            if( rest*den(i,k) .lt. -cor(i,k) ) then
               fac(i,k) = rest*den(i,k) / (-cor(i,k))
               flxe  (i,k,mp,ifgr,2) = fac(i,k) * flxe  (i,k,mp,ifgr,2)
               flxux (i,k,mp,ifgr,2) = fac(i,k) * flxux (i,k,mp,ifgr,2)
               flxuy (i,k,mp,ifgr,2) = fac(i,k) * flxuy (i,k,mp,ifgr,2)
               flxuz (i,k,mp,ifgr,2) = fac(i,k) * flxuz (i,k,mp,ifgr,2)
               cor(i,k) = fac(i,k) * cor(i,k)
            endif
            if(fac(i,k).ne.one)write(*,*) 'deyfix ',i,k,mp,icgr,fac(i,k)
            flxrho(i,k,mp,ifgr,2) = den(i,k)
            densty(icf(1)+i,jcf(mp),kcf(1)+k,icgr) = den(i,k) + cor(i,k)
         enddo
         enddo
C$OMP END DO
C$OMP END PARALLEL

c --- now chem

      do m = 1, qc
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,k),
C$OMP+            SHARED(icf,jcf,kcf,chem,den,cor,fac,flxche,
C$OMP+                   nxh,nzh,mp,msgn,dtx,icgr,ifgr,m)
         do i = 1, nxh
         do k = 1, nzh
            den(i,k) = chem(icf(1)+i,jcf(mp),kcf(1)+k,icgr,m)
            cor(i,k) = fac(i,k) * dtx * msgn * flxche(i,k,mp,ifgr,2,m)
            if(rest*den(i,k) .lt. -cor(i,k)) cor(i,k) = -rest * den(i,k)
            chem(icf(1)+i,jcf(mp),kcf(1)+k,icgr,m) = den(i,k) + cor(i,k)
         enddo
         enddo
C$OMP END PARALLEL DO
      enddo

c ---

 10   continue


      return
      end
c-----------------------------------------------------------------------
      subroutine dezfix( icgr, ifgr, dtx )

c --- fix coarse grid cells surrunding fine grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      real*8 den(q,q), cor(q,q), fac(q,q), rest
      parameter ( rest = 0.9D0 )

      !print *, 'grids, dezfix'


      do 10 mp = 2, 2      !  1 is not needed: reflecting boundry
         msgn = lthree-ltwo*mp
C$OMP PARALLEL DEFAULT(NONE), PRIVATE(j,i),
C$OMP+         SHARED(nxh,nyh,mp,msgn,dtx,icgr,ifgr,icf,jcf,kcf,
C$OMP+               densty,den,cor,fac,flxrho,flxe,flxux,flxuy,flxuz)
C$OMP DO
         do j = 1, nyh
         do i = 1, nxh
            den(i,j) = densty(icf(1)+i,jcf(1)+j,kcf(mp),icgr)
            cor(i,j) = dtx * msgn * flxrho(i,j,mp,ifgr,3)
            fac(i,j) = one
         enddo
         enddo
C$OMP END DO
C$OMP DO
         do j = 1, nyh
         do i = 1, nxh
            if( rest*den(i,j) .lt. -cor(i,j) ) then
               fac(i,j) = rest*den(i,j) / (-cor(i,j))
               flxe  (i,j,mp,ifgr,3) = fac(i,j) * flxe  (i,j,mp,ifgr,3)
               flxux (i,j,mp,ifgr,3) = fac(i,j) * flxux (i,j,mp,ifgr,3)
               flxuy (i,j,mp,ifgr,3) = fac(i,j) * flxuy (i,j,mp,ifgr,3)
               flxuz (i,j,mp,ifgr,3) = fac(i,j) * flxuz (i,j,mp,ifgr,3)
               cor(i,j) = fac(i,j) * cor(i,j)
            endif
            flxrho(i,j,mp,ifgr,3) = den(i,j)
            densty(icf(1)+i,jcf(1)+j,kcf(mp),icgr) = den(i,j) + cor(i,j)
         enddo
         enddo
C$OMP END DO
C$OMP END PARALLEL

c --- now chem

      do m = 1, qc
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,i),
C$OMP+            SHARED(icf,jcf,kcf,chem,den,cor,fac,flxche,
C$OMP+            nyh,nxh,mp,msgn,dtx,icgr,ifgr,m)
         do j = 1, nyh
         do i = 1, nxh
            den(i,j) = chem(icf(1)+i,jcf(1)+j,kcf(mp),icgr,m)
            cor(i,j) = fac(i,j) * dtx * msgn * flxche(i,j,mp,ifgr,3,m)
            if(rest*den(i,j) .lt. -cor(i,j)) cor(i,j) = -rest * den(i,j)
            chem(icf(1)+i,jcf(1)+j,kcf(mp),icgr,m) = den(i,j) + cor(i,j)
         enddo
         enddo
C$OMP END PARALLEL DO
      enddo

c---

 10   continue


      return
      end
c-----------------------------------------------------------------------
      subroutine ergfix( icgr )

c --- fix coarse grid cells surrunding fine grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'
      include   'ppms.cmn'

      real*8 ekp(q,q), ep(q,q), enup(q,q)

      !print *, 'grids, ergfix'

      do 10 mp= 1, 2
C$OMP PARALLEL DEFAULT(NONE), PRIVATE(j,k),
C$OMP+         SHARED(icgr,mp,nyh,nzh,small,kcf,jcf,icf,ekp,ep,enup,
C$OMP+                energy,velx,vely,velz,vxold,vyold,vzold,eold)
C$OMP DO
         do j = 1, nyh
         do k = 1, nzh
            ekp (j,k) =  half * (
     &                   velx(icf(mp),jcf(1)+j,kcf(1)+k,icgr)**2
     &                +  vely(icf(mp),jcf(1)+j,kcf(1)+k,icgr)**2
     &                +  velz(icf(mp),jcf(1)+j,kcf(1)+k,icgr)**2 )
            ep  (j,k) = energy(icf(mp),jcf(1)+j,kcf(1)+k,icgr)
            enup(j,k) = half * 0.25D0 * (
     &                  vxold(icf(mp)+1,jcf(1)+j+1,kcf(1)+k,icgr)**2
     &                + vyold(icf(mp)+1,jcf(1)+j+1,kcf(1)+k,icgr)**2
     &                + vzold(icf(mp)+1,jcf(1)+j+1,kcf(1)+k,icgr)**2
     &             +    vxold(icf(mp)-1,jcf(1)+j+1,kcf(1)+k,icgr)**2
     &                + vyold(icf(mp)-1,jcf(1)+j+1,kcf(1)+k,icgr)**2
     &                + vzold(icf(mp)-1,jcf(1)+j+1,kcf(1)+k,icgr)**2
     &             +    vxold(icf(mp)+1,jcf(1)+j-1,kcf(1)+k,icgr)**2
     &                + vyold(icf(mp)+1,jcf(1)+j-1,kcf(1)+k,icgr)**2
     &                + vzold(icf(mp)+1,jcf(1)+j-1,kcf(1)+k,icgr)**2
     &             +    vxold(icf(mp)-1,jcf(1)+j-1,kcf(1)+k,icgr)**2
     &                + vyold(icf(mp)-1,jcf(1)+j-1,kcf(1)+k,icgr)**2
     &                + vzold(icf(mp)-1,jcf(1)+j-1,kcf(1)+k,icgr)**2 )
         enddo
         enddo
C$OMP END DO
C$OMP DO
         do j = 1, nyh
         do k = 1, nzh
            if ( ep(j,k) .lt. (one+small)*ekp(j,k) ) then
               energy(icf(mp),jcf(1)+j,kcf(1)+k,icgr) = ekp(j,k) +
     &            max( ((eold(icf(mp)+1,jcf(1)+j+1,kcf(1)+k,icgr)+
     &                eold(icf(mp)-1,jcf(1)+j+1,kcf(1)+k,icgr)) +
     &               (eold(icf(mp)+1,jcf(1)+j-1,kcf(1)+k,icgr)+
     &                eold(icf(mp)-1,jcf(1)+j-1,kcf(1)+k,icgr)) )*0.25D0
     &                -enup(j,k), small*ekp(j,k) )
            endif
            if( ep(j,k) .lt. 0.95D0 *ekp(j,k) ) then
               write(*,'(a,4i4,1P,2E19.10)')
     &            'flufix1 energy:',mp,j,k,icgr, ekp(j,k),ep(j,k)
            endif
        enddo
        enddo
C$OMP END DO
C$OMP END PARALLEL
 10   continue


      do 20 mp= 1, 2
C$OMP PARALLEL DEFAULT(NONE), PRIVATE(k,i),
C$OMP+         SHARED(icgr,mp,nxh,nzh,small,kcf,jcf,icf,ekp,ep,enup,
C$OMP+                energy,velx,vely,velz,vxold,vyold,vzold,eold)
C$OMP DO
         do i = 1, nxh
         do k = 1, nzh
            ekp (i,k) =  half * (
     &                   velx(icf(1)+i,jcf(mp),kcf(1)+k,icgr)**2
     &                +  vely(icf(1)+i,jcf(mp),kcf(1)+k,icgr)**2
     &                +  velz(icf(1)+i,jcf(mp),kcf(1)+k,icgr)**2 )
            ep  (i,k) = energy(icf(1)+i,jcf(mp),kcf(1)+k,icgr)
            enup(i,k) = half * 0.25D0 * (
     &                  vxold(icf(1)+i+1,jcf(mp)+1,kcf(1)+k,icgr)**2
     &                + vyold(icf(1)+i+1,jcf(mp)+1,kcf(1)+k,icgr)**2
     &                + vzold(icf(1)+i+1,jcf(mp)+1,kcf(1)+k,icgr)**2
     &             +    vxold(icf(1)+i-1,jcf(mp)+1,kcf(1)+k,icgr)**2
     &                + vyold(icf(1)+i-1,jcf(mp)+1,kcf(1)+k,icgr)**2
     &                + vzold(icf(1)+i-1,jcf(mp)+1,kcf(1)+k,icgr)**2
     &             +    vxold(icf(1)+i+1,jcf(mp)-1,kcf(1)+k,icgr)**2
     &                + vyold(icf(1)+i+1,jcf(mp)-1,kcf(1)+k,icgr)**2
     &                + vzold(icf(1)+i+1,jcf(mp)-1,kcf(1)+k,icgr)**2
     &             +    vxold(icf(1)+i-1,jcf(mp)-1,kcf(1)+k,icgr)**2
     &                + vyold(icf(1)+i-1,jcf(mp)-1,kcf(1)+k,icgr)**2
     &                + vzold(icf(1)+i-1,jcf(mp)-1,kcf(1)+k,icgr)**2 )
         enddo
         enddo
C$OMP END DO
C$OMP DO
         do i = 1, nxh
         do k = 1, nzh
            if ( ep(i,k) .lt. (one+small)*ekp(i,k) ) then
               energy(icf(1)+i,jcf(mp),kcf(1)+k,icgr) = ekp(i,k) +
     &            max( ((eold(icf(1)+i+1,jcf(mp)+1,kcf(1)+k,icgr)+
     &                eold(icf(1)+i-1,jcf(mp)+1,kcf(1)+k,icgr)) +
     &               (eold(icf(1)+i+1,jcf(mp)-1,kcf(1)+k,icgr)+
     &                eold(icf(1)+i-1,jcf(mp)-1,kcf(1)+k,icgr)) )*0.25D0
     &                -enup(i,k), small*ekp(i,k) )
            endif
            if( ep(i,k) .lt. 0.95D0 *ekp(i,k) ) then
               write(*,'(a,4i4,1P,2E19.10)')
     &        'flufix2 energy:',i,mp,k,icgr, ekp(i,k),ep(i,k)
            endif
         enddo
         enddo
C$OMP END DO
C$OMP END PARALLEL
 20   continue


      do 30 mp= 2, 2    ! reflecting z=0 - boundary
C$OMP PARALLEL DEFAULT(NONE), PRIVATE(j,i),
C$OMP+         SHARED(icgr,mp,nxh,nyh,small,kcf,jcf,icf,ekp,ep,enup,
C$OMP+                energy,velx,vely,velz,vxold,vyold,vzold,eold)
C$OMP DO
         do j = 1, nyh
         do i = 1, nxh
            ekp (i,j) =  half * (
     &                   velx(icf(1)+i,jcf(1)+j,kcf(mp),icgr)**2
     &                +  vely(icf(1)+i,jcf(1)+j,kcf(mp),icgr)**2
     &                +  velz(icf(1)+i,jcf(1)+j,kcf(mp),icgr)**2 )
            ep  (i,j) = energy(icf(1)+i,jcf(1)+j,kcf(mp),icgr)
            enup(i,j) = half * 0.25D0 * (
     &                  vxold(icf(1)+i+1,jcf(1)+j+1,kcf(mp),icgr)**2
     &                + vyold(icf(1)+i+1,jcf(1)+j+1,kcf(mp),icgr)**2
     &                + vzold(icf(1)+i+1,jcf(1)+j+1,kcf(mp),icgr)**2
     &             +    vxold(icf(1)+i-1,jcf(1)+j+1,kcf(mp),icgr)**2
     &                + vyold(icf(1)+i-1,jcf(1)+j+1,kcf(mp),icgr)**2
     &                + vzold(icf(1)+i-1,jcf(1)+j+1,kcf(mp),icgr)**2
     &             +    vxold(icf(1)+i+1,jcf(1)+j-1,kcf(mp),icgr)**2
     &                + vyold(icf(1)+i+1,jcf(1)+j-1,kcf(mp),icgr)**2
     &                + vzold(icf(1)+i+1,jcf(1)+j-1,kcf(mp),icgr)**2
     &             +    vxold(icf(1)+i-1,jcf(1)+j-1,kcf(mp),icgr)**2
     &                + vyold(icf(1)+i-1,jcf(1)+j-1,kcf(mp),icgr)**2
     &                + vzold(icf(1)+i-1,jcf(1)+j-1,kcf(mp),icgr)**2 )
         enddo
         enddo
C$OMP END DO
C$OMP DO
         do j = 1, nyh
         do i = 1, nxh
            if ( ep(i,j) .lt. (one+small)*ekp(i,j) ) then
               energy(icf(1)+i,jcf(1)+j,kcf(mp),icgr) = ekp(i,j) +
     &         max( ((eold(icf(1)+i+1,jcf(1)+j+1,kcf(mp),icgr)+
     &                eold(icf(1)+i-1,jcf(1)+j+1,kcf(mp),icgr)) +
     &               (eold(icf(1)+i+1,jcf(1)+j-1,kcf(mp),icgr)+
     &                eold(icf(1)+i-1,jcf(1)+j-1,kcf(mp),icgr)) )*0.25D0
     &                -enup(i,j), small*ekp(i,j) )
            endif
            if( ep(i,j) .lt. 0.95D0 *ekp(i,j) ) then
               write(*,'(a,4i4,1P,2E19.10)')
     &        'flufix3 energy:',i,j,mp,icgr, ekp(i,j),ep(i,j)
            endif
        enddo
        enddo
C$OMP END DO
C$OMP END PARALLEL
 30   continue


      return
      end
c-----------------------------------------------------------------------
      subroutine rstzero

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'bounds.cmn'

      !print *, 'grids, rstzero'

      do igrd = 1, ngrd
         call fluzero(igrd)
         velmax(igrd) = one
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i), SHARED(gchi,igrd)
         do k = 1, qz
         do j = 1, qy
         do i = 1, qx
            gchi(i,j,k,igrd) = zero   !     just in case nconoff = 0
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO
      enddo

      blackhole = .false.
c      newwr     = .true.  !  set this variable for Neutrino output
      newwr     = .false.
      dij2dt = zero
      dij3dt = zero
      elot   = zero
      gmunb  = zero
      angunb = zero
      enelo   = zero
      enalo   = zero
      enxlo   = zero
      anelo   = zero
      analo   = zero
      anxlo   = zero

      sumbulk  = zero
      sumshear = zero
      sumneutr = zero
      sumshock = zero

c------- initialize background potentials, topmost level has no background

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i),
C$OMP+            SHARED(bakpot,baktpo,bakchi,bakrg)
      do j = 1, qy
      do k = 1, qz
      do i = 1, qx
         bakpot(i,j,k,1) = zero
         baktpo(i,j,k,1) = zero
         bakchi(i,j,k,1) = zero
         bakrg (i,j,k,1) = zero
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c     rest of potential is initialised in inigrav

      return
      end
c-----------------------------------------------------------------------
      subroutine flumem(ii,jj,igr,nsw,nswi,nswj)

      include 'qparam.cmn'
      include 'squants.cmn'
      include   'ppms.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      integer*4 nn(3)

      !print *, 'grids, flumem'

c --- save coarse fluxes at position of finer grid boundry onto
c     finer grid flxbnd

      ifgr = idxfgr(igr)
      nn(1) = nx
      nn(2) = ny
      nn(3) = nz

      if ( ifgr .ne. -lone ) then    ! so igr is not the finest

      ic1 = icofi(     lone, nswi, ifgr)
      ic2 = icofi( nn(nswi), nswi, ifgr)
      jc1 = icofi(     lone, nswj, ifgr)
      jc2 = icofi( nn(nswj), nswj, ifgr)

      if ( (ii.ge.ic1).and.(ii.le.ic2).and.
     &     (jj.ge.jc1).and.(jj.le.jc2) ) then

         i = ii - ic1 + 1
         j = jj - jc1 + 1

         kc1 = icofi(    lone, nsw, ifgr)       + 4  ! the rows of PPM
         kc2 = icofi( nn(nsw), nsw, ifgr) + 1   + 4  ! are shifted by 4

         flxrho(i,j,1,ifgr,nsw) = flxrho(i,j,1,ifgr,nsw) + rhoflx(kc1)
         flxe  (i,j,1,ifgr,nsw) = flxe  (i,j,1,ifgr,nsw) +   eflx(kc1)

         flxrho(i,j,2,ifgr,nsw) = flxrho(i,j,2,ifgr,nsw) + rhoflx(kc2)
         flxe  (i,j,2,ifgr,nsw) = flxe  (i,j,2,ifgr,nsw) +   eflx(kc2)

         do m = 1, qc
         flxche(i,j,1,ifgr,nsw,m)=flxche(i,j,1,ifgr,nsw,m)+cheflx(kc1,m)
         flxche(i,j,2,ifgr,nsw,m)=flxche(i,j,2,ifgr,nsw,m)+cheflx(kc2,m)
         enddo

         if    ( nsw .eq. lone ) then

          flxux (i,j,1,ifgr,nsw) = flxux (i,j,1,ifgr,nsw) +   uflx(kc1)
          flxuy (i,j,1,ifgr,nsw) = flxuy (i,j,1,ifgr,nsw) +  utflx(kc1)
          flxuz (i,j,1,ifgr,nsw) = flxuz (i,j,1,ifgr,nsw) + uttflx(kc1)

          flxux (i,j,2,ifgr,nsw) = flxux (i,j,2,ifgr,nsw) +   uflx(kc2)
          flxuy (i,j,2,ifgr,nsw) = flxuy (i,j,2,ifgr,nsw) +  utflx(kc2)
          flxuz (i,j,2,ifgr,nsw) = flxuz (i,j,2,ifgr,nsw) + uttflx(kc2)

         elseif( nsw .eq. ltwo ) then

          flxuy (i,j,1,ifgr,nsw) = flxuy (i,j,1,ifgr,nsw) +   uflx(kc1)
          flxux (i,j,1,ifgr,nsw) = flxux (i,j,1,ifgr,nsw) +  utflx(kc1)
          flxuz (i,j,1,ifgr,nsw) = flxuz (i,j,1,ifgr,nsw) + uttflx(kc1)

          flxuy (i,j,2,ifgr,nsw) = flxuy (i,j,2,ifgr,nsw) +   uflx(kc2)
          flxux (i,j,2,ifgr,nsw) = flxux (i,j,2,ifgr,nsw) +  utflx(kc2)
          flxuz (i,j,2,ifgr,nsw) = flxuz (i,j,2,ifgr,nsw) + uttflx(kc2)

         elseif( nsw .eq. lthree ) then

          flxuz (i,j,1,ifgr,nsw) = flxuz (i,j,1,ifgr,nsw) +   uflx(kc1)
          flxux (i,j,1,ifgr,nsw) = flxux (i,j,1,ifgr,nsw) +  utflx(kc1)
          flxuy (i,j,1,ifgr,nsw) = flxuy (i,j,1,ifgr,nsw) + uttflx(kc1)

          flxuz (i,j,2,ifgr,nsw) = flxuz (i,j,2,ifgr,nsw) +   uflx(kc2)
          flxux (i,j,2,ifgr,nsw) = flxux (i,j,2,ifgr,nsw) +  utflx(kc2)
          flxuy (i,j,2,ifgr,nsw) = flxuy (i,j,2,ifgr,nsw) + uttflx(kc2)

         endif

      endif
      endif

c --- add difference coarse boundry fluxes onto coarse grid flxbnd

c     nfine timesteps of fine grid per coarse grid
c     nfine**2 fine fluxes for each coarse flux
      rfin3 = dble(nfine**3)

c     if ( idxcgr(igr) .ne. -lone ) then   ! so igr is not coarsest

      i = (ii+lone)/nfine
      j = (jj+lone)/nfine

      flxrho(i,j,1,igr,nsw) = flxrho(i,j,1,igr,nsw) - rhoflx(5) / rfin3
      flxe  (i,j,1,igr,nsw) = flxe  (i,j,1,igr,nsw) -   eflx(5) / rfin3

      flxrho(i,j,2,igr,nsw) = flxrho(i,j,2,igr,nsw) - rhoflx(np5)/rfin3
      flxe  (i,j,2,igr,nsw) = flxe  (i,j,2,igr,nsw) -   eflx(np5)/rfin3

      do m = 1, qc
      flxche(i,j,1,igr,nsw,m)=flxche(i,j,1,igr,nsw,m)-cheflx(5,m)/rfin3
      flxche(i,j,2,igr,nsw,m)=flxche(i,j,2,igr,nsw,m)
     &                                             -cheflx(np5,m)/rfin3
      enddo

      if    ( nsw .eq. lone ) then

      flxux (i,j,1,igr,nsw) = flxux (i,j,1,igr,nsw) -   uflx(5) / rfin3
      flxuy (i,j,1,igr,nsw) = flxuy (i,j,1,igr,nsw) -  utflx(5) / rfin3
      flxuz (i,j,1,igr,nsw) = flxuz (i,j,1,igr,nsw) - uttflx(5) / rfin3

      flxux (i,j,2,igr,nsw) = flxux (i,j,2,igr,nsw) -   uflx(np5)/rfin3
      flxuy (i,j,2,igr,nsw) = flxuy (i,j,2,igr,nsw) -  utflx(np5)/rfin3
      flxuz (i,j,2,igr,nsw) = flxuz (i,j,2,igr,nsw) - uttflx(np5)/rfin3

      elseif( nsw .eq. ltwo ) then

      flxuy (i,j,1,igr,nsw) = flxuy (i,j,1,igr,nsw) -   uflx(5) / rfin3
      flxux (i,j,1,igr,nsw) = flxux (i,j,1,igr,nsw) -  utflx(5) / rfin3
      flxuz (i,j,1,igr,nsw) = flxuz (i,j,1,igr,nsw) - uttflx(5) / rfin3

      flxuy (i,j,2,igr,nsw) = flxuy (i,j,2,igr,nsw) -   uflx(np5)/rfin3
      flxux (i,j,2,igr,nsw) = flxux (i,j,2,igr,nsw) -  utflx(np5)/rfin3
      flxuz (i,j,2,igr,nsw) = flxuz (i,j,2,igr,nsw) - uttflx(np5)/rfin3

      elseif( nsw .eq. lthree ) then

      flxuz (i,j,1,igr,nsw) = flxuz (i,j,1,igr,nsw) -   uflx(5) / rfin3
      flxux (i,j,1,igr,nsw) = flxux (i,j,1,igr,nsw) -  utflx(5) / rfin3
      flxuy (i,j,1,igr,nsw) = flxuy (i,j,1,igr,nsw) - uttflx(5) / rfin3

      flxuz (i,j,2,igr,nsw) = flxuz (i,j,2,igr,nsw) -   uflx(np5)/rfin3
      flxux (i,j,2,igr,nsw) = flxux (i,j,2,igr,nsw) -  utflx(np5)/rfin3
      flxuy (i,j,2,igr,nsw) = flxuy (i,j,2,igr,nsw) - uttflx(np5)/rfin3

      endif
c     endif  ! if igr is not coarsest grid

      return
      end
c-----------------------------------------------------------------------
      subroutine flufix(ifgr,icgr)

c --- fix coarse grid cells surrunding fine grid

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'bounds.cmn'

      !print *, 'grids, flufix'

      do i = 1, q
         dtdx(i) = dtgr(icgr) / delx(icgr)
      enddo

      dtdxs = dtdx(1)

      icf(1) = icofi( lone, lone,   ifgr) - 1
      icf(2) = icofi(   qx, lone,   ifgr) + 1
      jcf(1) = icofi( lone, ltwo,   ifgr) - 1
      jcf(2) = icofi(   qy, ltwo,   ifgr) + 1
      kcf(1) = icofi( lone, lthree, ifgr) - 1
      kcf(2) = icofi(   qz, lthree, ifgr) + 1

      call dexfix( icgr, ifgr, dtdxs )
      call deyfix( icgr, ifgr, dtdxs )
      call dezfix( icgr, ifgr, dtdxs )

c     flxrho contains old density values

      call qtxfix( energy(1,1,1,icgr), flxe  (1,1,1,ifgr,1), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,1)        )
      call qtyfix( energy(1,1,1,icgr), flxe  (1,1,1,ifgr,2), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,2)        )
      call qtzfix( energy(1,1,1,icgr), flxe  (1,1,1,ifgr,3), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,3)        )

      call qtxfix( velx  (1,1,1,icgr), flxux (1,1,1,ifgr,1), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,1)        )
      call qtyfix( velx  (1,1,1,icgr), flxux (1,1,1,ifgr,2), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,2)        )
      call qtzfix( velx  (1,1,1,icgr), flxux (1,1,1,ifgr,3), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,3)        )

      call qtxfix( vely  (1,1,1,icgr), flxuy (1,1,1,ifgr,1), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,1)        )
      call qtyfix( vely  (1,1,1,icgr), flxuy (1,1,1,ifgr,2), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,2)        )
      call qtzfix( vely  (1,1,1,icgr), flxuy (1,1,1,ifgr,3), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,3)        )

      call qtxfix( velz  (1,1,1,icgr), flxuz (1,1,1,ifgr,1), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,1)        )
      call qtyfix( velz  (1,1,1,icgr), flxuz (1,1,1,ifgr,2), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,2)        )
      call qtzfix( velz  (1,1,1,icgr), flxuz (1,1,1,ifgr,3), dtdxs,
     &             densty(1,1,1,icgr), flxrho(1,1,1,ifgr,3)        )

      call ergfix(icgr)

      return
      end
