      program grb
c-----------------------------------------------------------------------
c    Version graviwaves backreaction, neutrino emission (Dec 1995)
c-----------------------------------------------------------------------
c
c                     -----------------------
c                     |  C H A R Y B D I S  |
c                     -----------------------
c
c-----------------------------------------------------------------------
c    Code for Hydrodynamic binARY star orBital Decay In Spiralls
c-----------------------------------------------------------------------
c
c    The code was developed by
c
c        B. A. Fryxell,   University of Chicago,  Chicago,   USA
c        E. Mueller,      MPI f. Astrophysik,     Garching,  FRG
c        C. N. Arnold,    ETA Systems Inc.,       St.Paul,   USA
c        M. Ruffert       MPI f. Astrophysik      Garching,  FRG
c
c-----------------------------------------------------------------------
c
c     CHARYBDIS contains the following features:
c
c     1) Hydrodynamics; this part originated from a direct eulerian
c        PPM-code as described in Colella and Woodward
c        (JCP, 54 (1984), 174) and implemented in PROMETHEUS.
c
c     2) An equation of state containing gas and radiation pressure;
c        The riemann solver is according to Colella and Glaz
c        ( JCP, 59 (1985), 264 ).
c
c     3) Self-gravity; a Fourier-transform algorithm is used as
c        described in Miller ...
c
c     4) Two point-masses, which change the gravitational potential
c        and move according to the leapfrog-algorithm.
c
c     5) The subroutine video that produces film compressed data;
c
c     6) A fixed nested grid structure inspired by
c        Berger and Colella (JCP, 82 (1989), 64).
c
c=======================================================================


      include 'qparam.cmn'
      include 'ppdisk.cmn'
      include 'squants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'
      include 'rays.cmn'
      include 'aquants.cmn'

cIBM      integer*4  mclock
      real*8 eigrb, ekgrb, tmpgrb
      real*4 etim(2), etime
      logical stoprun

      real timeleft,tlim
      cc     = 2.99792458D10
      g      = 6.6726D-8
      pi     = 4.0D0*atan(one)
      amu    = 1.6605D-24
      forthd = 4.D0 / 3.D0
      solmas = 1.989D33
      gamma  = 1.4D0
      gascon = 8.314510D7
      ! Mean molar mass (Diatomic hydrogen in this case)
      gasmol  =  2.01588D0
      ! Specific gas constant
      specr  = gascon / gasmol
      radcon = 7.56591D-15
      au     = 1.4959787D13
      !Stefan-Boltzmann constant
      stefan = 5.670374D-5
      !ergirr = specr * tirr / (gamma - one)
c      tlim = 84600.
c      timelimit for job (queue !)

      call rstzero
      call input
      call pgen (g,soparg**2)
      call readeos

      nxh = nx / 2
      nyh = ny / 2
      nzh = nz / 2

      dtc = 1.D10*dtmax
      dtg = 1.D10*dtmax

      if(irstrt .eq. lzero) then
c        construct new model

         time   = zero
         nstep  = lzero
         nbegin = lzero

         if(suffix .ne. 'aa') then
            write(*,*) 'check input file'
            stop' suffix '
         endif

         write(*,*)  ' **** constructing new initial model **** '
         call inigrid
         call inippm

         !print *, 'inippm done'
c         do i = 1, nx
c           do j = 1, ny
c             do k = 1, nz
c               ekgrb = half * (velx(i,j,k,1)**2 + vely(i,j,k,1)**2
c     &                         + velz(i,j,k,1)**2)
c               eigrb = energy(i,j,k,1) - ekgrb
c               tmpgrb = (gamma-one)*eigrb / specr
c               if(tmpgrb.gt.1E3) then
c                 print *, i,j,k, tmpgrb, ekgrb, eigrb
c                 stop'high temp grb'
c               endif
c             enddo
c           enddo
c         enddo

         if ( ifilm .ne. lzero ) call inivid
         call opeplt
         call restrt(lzero)
         call datout

      else
c        read restart file

         call opeplt
         call restrt(lone)
c         call justrestrt(lone)
c         stop 'after justrestrt'

         if ( ifilm .ne. lzero ) call inivid

         nbegin = nstep
         nend   = nstep + nend
      endif

      nout1 = lzero
      tout1 = zero
      nrst  = lzero
      trst  = zero

      ntst = max( nx * 4  / 100 , lone )
      vtst = 0.1D0

      do igrd = 1, ngrd
         call fillold(igrd)
      enddo

c-------  start of cycle loop

      do 100  ncycl = nbegin, nend-2, 2

c      if ( blackhole ) goto 500

      if     ( ncycl .eq. nbegin ) then
cc         call second(tleft)
cIBM         tleft = dble(mclock())
         tleft = etime(etim)
      elseif ( ncycl .eq. nbegin+2 ) then
cc         call second(timebuf)
cc         timebuf = 1.2D0*( timebuf - tleft)
         timebuf = dble(etime(etim)) - tleft
         write(*,'(a,1F8.1,a)') 'two timesteps took ', timebuf,' sec'
         write(*,'(a,1F8.1,a)') 'included video with ',timebuf2,' sec'
c         stop
cc      else
cc         call tremain(tleft)
cc         if ( tleft .lt. timebuf ) then
cc            write(*,*) '***** timelimit reached !! ******'
cc            write(*,*)
cc            goto 500
cc         endif
      endif

      if ( (nrst.ge.nrstrt) .or. (trst.ge.trstrt) )  call restrt(lzero)

c     give all grids the velocity of WD
c     call grdvel

c --- do one timestep of coarsest grid, including all finer ones
c --- the recursive scheduling is contained in shedul and ifinup

      if ( (ifilm.ne.lzero) .and. (mod(ncycl,ifilm ).eq.lzero) .and.
     &  ((nstep.lt.100).or.(time-tframe.gt.1.4*time/dble(nstep+1))))then

c         tleft2 = dble(mclock())
c         timebuf2 = 1.0D-2*(dble(mclock()) - tleft2)
         tleft2 = dble(etime(etim))
           call video(0)          !  for mono view
cs         call video(-1)    !  for stereo view
cs         call video(+1)    !  for stereo view
c         stop 'grb.f, after video'
         timebuf2 = dble(etime(etim)) - tleft2
c         write(*,'(a,1F8.1,a)') 'video took ',timebuf2,' sec'
c         call cutout(lone)
      endif

c      alo = log10(dble(nstep/2+1)) / alog10(two)
c      if ( abs(alo-dnint(alo)) .lt. 1.D-5 ) call freeze
c      call freeze

      do 300 iot = 1, 2

c      if ( relaxon .eq. lone ) call freeze

      do 300 node = 1, nodenr

         levl = ishedul(node)

         if ( (node.eq.lone) .and. (iot.eq.lone) .and.
     &        (mod(ncycl,nx*2).eq.lzero) ) then
c     &     (mod(ncycl,100      ).eq.lzero) ) then
c     &     (mod(ncycl,lone).eq.lzero) ) then
cc         if ( ncycl.ne.nbegin )  then
            call plotieee
cc         stop
         endif
c     if(node.eq.20) stop

c        check if grids have to be moved
cr       call regrid

c        move all matter on all grids successively

         if ( ifinup(node) )  then
            call evolve(levl,ltwo)
c            use finer grid values for coarser grids where possible
            call fineup(levl,lone)
         elseif ( (node.eq.lone) .and. (iot.eq.ltwo) ) then
            call evolve(levl,ltwo)
         else
            call evolve(levl,lone)
         endif
         if (node .eq .nodenr)  call pushbh


 300  continue

c --- 2 timesteps = 1 cycle
      itopgr = indxgr(lone,lone)
      time   = time   + two*dtgr(itopgr)
      nstep  = ncycl  + ltwo
      nout1  = nout1  + ltwo
      tout1  = tout1  + one*dtgr(itopgr)
      nrst   = nrst   + ltwo
      trst   = trst   + one*dtgr(itopgr)

      if (nout1 .ge. nout    .or.  tout1 .ge. tout)    call datout
      call tstep(lone)

      if ( time .gt. tmax ) then
         write(*,*) '+++++ time gt tmax +++++'
         goto 500
      else
        print *, ''
        print *, 'Timestep:', dtgr(itopgr)
        print *, 'Time:', time
        print *, 'maxvelx,y,z:', maxval(velx)/cc, maxval(vely)/cc,
     &                            maxval(velz)/cc
        ! To find max temperature in the disk
        print *, 'maxtemp:', maxval(temper)
        print *, ''
        do i = 1,nx
        do j = 1,ny
        do k = 1,nz
          if(temper(i,j,k,1).gt.10000) print *,'TEMP 1', i, j, k
        enddo
        enddo
        enddo

        do i = 1,nx
        do j = 1,ny
        do k = 1,nz
          if(temper(i,j,k,2).gt.10000) print *,'TEMP 2', i, j, k
        enddo
        enddo
        enddo

        do i = 1,nx
        do j = 1,ny
        do k = 1,nz
          if(temper(i,j,k,3).gt.10000) print *,'TEMP 3', i, j, k
        enddo
        enddo
        enddo

        do i = 1,nx
        do j = 1,ny
        do k = 1,nz
          if(temper(i,j,k,4).gt.10000) print *,'TEMP 4', i, j, k
        enddo
        enddo
        enddo
      endif
      if ( stoprun() ) goto 500

 100  continue
c------- end of cycle loop

 500  continue

      do igrd = 1, ngrd
         call calm(igrd)
      enddo

      if ( ifilm .ne. lzero ) call cloplt
      call restrt(lzero)

      if ( time .gt. tmax ) then
         call exit(1)
      endif

      print *, 'time', time

      stop' charyb'
      end
c-----------------------------------------------------------------------
      subroutine evolve(levl,icycl)
c --- evolve all matter on grids igrd at level levl
c --- 2 timesteps = 1 cycle

      include 'qparam.cmn'
      include 'grnest.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'

      do 300 inum = 1, indxgr(lzero,levl)


         igrd = indxgr(inum,levl)
c         write(*,*) 'evolve ',igrd,inum,levl,icycl
         call calm(igrd)
c         call velfudge(igrd)
         call fillold(igrd)
         call ppmstep(icycl,igrd)
         call accel   (igrd)         ! includes graviwaves backreaction

         ! Compute internal energy loss due to cooling
         if(sbyn .eq. lone) then
           call stamcool(igrd)
         elseif(bcool .gt. zero) then
           call betacool(igrd)
         else
           if(igrd .eq. lone) print *, 'Cooling off'
         endif

!         print *, 'maxvals evolve:', maxval(temper), maxval(velx),
!     &                               maxval(vely), maxval(velz)

!         call neutrino(igrd)         ! neutrino loss source terms      !Commented out - DRW
!         call entropy (igrd,icycl)   ! entropy shock source terms      !!!
c         stop'evolve'

  300 continue


      return
      end
c-----------------------------------------------------------------------
      subroutine ppmstep(icycl,igrd)
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'
      include 'ppdisk.cmn' ! eiprior needed for cooling

c     nxyzsw contains information on sweep direction; only
c     required if artificial viscosity and/or diffusion is applied
c     in subroutine hydro and/or in subroutine flaten
c        nxyzsw = 1  ===>  x-sweep
c        nxyzsw = 2  ===>  y-sweep
c        nxyzsw = 3  ===>  z-sweep

      if(sbyn.eq.one) then
        ! Save the specific internal energy before the hydrodynamic evolution.
        ! This is needed for the Stamatellos scheme.
        do i = 1,nx
        do j = 1,ny
        do k = 1,nz
          ekprior = half * ( velx(i,j,k,igrd)**2 +
     &                       vely(i,j,k,igrd)**2 +
     &                       velz(i,j,k,igrd)**2 )
          eiprior(i,j,k,igrd) = energy(i,j,k,igrd) - ekprior
        enddo
        enddo
        enddo
      endif

      do i = 1,nx
      do j = 1, ny
      do k = 1, nz
        if(velz(i,j,k,igrd).gt.half*cc) then
          print *, 'velz big start of ppmstep', igrd
          stop
        endif
      enddo
      enddo
      enddo

c-----
      if (icycl .eq. lone) then
c-----

c       call gasumas('ppms x',igrd)

      !print *, 'x-sweep 1', maxval(temper), maxval(velz)
c-------  x -  s w e e p

      nxyzsw = lone

      call boundry(lone,igrd,nxyzsw,icycl)

C$OMP PARALLEL DO PRIVATE(j,k), SHARED(ny,nz,igrd,nxyzsw)

      do j = 1, ny
      do k = 1, nz
         call getrwx (j, k, igrd)
         call hydrow (j, k, igrd, nxyzsw)
         call putrwx (j, k, igrd)
         call flupmem(j, k, igrd, nxyzsw,ltwo,lthree)
      enddo
      enddo

C$OMP END PARALLEL DO

c       call gasumas('ppms y',igrd)
      !print *, 'y-sweep 1', maxval(temper), maxval(velz)
c-------  y -  s w e e p  (only for 2-d and 3-d problems)

      nxyzsw = ltwo
      call boundry(lone,igrd,nxyzsw,icycl)

C$OMP PARALLEL DO PRIVATE(i,k), SHARED(nx,nz,igrd,nxyzsw)
      do i = 1, nx
      do k = 1, nz
         call getrwy (i, k, igrd)
         call hydrow (i, k, igrd, nxyzsw)
         call putrwy (i, k, igrd)
         call flupmem(i, k, igrd, nxyzsw,lone,lthree)
      enddo
      enddo

C$OMP END PARALLEL DO

c       call gasumas('ppms z',igrd)
      !print *, 'z-sweep 1', maxval(temper), maxval(velz)
c-------  z -  s w e e p  (only for 3-d problems)

      nxyzsw = lthree
      call boundry(lone,igrd,nxyzsw,icycl)

C$OMP PARALLEL DO PRIVATE(j,i), SHARED(ny,nx,igrd,nxyzsw)
      do j = 1, ny
      do i = 1, nx
         call getrwz (i, j, igrd)
         call hydrow (i, j, igrd, nxyzsw)
         call putrwz (i, j, igrd)
         call flupmem(i, j, igrd, nxyzsw,lone,ltwo)
      enddo
      enddo

      do i = 1,nx
      do j = 1, ny
      do k = 1, nz
        if(velz(i,j,k,igrd).gt.half*cc) then
          print *, 'velz big z-sweep ppmstep', igrd
          stop
        endif
      enddo
      enddo
      enddo

C$OMP END PARALLEL DO

c-----
      elseif (icycl .eq. ltwo) then
c-----
c --- begin of second half of cycle with sweep directions reversed

c       call gasumas('ppms z',igrd)
      !print *, 'z-sweep 2', maxval(temper), maxval(velz)
c-------  z -  s w e e p  (only for 3-d problems)

      nxyzsw = lthree

      call boundry(lone,igrd,nxyzsw,icycl)

C$OMP PARALLEL DO PRIVATE(j,i), SHARED(ny,nx,igrd,nxyzsw)
      do j = 1, ny
      do i = 1, nx
         call getrwz (i, j, igrd)
         call hydrow (i, j, igrd, nxyzsw)
         call putrwz (i, j, igrd)
         call flupmem(i, j, igrd, nxyzsw,lone,ltwo)
      enddo
      enddo
C$OMP END PARALLEL DO

c       call gasumas('ppms y',igrd)
      !print *, 'y-sweep 2', maxval(temper), maxval(velz)
c-------  y -  s w e e p  (only for 2-d and 3-d problems)

      nxyzsw = ltwo

      call boundry(lone,igrd,nxyzsw,icycl)

C$OMP PARALLEL DO PRIVATE(i,k), SHARED(nx,nz,igrd,nxyzsw)
      do i = 1, nx
      do k = 1, nz
         call getrwy (i, k, igrd)
         call hydrow (i, k, igrd, nxyzsw)
         call putrwy (i, k, igrd)
         call flupmem(i, k, igrd, nxyzsw,lone,lthree)
      enddo
      enddo

C$OMP END PARALLEL DO

c       call gasumas('ppms x',igrd)
      !print *, 'x-sweep 2', maxval(temper), maxval(velz)
c-------  x -  s w e e p

      nxyzsw = lone

      call boundry(lone,igrd,nxyzsw,icycl)

C$OMP PARALLEL DO PRIVATE(j,k), SHARED(ny,nz,igrd,nxyzsw)
      do j = 1, ny
      do k = 1, nz
         call getrwx (j, k, igrd)
         call hydrow (j, k, igrd, nxyzsw)
         call putrwx (j, k, igrd)
         call flupmem(j, k, igrd, nxyzsw,ltwo,lthree)
      enddo
      enddo
C$OMP END PARALLEL DO

c-----
      endif
c-----

c       call gasumas('ppms e',igrd)

      return
      end
c-----------------------------------------------------------------------
      subroutine gravty(igrd)
c     find gravitational potential for grid igrd; output is gpot(core);
c     bakpot (..., igrd) must be given.
c     produces bakpot for next finer grids.

      include 'qparam.cmn'
      include 'grnest.cmn'

c --- calculate background potentials for
c     all grids one level finer and covered by igrd

      !print *, 'grb, gravty'

      do ifgr = 1, ngrd
         if ( idxcgr(ifgr) .eq. igrd ) then
            call detoro(igrd,lone)
            call rozero(ifgr)
            call potent
            call adbakf(igrd,ifgr,lone)
            call tripo(igrd,ifgr,lone)
         endif
      enddo

c --- find full potential for this grid igrd

      call detoro(igrd,lone)
      call potent
      call adbakcf(igrd,lone)

c      if ( relaxon .eq. lone ) call pseudopot(igrd)

      return
      end
c-----------------------------------------------------------------------
      subroutine dtgravty(igrd)
c     find time derivative of pot

      include 'qparam.cmn'
      include 'grnest.cmn'

c --- calculate background potentials for
c     all grids one level finer and covered by igrd

      !print *, 'grb, dtgravty'

      do ifgr = 1, ngrd
         if ( idxcgr(ifgr) .eq. igrd ) then
            call detoro(igrd,ltwo)
            call rozero(ifgr)
            call potent
            call adbakf(igrd,ifgr,ltwo)
            call tripo(igrd,ifgr,ltwo)
         endif
      enddo

c --- find full potential for this grid igrd

      call detoro(igrd,ltwo)
      call potent
      call adbakcf(igrd,ltwo)

      return
      end
c-----------------------------------------------------------------------
      subroutine chigravty(igrd)
c     find potential including energy corrections for deep potential

      include 'qparam.cmn'
      include 'grnest.cmn'

c --- calculate background potentials for
c     all grids one level finer and covered by igrd

      !print *, 'grb, chigravty'

      do ifgr = 1, ngrd
         if ( idxcgr(ifgr) .eq. igrd ) then
            call detoro(igrd,lthree)
            call rozero(ifgr)
            call potent
            call adbakf(igrd,ifgr,lthree)
            call tripo(igrd,ifgr,lthree)
         endif
      enddo

c --- find full potential for this grid igrd

      call detoro(igrd,lthree)
      call potent
      call adbakcf(igrd,lthree)

      return
      end
c-----------------------------------------------------------------------
      subroutine rgpot(igrd)
c     find quantity pot R in equations

      include 'qparam.cmn'
      include 'grnest.cmn'

c --- calculate background potentials for
c     all grids one level finer and covered by igrd

      !print *, 'grb, rgpot'

      do ifgr = 1, ngrd
         if ( idxcgr(ifgr) .eq. igrd ) then
            call der2ro(igrd)
            call rozero(ifgr)
            call potent
            call adbakr(igrd,ifgr)
            call tripo(igrd,ifgr,ltwo+ltwo)
         endif
      enddo

c --- find full potential for this grid igrd

      call der2ro(igrd)
      call potent
      call adbakcr(igrd)

      return
      end
c-----------------------------------------------------------------------
      subroutine calm(igrd)

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'eos.cmn'

      real*8  rh(qx), tm(qx), ei(qx), pr(qx), gc(qx), ge(qx), yea(qx),
     &        dtm(qx), ek(qx), tgr, rhl(qx), velfac, ent(qx), ye(qx)
      integer*4 inum, inx(qx), igrd

c -- calm down matter 'outside'

c      denmax = zero

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+            PRIVATE(j,k,i,ek,velfac,inum,inx,rh,ye,tm,ei,pr,gc,ge),
C$OMP+             SHARED(densty,rhoin,velx,vely,velz,energy,small,
C$OMP+                    chein,tin,chem,igrd)

      !print *, 'grb, calm'
      do  k = 1, qz
      do  j = 1, qy

         do  i = 1, qx
c            denmax = max ( denmax, densty(i,j,k,igrd) )
         if ( densty(i,j,k,igrd) .lt. 4.D0*rhoin ) then
           ek(i) = half * (
     &      velx(i,j,k,igrd)**2+vely(i,j,k,igrd)**2+velz(i,j,k,igrd)**2)
           energy(i,j,k,igrd) = energy(i,j,k,igrd) - ek(i)
           !velfac = sqrt( 0.5D14 / ek(i) )
           velfac = 1.0                                                 ! DRW - Dummy
           velx(i,j,k,igrd) = velx(i,j,k,igrd) * velfac
           vely(i,j,k,igrd) = vely(i,j,k,igrd) * velfac
           velz(i,j,k,igrd) = velz(i,j,k,igrd) * velfac
           ek(i) = half * (
     &      velx(i,j,k,igrd)**2+vely(i,j,k,igrd)**2+velz(i,j,k,igrd)**2)
           energy(i,j,k,igrd) = energy(i,j,k,igrd) + (one+small)*ek(i)
         endif
         enddo

         inum = lzero
         do  i = 1, qx
            if ( densty(i,j,k,igrd) .lt. ten*rhoin ) then
               inum = inum + lone
               inx(inum) = i
               rh(inum)  = densty(i,j,k,igrd)
               ye(inum)  = chein(lone) / rhoin
               tm(inum)  = tin
            endif
         enddo

         if ( inum .ne. lzero ) then

            call eos (rh, tm, ei, pr, ye, gc, ge, inum, lone)

            do  i = 1, inum
               ek(i) = half * ( velx(inx(i),j,k,igrd)**2 +
     &             vely(inx(i),j,k,igrd)**2 + velz(inx(i),j,k,igrd)**2 )
               energy(inx(i),j,k,igrd) = ei(i)/rh(i) + (one+small)*ek(i)
               chem(inx(i),j,k,igrd,lone) = ye(i) * rh(i)
            enddo

         endif

      enddo
      enddo

C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine freeze
c     slow and cool matter for relaxation test runs

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'eos.cmn'

      real*8  rh(qx), tm(qx), ei(qx), pr(qx), gc(qx), ge(qx), yea(qx),
     &        dtm(qx), ek(qx), tgr, rhl(qx), velfac, ent(qx), ye(qx)
      integer*4 inum, inx(qx), igrd

      !print *, 'grb, freeze'

      do igrd = 1, ngrd

c -- slow matter

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+            PRIVATE(j,k,i,rh,ye,tm,ei,ek,pr,gc,ge,ent),
C$OMP+             SHARED(densty,velx,vely,velz,energy,temper,tmpent,
C$OMP+                    chem,chein,rhoin,tin,igrd,small)
      do  j = 1, qy
      do  k = 1, qz

         do  i = 1, qx
         if ( densty(i,j,k,igrd) .ge. 4.D0*rhoin ) then
           stop'Big old density, grb, 1'
           ek(i) = half * (
     &      velx(i,j,k,igrd)**2+vely(i,j,k,igrd)**2+velz(i,j,k,igrd)**2)
           energy(i,j,k,igrd) = energy(i,j,k,igrd) - ek(i)
           velx(i,j,k,igrd) = velx(i,j,k,igrd) / two
           vely(i,j,k,igrd) = vely(i,j,k,igrd) / two
           velz(i,j,k,igrd) = velz(i,j,k,igrd) / two
           energy(i,j,k,igrd) = energy(i,j,k,igrd) + half * (
     &      velx(i,j,k,igrd)**2+vely(i,j,k,igrd)**2+velz(i,j,k,igrd)**2)
         endif
         enddo


c --- freeze matter
         if (.false.) then

         do  i = 1, qx
               rh(i)  = densty(i,j,k,igrd)
               ye(i)  = chein(lone) / rhoin
               tm(i)  = max( tin, tmpent(i,j,k,igrd)/two )
         enddo

         call eos (rh, tm, ei, pr, ye, gc, ge, qx, lone)
         call eosent( rh, tm, ye, ent, qx, lone ) !  tempr -> entropy

         do  i = 1, qx
         if ( densty(i,j,k,igrd) .ge. 4.D0*rhoin ) then
            stop'Big old density, grb, two'
            energy(i,j,k,igrd) = ei(i)/rh(i) + (one+small)*half * (
     &      velx(i,j,k,igrd)**2+vely(i,j,k,igrd)**2+velz(i,j,k,igrd)**2)
            chem(i,j,k,igrd,lone) = ye(i) * rh(i)
            chem(i,j,k,igrd,ltwo) = ent(i) * rh(i)
            temper(i,j,k,igrd) = tm(i)
            tmpent(i,j,k,igrd) = tm(i)
         endif
         enddo

      endif

      enddo
      enddo
C$OMP END PARALLEL DO

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pseudopot(igrd)
c     add centrifugal force pseudo potential for rotating frame of reference

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      integer*4 igrd

      !print *, 'grb, pseudopot'

      d1x2 = delx(indxgr(1,1))**2

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+            PRIVATE(j,k,i,tmy1,rrr),
C$OMP+             SHARED(ny,nz,nx,igrd,ycm,xcm,rhoin,omega,d1x2,
C$OMP+                    densty,gpot)
      do j=1,ny
         tmy1  = (topos(dble(j),ltwo,igrd) - ycm)**2
      do k=1,nz
      do i=1,nx
         if ( densty(i,j,k,igrd) .ge. 4.D0*rhoin ) then
            stop'Big old density, grb, three'
c        matter inside NS as taken in subroutine calm
            rrr  = (tmy1 + (topos(dble(i),lone,igrd)-xcm)**2) * d1x2
            gpot(i,j,k,igrd) = gpot(i,j,k,igrd)! - half* omega**2 * rrr !DRW - Removing all omegas.
         endif
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine wakeup
c     give BH and NS orbital velocities

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      real*8 erg(qx)

      !print *, 'grb, wakeup'

      itopgr=  indxgr(lone,lone)
      d1x = delx(itopgr)

      pmass2 = gasmas

      vrns = zero
      dis1 = bidist * pmass2 / (pmass1+pmass2)
      dis2 = bidist * pmass1 / (pmass1+pmass2)

      if ( nspin .eq. 4 ) then   !colli
         velx1 = - vrns * dis1 / bidist
     &        - sqrt( two * g * pmass2**2 / bidist / (pmass1+pmass2) )
         vely1 = zero
         velx2 = + vrns * dis2 / bidist
     &        + sqrt( two * g * pmass1**2 / bidist / (pmass1+pmass2) )
         vely2 = zero
      else
         velx1 = +vrns  * dis1 / bidist
         vely1 = +omega * dis1
         velx2 = -vrns  * dis2 / bidist
         vely2 = -omega * dis2
      endif

      vnsx = velx2
      vnsy = vely2

      pmass2 = zero    ! not needed any more

      do 301 igrd = 1, ngrd

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,tmy1,tmx1,erg),
C$OMP+            SHARED(nx,ny,nz,igrd,cenx,ceny,rhoin,d1x,omega,
C$OMP+                   velx,vely,energy,densty,bhdens)

         do j=1,ny
            tmy1  = topos(dble(j),ltwo,igrd) - cenx
         do k=1,nz
         do i=1,nx
c            if ( velx(i,j,k,igrd) .ne. zero ) then ! inside NS
           if ( (bhdens(i,j,k,igrd) .eq. zero) .and.
     &          (densty(i,j,k,igrd) .ge. 4.D0*rhoin)    )  then
c           matter inside NS as taken in subroutine calm

               tmx1  = topos(dble(i),lone,igrd) - cenx

               erg(i)= half*( velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2)

               velx(i,j,k,igrd) = velx(i,j,k,igrd)! - omega*tmy1*d1x    !DRW - Commented out omega part.

               vely(i,j,k,igrd) = vely(i,j,k,igrd)! + omega*tmx1*d1x    !!

               energy(i,j,k,igrd) = energy(i,j,k,igrd)  - erg(i)
     &               + half*( velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2)

            endif

         enddo
         enddo
         enddo

C$OMP END PARALLEL DO

  301 continue


      vlox1 = velx1
      vloy1 = vely1
      vlox2 = velx2
      vloy2 = vely2


      return
      end
c-----------------------------------------------------------------------
      subroutine nospike( rhw, tmw, tme, eiw, prw, yew, gcw, gew, np )
c --- reduce temperature spikes
C     PARALLEL OK

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'eos.cmn'

      real*8  rhw(np), tmw(np), tme(np), eiw(np), prw(np), yew(np),
     &        gcw(np), gew(np)
      real*8  rh(q), tm(q), ei(q), pr(q), gc(q), ge(q), ye(q)
      real*8  dtm(q), tgr, tmi, tgr1, tgr0, tave(q)
      integer*4 inum, inx(q), np, i

      !print *, 'grb, nospike'

      if ( np .gt. q ) stop ' nospike'

      tgr = one    ! 1MeV per zone as temperature gradient
      tmi = smltme

      inum = lzero

      do  i = 1, np-1
         dtm(i) =    tmw(i+1) - tmw(i)
         tave(i) = ( tmw(i+1) + tmw(i) )
      enddo


      do  i = 2, np-1
         if ( (tmw(i).lt. 9.D0) .or. (rhw(i) .lt. 1.5D14) ) then
            if ( (dtm(i) .lt. -tgr) .and. (dtm(i-1) .gt.  tgr) ) then
               inum = inum + lone              !   spike up
               inx(inum)= i
               rh(inum) = rhw(i)
               ye(inum) = yew(i)
               tm(inum) = max( half*(tmw(i-1)+tmw(i+1)), tmi,
     &                            min( tme(i), tmw(i) )        )
            endif
            tgr1 = tgr / max( min(tgr,tave(i)), 0.1D0 )
            tgr0 = tgr / max( min(tgr,tave(i-1)), 0.1D0 )
            if ( (dtm(i) .gt.  tgr1) .and. (dtm(i-1) .lt. -tgr0) ) then
               inum = inum + lone              !   spike down
               inx(inum)= i
               rh(inum) = rhw(i)
               ye(inum) = yew(i)
               tm(inum) = max( half*(tmw(i-1)+tmw(i+1)), tmi )
            endif
         endif
      enddo


      if ( inum .ne. lzero ) then

         call eos (rh, tm, ei, pr, ye, gc, ge, inum, lone)

         do  i = 1, inum
            rhw(inx(i)) = rh(i)
            tmw(inx(i)) = tm(i)
            eiw(inx(i)) = ei(i)
            prw(inx(i)) = pr(i)
            yew(inx(i)) = ye(i)
            gcw(inx(i)) = gc(i)
            gew(inx(i)) = ge(i)
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nospike3d(igrd)
c     call nospike in all 3 dimensions

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'eos.cmn'

      real*8  rhw(q), tmw(q), tme(q),eiw(q),prw(q),yew(q),gcw(q),gew(q)
      integer*4 igrd

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i),
C$OMP+           SHARED(ny,nx,nz,igrd,energy,densty,velx,vely,velz,chem)

      !print *, 'grb, nospike3d'

      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         energy(i,j,k,igrd) =
     &         densty(i,j,k,igrd) * ( energy(i,j,k,igrd) - half * (
     &     velx(i,j,k,igrd)**2+vely(i,j,k,igrd)**2+velz(i,j,k,igrd)**2))
         chem(i,j,k,igrd,1) = chem(i,j,k,igrd,lone) / densty(i,j,k,igrd)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c --- x-direction

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+     PRIVATE(j,k,i,rhw,tmw,tme,prw,eiw,yew,gcw,gew),
C$OMP+     SHARED(densty,temper,tmpent,press,energy,chem,nx,ny,nz,igrd)
      do j = 1, ny
      do k = 1, nz

      do i = 1, nx
         rhw(i) = densty(i,j,k,igrd)
         tmw(i) = temper(i,j,k,igrd)
         tme(i) = tmpent(i,j,k,igrd)
         prw(i) = press (i,j,k,igrd)
         eiw(i) = energy(i,j,k,igrd)
         yew(i) = chem(i,j,k,igrd,lone)
      enddo

      call nospike( rhw, tmw, tme, eiw, prw, yew, gcw, gew, nx )

      do i = 1, nx
         densty(i,j,k,igrd) = rhw(i)
         temper(i,j,k,igrd) = tmw(i)
         press (i,j,k,igrd) = prw(i)
         energy(i,j,k,igrd) = eiw(i)
         chem(i,j,k,igrd,lone) = yew(i)
      enddo

      enddo
      enddo
C$OMP END PARALLEL DO

c --- y-direction

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+     PRIVATE(i,k,j,rhw,tmw,tme,prw,eiw,yew,gcw,gew),
C$OMP+     SHARED(densty,temper,tmpent,press,energy,chem,nx,ny,nz,igrd)
      do i = 1, nx
      do k = 1, nz

      do j = 1, ny
         rhw(j) = densty(i,j,k,igrd)
         tmw(j) = temper(i,j,k,igrd)
         tme(j) = tmpent(i,j,k,igrd)
         prw(j) = press (i,j,k,igrd)
         eiw(j) = energy(i,j,k,igrd)
         yew(j) = chem(i,j,k,igrd,lone)
      enddo

      call nospike( rhw, tmw, tme, eiw, prw, yew, gcw, gew, ny )

      do j = 1, ny
         densty(i,j,k,igrd) = rhw(j)
         temper(i,j,k,igrd) = tmw(j)
         press (i,j,k,igrd) = prw(j)
         energy(i,j,k,igrd) = eiw(j)
         chem(i,j,k,igrd,lone) = yew(j)
      enddo

      enddo
      enddo
C$OMP END PARALLEL DO

c --- z-direction

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+     PRIVATE(j,i,k,rhw,tmw,tme,prw,eiw,yew,gcw,gew),
C$OMP+     SHARED(densty,temper,tmpent,press,energy,chem,nx,ny,nz,igrd)
      do j = 1, ny
      do i = 1, nx

      do k = 1, nz
         rhw(k) = densty(i,j,k,igrd)
         tmw(k) = temper(i,j,k,igrd)
         tme(k) = tmpent(i,j,k,igrd)
         prw(k) = press (i,j,k,igrd)
         eiw(k) = energy(i,j,k,igrd)
         yew(k) = chem(i,j,k,igrd,lone)
      enddo

      call nospike( rhw, tmw, tme, eiw, prw, yew, gcw, gew, nz )

      do k = 1, nz
         densty(i,j,k,igrd) = rhw(k)
         temper(i,j,k,igrd) = tmw(k)
         press (i,j,k,igrd) = prw(k)
         energy(i,j,k,igrd) = eiw(k)
         chem(i,j,k,igrd,lone) = yew(k)
      enddo

      enddo
      enddo
C$OMP END PARALLEL DO

c ---

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,i),
C$OMP+        SHARED(nx,ny,nz,igrd,energy,densty,velx,vely,velz,chem)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         energy(i,j,k,igrd) =
     &     energy(i,j,k,igrd) / densty(i,j,k,igrd)  + half * (
     &     velx(i,j,k,igrd)**2+vely(i,j,k,igrd)**2+velz(i,j,k,igrd)**2)
         chem(i,j,k,igrd,lone)=chem(i,j,k,igrd,lone)*densty(i,j,k,igrd)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine lightning(u,ut,utt,e,ek,nn)
c --- slow down matter below the speed of light
C     PARALLEL OK

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'

      real*8  u(nn), ut(nn), utt(nn), e(nn), ek(nn),
     &        velfac, cfudge, ekinmax
      integer*4 nn

      !print *, 'grb, lightning'

      cfudge  = two
      ekinmax = half * (cfudge*cc)**2

c      do  i = 1, nn
c         if ( ek(i) .gt. ekinmax )
c     &      write(*,*) 'lightning ',i,sqrt(two*ek(i))/cc
c      enddo

      do  i = 1, nn
         if ( ek(i) .gt. ekinmax ) then
            e(i)   = e(i) - ek(i)
            velfac = sqrt( ekinmax / ek(i) )
            u  (i) = u  (i) * velfac
            ut (i) = ut (i) * velfac
            utt(i) = utt(i) * velfac
            ek (i) = half * ( u(i)**2 + ut(i)**2 + utt(i)**2 )
            e  (i) = e(i) + ek(i)
         endif
      enddo

      return
      end
c------------------------------------------------------------------------
      subroutine sucks(igrd,avrdo)
c --  sum mass and momentum within 1 Rs into point mass pmass1
c     includes:  1) put density distrib of bh-delta-fkt on grid igrd
c                2) unify this grid with other grids for consistent bh-mass
c                3) fill interior of bh with thermodyn. quantities
c                4) fill density of gas to position of bh for gravi waves only
c     avrdo=0: only for init

      include 'qparam.cmn'
      include 'ppdisk.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      integer*4 igrd, avrdo
      real*8    rs, psumas, angmsu, baksum, disx, disy, disq, disqz,
     &          bhmas(ngrd), psmomx, psmomy, rhom, dmass1, disqy,
     &          blhomas(ngrd),angacc, angcc1

      rs   = two * g * pmass1   / cc**2
      itopgr=  indxgr(lone,lone)
      dxt  = delx(itopgr)

      if ( akerr .eq. zero ) then
         rota = zero
      elseif ( abs(akerr) .ge. one ) then
         rota = abs(angm1) * cc / (g * pmass1**2)
      else
         rota = abs(akerr)
      endif

      !print *, 'O N E'

      if ( (rota.lt.zero) .or. (rota.gt.one) ) then
         write(*,*) 'rota:',rota,'  akerr:',akerr
         stop
      endif

      !print *, 'T W O'

      thrd = one/3D0
      z1 = one  +  (one-rota**2)**thrd
     &                    *( (one+rota)**thrd + (one-rota)**thrd )
      z2 = sqrt( 3D0*rota**2 + z1**2 )
      rin = 3D0 + z2 - sign( sqrt( (3D0-z1)*(3D0+z1+two*z2) ), angm1 )
c      rin = 3D0 + z2 - sqrt( (3D0-z1)*(3D0+z1+two*z2) ) ! prograde
c      rin = 3D0 + z2 + sqrt( (3D0-z1)*(3D0+z1+two*z2) ) ! retrograde
      r1 = one + sqrt(one-rota**2)
      rbhkerr = ( rin + r1 ) / two   *rs/two

c DRW - I will put this in for now. ====================================

      !print *, 'r_in', r_in
      rbhkerr = r_in

c ======================================================================

c                        average      rin,r1 in units of rs/2

c      frsq = (two*rs / dxt)**2    ! hydro-size of BH is two*rs !
      frsq = (rbhkerr / dxt)**2    ! hydro-size of BH is rbhkerr
      dxi  = delx(igrd)
      d3x2 = two * dxi**3
      frhoin = 5.D0 * rhoin

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,densty,ro,igrd)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         ro(i,j,k) = densty(i,j,k,igrd)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c     zero area covered by finer grids
      ifgr = idxfgr(igrd)
      if( ifgr .ne. -lone ) call rozero(ifgr)

      psumas = zero
      angmsu = zero
      psmomx = zero
      psmomy = zero
      angacc = zero
      psueint = zero

c --- E.T. addition
c      itopgr=indxgr(lone,lone)
c      d1x   =delx(itopgr)
      rkf  =one / (dble(nfine)**(igrd-1))
      xopos=topos(zero,lone,igrd)
      yopos=topos(zero,ltwo,igrd)
      topx1=posx1
      topy1=posy1
      cmcax = (cmgax*gasmas + topx1*pmass1)/(gasmas+pmass1)
      cmcay = (cmgay*gasmas + topy1*pmass1)/(gasmas+pmass1)
      vcmcx = (vcmgx*gasmas + velx1*pmass1)/(gasmas+pmass1)
      vcmcy = (vcmgy*gasmas + vely1*pmass1)/(gasmas+pmass1)
c      write(*,*) 'sucks : ',cmcax,cmcay,vcmcx,vcmcy
c --- end addition
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,tpy,tpx,disy,disx,
C$OMP+                 disqy,disqx,disqz,dcmx,dcmy,disq,rhom),
C$OMP+     SHARED(nx,ny,nz,igrd,cmcay,cmcax,posy1,posx1,frsq,ro,frhoin,
C$OMP+            velx,vely,velz,vcmcy,vcmcx,energy,bhdens,velx1,vely1),
C$OMP+     REDUCTION(+:psumas,angacc,angmsu,psmomx,psmomy,psueint)

      do j = 1, ny
         tpy  = topos(dble(j),ltwo,igrd)
         dcmy = tpy - cmcay
         disy = tpy - posy1
         disqy = disy**2
      do k = 1, nz
         disqz = ( topos(dble(k),lthree,igrd) - half )**2 + disqy
      do i = 1, nx
         tpx  = topos(dble(i),lone,igrd)
         dcmx = tpx - cmcax
         disx = tpx - posx1
         disq  = disx**2  + disqz
         if ( disq .lt. frsq ) then
            rhom = max(ro(i,j,k)-frhoin,zero)
            psumas = psumas + rhom
            angacc = angacc + rhom *
     &           (    dcmx * (vely(i,j,k,igrd)-vcmcy)
     &              - dcmy * (velx(i,j,k,igrd)-vcmcx) )
            angmsu = angmsu + rhom *
c     &           (    disx * (vely(i,j,k,igrd)-vely1)      ! right?
c     &              - disy * (velx(i,j,k,igrd)-velx1)  )   ! right?
     &           (    disx * vely(i,j,k,igrd)             ! wrong?
     &              - disy * velx(i,j,k,igrd)  )          ! wrong?
            psmomx = psmomx + rhom * velx(i,j,k,igrd)
            psmomy = psmomy + rhom * vely(i,j,k,igrd)
            psueint = psueint + rhom * ( max( zero,
     &            energy(i,j,k,igrd) - half * (velx(i,j,k,igrd)**2 +
     &                   vely(i,j,k,igrd)**2 + velz(i,j,k,igrd)**2)  ))
c            bhdens(i,j,k,igrd) = (one - disq/frsq)**2
c try to steepen density profile here: **3 instead of **2
            bhdens(i,j,k,igrd)=max(exp(-4D0*(disq/frsq)**2)-0.05,1D-20)
         else
            bhdens(i,j,k,igrd) = zero
         endif
      enddo
      enddo
      enddo

C$OMP END PARALLEL DO

      dmass1= psumas*d3x2
      velx1 = ( velx1*pmass1 + psmomx*d3x2 ) / (pmass1+dmass1)
      vely1 = ( vely1*pmass1 + psmomy*d3x2 ) / (pmass1+dmass1)
      pmass1 = pmass1 + dmass1
      gasmas = gasmas - dmass1
      angm1  = angm1  + angmsu*d3x2 * dxt
      angcc1 = angcc1 + angacc*d3x2 * dxt
      eint1 = eint1 + psueint*d3x2
c      write(*,*) 'angm1',angm1
c      write(*,*) 'angcc1',angcc1

      if ( avrdo .eq. 1 ) then
c --- scan all levels, starting with finest
         do 300 levl = mlev, 2, -1
c --- take all grids on level levl
         do 300 inum = 1, indxgr(lzero,levl)
            ifgr = indxgr(inum,levl)
c --- take all grids one level coarser
         do 300 inuf = 1, indxgr(lzero,levl-lone)
            icgr = indxgr(inuf,levl-lone)

            if ( icgr .eq. idxcgr(ifgr) ) then
            call avrqty ( bhdens(1,1,1,ifgr), bhdens(1,1,1,icgr), ifgr,
     &                    densty(1,1,1,ifgr), densty(1,1,1,icgr), lzero)
            endif
  300    continue

         lgrd = itopgr
         bhmas(lgrd) = zero
         bhmas(igrd) = zero
         blhomas(igrd)= zero

         bhma = zero
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,bhdens,lgrd), REDUCTION(+:bhma)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            bhma = bhma + bhdens(i,j,k,lgrd)
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO
         bhmas(lgrd) = bhmas(lgrd) + bhma
         bhdfac = pmass1 / bhmas(lgrd) / (two*dxt**3)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,bhdens,velx,vely,velz,
C$OMP+                   energy,press,temper,chem,densty,velx1,vely1,
C$OMP+                   smalle,smallp,tin,bhdfac)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            if ( bhdens(i,j,k,igrd) .ne. zero ) then
               velx  (i,j,k,igrd) = velx1
               vely  (i,j,k,igrd) = vely1
               velz  (i,j,k,igrd) = zero
               energy(i,j,k,igrd) = smalle + half *(velx1**2 + vely1**2)
               press (i,j,k,igrd) = smallp
               temper(i,j,k,igrd) = tin
               chem  (i,j,k,igrd,lone) = 1.D-3*densty(i,j,k,igrd) ! Ye * rho
               chem  (i,j,k,igrd,ltwo) = 1.D-2*densty(i,j,k,igrd) ! ent* rho
               densty(i,j,k,igrd) = bhdens(i,j,k,igrd)*bhdfac
c               blhomas(igrd) = blhomas(igrd) + densty(i,j,k,igrd) ! test
            endif
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

      endif   ! of if avrdo

      return
      end
c------------------------------------------------------------------------
      subroutine kickbh(igrd,kif)
c --- update velocity of BH as average of changed values by gravi waves backr.
c --- reset density within bh horizon to small values.

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      integer*4 igrd, kif
      real*8    pvx, pvy, pma, dm, d3x2

      !print *, 'grb, kickbh'

      frhoin = 5.D0 * rhoin    ! factor 5. because of subroutine calm

      if ( kif .eq. lone ) then

         pma = zero
         pvx = zero
         pvy = zero

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,densty,ro)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            ro(i,j,k) = densty(i,j,k,igrd)
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

c     zero area covered by finer grids
         ifgr = idxfgr(igrd)
         if( ifgr .ne. -lone ) call rozero(ifgr)

         d3x2 = two * delx(igrd)**3

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,dm),
C$OMP+            SHARED(nx,ny,nz,igrd,bhdens,velx,vely,ro,d3x2,frhoin),
C$OMP+            REDUCTION(+:pma,pvx,pvy)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            if ( bhdens(i,j,k,igrd) .ne. zero ) then
               dm =  d3x2 * max(ro(i,j,k)-frhoin,zero)
               pma = pma +                    dm
               pvx = pvx + velx(i,j,k,igrd) * dm
               pvy = pvy + vely(i,j,k,igrd) * dm
            endif
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO

         velx1 = velx1 + ( pvx - velx1*pma ) / pmass1
         vely1 = vely1 + ( pvy - vely1*pma ) / pmass1

      endif

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,bhdens,densty,frhoin)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         if ( bhdens(i,j,k,igrd) .ne. zero )
     &      densty(i,j,k,igrd) = frhoin
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO


      return
      end
c------------------------------------------------------------------------
      subroutine PaWi(igrd)                                             !DRW - Does this ever get called?
c --- change interaction between gas and BH to Paczynski-Wiita potential
c --- update velocity of BH and velocity of gas accordingly

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      integer*4 igrd
      real*8    diffpawi, dpawimax, d3x2, rs, pawidx, pawidy, pawidz
      real*8    dis, disq, disx, disy, disz, disqy, disqz, dxt, ek(q)
      real*8    rota, z1, z2, beta, beta2, rin, r1


      d3x2 = two * delx(igrd)**3
      dxt = delx(indxgr(lone,lone))

      gmt = g * pmass1 * dtgr(igrd)
      rs   = two * g * pmass1 / cc**2

      if ( akerr .eq. zero ) then
         rota = zero
      elseif ( abs(akerr) .ge. one ) then
         rota = abs(angm1) * cc / (g * pmass1**2)
      else
         rota = abs(akerr)
      endif

      print *, 'PaWi???'

      if ( (rota.lt.zero) .or. (rota.gt.one) ) then
         write(*,*) 'rota:',rota,'  akerr:',akerr
         stop
      endif

      thrd = one/3D0
      z1 = one  +  (one-rota**2)**thrd
     &                    *( (one+rota)**thrd + (one-rota)**thrd )
      z2 = sqrt( 3D0*rota**2 + z1**2 )
      rin = 3D0 + z2 - sign( sqrt( (3D0-z1)*(3D0+z1+two*z2) ), angm1 )
c      rin = 3D0 + z2 - sqrt( (3D0-z1)*(3D0+z1+two*z2) )    ! prograde
c      rin = 3D0 + z2 + sqrt( (3D0-z1)*(3D0+z1+two*z2) )    ! retrograde
      r1 = one + sqrt(one-rota**2)
      rbhkerr = ( rin + r1 ) / two   *rs/two
c                        average      rin,r1 in units of rs/2

      beta = rin / r1  - one
      beta2 = two - beta
      r1 = r1/two * rs       ! change to  Schwarzschild Radius at rs

      disl = 1.8D0/two * rbhkerr
      dpawimax = one/(disl-r1)**2 - one/disl**2
      dpfac = one/(disl-r1)**2

      pvx = zero
      pvy = zero

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,densty,ro)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         ro(i,j,k) = densty(i,j,k,igrd)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

c     zero area covered by finer grids
      ifgr = idxfgr(igrd)
      if( ifgr .ne. -lone ) call rozero(ifgr)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,disx,disy,disz,
C$OMP+           disqy,disqz,disq,dis,ek,diffpawi,pawidx,pawidy,pawidz),
C$OMP+     SHARED(nx,ny,nz,igrd,posy1,posx1,dxt,bhdens,r1,ro,beta,beta2,
C$OMP+             small,velx,vely,velz,energy,dpawimax,dpfac,disl,gmt),
C$OMP+     REDUCTION(-:pvx,pvy)
      do j = 1, ny
         disy  = ( posy1 - topos(dble(j),ltwo,igrd) ) *dxt
         disqy = disy**2
      do k = 1, nz
         disz  = ( half  - topos(dble(k),lthree,igrd))*dxt
         disqz = disz**2  + disqy
      do i = 1, nx
         disx  = ( posx1 - topos(dble(i),lone,igrd) ) *dxt
         disq  = disx**2  + disqz
         if ( bhdens(i,j,k,igrd) .eq. zero ) then   !OUTside BH
            disq = max( disq, (0.75*r1)**2 )
            dis = sqrt(disq)
            ek(i) =  half *( velx(i,j,k,igrd)**2 + vely(i,j,k,igrd)**2 +
     &                       velz(i,j,k,igrd)**2 )
            energy(i,j,k,igrd) = energy(i,j,k,igrd) - ek(i)
            diffpawi = one / (dis**beta2 * (dis-r1)**beta)  -  one/disq
            if ( diffpawi .gt. dpawimax )
     &         diffpawi = dpfac*(disl/dis)**2.2 - one/disq
            diffpawi = gmt * diffpawi / dis
            pawidx = diffpawi * disx
            pawidy = diffpawi * disy
            pawidz = diffpawi * disz
            velx(i,j,k,igrd) = velx(i,j,k,igrd)  +  pawidx
            vely(i,j,k,igrd) = vely(i,j,k,igrd)  +  pawidy
            velz(i,j,k,igrd) = velz(i,j,k,igrd)  +  pawidz
            pvx = pvx - pawidx * ro(i,j,k)
            pvy = pvy - pawidy * ro(i,j,k)
            ek(i) =  half *( velx(i,j,k,igrd)**2
     &                     + vely(i,j,k,igrd)**2 + velz(i,j,k,igrd)**2 )
            energy(i,j,k,igrd) =  ek(i) +
     &                 max( small*ek(i), energy(i,j,k,igrd) )
         endif
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      velx1 = velx1 + pvx / pmass1 * d3x2
      vely1 = vely1 + pvy / pmass1 * d3x2

C$OMP PARALLEL DO DEFAULT (NONE), PRIVATE(i,j,k), SHARED(nx,ny,nz,
C$OMP+       igrd,bhdens,velx1,vely1,energy,smalle,velx,vely,velz)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         if ( bhdens(i,j,k,igrd) .ne. zero ) then   ! INside BH
            velx  (i,j,k,igrd) = velx1
            vely  (i,j,k,igrd) = vely1
            velz  (i,j,k,igrd) = zero
            energy(i,j,k,igrd) = smalle + half*( velx1**2 + vely1**2 )
         endif
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c------------------------------------------------------------------------
      subroutine pushbh
c --- advance black hole position, time-centered algorithim

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      !print *, 'grb, pushbh'

      itopgr=  indxgr(lone,lone)

cc add acceleration from centrifugal force in rotating frame of reference
cc to change velocities
c      if ( relaxon .eq. lone ) then
cc         velx1 = velx1 + dtgr(itopgr)*omega**2*(posx1-xcm)*delx(itopgr)
cc         vely1 = vely1 + dtgr(itopgr)*omega**2*(posy1-ycm)*delx(itopgr)
cc         something wrong, so brute force:
c          velx1 = zero
c          vely1 = zero
c      endif

      posx1 = posx1 + half*(velx1+vlox1)*dtgr(itopgr)/delx(itopgr)
      posy1 = posy1 + half*(vely1+vloy1)*dtgr(itopgr)/delx(itopgr)

      if ( (posx1.lt.1) .or. (posx1.gt.qx) .or.
     &     (posy1.lt.1) .or. (posy1.gt.qy)       ) then
         call datout
         write(*,*) '*** BH beyond grids! ***'
         call exit(1)
         stop'*** stop in pushbh ***'
      endif

      vlox1 = velx1
      vloy1 = vely1

      do lev = mlev, 2, -1
         call fineup(lev,lzero)
      enddo

      return
      end
c------------------------------------------------------------------------
      subroutine gravtyb(igrd)
c --- set correct density for black hole, only for calculating pot.

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      integer*4 igrd

      !print *, 'grb, gravtyb'

      call bhdenup(igrd)
      call gravty(igrd)
      call bhdendown(igrd)

      return
      end
c------------------------------------------------------------------------
      subroutine bhdenup(igrd)
c --- set correct density for black hole, only for calculating pot.
      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      integer*4 igrd

C$OMP PARALLEL DO DEFAULT(NONE), private(i,j,k),
C$OMP+            SHARED(nx,ny,nz,bhdens,densty,bhdfac,igrd)

      !print *, 'grb, bhdenup'

      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         if ( bhdens(i,j,k,igrd) .ne. zero )
     &      densty(i,j,k,igrd) = bhdens(i,j,k,igrd)*bhdfac
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end

c------------------------------------------------------------------------
      subroutine bhdendown(igrd)
c --- reset density back down

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'compct.cmn'

      integer*4 igrd

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+         SHARED(nx,ny,nz,igrd,bhdens,densty,rhoin)

      !print *, 'grb, bhdendown'

      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         if ( bhdens(i,j,k,igrd) .ne. zero )
     &      densty(i,j,k,igrd) = 5.D0*rhoin  ! 5 because of calm
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c-----------------------------------------------------------------------
      subroutine velfudge(igrd)
c --- for 64 grid models push binary back together by amount
c     lost through low numerical resolution: 0.5km/ms
c --- instead change kinetic energy: 1% per ms

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'eos.cmn'
      real*8  rh(qx), tm(qx), ei(qx), pr(qx), gc(qx), ge(qx), yea(qx),
     &        dtm(qx), ek(qx), tgr, rhl(qx), ent(qx),
     &        ergkinfac, velfac
      integer*4 inum, inx(qx), igrd

      !print *, 'grb, velfudge'

      if ( qx .ne. 64 ) return
      if ( gasmas .lt. 1.62E0*solmas ) return

      ergkinfac = 1. + 0.01D0 / 1.D-3 * dtgr(igrd)
      velfac = sqrt(ergkinfac)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,ek),
C$OMP+            SHARED(igrd,densty,velx,vely,velz,energy,velfac,rhoin)
      do  j = 1, qy
      do  k = 1, qz

         do  i = 1, qx
         if ( densty(i,j,k,igrd) .gt. 4.D0*rhoin ) then
           ek(i) = half * (
     &      velx(i,j,k,igrd)**2+vely(i,j,k,igrd)**2+velz(i,j,k,igrd)**2)
           energy(i,j,k,igrd) = energy(i,j,k,igrd) - ek(i)
           velx(i,j,k,igrd) = velx(i,j,k,igrd) * velfac
           vely(i,j,k,igrd) = vely(i,j,k,igrd) * velfac
           velz(i,j,k,igrd) = velz(i,j,k,igrd) * velfac
           energy(i,j,k,igrd) = energy(i,j,k,igrd) + half * (
     &      velx(i,j,k,igrd)**2+vely(i,j,k,igrd)**2+velz(i,j,k,igrd)**2)
         endif
         enddo

      enddo
      enddo
C$OMP END PARALLEL DO

      velx1 = velx1 * velfac
      vely1 = vely1 * velfac

      return
      end
c-----------------------------------------------------
      subroutine gasumaso(wri)

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'eos.cmn'

      character*6 wri

      !print *, 'grb, gasumaso'

      itopgr=  indxgr(1,1)
      gasmas = zero
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,j,i), REDUCTION(+:gasmas),
C$OMP+            SHARED(nx,ny,nz,densty,itopgr)
      do 244 j=1,ny
      do 244 k=1,nz
      do 244 i=1,nx
         gasmas = gasmas + densty(i,j,k,itopgr)
  244 continue
C$OMP END PARALLEL DO
      gasmas = two*gasmas * delx(itopgr)**3
      write(*,'(2X,A,A,1F16.12,1P,1E20.10)') wri,' gasumas :',
     &             gasmas/solmas ! , densty(9,14,5,3)


      return
      end
c-----------------------------------------------------
      subroutine gasumas(wri,igr)

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'grnest.cmn'
      include 'eos.cmn'

      character*6 wri
      integer*4 igr

      !print *, 'grb, gasumas'

      gasmas = zero
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,j,i), REDUCTION(+:gasmas),
C$OMP+            SHARED(nx,ny,nz,densty,energy,igr,velx,vely,velz)
      do 244 j=1,ny
      do 244 k=1,nz
      do 244 i=1,nx
c         gasmas = gasmas + densty(i,j,k,igr)
c         gasmas = gasmas + energy(i,j,k,igr)*densty(i,j,k,igr)
         gasmas = gasmas + densty(i,j,k,igr) * (
     &        velx(i,j,k,igr)**2 + vely(i,j,k,igr)**2 +
     &        velz(i,j,k,igr)**2  )
  244 continue
C$OMP END PARALLEL DO
      gasmas = two*gasmas * delx(igr)**3
c      write(*,'(2X,A,A,1F16.12,1P,1E20.10)') wri,' gasumas :',
c     &             gasmas/solmas ! , densty(9,14,5,3)
      write(*,'(2X,A,A,1E20.10,1P,1E20.10)') wri,' gasumas :',
     &             gasmas/solmas ! , densty(9,14,5,3)


      return
      end
