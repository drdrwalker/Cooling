c-----------------------------------------------------------------------
      subroutine input

c     read input parameters from file inidata.charyb

      include 'qparam.cmn'
      include 'ppdisk.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'rays.cmn'
      include 'compct.cmn'

      character*8  label1

      character*72 text
      character*44 txt
      character*5  txtxt
      data         txtxt /'...  '/

      home = '/home/DiskSimulations/Cooling/'
      homedata = '/home/DiskSimulations/Cooling/'

      !print *, 'inout, input'

      if ( qx.eq.32 ) then
         open(4,file=home//'inidata_merge_32'
     &        ,form='formatted')
      elseif ( qx.eq.64 ) then
         open(4,file=home//'inidata_merge_64'
     &        ,form='formatted')
      elseif ( qx.eq.128 ) then
         open(4,file=home//'inidata_merge_128'
     &        ,form='formatted')
      elseif ( qx.eq.256 ) then
         open(4,file=home//'inidata_merge_256'
     &        ,form='formatted')
      else
         write(*,*) 'input, inout,f,  wrong qx :',qx
         stop
      endif

      do 5 n=1,4
        read (4,*) text
        write(*,*) text
5       continue
      print *, ' '

      read (4,*) txt, label, nend
      write(*,'(a,a,a,1I6)') txt, txtxt, label, nend
      label1 = 'nend'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, tmax
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, tmax
      label1 = 'tmax'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, irstrt
      write(*,'(a,a,a,1I6)') txt, txtxt, label, irstrt
      label1 = 'irstrt'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, basenm
      write(*,'(a,a,a,a)') txt, txtxt, label, basenm
      label1 = 'basenm'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, suffix
      label1 = 'suffix'
      if(label .ne. label1) go to 10
      if ( irstrt .ne. lzero ) then
         open(1,file=basenm//'_next',form='formatted')
         read (1,'(a)',end=11) suffix
         close(1)
c         close(1,status='delete')
      endif
      write(*,'(a,a,a,a)') txt, txtxt, label, suffix

      read (4,*) txt, label, nrstrt
      write(*,'(a,a,a,1I6)') txt, txtxt, label, nrstrt
      label1 = 'nrstrt'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, trstrt
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, trstrt
      label1 = 'trstrt'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, nout
      write(*,'(a,a,a,1I6)') txt, txtxt, label, nout
      label1 = 'nout'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, tout
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, tout
      label1 = 'tout'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, dtini
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, dtini
      label1 = 'dtini'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, dtmin
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, dtmin
      label1 = 'dtmin'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, dtmax
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, dtmax
      label1 = 'dtmax'
      if(label .ne. label1) go to 10

c     read (4,*) txt, label, cfl
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, cfl
c     label1 = 'cfl'
c     if(label .ne. label1) go to 10
c' Courant factor ............................' 'cfl'       0.8
      cfl = 0.8D0

c     read (4,*) txt, label, cgr
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, cgr
c     label1 = 'cgr'
c     if(label .ne. label1) go to 10
c' Gravitational timestep factor .............' 'cgr'       0.7
      cgr = 0.6D0

c     read (4,*) txt, label, nriem
c     write(*,'(a,a,a,1I6)') txt, txtxt, label, nriem
c     label1 = 'nriem'
c     if(label .ne. label1) go to 10
c' Number of iterations in Riemann solver ....' 'nriem'     4
      nriem = 5

c     read (4,*) txt, label, cvisc
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, cvisc
c     label1 = 'cvisc'
c     if(label .ne. label1) go to 10
c' Artificial viscosity constant .............' 'cvisc'     0.1
      cvisc = zero

c     read (4,*) txt, label, small
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, small
c     label1 = 'small'
c     if(label .ne. label1) go to 10
c' Cut-off value .............................' 'small'     1.d-10
      small = 1.D-20

c     read (4,*) txt, label, smlrho
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, smlrho
c     label1 = 'smlrho'
c     if(label .ne. label1) go to 10
c' Cut-off value for density .................' 'smlrho'    1.d-35
      !smlrho = 1.D5
      smlrho = 1.D-30

c     read (4,*) txt, label, smallp
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, smallp
c     label1 = 'smallp'
c     if(label .ne. label1) go to 10
c' Cut-off value for pressure ................' 'smallp'    1.d-25
      smallp = 1.D-25

c     read (4,*) txt, label, smalle
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, smalle
c     label1 = 'smalle'
c     if(label .ne. label1) go to 10
c' Cut-off value for energy ..................' 'smalle'    1.d-10
      smalle = 1.D-20

c     read (4,*) txt, label, smallu
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, smallu
c     label1 = 'smallu'
c     if(label .ne. label1) go to 10
c' Cut-off value for velocity ................' 'smallu'    1.d-05
      smallu = 1.D-21
c      smallu = sqrt(smalle / 6) !This was chosen to avoid supercritical
                                !flows in small value regions.

c     read (4,*) txt, label, smlche(1)
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, smlche
c     label1 = 'smlche'
c     if(label .ne. label1) go to 10
c' Cut-off value for chemical composition ....' 'smlche'     1.d-10
      do m = 1, qc
         smlche(m) = 1.D-10
      enddo

c     read (4,*) txt, label, igodu
c     write(*,'(a,a,a,1I6)') txt, txtxt, label, igodu
c     label1 = 'igodu'
c     if(label .ne. label1) go to 10
c' Use Godunov method, if  IGODU = 1 .........' 'igodu'     0
      igodu = 0

c     read (4,*) txt, label, nsdim
c     write(*,'(a,a,a,1I6)') txt, txtxt, label, nsdim
c     label1 = 'nsdim'
c     if(label .ne. label1) go to 10
      nsdim = 3

      read (4,*) txt, label, gridlx
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, gridlx
      label1 = 'gridlx'
      if(label .ne. label1) go to 10

      ! Outer disk radius
      r_out = gridlx / 2

c     read (4,*) txt, label, gridly
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, gridly
c     label1 = 'gridly'
c     if(label .ne. label1) go to 10
      gridly = gridlx

c     read (4,*) txt, label, gridlz
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, gridlz
c     label1 = 'gridlz'
c     if(label .ne. label1) go to 10
      gridlz = gridlx*dble(qz)/dble(qx)

      read (4,*) txt, label, nx
      nx = qx
      write(*,'(a,a,a,1I6)') txt, txtxt, label, nx
      label1 = 'nx'
      if(label .ne. label1) go to 10
      if (nx .gt. qx)   stop'** nx **'

      read (4,*) txt, label, rifrac
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, rifrac
      label1 = 'rifrac'
      if(label .ne. label1) go to 10

      ! Inner disk radius
      r_in  = rifrac * r_out

c     read (4,*) txt, label, ny
c     write(*,'(a,a,a,1I6)') txt, txtxt, label, ny
c     label1 = 'ny'
c     if(label .ne. label1) go to 10
c     if (ny .gt. qy)   stop'** ny **'
      ny = qy

c     read (4,*) txt, label, nz
c     write(*,'(a,a,a,1I6)') txt, txtxt, label, nz
c     label1 = 'nz'
c     if(label .ne. label1) go to 10
c     if (nz .gt. qz)   stop'** nz **'
      nz = qz

      read (4,*) txt, label, itstp
      write(*,'(a,a,a,1I6)') txt, txtxt, label, itstp
      label1 = 'itstp'
      if(label .ne. label1) go to 10

c     read (4,*) txt, label, epsiln
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, epsiln
c     label1 = 'epsiln'
c     if(label .ne. label1) go to 10
c' PPM parameter used to detect shocked zones ' 'epsiln'    0.33
      epsiln = 0.33D0

c     read (4,*) txt, label, omg1
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, omg1
c     label1 = 'omg1'
c     if(label .ne. label1) go to 10
c' PPM dissipation parameter omega1 ..........' 'omg1'      0.75
      omg1 = 0.75D0

c     read (4,*) txt, label, omg2
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, omg2
c     label1 = 'omg2'
c     if(label .ne. label1) go to 10
c' PPM dissipation parameter omega2 ..........' 'omg2'      10.0
      omg2 = 10.0D0

      read (4,*) txt, label, bidist
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, bidist
      label1 = 'bidist'
      if(label .ne. label1) go to 10
c  Distance of White Dwarf from center (in cm)' 'bidist'    1.D+0
c      bidist = 1.D0

      read (4,*) txt, label, pmass1
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, pmass1
      label1 = 'pmass1'
      if(label .ne. label1) go to 10

      pmass1 = pmass1 * solmas

      read (4,*) txt, label, pmass2
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, pmass2
      label1 = 'pmass2'
      if(label .ne. label1) go to 10

      pmass2 = pmass2 * pmass1

      read (4,*) txt, label, npawi
      write(*,'(a,a,a,1I6)') txt, txtxt, label, npawi
      label1 = 'npawi'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, nspin
      write(*,'(a,a,a,1I6)') txt, txtxt, label, nspin
      label1 = 'nspin'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, akerr
      write(*,'(a,a,a,1P,1E16.9,0P)') txt, txtxt, label, akerr
      label1 = 'akerr'
      if(label .ne. label1) go to 10

c     read (4,*) txt, label, wdrad
c     write(*,'(a,a,a,1P,1E16.9, 0P)') txt, txtxt, label, wdrad
c     label1 = 'wdrad'
c     if(label .ne. label1) go to 10
      wdrad = 0.1D0

c     read (4,*) txt, label, soparp
c     write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, soparp
c     label1 = 'soparp'
c     if(label .ne. label1) go to 10
      soparp = 0.5D0

      read (4,*) txt, label, soparg
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, soparg
      label1 = 'soparg'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, mgonoff
      if ( relaxon .eq. lone ) mgonoff = 0
      write(*,'(a,a,a,1I6)') txt, txtxt, label, mgonoff
      label1 = 'mgonoff'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, mconoff
      write(*,'(a,a,a,1I6)') txt, txtxt, label, mconoff
      label1 = 'mconoff'
      if(label .ne. label1) go to 10

      corrmax = 0.99D0

c     read (4,*) txt, label, ifilm
c     write(*,'(a,a,a,1I6)') txt, txtxt, label, ifilm
c     label1 = 'ifilm'
c     if(label .ne. label1) go to 10
      ifilm = 4
      if ( nx .eq. 32 ) ifilm = 2
      if ( nx .eq. 64 ) ifilm = 4
      if ( nx .eq. 128) ifilm = 8
      if ( nx .eq. 256) ifilm = 8

      read (4,*) txt, label, dymax
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, dymax
      label1 = 'dymax'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, tpmax
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, tpmax
      label1 = 'tpmax'
      if(label .ne. label1) go to 10

      read (4,*) txt, label, bcool
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, bcool
      label1 = 'bcool'
      if(label .ne. label1) go to 10

      ! Are we using the Stamatellos method? One for yes, zero for no.
      sbyn = zero                                                        !DRW - Stefan-Bolztmann yes/no

      read (4,*) txt, label, tirr
      write(*,'(a,a,a,1P,1E17.10,0P)') txt, txtxt, label, tirr
      label1 = 'tirr'
      if(label .ne. label1) go to 10

      ! In the Stamatellos method we need tirr**4                       ! DRW
      tirr4 = tirr**4

      ! Irradiation energy                                              ! DRW
      ergirr = specr * tirr / (gamma-one)

      ! Opacity parameters                                              ! DRW
      ! Neater to put them in arrays maybe, but there is only a few and
      ! they vary in size significantly. Not likely to need beyond the
      ! first four.
      ! kappa parameter
      kap1 = 2.D-4
      kap2 = 2.D16
      kap3 = 0.1D0
      kap4 = 2.D81
      kap5 = 1.D-8
      kap6 = 1.D-36
      kap7 = 1.5D20
      kap8 = 0.348D0
      ! Exponent on the density. For regions 1, 2, 3 and 8 this
      ! parameter is 0 and we do not need to include it.
      rhoex4 = 1.D0
      rhoex5 = 2/3
      rhoex6 = 1/3
      rhoex7 = 1.D0
      ! Exponent on the temperature. Not included is temex8 which is
      ! just 0 anyway, and this is a region that should not be reached
      ! anyway.
      temex1 = 2.D0
      temex2 = -7.D0
      temex3 = 0.5D0
      temex4 = -24.D0
      temex5 = 3.D0
      temex6 = 10.D0
      temex7 = -2.5D0
      ! First two transition temperatures are fixed.
      ttran1 = 166.8101D0
      ttran2 = 202.6768D0

c We either solve the cooling eqn using Euler method or CN. I don't
c think I will be likely to change this so will leave it like this.
      cnyn = zero

      close(4)

      open(1,file='./'//basenm//'_cont',form='formatted')
      write(1,'(1i1)') lone
      close(1)

      write(*,*)
      return


  10  continue
      close(1)
      print *, ' '
      write(*,1000)
 1000 format(' incorrect input deck')
      write(*,1001) label,label1
 1001 format(' label = ',a6,'  expected label = ',a6)
      stop

  11  continue
      write(*,*) basenm//'_next FILE is empty'
      call exit(1)
      stop

      end


c-----------------------------------------------------------------------
      subroutine restrt(ird)

c     ird=1 : read , =2 : write , data for deferred calculation

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'
      include   'rays.cmn'
      real*4  dty(qx,qy)
      character*7 testfil


c      include   'ppms.cmn'

c      testfil= basenm // 'T' // suffix
c ---

c      open(33,file=testfil,form='unformatted')

      !print *, 'inout, restrt'

      if(ird.eq.1) then         ! read

         open(2,file=rstfil,form='unformatted')
         read(2)   ! dummy read to skip first record
         read(2)   densty
         read(2)   velx
         read(2)   vely
         read(2)   velz
         read(2)   energy
         read(2)   chem
         read(2)   temper
         read(2)   tmpent
         read(2)   time, deltax,
     &             bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,
     &             rhoin, uin, utin, uttin, pin, ein, gin, tin,
     &             gamein, gamcin, chein,
     &             nstep, igeom, dto, cenx, ceny, cenz, iframe,
     &             pmass1, posx1, posy1, velx1, vely1, vlox1, vloy1,
     &             pmass2, posx2, posy2, velx2, vely2, vlox2, vloy2,
     &             enuloss, eneloss, enaloss, enxloss, bhdfac, angm1,
     &             indxgr, idxcgr, idxfgr, levlgr, norgin, angm2, eint1,
     &             ishedul, ifinup, nodenr, delx, dtgr, maxlev, dtgg,
     &             cmgax, cmgay, signif, dij3dt, egquel, gmunb, angunb,
     &             angcc1,  thestar, phistar, iendf, tlzgrav, omega
         close(2)

         if ( suffix .eq. 'aa' ) then
            nstep = 0
            time = zero
         endif

         write(*,'(a,a,a,i5,a,f8.3)') 'restarting from file ',rstfil,
     &      '   nstep=',nstep,'   time=',time*1.D3
         write(*,*)

         call filnam

         itopgr=  indxgr(1,1)
         gasmas = zero
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,j,i),
C$OMP+            SHARED(nx,ny,nz,densty,itopgr), REDUCTION(+:gasmas)
         do 244 j=1,ny
         do 244 k=1,nz
         do 244 i=1,nx
            gasmas = gasmas + densty(i,j,k,itopgr)
  244    continue
C$OMP END PARALLEL DO
         gasmas = two*gasmas * delx(itopgr)**3
         write(*,*) 'restrt, gasmas :', gasmas/solmas
c --- first without fineup
         do igrd = ngrd, 1, -1
            call sucks(igrd,lzero)
         enddo
c     now with fineup
         do igrd = ngrd, 1, -1
            call sucks(igrd,lone)
         enddo

         if ( relaxon .eq. ltwo ) call wakeup

c         angm1 = sign( angm1, akerr )

      else                      ! write

         if ( relaxon .eq. lone ) then
            time  = zero
            nstep = 0
         endif

         open(2,file=rstfil,form='unformatted')
         write(2)  qx, qy, qz, qb, qc, q, ngrd, nfine, mlev
         write(2)  densty
         write(2)  velx
         write(2)  vely
         write(2)  velz
         write(2)  energy
         write(2)  chem
         write(2)  temper
         write(2)  tmpent
         write(2)  time, deltax,
     &             bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,
     &             rhoin, uin, utin, uttin, pin, ein, gin, tin,
     &             gamein, gamcin, chein,
     &             nstep, igeom, dto, cenx, ceny, cenz, iframe,
     &             pmass1, posx1, posy1, velx1, vely1, vlox1, vloy1,
     &             pmass2, posx2, posy2, velx2, vely2, vlox2, vloy2,
     &             enuloss, eneloss, enaloss, enxloss, bhdfac, angm1,
     &             indxgr, idxcgr, idxfgr, levlgr, norgin, angm2, eint1,
     &             ishedul, ifinup, nodenr, delx, dtgr, maxlev, dtgg,
     &             cmgax, cmgay, signif, dij3dt, egquel, gmunb, angunb,
     &             angcc1, thestar, phistar, iendf, tlzgrav, omega
         close(2)

         open(1,file=basenm//'_next',form='formatted')
         write(1,'(a)') suffix
         close(1,status='keep')

         call filnam

         nrst = 0
         trst = zero

c         call pltout

      endif                     ! of read or write

      call inigrav

c --- density contains now BH for printing out
c --- print out gpot and density, when they are fresh
c      write(33) time
c      do 56 k=1,2
c       do 51 nngr = 1, ngrd
c
c         do 55 j=1,qy
c         do 55 i=1,qx
c           dty(i,j) = densty(i,j,k,nngr)
c   55    continue
c         write(33) dty


c         do 52 j=1,qy
c         do 52 i=1,qx
c           dty(i,j) = sign(
c     &        max(1.D-30,abs(gpot(i,j,k,nngr))),gpot(i,j,k,nngr))
c   52    continue
c         write(33) dty

c   51  continue
c   56  continue
c       close(33)

c --- set back density to small values for Hydro
c --- do this on all grids !
       do 172 lev = 1, mlev
      do 172 num = 1, indxgr(lzero,lev)
         igrd = indxgr(num,lev)
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,igrd,bhdens,densty,rhoin)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         if ( bhdens(i,j,k,igrd) .ne. zero )
     &      densty(i,j,k,igrd) = 5.D0*rhoin  ! 5 because of calm
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
  172 continue
      polx2 = posx2
      poly2 = posy2

      write(*,41) '  ms       posx1      posy1     velx1     vely1  '//
     &           '  bidist       eglum       getherm   '//
     &           '  egquel    rhmax   mx  my  mz g  gmunb  '
      write(*,41) ' nstep    cmcax    cmcay     vcmcx       vcmcy   '//
     &           '    lzgas       gekin      geintern    epotot '//
     &           '       ep12       epot1  '
      write(*,41) '          velxns     velyns     hplus      hcross'//
     &           '    ergboun   gasmas    pmas   enuloss   eneloss  '//
     &           '    enaloss      enxloss   ergtot'
      write(*,41) '          d2dtxx         d2dtyy        d2dtzz  '//
     &           '     d2dtxy        d1mass '//
     &           '     d2mass     d3mass       netot        ltot '
      write(*,41) '          avene   avena   avenx    tlzgrav     '//
     &           ' sumbulk      sumshear    sumneutr    sumshock     '//
     &           ' dkmass       dcmass     pmass1'
      write(*,41) '          posx2      posy2     cyl1mas   cyl2mas  '//
     &           '  cyl3mas    con1mas    con2mas    con3mas   '//
     &           '    angm1       eint1      limc'
      write(*,41) '          dij3dt1       dij3dt2      dij3dt3   '//
     &           '     dij3dt5       dij3dt6      dij3dt9   '
      write(*,41) '          vlmax   mx my mz mg tmpmax itm jtm ktm '//
     &           ' itmgr  ekin1     angcc1     egprl   bidiav'
      write(*,41) ' '

   41 format (A)

c --- 1ms: time in millisecs,2pos(x/y)1: BH position in zones on top grid
c --- 3vel(x/y)1: BH velocity in cm/sec,4bidist: distance density max of
c --- compact objects in km,5eglum: gw lum in erg/s calc by quadrupole formula,
c --- 6getherm: thermal energy in erg,egquel: gw energy in erg as sum over all
c --- cells(local !),7rhmax: density max in g/cm^3, 8,9,10,11coords of density max, grid
c --- on which dens max lies,12gmunb: unbound gas mass
c --- 13nstep: timestep(two make one cycle), 14,15cmca(x/y): coord of total cm on top
c --- grid,vcmc(x/y): velocity of cm, lzgas: angular mom of gas relative to cm
c --- gekin: kinetic energy of gas,geintern: internal energy of gas as difference
c --- of total (energy in hydro) and kinetic,epotot:complete potential energy on grid
c --- ,ep12:potential energy between BH and NS,epot1: self potential energy of BH,
c --- vlmax: maximum velocity, coords and grid of max vel, tmpmax: max temperature,
c --- vel(x/y)ns: velocity of NS, hplus/hcross: amplitudes of gw, ergboun:??
c --- gasmas: mass of gas in solmas,pmas: mass of NS as 30km sphere,
c --- enuloss: energy loss by all neutrinos,eneloss: lum by enu,enaloss: lum by
c --- e.anti.nu,enxloss: lum by rest species, ergtot: total energy on grid
c --- d2dtxx:second time derivative of quadrupol tensor,d(1/2/3)mass: mass up to
c --- a certain density,netot: electron number,ltot: total angular momentum,
c --- avene/a/x: average neutrino energy, tlzgrav: total ang mom loss due to gw
c --- sumbulk: entropy source,bulk viscosity?,sumshear: entr. source shear viscosity?,
c --- sumneutr:entr source,neutrino term,sumshock: entropy generation in shocks,
c --- dkmass: mass with ang mom over certain limit->stays in disc,dcmass: like dkmass
c --- different limit, pmass1: BH mass,pos(x/y)2: coords of NS (cm),cyl(1/2/3)mas:
c --- mass within a cylinder over the orbital plane,con(1/2/3)mas: mass within a cone,
c --- angm1: self ang mom of BH due to accretion (mor-Style), dij3dti: third time
c --- derivative of quadrupol tensor, limc: orbital ang mom of BH relative to the total mass
c --- cm (values after acretion used!),itm/jtm/ktm/itmgr: coords of temp max and grid,
c --- ekin1: kinetic energy of BH,angcc1: ang mom of accreted gas relative to tot mass cm,
c --- egprl: ang mom of mass that has left grid
      return
      end
c-----------------------------------------------------------------------
      subroutine datout

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'

c  use the following arrays instead of defining new ones.
c  implies that these arrays are not used otherwise at the moment.


      real*8 pot(q), vzsq(q), vsq(q), ecm(q), ekcm(q), rho(q), ei(q),
     &     ewd(q), ekwd(q), vgx(q), vgy(q), eth(q), tmpm(q), v(q)
      equivalence  (pot,u), (vzsq,ut), (vsq,utt), (ecm,enu),
     &             (ekcm,unu), (vgx,utl), (vgy,uttl),
     &             (ewd,vl), (ekwd,vr)

      real*8 dsq, dsqz,gco2, d1mass, d2mass, d3mass, dkmass, dcmass
      real*8 limp, limc, lzgas, ltot, lzgtmp, rnetot, visfac,
     &       jstark, jstarc, rs3sq, dxysq, cyl1mas, cyl2mas, cyl3mas,
     &       con1mas, con2mas, con3mas, dyma,fakemass


c     logical overlap
c     logical overx1, overx2, overy1, overy2

      include 'eos.cmn'

      real*8  gcc, gasm(ngrd)
      integer*4 iye(qx)

c--------

c      do 160 lev = mlev, 2, -1
c  160    call fineup(lev,lzero)

      !print *, 'inout, datout'

      itopgr=indxgr(lone,lone)

      d1x = delx(itopgr)
      d3x = delx(itopgr)**3
      tms = time *1000.D0

      topso= soparp / dble( 2 ** (maxlev(lone)-lone) )
c ---
c      fakemass = 1.63D0 * solmas
      rsbh  = two*g*pmass1/cc**2  / d1x
      rsbh3 = 3.D0 * rsbh
      twothrd = two/3.D0

c  black hole limit radius

      rfac = 3.D0
      bhfac = two * g / cc**2
      gasmas = 3.28D0 * solmas  ! only for next statement  !
      radbhl= bhfac * gasmas  * rfac
      rbhl1 = radbhl /ten      / d1x
      rbhl2 = radbhl /ten*2.D0 / d1x
      rbhl3 = radbhl /ten*3.D0 / d1x
      rbhl4 = radbhl /ten*4.D0 / d1x
      rbhl5 = radbhl /ten*5.D0 / d1x
      rbhl6 = radbhl /ten*6.D0 / d1x
      rbhl7 = radbhl /ten*7.D0 / d1x
      rbhl8 = radbhl /ten*8.D0 / d1x
      rbhl9 = radbhl /ten*9.D0 / d1x
      rbhl0 = radbhl /ten*ten/ d1x
      gbhl1 = zero
      gbhl2 = zero
      gbhl3 = zero
      gbhl4 = zero
      gbhl5 = zero
      gbhl6 = zero
      gbhl7 = zero
      gbhl8 = zero
      gbhl9 = zero
      gbhl0 = zero
      gmunb = zero
      angunb= zero

c     iwdgr = indxgr(ltwo,maxlev(ltwo))
      iwdgr = indxgr(lone,maxlev(lone))
      topx1 = posx1
      topy1 = posy1
      topx2 = posx2
      topy2 = posy2
c      topx1 = topos(posx1,lone,indxgr(lone,maxlev(lone)))
c      topy1 = topos(posy1,ltwo,indxgr(lone,maxlev(lone)))
c      topx2 = topos(posx2,lone,iwdgr)
c      topy2 = topos(posy2,ltwo,iwdgr)

c  total gas mass on grid, gas mass around WD,
c  position and velocity of gas centre of mass CM

      gama2 = zero
      gasmas= zero
      d1mass= zero
      d2mass= zero
      d3mass= zero
      cyl1mas=zero
      cyl2mas=zero
      cyl3mas=zero
      con1mas=zero
      con2mas=zero
      con3mas=zero
      denpo = zero
      cmgax = zero
      cmgay = zero
      vcmgx = zero
      vcmgy = zero
      rnetot= zero
      vlmax = zero
      rhmax = zero
      tmpmax= 0.01D0
      rnxh= dble(nx/2) + half

c--gridloop start
      do 500 levl = 1, mlev          ! loop over all levels
      rkf  = one / ( dble(nfine)**(levl-lone) )
      rkf3 = rkf**3
      do 500 num = 1, indxgr(lzero,levl) ! loop over all grids on levl
         lgrd = indxgr(num,levl)! lgrd=1 here --> one grid per level
         xopos= topos( zero, lone, lgrd)
         yopos= topos( zero, ltwo, lgrd)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,densty,chem,ro,po,lgrd)
         do j = 1, ny
         do k = 1, nz
         do i = 1, nx
            ro(i,j,k) = densty(i,j,k,lgrd)
            po(i,j,k) = chem(i,j,k,lgrd,lone) / ro(i,j,k)
         enddo
         enddo
         enddo
C$OMP END PARALLEL DO
c         po actually is Ye !!

c        zero area covered by finer grids
         ifgr = idxfgr(lgrd)
         if( ifgr .ne. -lone ) call rozero(ifgr)

c--gridloop start

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(nx,ny,nz,yopos,xopos,rkf,
C$OMP+                   ro,rkf3,rsbh,rsbh3,rnxh,topy1,topx1,twothrd),
C$OMP+            PRIVATE(i,j,k,dsq,dsqy,dsqz,dzc,rd,cyd,rho),
C$OMP+            REDUCTION(+:cyl1mas,con1mas,cyl2mas,con2mas,
C$OMP+                        cyl3mas,con3mas)
      do 201 j=1,ny
         dsq = (yopos+dble(j)*rkf-rnxh )**2
         dsqy= (yopos+dble(j)*rkf-topy1)**2
      do 201 k=1,nz
         dsqz= (     (dble(k)-half)*rkf )**2 + dsq
         dzc =       (dble(k-lone)-half)*rkf
      do 201 i=1,nx
         rd  = (xopos+dble(i)*rkf-rnxh )**2 + dsqz
         cyd = (xopos+dble(i)*rkf-topx1)**2 + dsqy
         rho(i) = ro(i,j,k) * rkf3

c        if (rd.le.rbhl1**2) gbhl1 = gbhl1 + rho(i)
c        if (rd.le.rbhl2**2) gbhl2 = gbhl2 + rho(i)
c        if (rd.le.rbhl3**2) gbhl3 = gbhl3 + rho(i)
c        if (rd.le.rbhl4**2) gbhl4 = gbhl4 + rho(i)
c        if (rd.le.rbhl5**2) gbhl5 = gbhl5 + rho(i)
c        if (rd.le.rbhl6**2) gbhl6 = gbhl6 + rho(i)
c        if (rd.le.rbhl7**2) gbhl7 = gbhl7 + rho(i)
c        if (rd.le.rbhl8**2) gbhl8 = gbhl8 + rho(i)
c        if (rd.le.rbhl9**2) gbhl9 = gbhl9 + rho(i)
c        if (rd.le.rbhl0**2) gbhl0 = gbhl0 + rho(i)

         if ( dzc .gt. rsbh ) then
            if ( cyd .le. rsbh3**2 ) cyl1mas = cyl1mas + rho(i)
            if ( dzc**2 .gt. cyd*( (sqrt(cyd)/rsbh3)**twothrd-one) )
     &            con1mas = con1mas + rho(i)
         endif
         if ( dzc .gt. two*rsbh ) then
            if ( cyd .le. rsbh3**2 ) cyl2mas = cyl2mas + rho(i)
            if ( dzc**2 .gt. cyd*( (sqrt(cyd)/rsbh3)**twothrd-one) )
     &            con2mas = con2mas + rho(i)
         endif
         if ( dzc .gt. 3.D0*rsbh ) then
            if ( cyd .le. rsbh3**2 ) cyl3mas = cyl3mas + rho(i)
            if ( dzc**2 .gt. cyd*( (sqrt(cyd)/rsbh3)**twothrd-one) )
     &            con3mas = con3mas + rho(i)
         endif

 201  continue
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,vv,tmpmm,rone,rho),
C$OMP+           SHARED(nx,ny,nz,lgrd,ro,rkf3,xopos,yopos,rkf,
C$OMP+                  velx,vely,velz,po,densty,smlrho,tmpent),
C$OMP+           REDUCTION(+:cmgax,cmgay,vcmgx,vcmgy,gasmas,rnetot,
C$OMP+                      d1mass,d2mass,d3mass),
C$OMP+           REDUCTION(max:vlmax,tmpmax,rhmax)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         rho(i) = ro(i,j,k)*rkf3

         cmgax  = cmgax  + rho(i)*(xopos+dble(i)*rkf)
         cmgay  = cmgay  + rho(i)*(yopos+dble(j)*rkf)
         vcmgx  = vcmgx  + rho(i)*velx(i,j,k,lgrd)
         vcmgy  = vcmgy  + rho(i)*vely(i,j,k,lgrd)
         gasmas = gasmas + rho(i)
         rnetot = rnetot + rho(i)*po(i,j,k)

         if ( densty(i,j,k,lgrd) .lt. 1.D10 ) d1mass = d1mass + rho(i)
         if ( densty(i,j,k,lgrd) .lt. 1.D11 ) d2mass = d2mass + rho(i)
         if ( densty(i,j,k,lgrd) .lt. 1.D12 ) d3mass = d3mass + rho(i)

         rone = rho(i) / max(rho(i),smlrho)
         tmpmm= tmpent(i,j,k,lgrd) * rone
         vv   = sqrt( velx(i,j,k,lgrd)**2 + vely(i,j,k,lgrd)**2  +
     &                velz(i,j,k,lgrd)**2 )  * rone

         vlmax = max( vlmax, vv )
         rhmax = max( rhmax, ro(i,j,k) )
         tmpmax= max( tmpmax, tmpmm )
      enddo
      enddo
      enddo

!      print *, "rnetot1:", rnetot                                      !DRW - Dealing with a pesky NaN - Not a problem any more?

C$OMP END PARALLEL DO

c the following loop has a race condition, but harmless, because
c only at most one CPU actually fills memory location
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,rhmax,ro,imaxx,imaxy,imaxz,imaxg,lgrd)
      do j = 1, ny
      do k = 1, nz
      do i = 1, nx
         if ( rhmax .eq. ro(i,j,k) ) then
            imaxx = i
            imaxy = j
            imaxz = k
            imaxg = lgrd
         endif
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      imavx = 0
      imavy = 0
      imavz = 0
      imavg = 0
c      imaxx = 0
c      imaxy = 0
c      imaxz = 0
c      imaxg = 0
      itm   = 0
      jtm   = 0
      ktm   = 0
      itmgr = 0

c      if (.true.) goto 500

c  find max of velocity

c      ivlmax = isamax(nx, v(1), lone)
c      vtmax = v(ivlmax)
c      if ( vlmax .lt. vtmax) then
c         vlmax = vtmax
c         imavx = ivlmax
c         imavy = j
c         imavz = k
c         imavg = lgrd
c      endif

c      do i = 1, nx
c      if ( abs(v(i)) .gt. vlmax ) then
c         vlmax = v(i)
c         imavx = i
c         imavy = j
c         imavz = k
c         imavg = lgrd
c      endif
c      enddo

c  find max of density

c      irhmax = isamax(nx, ro(1,j,k), lone)
c      rtmax = ro(irhmax,j,k)
c      if ( rhmax .lt. rtmax) then
c         rhmax = rtmax
c         imaxx = irhmax
c         imaxy = j
c         imaxz = k
c         imaxg = lgrd
c      endif

c      do i = 1, nx
c      rtmax = rho(i) / rkf3
c      if ( rtmax .gt. rhmax ) then
c         rhmax = rtmax
c         imaxx = i
c         imaxy = j
c         imaxz = k
c         imaxg = lgrd
c      endif
c      enddo

c  find max of temperature

c      itmm = ismax( nx, tmpm(1), lone )
c      if ( tmpmax .le. tmpm(itmm) )  then
c         itm  = itmm
c         jtm  = j
c         ktm  = k
c         itmgr = lgrd
c         tmpmax = tmpm(itmm)
c      endif

c      do i = 1, nx
c      if ( tmpm(i) .gt. tmpmax ) then
c         tmpmax= tmpm(i)
c         itm   = i
c         jtm   = j
c         ktm   = k
c         itmgr = lgrd
c      endif
c      enddo


 500  continue
c--gridloop end
c      write(*,*)'cmgax :',cmgax
c      write(*,*)'cmgax/gasmas :',cmgax/gasmas

      cmgax = cmgax / gasmas                                            !DRW - gasmas is 0.0, causung NaN cmgax and cmgay
      cmgay = cmgay / gasmas
      vcmgx = vcmgx / gasmas
      vcmgy = vcmgy / gasmas
      gasmas= two * gasmas * d3x
      d1mass= two * d1mass * d3x
      d2mass= two * d2mass * d3x
      d3mass= two * d3mass * d3x
      rnetot= two * rnetot * d3x
      gbhl1 = two * gbhl1  * d3x   * bhfac / d1x
      gbhl2 = two * gbhl2  * d3x   * bhfac / d1x
      gbhl3 = two * gbhl3  * d3x   * bhfac / d1x
      gbhl4 = two * gbhl4  * d3x   * bhfac / d1x
      gbhl5 = two * gbhl5  * d3x   * bhfac / d1x
      gbhl6 = two * gbhl6  * d3x   * bhfac / d1x
      gbhl7 = two * gbhl7  * d3x   * bhfac / d1x
      gbhl8 = two * gbhl8  * d3x   * bhfac / d1x
      gbhl9 = two * gbhl9  * d3x   * bhfac / d1x
      gbhl0 = two * gbhl0  * d3x   * bhfac / d1x
c     denpo = denpo / dble(nx*ny*nz)
      cyl1mas= two * cyl1mas * d3x
      cyl2mas= two * cyl2mas * d3x
      cyl3mas= two * cyl3mas * d3x
      con1mas= two * con1mas * d3x
      con2mas= two * con2mas * d3x
      con3mas= two * con3mas * d3x

c  check black hole condition

      if ( ( gbhl1 .gt. rbhl1 ) .or. ( gbhl2 .gt. rbhl2 ) .or.
     &     ( gbhl3 .gt. rbhl3 ) .or. ( gbhl4 .gt. rbhl4 ) .or.
     &     ( gbhl5 .gt. rbhl5 ) .or. ( gbhl6 .gt. rbhl6 ) .or.
     &     ( gbhl7 .gt. rbhl7 ) .or. ( gbhl8 .gt. rbhl8 ) .or.
     &     ( gbhl9 .gt. rbhl9 ) .or. ( gbhl0 .gt. rbhl0 )     ) then
        blackhole = .true.
        write(*,*) 'inout:  black hole!'
        write(*,'(a,10F13.8)') 'rad:  ', rbhl1, rbhl2, rbhl3, rbhl4,
     &       rbhl5, rbhl6, rbhl7, rbhl8, rbhl9, rbhl0
        write(*,'(a,10F13.8)') 'mass: ', gbhl1, gbhl2, gbhl3, gbhl4,
     &       gbhl5, gbhl6, gbhl7, gbhl8, gbhl9, gbhl0
      endif
c        stop

c  total and reduced masses
c caution:pmass2 is zero !!
      tgamas= gasmas + pmass1
      totmas= tgamas + pmass2
      redmac= pmass1 * gasmas / tgamas
      redmat= pmass2 * tgamas / totmas


c  CM position and velocity of
c       a) gas mass and black hole  b) everything

      cmcax = (cmgax*gasmas + topx1*pmass1) / tgamas
      cmcay = (cmgay*gasmas + topy1*pmass1) / tgamas
c b)
      cmx   = (cmcax*tgamas + topx2*pmass2) / totmas
      cmy   = (cmcay*tgamas + topy2*pmass2) / totmas
c      write(*,*)'pmass2,totmas :',pmass2,totmas
c      write(*,*)'cmx,cmcax :',cmx,cmcax
      vcmcx = (vcmgx*gasmas + velx1*pmass1) / tgamas                    !DRW - vcmcx becomes NaN here
      vcmcy = (vcmgy*gasmas + vely1*pmass1) / tgamas                    !DRW - vcmcy becomes NaN here
c b)
      vcmx  = (vcmcx*tgamas + velx2*pmass2) / totmas
      vcmy  = (vcmcy*tgamas + vely2*pmass2) / totmas


c  kinetic, thermal, total energies, relative to WD and CM;
c  angular momenta

      gco2 = -g * pmass2 / d1x

      gepot = zero
      emixpot=zero
      gekin = zero
      geintern= zero
      getherm = zero
      geges = zero
      lzgas = zero
      bmcm  = zero
      bmwd  = zero
c     getot = zero
      dkmass= zero
      dcmass= zero

      jstark = g * totmas / cc * sqrt(6.D0)
      jstarc = g * totmas / cc * 6.D0
      rs3sq  =(g * totmas / cc/cc * 6.D0)**2

      posx2 = zero
      posy2 = zero
      velx2 = zero
      vely2 = zero
      pmas  = zero
      rnssq = (20.D5 / d1x)**2
      rmtx = topos(dble(imaxx), lone, imaxg)
      rmty = topos(dble(imaxy), ltwo, imaxg)

c--gridloop start
      do 502 levl = 1, mlev                ! loop over all levels
      rkf  = one / ( dble(nfine)**(levl-1) )
      rkf3 = rkf**3
      topz2 = half*rkf
      topz1 = topz2
      rmtz  = half*rkf

      do 501 num = 1, indxgr(lzero,levl)       ! loop over all grids on levl
         lgrd = indxgr(num,levl)
         xopos = topos( zero, lone, lgrd)
         yopos = topos( zero, ltwo, lgrd)
c         zopos = topos( zero, lthree, igrd)

c        call addwdpo(lgrd)
c         call bhdenup(lgrd)
cp         pmin = zero


C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k),
C$OMP+            SHARED(nx,ny,nz,lgrd,ro,po,densty,chem)
         do 112 j = 1, ny
         do 112 k = 1, nz
         do 112 i = 1, nx
            ro(i,j,k) = densty(i,j,k,lgrd)
            po(i,j,k) = chem(i,j,k,lgrd,lone) / ro(i,j,k)
  112    continue
c         po actually is Ye !!
C$OMP END PARALLEL DO

c        zero area covered by finer grids
         ifgr = idxfgr(lgrd)
         if( ifgr .ne. -lone ) call rozero(ifgr)

c        zero area covered by WD on finest grid
c        if ( lgrd.eq.iwdgr ) then
c          do 113 k = 1, nz
c          do 113 j = 1, ny
c          do 113 i = 1, nx
c113         ro(i,j,k) = ro(i,j,k) * dble(inwd(i,j,k))
c        endif
c--gridloop start


C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,k,dsqy,edsqy,dysq,dymasq,
C$OMP+                   dsq,edsq,rds,erds,dxysq,dmasq,pot,vgx,vgy,ei,v,
C$OMP+                   vzsq,vsq,lzgtmp,ecm,ekcm,ewd,ekwd,dzmasq,rho),
C$OMP+            SHARED(nx,ny,nz,rkf,cmy,rmty,cmx,rmtx,rmtz,pmass1,
C$OMP+                   yopos,xopos,zopos,d1x,vcmx,vcmy,jstark,jstarc,
C$OMP+                   topy2,topy1,topx2,topx1,topz2,topz1,energy,
C$OMP+                   lgrd,gco2,rnssq,smlrho,ro,gpot,bhdens,rs3sq,
C$OMP+                   velx,vely,velz,rkf3,rsbh,g),
C$OMP+            REDUCTION(+:gepot,lzgas,dkmass,dcmass,geintern),
C$OMP+            REDUCTION(+:gekin,pmas,posx2,posy2,velx2,vely2),
C$OMP+            REDUCTION(-:emixpot)

      do j = 1, ny
         dsqy   =   (yopos+dble(j)*rkf-topy2)**2
         edsqy  =   (yopos+dble(j)*rkf-topy1)**2
         dysq   =   (yopos+dble(j)*rkf-cmy  )**2
         dymasq =   (yopos+dble(j)*rkf-rmty )**2
      do k = 1, nz
         dsq    =   (      dble(k)*rkf-topz2)**2 + dsqy
         edsq   =   (      dble(k)*rkf-topz1)**2 + edsqy
         dzmasq =   (      dble(k)*rkf-rmtz )**2 + dymasq
      do i = 1, nx
         rds = sqrt((xopos+dble(i)*rkf-topx2)**2 + dsq)
         erds= sqrt((xopos+dble(i)*rkf-topx1)**2 + edsq)
         dxysq =  ( (xopos+dble(i)*rkf-cmx  )**2 + dysq ) * d1x*d1x
         dmasq =    (xopos+dble(i)*rkf-rmtx )**2 + dzmasq

         rho(i) = ro  (i,j,k)*rkf3
         pot(i) = gpot(i,j,k,lgrd)
         vgx(i) = velx(i,j,k,lgrd)
         vgy(i) = vely(i,j,k,lgrd)

         vzsq(i)= velz(i,j,k,lgrd)**2
         vsq(i) = vgx(i)**2 + vgy(i)**2 + vzsq(i)

c  total potential energy of gas

         gepot = gepot + rho(i)*pot(i)
cp       pmin = min(pmin,pot(i))
         if (erds .gt. rsbh) then
            emixpot = emixpot - rho(i)*g*pmass1/erds/d1x
         endif
c  angular momentum of gas about the centre of mass
c  exclude BH area, accreted ang mom is calc in subroutine sucks
         if (bhdens(i,j,k,lgrd) .eq. zero) then
            lzgtmp = ( (xopos+dble(i)*rkf-cmx) * (vgy(i)-vcmy) -
     &            (yopos+dble(j)*rkf-cmy) * (vgx(i)-vcmx)   ) * d1x

            lzgas = lzgas + lzgtmp * rho(i)
         endif

         if ( dxysq .gt. rs3sq ) then
            if ( lzgtmp .gt. jstark ) dkmass = dkmass + rho(i)
            if ( lzgtmp .gt. jstarc ) dcmass = dcmass + rho(i)
         endif

c  kinetic energy of gas relative to CM and WD

         ei(i)= ( energy(i,j,k,lgrd) - half * vsq(i) )

         ekcm(i)=half*( (vgx(i)-vcmx )**2 + (vgy(i)-vcmy )**2 + vzsq(i))
         ekwd(i)=half*( (vgx(i)-velx2)**2 + (vgy(i)-vely2)**2 + vzsq(i))
         ecm(i) = ei(i) + ekcm(i) + pot(i)
         ewd(i) = ei(i) + ekwd(i) + gco2/rds

c  sum of thermal and total (kinetic relative to CM) energy

c        getot  = getot  + rho(i) * energy(i,j,k,lgrd)
         geintern = geintern + rho(i) * ei(i)
         gekin    = gekin    + rho(i) * ekcm(i)
cc       geges = geges + rho(i) * energy(i,j,k,lgrd)

         v(i) = sqrt(vsq(i)) * rho(i)/(rho(i)+smlrho)

         if ( dmasq .lt. rnssq ) then                                   !DRW - It seems that this condition is ~always~ false
c         write(*,*) 'adding pmas',pmas
            pmas  = pmas  + rho(i)
            posx2 = posx2 + rho(i) * (xopos+dble(i)*rkf)
            posy2 = posy2 + rho(i) * (yopos+dble(j)*rkf)
            velx2 = velx2 + rho(i) * vgx(i)
            vely2 = vely2 + rho(i) * vgy(i)
         endif

c DRW - The above condition is never met. So I will update these        !DRW.
c quantities anyway.

         pmas  = pmas  + rho(i)
         posx2 = posx2 + rho(i) * (xopos+dble(i)*rkf)
         posy2 = posy2 + rho(i) * (yopos+dble(j)*rkf)
         velx2 = velx2 + rho(i) * vgx(i)
         vely2 = vely2 + rho(i) * vgy(i)

c Back to MRR                                                           !!!!

      enddo
      enddo
      enddo

C$OMP END PARALLEL DO

      if ( lgrd .eq. itopgr ) then

      k = nz
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,j,rho),
C$OMP+            SHARED(nx,ny,ro,rkf3,lgrd,velz,energy,gpot,d1x,k,
C$OMP+                   xopos,yopos,rkf,cmx,cmy,velx,vely,vcmx,vcmy)
C$OMP+            REDUCTION(+:gmunb,angunb)
      do j = 1, ny
      do i = 1, nx
         rho(i) = ro(i,j,k) * rkf3

c         iye(i) = nint( max( one, min( dble(nye-1),
c     &    one+dble(nye-1)*(po(i,j,k)-yei(1)) / (yei(nye)-yei(1))))  )
c                          po actually is Ye !!

      if ( ( velz  (i,j,k,lgrd) .gt. zero) .and.
     &     ( energy(i,j,k,lgrd) ! - eco(iye(i))
     &            + gpot(i,j,k,lgrd).gt.zero) ) then
        gmunb = gmunb + rho(i)*velz(i,j,k,lgrd)
      endif
      if   ( velz  (i,j,k,lgrd) .gt. zero) then
        angunb=angunb + rho(i)*velz(i,j,k,lgrd) *
     &   ( (xopos+dble(i)*rkf-cmx) * (vely(i,j,k,lgrd)-vcmy) -
     &     (yopos+dble(j)*rkf-cmy) * (velx(i,j,k,lgrd)-vcmx)  ) * d1x
      endif

      enddo
      enddo

C$OMP END PARALLEL DO

      j = ny!-1
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,k,rho),
C$OMP+            SHARED(nx,nz,ro,rkf3,lgrd,vely,energy,gpot,d1x,j,
C$OMP+                   xopos,yopos,rkf,cmx,cmy,velx,vcmx,vcmy)
C$OMP+            REDUCTION(+:gmunb,angunb)
      do i = 1, nx
      do k = 1, nz
         rho(i) = ro(i,j,k) * rkf3

c         iye(i) = nint( max( one, min( dble(nye-1),
c     &    one+dble(nye-1)*(po(i,j,k)-yei(1)) / (yei(nye)-yei(1))))  )

      if ( ( vely  (i,j,k,lgrd) .gt. zero) .and.
     &     ( energy(i,j,k,lgrd) ! - eco(iye(i))
     &            + gpot(i,j,k,lgrd).gt.zero) ) then
        gmunb = gmunb + rho(i)*vely(i,j,k,lgrd)
      endif
      if   ( vely  (i,j,k,lgrd) .gt. zero) then
        angunb=angunb + rho(i)*vely(i,j,k,lgrd) *
     &   ( (xopos+dble(i)*rkf-cmx) * (vely(i,j,k,lgrd)-vcmy) -
     &     (yopos+dble(j)*rkf-cmy) * (velx(i,j,k,lgrd)-vcmx)  ) * d1x
      endif

      enddo
      enddo
C$OMP END PARALLEL DO

      j = 1
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(i,k,rho),
C$OMP+            SHARED(nx,nz,ro,rkf3,lgrd,vely,energy,gpot,d1x,j,
C$OMP+                   xopos,yopos,rkf,cmx,cmy,velx,vcmx,vcmy)
C$OMP+            REDUCTION(+:gmunb,angunb)
      do i = 1, nx
      do k = 1, nz
         rho(i) = ro(i,j,k) * rkf3

c         iye(i) = nint( max( one, min( dble(nye-1),
c     &    one+dble(nye-1)*(po(i,j,k)-yei(1)) / (yei(nye)-yei(1))))  )

      if ( ( vely  (i,j,k,lgrd) .lt. zero) .and.
     &     ( energy(i,j,k,lgrd) ! - eco(iye(i))
     &            + gpot(i,j,k,lgrd).gt.zero) ) then
        gmunb = gmunb + rho(i)*abs(vely(i,j,k,lgrd))
      endif
      if   ( vely  (i,j,k,lgrd) .lt. zero) then
        angunb=angunb + rho(i)*abs(vely(i,j,k,lgrd)) *
     &   ( (xopos+dble(i)*rkf-cmx) * (vely(i,j,k,lgrd)-vcmy) -
     &     (yopos+dble(j)*rkf-cmy) * (velx(i,j,k,lgrd)-vcmx)  ) * d1x
      endif

      enddo
      enddo
C$OMP END PARALLEL DO

      i = nx!-1
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,rho),
C$OMP+            SHARED(ny,nz,ro,rkf3,lgrd,vely,energy,gpot,d1x,i,
C$OMP+                   xopos,yopos,rkf,cmx,cmy,velx,vcmx,vcmy)
C$OMP+            REDUCTION(+:gmunb,angunb)
      do j = 1, ny
      do k = 1, nz
         rho(i) = ro(i,j,k) * rkf3

c         iye(i) = nint( max( one, min( dble(nye-1),
c     &    one+dble(nye-1)*(po(i,j,k)-yei(1)) / (yei(nye)-yei(1))))  )

      if ( ( velx  (i,j,k,lgrd) .gt. zero) .and.
     &     ( energy(i,j,k,lgrd) ! - eco(iye(i))
     &                + gpot(i,j,k,lgrd).gt.zero ) ) then
        gmunb = gmunb + rho(i)*velx(i,j,k,lgrd)
      endif
      if   ( velx  (i,j,k,lgrd) .gt. zero) then
        angunb=angunb + rho(i)*velx(i,j,k,lgrd) *
     &   ( (xopos+dble(i)*rkf-cmx) * (vely(i,j,k,lgrd)-vcmy) -
     &     (yopos+dble(j)*rkf-cmy) * (velx(i,j,k,lgrd)-vcmx)  ) * d1x
      endif

      enddo
      enddo
C$OMP END PARALLEL DO

      i = 1
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(j,k,rho),
C$OMP+            SHARED(ny,nz,ro,rkf3,lgrd,vely,energy,gpot,d1x,i,
C$OMP+                   xopos,yopos,rkf,cmx,cmy,velx,vcmx,vcmy),
C$OMP+            REDUCTION(+:gmunb,angunb)
      do j = 1, ny
      do k = 1, nz
         rho(i) = ro(i,j,k) * rkf3

c         iye(i) = nint( max( one, min( dble(nye-1),
c     &    one+dble(nye-1)*(po(i,j,k)-yei(1)) / (yei(nye)-yei(1))))  )

      if ( ( velx  (i,j,k,lgrd) .lt. zero) .and.
     &     ( energy(i,j,k,lgrd) ! - eco(iye(i))
     &                + gpot(i,j,k,lgrd).gt.zero ) ) then
        gmunb = gmunb + rho(i)*abs(velx(i,j,k,lgrd))
      endif
      if   ( velx  (i,j,k,lgrd) .lt. zero) then
        angunb=angunb + rho(i)*abs(velx(i,j,k,lgrd)) *
     &   ( (xopos+dble(i)*rkf-cmx) * (vely(i,j,k,lgrd)-vcmy) -
     &     (yopos+dble(j)*rkf-cmy) * (velx(i,j,k,lgrd)-vcmx)  ) * d1x
      endif

      enddo
      enddo
C$OMP END PARALLEL DO

      endif ! of if lgrd eq itopgr

c  find thermal energy

c      do i = 1, nx
c        rho(i) = max( ro(i,j,k), smlrho )
c        ye(i)  = po(i,j,k)
c        ei(i)  = rho(i) * ei(i)
c        tmp(i) = ten**tel(lone)
c      enddo

c      do i = 1, nx
c        getherm = getherm + ei(i) * rkf3
c      enddo

c      call eos(rho(1),tmp(1),ei(1),p(1),ye(1),gamc(1),game(1),nx,lone)
c *** watch out: eos changes also ei and ye where necessary !!!

c      do i = 1, nx
c        getherm = getherm - ei(i) * rkf3
c      enddo

c  how much mass is bound to WD and to CM?

c      do 205 i=1,nx
c          if( ecm(i) .lt. zero )      bmcm = bmcm + rho(i)
c          if( ewd(i) .lt. -6.67D14 ) bmwd = bmwd + rho(i)
c 205  continue


c  204 continue

  501 continue
  502 continue
c--gridloop end

cp      write(*,*) 'pmin  ',pmin
c      write(*,*) 'NS mass :' ,pmas/solmas
c     getot   = two * getot * d3x
      gepot   = two * gepot * d3x
      emixpot = two * emixpot *d3x
      gekin   = two * gekin * d3x
      getherm = two * getherm * d3x
      geintern= two * geintern* d3x
      geges   = gekin + geintern
      lzgas   = two * lzgas * d3x
c      bmcm    = two * bmcm  * d3x / gasmas
c      bmwd    = two * bmwd  * d3x / gasmas
      dkmass  = two * dkmass * d3x
      dcmass  = two * dcmass * d3x
c      write(*,*) 'Mixpot:',emixpot
c  position of NS relative to gas center of mass
c      bdx   = topos( dble(imaxx), lone, imaxg) - cmgax
c      bdy   = topos( dble(imaxy), ltwo, imaxg) - cmgay
c      angle = 180.D0 / pi * atan2(bdy,bdx)
c      didi  = bdx**2 + bdy**2
c      bidist= sqrt(didi) * d1x * two
c      velxns =      velx(imaxx,imaxy,imaxz,imaxg) - vcmgx
c      velyns =      vely(imaxx,imaxy,imaxz,imaxg) - vcmgy
c      velns = sqrt( velxns**2 + velyns**2 )

c      posx2 = cmgax
c      posy2 = cmgay
c      velx2 = vcmgx
c      vely2 = vcmgy

      posx2 = posx2 / pmas                                              !DRW - pmas = zero causes posx2, posy2 to be NaN
      posy2 = posy2 / pmas
      velx2 = velx2 / pmas
      vely2 = vely2 / pmas
      pmas = two * pmas * d3x
      bidist = sqrt( (posx1-posx2)**2 + (posy1-posy2)**2 ) * d1x
      bidiav = sqrt( (posx1-cmgax)**2 + (posy1-cmgay)**2 ) * d1x
c      write(*,*)'bidist in datout :', bidist
      velxns = velx2
      velyns = vely2
      velns = sqrt( velxns**2 + velyns**2 )

c  angular momenta: orbits, total
c --- ang. mom. orbit of BH, (spin angm1 in "sucks")
      vrex = velx1 - vcmx
      vrey = vely1 - vcmy

c      limc = redmac * d1x *
c     &             ( (topx1-cmx)*vrey - (topy1-cmy)*vrex )

      limc=pmass1*d1x*((topx1-cmx)*vrey-(topy1-cmy)*vrex )

c      limc = 0.d0
      vrex = velx2 - vcmx
      vrey = vely2 - vcmy
c      limp = redmat * d1x *
c     &             ( (topx2-cmx)*vrey - (topy2-cmy)*vrex )
      limp = 0.d0
c      ltot = limc + limp + lzgas + angm1 + angm2
c      write(*,*) 'angcc1:datout',angcc1
      ltot = limc + lzgas + angm1
c  energies: potential, kinetic, total
      ep12 = emixpot
c      ep12 = -g * pmass1 * gasmas / d1x /
c     &      sqrt( (topx1-topx2)**2 + (topy1-topy2)**2 + topso**2 )
c      write(*,*)'ep12 with NS pos',ep12
c      ep12 = -g * pmass1 * gasmas / d1x /
c     &      sqrt( (topx1-cmgax)**2 + (topy1-cmgay)**2 + topso**2 )
c      write(*,*)'ep12 with CM pos',ep12
c     ep12 = -g * pmass1 * pmass2 / d1x /
c    &      sqrt( (topx1-topx2)**2 + (topy1-topy2)**2 + topso**2 )
c     epot1  = pmass1 * potint( posx1, posy1, indxgr(1,maxlev(1)) )
c     epot2  = pmass2 * potint( posx2, posy2, indxgr(2,maxlev(2)) )
c     epot2  = pmass2 * gpot(16,16,1,9)
      epot1  = -25D0/28D0 * g * pmass1**2 / (two*two*g*pmass1/cc**2)
c?      epot1  = -21D0/22D0 * g * pmass1**2 / (two*two*g*pmass1/cc**2)
      epot2  = zero                         ! BH is 2Rs in radius!
      epotot = half*gepot + epot1 + epot2 + ep12
cp      epotot = half*gepot
c     write(*,*) 'datout   ', epot2/pmass2 , gpot(16,16,1,9)

      ekin1 = half * pmass1 * ( (velx1 - vcmx)**2 + (vely1 - vcmy)**2 )
      ekin2 = half * pmass2 * ( (velx2 - vcmx)**2 + (vely2 - vcmy)**2 )

      tmpmaw  = log10(tmpmax)
      rhmaxw  = log10(rhmax)
c     denpo  = log10(densty( nint(topx2), nint(topy2), 1, itopgr ) )
      denpo  = log10(densty( nint(cmgax), nint(cmgay), 1, itopgr ) )
c     denpo  = log10
c    &    (densty( nint(posx2), nint(posy2), lone, iwdgr ) )
c     tempow  = log10(tempo)

c      tempow = one
c      itm = lzero
c      jtm = lzero
c      ktm = lzero
c      itmgr = lzero

      dlzgrav =  -0.4D0 * g / cc**5  *  (
     &         dij2dt(1,1) * dij3dt(1,2)  -  dij2dt(1,2) * dij3dt(1,1)
     &       + dij2dt(1,2) * dij3dt(2,2)  -  dij2dt(2,2) * dij3dt(1,2)
     &       + dij2dt(1,3) * dij3dt(2,3)  -  dij2dt(2,3) * dij3dt(1,3) )
      tlzgrav = tlzgrav + dlzgrav * dtgr(itopgr)*two   ! two timesteps

      eglum = g / 5.D0 / cc**5 *
     &   (   dij3dt(1,1)**2 + dij3dt(1,2)**2 + dij3dt(1,3)**2
     &     + dij3dt(2,1)**2 + dij3dt(2,2)**2 + dij3dt(2,3)**2
     &     + dij3dt(3,1)**2 + dij3dt(3,2)**2 + dij3dt(3,3)**2  )

      gcc = g / cc**4
      do j = 1, 3
      do i = 1, 3
         dij2dt(i,j) = gcc * dij2dt(i,j)
      enddo
      enddo
      hplus = dij2dt(1,1) - dij2dt(2,2)
      hcross= two * dij2dt(1,2)
      egpr =  2.D3 * gmunb *d1x*d1x * dtgr(itopgr) !factor 1.E3 for printout
      egprl=  two  * angunb*d1x*d1x * dtgr(itopgr)                      !DRW - egprl becomes NaN here because angunb is NaN
      ergtot= geges + ekin1 + epotot + enuloss + egquel + eint1
c     getot = getot + epotot + ekin1 + ekin2

      visfac = one/gasmas/two/dtgr(itopgr)  ! two timesteps are added in
                                            ! subr. entropy

c  print out values                                                     !12345

      write (6,354) tms, topx1, topy1, velx1, vely1,
     &         bidist, eglum/1.D50,
     &         getherm/1.D50, egquel/1.D50,
     &         rhmaxw, imaxx, imaxy, imaxz, imaxg, egpr/solmas
  354 format ( 1F7.4, 2F11.6, 1P, 3E11.3, 3E12.4,
     &             0P,1(F7.3, 2I4,1I3, 1I2), 1F8.5 )
      write (6,355) nstep, cmcax, cmcay, vcmcx, vcmcy,
     &         lzgas/1.D50, gekin/1.D50, geintern/1.D50, epotot/1.D50,
     &         ep12/1.D50,  epot1/1.D50
  355 format (I6,' ', 2F10.5, 1P, 8E12.4 )
      write (6,356) velxns, velyns, hplus, hcross, ergboun/1.D50,
     &              gasmas/solmas,pmas/solmas,
     &              enuloss/1.D50, eneloss/1.D50, enaloss/1.D50,
     &              enxloss/1.D50, ergtot/1.D50
  356 format (7X, 1P, 5E11.3, 0P, 2F8.4, 1P, 4E12.4, 1E13.5 )
      write (6,357) dij2dt(1,1), dij2dt(2,2), dij2dt(3,3),
     &              dij2dt(1,2), d1mass/solmas, d2mass/solmas,
     &              d3mass/solmas, rnetot/amu/1.D50, ltot/1.D50
  357 format (7X, 1P, 4E14.6, 3E12.4, 2E13.5 )
      write (6,358) enelo/(anelo+1.D-30), enalo/(analo+1.D-30),
     &              enxlo/(anxlo+1.D-30), tlzgrav/1.D50,
     &              sumbulk*visfac,  sumshear*visfac,
     &              sumneutr*visfac, sumshock*visfac,
     &              dkmass/solmas, dcmass/solmas, pmass1/solmas
  358 format (7X, 0P, 3F8.3, 1P, 7E13.4, 0P, 1F8.4)
      write (6,359) topx2, topy2,
     &              cyl1mas/solmas, cyl2mas/solmas, cyl3mas/solmas,
     &              con1mas/solmas, con2mas/solmas, con3mas/solmas,
     &              angm1/1.D50, eint1/1.D50, limc/1.D50
  359 format (7X, 0P, 2F11.6, 1P, 6E11.3, 3E13.5)
      write (6,361)
     &            dij3dt(1,1)/1.D50,dij3dt(1,2)/1.D50,dij3dt(1,3)/1.D50,
     &            dij3dt(2,2)/1.D50,dij3dt(2,3)/1.D50,dij3dt(3,1)/1.D50
  361 format (7X, 1P, 6E14.6)
      write (6,360) vlmax, imavx, imavy, imavz, imavg, tmpmaw,
     &              itm, jtm, ktm, itmgr, ekin1/1.D50,
     &              angcc1/1.D50, egprl/1.D50, bidiav
  360 format (7X, 1P, 1E10.2, 0P, 2I4, 1I3, 1I2,
     &               1F7.3, 2I4, 1I3, 1I2, 1P, 4E12.4, 1E11.3 )
      write(*,*) ' '

c DRW - Going to use Python to find the NaNs ===========================

      open(unit=420,file='names.txt')

      write(420,*) 'ms ','posx1 ','posy1 ','velx1 ','vely1 ','bidist '
      write(420,*) 'eglum ','getherm ','egquel ','rhmax ','mxx ','mxy '
      write(420,*) 'mxz ','g ','gmunb ','nstep ','cmcax ','cmcay '
      write(420,*) 'vcmcx ','vcmcy ','lzgas ','gekin ','geintern '
      write(420,*) 'epotot ','ep12 ','epot1 ','velxns ','velyns '
      write(420,*) 'hplus ','hcross ','ergboun ','gasmas ','pmas '
      write(420,*) 'enuloss ','eneloss ','enaloss ','enxloss ','ergtot '
      write(420,*) 'd2dtxx ','d2dtyy ','d2dtzz ','d2dtxy ','d1mass '
      write(420,*) 'd2mass ','d3mass ','netot ','ltot ','avene '
      write(420,*) 'avena ','avenx ','tlzgrav ','sumbulk ','sumshear '
      write(420,*) 'sumneutr ','sumshock ','dkmass ','dcmass ','pmass1 '
      write(420,*) 'posx2 ','posy2 ','cyl1mas ','cyl2mas ','cyl3mas '
      write(420,*) 'con1mas ','con2mas ','con3mas ','angm1 ','eint1 '
      write(420,*) 'limc ','dij3dt1 ','dij3dt2 ','dij3dt3 ','dij3dt5 '
      write(420,*) 'dij3dt6 ','dij3dt9 ','vlmax ','mvx ','mvy ','mvz '
      write(420,*) 'mg ','tmpmax ','itm ','jtm ','ktm ','itmgr '
      write(420,*) 'ekin1 ','angcc1 ','egprl ','bidiav'

      close(420)

      open(unit=421,file='values.txt')

      write(421,*) tms,topx1,topy1,velx1,vely1,bidist,eglum/1.D50
      write(421,*) getherm/1.D50,egquel/1.D50,rhmaxw,imaxx,imaxy,imaxz
      write(421,*) imaxg,egpr/solmas,nstep,cmcax,cmcay,vcmcx,vcmcy
      write(421,*) lzgas/1.D50,gekin/1.D50,geintern/1.D50,epotot/1.D50
      write(421,*) ep12/1.D50,epot1/1.D50,velxns,velyns,hplus,hcross
      write(421,*) ergboun/1.D50,gasmas/solmas,pmas/solmas,enuloss/1.D50
      write(421,*) eneloss/1.D50,enaloss/1.D50,enxloss/1.D50
      write(421,*) ergtot/1.D50,dij2dt(1,1),dij2dt(2,2),dij2dt(3,3)
      write(421,*) dij2dt(1,2),d1mass/solmas,d2mass/solmas,d3mass/solmas
      write(421,*) rnetot/amu/1.D50,ltot/1.D50,enelo/(anelo+1.D-30)
      write(421,*) enalo/(analo+1.D-30),enxlo/(anxlo+1.D-30)
      write(421,*) tlzgrav/1.D50,sumbulk*visfac,sumshear*visfac
      write(421,*) sumneutr*visfac,sumshock*visfac,dkmass/solmas
      write(421,*) dcmass/solmas,pmass1/solmas,topx2,topy2
      write(421,*) cyl1mas/solmas,cyl2mas/solmas,cyl3mas/solmas
      write(421,*) con1mas/solmas,con2mas/solmas,con3mas/solmas
      write(421,*) angm1/1.D50,eint1/1.D50,limc/1.D50,dij3dt(1,1)/1.D50
      write(421,*) dij3dt(1,2)/1.D50,dij3dt(1,3)/1.D50,dij3dt(2,2)/1.D50
      write(421,*) dij3dt(2,3)/1.D50,dij3dt(3,1)/1.D50,vlmax,imavx,imavy
      write(421,*) imavz,imavg,tmpmaw,itm,jtm,ktm,itmgr,ekin1/1.D50
      write(421,*) angcc1/1.D50,egprl/1.D50,bidiav

      close(421)

c=======================================================================

      nout1 = lzero
      tout1 = zero

      enelo = zero
      enalo = zero
      enxlo = zero
      anelo = zero
      analo = zero
      anxlo = zero

c      stop
      return
      end
c-----------------------------------------------------------------------
      subroutine opeplt
c     open files pertinent to plots

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include   'rays.cmn'
      character*7 denfil

      !print *, 'inout, opeplt'

      outfil = basenm // 'O' // suffix
      rstfil = basenm // 'R' // suffix
      vidfil = basenm // 'V' // suffix
c      rgbcha = basenm // 'B' // suffix
      denfil = basenm // 'C' // suffix
      open (44, file=denfil, form='unformatted')
      open(8,file=outfil,form='unformatted')

      return
      end
c-----------------------------------------------------------------------
      subroutine cloplt
c     close files pertinent to plots

      include 'qparam.cmn'
      include 'squants.cmn'

      include   'rays.cmn'

      !print *, 'inout, cloplt'

c      close(3)
      close(8)
      close(44)
c      close(66)

      return
      end
c-----------------------------------------------------------------------
      subroutine cutout(igr)
c     output planes for 2D-film

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'

      real*4         outden(qx,qy), outtem(qx,qy), outye (qx,qy),
     &               outvx (qx,qy), outvy (qx,qy), outter(qx,qy),
     &               outvz (qx,qy), outerg(qx,qy), outent(qx,qy),
     &               outdeh(qx,qz), outteh(qx,qz), outyeh(qx,qz),
     &               outtrh(qx,qz), outenh(qx,qz),
     &               outtime,  outposx2, outposy2
      integer*4      outigrd,  outnorgin(ngrd,3)

      common/outdat/ outden       , outtem       , outye        ,
     &               outvx        , outvy        , outter       ,
     &               outvz        , outerg       , outent       ,
     &               outdeh       , outteh       , outyeh       ,
     &               outtrh       , outenh


      !print *, 'inout, cutout'

      do igrd = 1, ngrd
      do i = 1, 3
         outnorgin(igrd,i) = norgin(igrd,i)
      enddo
      enddo


      do 200 igrd = igr, ngrd

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(iy,ix),
C$OMP+            SHARED(nx,ny,densty,temper,tmpent,chem,igrd,
C$OMP+                   outden,outtem,outter,outye,outent)
      do iy=1,ny
      do ix=1,nx
         outden(ix,iy) =
     &          max( 1.D-30, min( 1.D+30, densty(ix,iy,1,igrd) ) )
         outtem(ix,iy) =
     &          max( 1.D-30, min( 1.D+30, temper(ix,iy,1,igrd) ) )
         outter(ix,iy) =
     &          max( 1.D-30, min( 1.D+30, tmpent(ix,iy,1,igrd) ) )
         outye (ix,iy) = max( 1.D-30, min( 1.D+30,
     &         chem(ix,iy,1,igrd,1) / densty(ix,iy,1,igrd) ))
         outent(ix,iy) = max( 1.D-30, min( 1.D+30,
     &         chem(ix,iy,1,igrd,2) / densty(ix,iy,1,igrd) ))
c        outerg(ix,iy) = energy(ix,iy,1,igrd)
c        outvx (ix,iy) = velx  (ix,iy,1,igrd)
c        outvy (ix,iy) = vely  (ix,iy,1,igrd)
c        outvz (ix,iy) = velz  (ix,iy,1,igrd)
      enddo
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(iy,ix),
C$OMP+            SHARED(nx,ny,bhdens,igrd,outden,outter,outye)
      do iy=1,ny
      do ix=1,nx
         if (bhdens(ix,iy,1,igrd) .ne. zero) then
           outden(ix,iy) = 1.D16
           outter(ix,iy) = one
c           outent(ix,iy) = two
           outye (ix,iy) = one
         endif
      enddo
      enddo
C$OMP END PARALLEL DO

c         outye = zero
c
c      do iy=1+1,ny-1
c      do ix=1+1,nx-1
c         outye (ix,iy) = velx(ix+1,iy,1,igrd)-velx(ix-1,iy,1,igrd)
c     &                  +vely(ix,iy+1,1,igrd)-vely(ix,iy-1,1,igrd)
c      enddo
c      enddo

      outtime  = time
      outposx2 = posx2
      outposy2 = posy2
      outigrd  = igrd

      write(8)  outtime, outposx2, outposy2, outigrd, outnorgin
      write(8)  outden
      write(8)  outtem
      write(8)  outter
      write(8)  outye
      write(8)  outent

      nyh = ny/ltwo

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(iz,ix),
C$OMP+            SHARED(nx,nz,densty,temper,tmpent,chem,igrd,nyh,
C$OMP+                   outdeh,outteh,outtrh,outyeh,outenh)
      do 220 ix=1,nx
      do 220 iz=1,nz
         outdeh(ix,iz) =
     &          max( 1.D-30, min( 1.D+30, densty(ix,nyh,iz,igrd) ) )
         outteh(ix,iz) =
     &          max( 1.D-30, min( 1.D+30, temper(ix,nyh,iz,igrd) ) )
         outtrh(ix,iz) =
     &          max( 1.D-30, min( 1.D+30, tmpent(ix,nyh,iz,igrd) ) )
         outyeh(ix,iz) = max( 1.D-30, min( 1.D+30,
     &         chem(ix,nyh,iz,igrd,1) / densty(ix,nyh,iz,igrd) ))
         outenh(ix,iz) = max( 1.D-30, min( 1.D+30,
     &         chem(ix,nyh,iz,igrd,2) / densty(ix,nyh,iz,igrd) ))
 220     continue
C$OMP END PARALLEL DO

      write(8)  outdeh
      write(8)  outteh
      write(8)  outtrh
      write(8)  outyeh
      write(8)  outenh

 200  continue


      return
      end


c-----------------------------------------------------------------------
      subroutine filnam

c     construct new filenames for output, plot, and restart files

      include 'qparam.cmn'
      include 'squants.cmn'

      character*1 sf1,sf2,char

      !print *, 'inout, filnam'

      if (suffix(2:2) .eq. 'z' .or. suffix(2:2) .eq. 'Z') then
         sf1 = suffix(1:1)
         isf1 = ichar(sf1)
         sf2 = suffix(2:2)
         isf2 = ichar(sf2)

         isf1 = isf1 + 1
         isf2 = isf2 - 25
         suffix(1:1) = char(isf1)
         suffix(2:2) = char(isf2)
      else
         sf2 = suffix(2:2)
         isf2 = ichar(sf2)
         isf2 = isf2 + 1
         suffix(2:2) = char(isf2)
      endif

      rstfil = basenm // 'R' // suffix
      outfil = basenm // 'O' // suffix
      vidfil = basenm // 'V' // suffix
      rgbcha = basenm // 'B' // suffix
      neufil = basenm // 'N' // suffix

      return
      end
c-----------------------------------------------------------------------
      logical function stoprun()

      include 'qparam.cmn'
      include 'squants.cmn'

      integer*4 in

      open(1,file='./'//basenm//'_cont',form='formatted',err=100)
      read(1,*,err=100,end=100) in
      close(1)

      if ( in .ne. lzero ) then
         stoprun = .false.
         return
      endif

  100 continue

      write(*,*) ' stop this run! '
      stoprun = .true.

      return
      end

c----------------------------------------------------------------------
      subroutine readeos

c Read tables for Shen EoS from file

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'eos.cmn'
      include 'ppdisk.cmn'
      include 'compct.cmn'

      integer*4 i,j,k,m
      integer*4 ierr
      real*8 eave(nroent)
      real*4 rolin(nro), telin(nte), yeiin(nye),
     &       prelin(nro,nte,nye), erglin(nro,nte,nye),
     &       gmcin (nro,nte,nye),  tcnin(nro,nte,nye),
     &       tcpin (nro,nte,nye)

      integer*4 idum1, idum2, idum3
      real*4 rolein(nroent), telein(nteent), yeiein(nyeent),
     &       entlin(nroent,nteent,nyeent)

      real*4 mup(nro,nte,nye), mun(nro,nte,nye), dummy(nro,nte,nye)

c Force the following variables to be stored on the heap as opposed to
c the stack so that we don't overflow the stack when running openMP

c      save rolin, telin, yeiin, prelin, erglin, gmcin, tcnin, tcpin
c      save mup, mun, dummy

      common /tmp2/ rolin, telin, yeiin, prelin, erglin, gmcin, tcnin,
     &              tcpin, mup, mun, dummy
C$OMP THREADPRIVATE (/tmp2/)

c ----- Read Shen EoS data file -----

      open(45, file=homedata//'shen_eos', form='unformatted', err=1000)
      read(45, err=1001) idum1, idum2, idum3, rolin, telin, yeiin,
     &     prelin, erglin, dummy, mun, mup, dummy, gmcin, entlin
      close(45)

c Throw an error and quit if the header doesn't match with the
c array sizes in the hydro code

      !print *, 'inout, readeos'

      if(idum1.ne.nro .or. idum2.ne.nte .or. idum3.ne.nye) then
         write(*,*) ' '
         write(*,*) 'ERROR: array sizes did not match with eos header'
         write(*,*) 'nro, nte, nye are:'
         write(*,*) nro, nte, nye
         write(*,*) 'header of "shen_eos" is:'
         write(*,*) idum1, idum2, idum3
         stop
      end if

c Convert tabulated values to double precision (and convert to the
c correct cgs units if necessary (1 erg = 1.60217733e6 MeV)

      do i = 1, nro
         rol(i) = dble(rolin(i))
         rolent(i) = rol(i)
      enddo

      do i = 1, nte
         tel(i) = dble(telin(i))
         telent(i) = tel(i)
      enddo

      do i = 1, nye
         yei(i) = dble(yeiin(i))
         yeient(i) = yei(i)
      enddo

      do k = 1, nye
         do j = 1, nte
            do i = 1, nro
               prel(i,j,k) = dble(prelin(i,j,k))+dlog10(1.60217733d-6)
               ergl(i,j,k) = dble(erglin(i,j,k))+dlog10(1.60217733d-6)
               gmc (i,j,k) = dble(gmcin (i,j,k))
               entl(i,j,k) = dble(log10(entlin(i,j,k)))
            enddo
         enddo
      enddo

c Convert chemical potentials to degeneracy parameters

      do k = 1, nye
         do j = 1, nte
            do i = 1, nro
               tcp(i,j,k) = dble(mup(i,j,k)+1.29) / (10.d0**tel(j))
               tcn(i,j,k) = dble(mun(i,j,k)+1.29) / (10.d0**tel(j))
            end do
         end do
      end do

c Output to screen the extrema of the indices for the Eos

      write(*,*) '===================================================='
      write(*,'(A,2F9.4, / A,2F9.4, / A,2F9.5)')
     &     'readeos  log min/max of density:', rol(1),rol(nro),
     &     '   log min/max of temperature:  ',tel(1),tel(nte),
     &     '   min/max of ye:               ',yei(1),yei(nye)
      write(*,*) '===================================================='

c Compute specific energies

      do i = 1, nye
         eco(i) = ten**ergl(1,1,i) / ten**rol(1)
      enddo
c-------------------------------------------------

c Define minimum values of thermodynamic quantities

      !smlrho = ten**rol(1)                                             !DRW - Commented out
      smlrho = 1D-25
      !smallp = ten**prel(1,1,1)                                        !DRW commented out
      smallp = 1D-17
c NB! smalle is internal specific energy
      !smalle = (ten**ergl(1,1,1))/smlrho                               !DRW - Commented out
      smalle = 1D-10
      smltem = 1D-18
      smltme = ten**tel(1)

      write(*,*) ' '
      write(*,*) '***** small values *****'
      write(*,'(a)') ' smlrho      smallp      smalle      smltme'
      write(*,'(1p,4e12.5)') smlrho, smallp, smalle, smltme
      write(*,*) ' '
      print *, '***** Disk and star properties *****'
      print *, 'Inner radius, outer radius:', r_in/au, r_out/au
      print *, 'Star mass, disk mass:', pmass1/solmas, pmass2/solmas
      print *, ''

      return

c File I/O errors *****************************************************

 1000 write(*,*) "Couldn't find shen_eos"
      stop
 1001 write(*,*) "Error reading shen_eos"
      stop

      end

c-----------------------------------------------------------------------
      subroutine justrestrt(ird)


      include 'qparam.cmn'
      include 'squants.cmn'
      include 'aquants.cmn'
      include 'compct.cmn'
      include 'grnest.cmn'
      include   'rays.cmn'
      real*4  dty(qx,qy)
      real*8  yeout( qx, qy, qz, ngrd )
      character*7 testfil

      !print *, 'inout, justrestrt'



         open(2,file=rstfil,form='unformatted')
         read(2)   ! dummy read to skip first record
         read(2)   densty
         read(2)   velx
         read(2)   vely
         read(2)   velz
         read(2)   energy
         read(2)   chem
         read(2)   temper
         read(2)   tmpent
         read(2)   time, deltax,
     &             bndmnx, bndmxx, bndmny, bndmxy, bndmnz, bndmxz,
     &             rhoin, uin, utin, uttin, pin, ein, gin, tin,
     &             gamein, gamcin, chein,
     &             nstep, igeom, dto, cenx, ceny, cenz, iframe,
     &             pmass1, posx1, posy1, velx1, vely1, vlox1, vloy1,
     &             pmass2, posx2, posy2, velx2, vely2, vlox2, vloy2,
     &             enuloss, eneloss, enaloss, enxloss, bhdfac, angm1,
     &             indxgr, idxcgr, idxfgr, levlgr, norgin, angm2, eint1,
     &             ishedul, ifinup, nodenr, delx, dtgr, maxlev, dtgg,
     &             cmgax, cmgay, signif, dij3dt, egquel, gmunb, angunb,
     &             angcc1,  thestar, phistar, iendf, tlzgrav
         close(2)


         write(*,'(a,a,a,i5,a,f8.3)') 'restarting from file ',rstfil,
     &      '   nstep=',nstep,'   time=',time*1.D3
         write(*,*)

         write(*,*)    pmass1, posx1, posy1, velx1, vely1, angm1
c        write(*,*)    pmass2, posx2, posy2, velx2, vely2, angm2
         stop

         itopgr=  indxgr(1,1)
         gasmas = zero
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(k,j,i),
C$OMP+            SHARED(nx,ny,nz,densty,itopgr), REDUCTION(+:gasmas)
         do 244 j=1,ny
         do 244 k=1,nz
         do 244 i=1,nx
            gasmas = gasmas + densty(i,j,k,itopgr)
  244    continue
C$OMP END PARALLEL DO
         gasmas = two*gasmas * delx(itopgr)**3
         write(*,*) 'restrt, gasmas :', gasmas/solmas
c --- first without fineup
         do igrd = ngrd, 1, -1
            call sucks(igrd,lzero)
         enddo
c     now with fineup
         do igrd = ngrd, 1, -1
            call sucks(igrd,lone)
         enddo


         do 11 igrd = 1, ngrd
         do 11 j = 1, ny
         do 11 k = 1, nz
         do 11 i = 1, nx
   11       yeout(i,j,k,igrd)=chem(i,j,k,igrd,lone)/densty(i,j,k,igrd)


         open(2,file='outdata',form='unformatted')
         write(2)  qx, qy, qz
         write(2)  densty
         write(2)  velx
         write(2)  vely
         write(2)  velz
         write(2)  energy
         write(2)  yeout
         write(2)  temper
         write(2)  tmpent
         close(2)




      return
      end
