      subroutine eos (rho, tem, e, p, yee, gamc, game, nos, input)
c
c     E Q U A T I O N   O F   S T A T E
c     Adiabatic, ideal gas eos
c     Danny Rhys Walker
c
c     The state variables are calculated at nos grid points
c
c     rho   :    Density
c     tem   :    Temperature
c     e     :    Energy density
c     p     :    Pressure
c     gamc  :    d(ln p) / d(ln rho)  at fixed entropy
c
c     nos   :    Actual number of grid points
c
c     State variables are in   cgs-units .
c
c     input :   1    <==>  input  :  RHO
c             other  <==>  input  :  P
c
c     output :  e, p, gamc, game   if  input = 1
c     output :  t, p, gamc, game   if  input = 0
c
c----------------------------------------------------------------------

      include 'qparam.cmn'
      include 'ppdisk.cmn'
      include 'squants.cmn'

      real*8  rho (nos), tem (nos), e(nos), p(nos), yee(nos),
     &        gamc(nos), game(nos)

      integer input, nos

      ! The quantities used are rho, te, e and p. The others are relics.

c =============================

      if (input .eq. lone) then
c   D e n s i t y,  (y e)  and  t e m p e r a t u r e  given
c   --------------------------------------------------------
        do i = 1, nos
          if(rho(i).eq.zero) then ! Make sure p is not zero! rho(i) should never be zero though.
            print *, 'Zero rho(i), eos 1', i
            p(i) = smallp
          else
            p(i) = specr * rho(i) * tem(i) !Ideal gas law
          endif
          e(i)    =  p(i) / (gamma - one)
          gamc(i) =  gamma
          game(i) =  gamma
        enddo
c =============================
      else    ! D e n s i t y,  (y e)  and  e n e r g y  given
c             ------------------------------------------------
        do i = 1, nos
          if(rho(i).eq.zero) then ! Make sure p is not zero.
            print *, 'Zero rho(i), eos 2', i
            p(i) = smallp
          else
            p(i) = e(i) * (gamma - one)
          endif
          tem(i)  =  p(i) / (specr * rho(i))
          if( tem(i).gt.1000 ) then
            print *, 'eos big tmp'
            print *, i, rho(i), p(i), smlrho
            stop
          endif
          gamc(i) =  gamma
          game(i) =  gamma
        enddo
c ==============================
      endif
c ==============================
      return
      end
c----------------------------------------------------------------------

      subroutine eosent (rho, tem, yee, ent, nos, mode)
C PARALLEL OK
c       Entropy in    E q u a t i o n   o f   s t a t e
c
c      Dummy!
c
c      The state variables are calculated at nos grid points
c
c      rho   :    Density           (g / cm^3)
c      tem   :    Temperature       (MeV)
c      yee   :    Electron Fraction
c      entr  :    Entropy           (kB per nucleon)
c
c      nos   :    Actual number of grid points
c
c      mode :   1  <==>  input :  T,   ye ,rho        output: entropy
c               0  <==>  input :  entropy, ye, rho    output: T
c
c----------------------------------------------------------------------

      include 'qparam.cmn'
      include 'squants.cmn'
      include 'eos.cmn'

      integer*4 mode, nos
      real*8  rho(nos),  tem(nos), yee(nos), ent(nos)

      integer*4 qxy
      parameter( qxy = q )

      real*8 rho0(qxy), tem0(qxy), ye0(qxy),
     &       ebl(qxy), ebh(qxy), elg(qxy), ebm(qxy),
     &       wr1(qxy), wr0(qxy), we1(qxy), we0(qxy),
     &       wy1(qxy), wy0(qxy), wt1(qxy), wt0(qxy)

      integer*4 ittl(qxy), itth(qxy), irho0(qxy), irho1(qxy),
     &        iye0(qxy), iye1(qxy), item0(qxy), item1(qxy),
     &        ittm(qxy), isect(qxy), sumisec

      common/tmp5/   rho0, tem0, ye0, ebl, ebh, elg, ebm,
     &       wr1, wr0, we1, we0, wy1, wy0, wt1, wt0,
     &       irho0, irho1, iye0, iye1, item0, item1
     &       ittl, itth, ittm, isect, sumisec
C$OMP THREADPRIVATE (/tmp5/)

      if ( nos .gt. qxy ) stop ' eosent'

c----------------------------------------------------------------------
      if (mode .eq. lone)     then
c-----  t e m p e r a t u r e   given
c----------------------------------------------------------------------

      do i = 1, nos
         ent(i) = one  ! dummy
      enddo

c----------------------------------------------------------------------
      else   ! --   e n t r o p y        given
c----------------------------------------------------------------------

         stop'eos'

c   -------------------------
      endif
c   -------------------------
      return
      end
c----------------------------------------------------------------------

c=======================================================================
c                          O P A C I T Y                               =
c=======================================================================

      subroutine eosopac(igrd,opacity)
c
c     O P A C I T Y   E Q U A T I O N   O F   S T A T E
c     Based on Bell & Lin (1994)
c     Danny Rhys Walker
c
c     This subroutine calculates the opacity at every grid point on the
c     current grid.
c
c     Requirements: Parameters kap1,2... temex1,2... rhoex4... as well
c                   as the first two transition temperatures ttran1,
c                   ttran2. We get these from ppdisk.cmn.
c
c     Input: igrd.
c
c     Output: opacity.
c
c-----------------------------------------------------------------------

      include 'qparam.cmn'  ! qx, qy, ngrd
      include 'ppdisk.cmn'  ! Opacity parameters
      include 'aquants.cmn' ! densty, temper
      include 'squants.cmn' ! Grid parameters

      real*8 rho, tem
      real*8 opacity(qx,qy,qz,ngrd)

      do i = 1,nx
        do j = 1,ny
          do k = 1,nz
            tem = temper(i,j,k,igrd)
            rho = densty(i,j,k,igrd)

            ! The first two transition temperatures are fixed
            if(tem .lt. ttran1) then
              ! rhoex1 is 0
              opacity(i,j,k,igrd) = kap1*(tem**temex1)

            elseif(tem .lt. ttran2) then
              ! rhoex2 is 0
              opacity(i,j,k,igrd) = kap2*(tem**temex2)

            else
              ! Transition temperatures from this point depend on density
              ! and need to be calculated
              ttran3 = ( (kap3/kap4) * rho**(-rhoex4) )
     &                  **(1/(temex4-temex3))
              ttran4 = ( (kap4/kap5) * rho**(rhoex4-rhoex5) )
     &                  **(1/(temex5-temex4))
              ttran5 = ( (kap5/kap6) * rho**(rhoex5-rhoex6) )
     &                  **(1/(temex6-temex5))
              ttran6 = ( (kap6/kap7) * rho**(rhoex6-rhoex7) )
     &                  **(1/(temex7-temex6))
              ttran7 = ( (kap7/kap8) * rho**(rhoex7-rhoex8) )
     &                  **(1/(temex8-temex7))

              if(tem .lt. ttran3) then
                ! rhoex3 is 0
                opacity(i,j,k,igrd) = kap3*(tem**temex3)

              elseif(tem .lt. ttran4) then
                opacity(i,j,k,igrd) = kap4*(rho**rhoex4)*(tem**temex4)

              elseif(tem .lt. ttran5) then
                opacity(i,j,k,igrd) = kap5*(rho**rhoex5)*(tem**temex5)

              elseif(tem .lt. ttran6) then
                opacity(i,j,k,igrd) = kap6*(rho**rhoex6)*(tem**temex6)

              elseif(tem .lt. ttran7) then
                opacity(i,j,k,igrd) = kap7*(rho**rhoex7)*(tem**temex7)

              else
                ! temex8 and rhoex8 are 0
                opacity(i,j,k,igrd) = kap8

              endif
            endif
          enddo
        enddo
      enddo

      return
      end
