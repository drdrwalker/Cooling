      subroutine eosiso (rho, tem, e, p, yee, gamc, game, nos, input)

c
c               E Q U A T I O N   O F   S T A T E
c               Isothermal protoplanetary disk
c               AUTHOR: Danny Rhys Walker
c
c      The state variables are calculated at nos grid points
c
c      rho   :    Density
c      tem   :    Temperature
c      e     :    Energy density
c      p     :    Pressure
c      gamc  :    d(ln p) / d(ln rho)  at fixed entropy
c
c      nos   :    Actual number of grid points
c
c      State variables are in   cgs-units .
c
c      input :   1    <==>  input  :  RHO
c              other  <==>  input  :  P
c
c          output :  e,   p, gamc, game, tem   if  input = 1
c          output :  rho, e, gamc, game, tem   if  input not 1
c
c----------------------------------------------------------------------

      include 'qparam.cmn'
      include 'squants.cmn'

      real*8  rho (nos), tem (nos), e(nos), p(nos), yee(nos),
     &        gamc(nos), game(nos)

      real*8 boltz,         !Boltzmann constant, cgs
     &       diskT,         !Disk temp, Kelvin
     &       muMol,         !Mean molecular weight, cgs (I THINK)
     &       mHyd,          !Mass of hydrogen atom, cgs
     &       cSound2,       !Sound speed squared

      cSound2 = (boltz*diskT) / (muMol * mHyd)

      integer input, nos

c   ---------------------------
      if (input .eq. 1) then
c                 D e n s i t y   given
c   ---------------------------
        do i = 1, nos
          p(i)    =  rho(i) * cSound2
          e(i)    =  p(i) / (gamma-one) !Perfect gas equation for internal energy
          gamc(i) =  gamma
          game(i) =  gamma
          tem(i)  =  diskT !Temperature of the isothermal disk
        enddo
c   --------------------------
      else    !    P r e s s u r e    given
c   --------------------------
        do i = 1, nos
          rho(i)  =  p(i) / cSound2
          gamc(i) =  gamma
          game(i) =  gamma
          e(i)    =  p(i) / (gamma-one) !Perfect gas equation for internal energy
          tem(i)  =  diskT !!Temperature of the isothermal disk
        enddo
c   -------------------------
      endif
c   -------------------------
      return
      end
c----------------------------------------------------------------------
