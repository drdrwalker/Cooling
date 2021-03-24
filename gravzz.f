C-----------------------------------------------------------------------

      SUBROUTINE POTENT

C  POTENTIAL SOLVER.  CHANGE TO BUFFER IN/OUT, 1 DEC 1981

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      COMMON /ROPO/ RO(qx,qx,qz), PO(qx,qx,qz)
      REAL*8 LOW(qxy,qx), UP(qxy,qx), FAKE(qx1*qx1*q2x-qxy*q2x)
      REAL*8 HOLD(qx1,qx1,qx1), fort(qx*qx*qx,3)
      COMMON/RAYPOT/ LOW, UP, FAKE, fort, dummy
      REAL*8 dummy(qx1*qx1*qx1)
      REAL*8 RHO(qx,qx,qx),POT(qx2,qx2,qx2)
      EQUIVALENCE (LOW(1,1),RHO(1,1,1)),(UP(1,1),POT(1,1,1))
      real*8 rhol(qx,qx,q2x)
      equivalence (rhol, rho)

      COMMON /CPGEN/ HOLD

      !print *, 'gravzz, POTENT'


      CALL COPYHT(RO,RHO)
      CALL MIRROR(ltwo,RHO)

C  START CALCULATION.

      CALL CLEAR(UP)
      CALL FWDFT(ltwo,LOW)   !2
      CALL PTRANS2D
      CALL XSAVE(UP,lone)
      CALL CLEAR(UP)

      CALL FWDFT(ltwo,LOW)   !2
      CALL MIRROR(lone,LOW)
      CALL PARATRANSP
      CALL XSAVE(UP,ltwo)
      CALL CLEAR(UP)

      CALL FWDFT(lone,LOW)   !1

      DO 10 M=1,q2x
         MR=MIN(M,q2x+2-M)
      DO 10 L=1,qx
      DO 10 K=1,qx
         rhol(K,L,M)=rhol(K,L,M)*HOLD(K,L,MR)
   10 CONTINUE

C  MULTIPLY.  3L,1L,2

      CALL REVFT(lone,LOW)   !1
      CALL XSAVE(LOW,lthree)

      CALL RESTOR(LOW,ltwo)
      CALL CLEAR(UP)
      CALL FWDFT(lone,LOW)   !1

      DO 20 M=1,q2x
         MR=MIN(M,q2x+2-M)
      DO 20 L=1,qx
      DO 20 K=1,qx
         KR=qx2-K
         rhol(K,L,M)=rhol(K,L,M)*HOLD(KR,L,MR)
   20 CONTINUE

C  MULTIPLY.  3L,1U,2

      CALL REVFT(lone,LOW)   !1
      CALL COPY(LOW,UP)

      CALL RESTOR(LOW,lthree)
      CALL PARARTRANS
      CALL REVFT(ltwo,LOW)   !2
      CALL XSAVE(LOW,ltwo)

      CALL RESTOR(LOW,lone)
      CALL CLEAR(UP)
      CALL FWDFT(ltwo,LOW)   !2
      CALL MIRROR(lone,LOW)
      CALL PARATRANSP
      CALL XSAVE(UP,lone)
      CALL CLEAR(UP)

      CALL FWDFT(lone,LOW)   !1

      DO 30 M=1,q2x
         MR=MIN(M,q2x+2-M)
      DO 30 L=1,qx
         LR=qx2-L
      DO 30 K=1,qx
         rhol(K,L,M)=rhol(K,L,M)*HOLD(K,LR,MR)
   30 CONTINUE

C  MULTIPLY.  3U,1L,2

      CALL REVFT(lone,LOW)   !1
      CALL XSAVE(LOW,lthree)

      CALL RESTOR(LOW,lone)
      CALL CLEAR(UP)
      CALL FWDFT(lone,LOW)   !1

      DO 40 M=1,q2x
         MR=MIN(M,q2x+2-M)
      DO 40 L=1,qx
         LR=qx2-L
      DO 40 K=1,qx
         KR=qx2-K
         rhol(K,L,M)=rhol(K,L,M)*HOLD(KR,LR,MR)
   40 CONTINUE

C  MULTIPLY.  3U,1U,2

      CALL REVFT(lone,LOW)   !1

      CALL COPY(LOW,UP)
      CALL RESTOR(LOW,lthree)
      CALL PARARTRANS
      CALL REVFT(ltwo,LOW)   !2

      CALL COPY(LOW,UP)
      CALL RESTOR(LOW,ltwo)
      CALL PTRANS2D
      CALL REVFT(ltwo,LOW)   !2
      CALL CLEAR(UP)

      CALL RCPYHT(LOW,PO)


      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FWDFT(NH,A)
      include 'qparam.cmn'
      INTEGER*4 NH, qxynh
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      REAL*8 A(qxy,q2x), TEMP(qxy)

      qxynh = qxy/NH
C  BIT-REVERSED TRANSPOSITION

      !print *, 'gravzz, FWDFT'

      j = lone
      do 30 i = 1, q2x
         if (j .gt. i) then
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(KK), SHARED(A,TEMP),
C$OMP+            FIRSTPRIVATE(qxynh,j,i)
            DO KK=1,qxynh
               TEMP(KK)=A(KK,j)
               A(KK,j)=A(KK,i)
               A(KK,i)=TEMP(KK)
            enddo
C$OMP END PARALLEL DO
         endif
         m = q2x / ltwo
 21      if ( (m .ge. ltwo) .and. (j .gt. m) ) then
             j = j - m
             m = m / ltwo
         goto 21
         endif
         j = j + m
 30   continue


      M=lone

  400 CONTINUE

C  CALL RADIX 2 REAL*8 KERNEL

      CALL R2RFTK(qxy,NH,A,M,lone)

      M=M*ltwo
      IF (M.LT.q2x) GO TO 400

      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE REVFT(NH,A)
      include 'qparam.cmn'
      INTEGER*4 NH, qxynh
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      REAL*8 A(qxy,q2x), TEMP(qxy)

      qxynh = qxy/NH

C  COMPLEX CONJUGATE DATA

      !print *, 'gravzz, REVFT'

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(J,KK), SHARED(A),
C$OMP+            FIRSTPRIVATE(qxynh)
      DO 200 J=qx2,q2x
         DO 200 KK=1,qxynh
  200       A(KK,J)=-A(KK,J)
C$OMP END PARALLEL DO

      M=q2x

  300 CONTINUE
      M=M/ltwo

C  CALL RADIX 2 HERMITE KERNEL

      CALL R2HFTK(qxy,NH,A,M,lone)
      IF (M.GT.lone) GO TO 300

C  BIT-REVERSED TRANSPOSITION

      j = lone
      do 30 i = 1, q2x
         if (j .gt. i) then
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(KK), SHARED(A,TEMP),
C$OMP+            FIRSTPRIVATE(qxynh,j,i)
            DO 20 KK=1,qxynh
               TEMP(KK)=A(KK,j)
               A(KK,j)=A(KK,i)
               A(KK,i)=TEMP(KK)
 20         CONTINUE
C$OMP END PARALLEL DO
         endif
         m = q2x / ltwo
 21      if ( (m .ge. ltwo) .and. (j .gt. m) ) then
             j = j - m
             m = m / ltwo
         goto 21
         endif
         j = j + m
 30   continue


      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE R2RFTK(QN,NH,A,M,L0)

C  RADIX 2 REAL*8 KERNEL, FROM SANDES ROUTINES.

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      INTEGER*4 QN,NH,M,L0, QNNH
      REAL*8 A(QN,q2x)

      INTEGER*4 J,J0,J1,K,K0,K1,KK,L1,MOVER2,M2
      INTEGER*4 K00,K01,K10,K11
      REAL*8 ANGLE,C,IS(qxy),IU(qxy),I1,RS(qxy),RU(qxy),R1,S,TWOPI
      REAL*8 T00,T01,T10,T11

      DATA TWOPI /6.2831 85307 17958 64769/

      QNNH = QN/NH
      if ( QN .gt. qxy ) stop ' R2RFTK'

      L1=L0+M
      MOVER2=(M-1)/ltwo
      M2=M*ltwo

      !print *, 'gravzz, R2RFTK'

      DO 100 K=1,q2x,M2
         K0=K+L0-1
         K1=K+L1-1
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(KK), SHARED(A,RS,RU),
C$OMP+            FIRSTPRIVATE(QNNH,K0,K1)
         DO KK=1,QNNH
            RS(KK)=A(KK,K0)
            RU(KK)=A(KK,K1)
            A(KK,K0)=RS(KK)+RU(KK)
            A(KK,K1)=RS(KK)-RU(KK)
         enddo
C$OMP END PARALLEL DO
  100 CONTINUE

      IF (MOVER2 .LT. lone) GO TO 400

      DO 300 J=1,MOVER2
         J0=J+lone
         J1=M-J*ltwo
         ANGLE=TWOPI*dble(J)/dble(M2)
         C=COS(ANGLE)
         S=SIN(ANGLE)

         DO 300 K0=J0,q2x,M2
            K1=K0+J1
            K10=L1+K0-lone
            K11=L1+K1-lone
            K00=L0+K0-lone
            K01=L0+K1-lone
C$OMP PARALLEL DO DEFAULT(NONE), SHARED(A,RS,IS,RU,IU),
C$OMP+    FIRSTPRIVATE(QNNH,K00,K01,K10,K11,S,C), PRIVATE(KK,I1,R1)
            DO KK=1,QNNH
               R1=A(KK,K10)*C+A(KK,K11)*S
               I1=A(KK,K11)*C-A(KK,K10)*S
               RS(KK)=A(KK,K00)+R1
               IS(KK)=A(KK,K01)+I1
               RU(KK)=A(KK,K00)-R1
               IU(KK)=A(KK,K01)-I1
               A(KK,K00)=RS(KK)
               A(KK,K01)=RU(KK)
               A(KK,K11)=IS(KK)
               A(KK,K10)=-IU(KK)
            enddo
C$OMP END PARALLEL DO
  300 CONTINUE

  400 CONTINUE

      IF (M.NE.(M/2)*ltwo) GO TO 600
      J=M/ltwo+lone

      DO 500 K=J,q2x,M2
         K1=L1+K-lone
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(KK), SHARED(A),
C$OMP+            FIRSTPRIVATE(QNNH,K1)
         do KK=1,QNNH
            A(KK,K1)=-A(KK,K1)
         enddo
C$OMP END PARALLEL DO
  500 CONTINUE
  600 CONTINUE

      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE R2HFTK(QN,NH,A,M,L0)

C  RADIX 2 HERMITE KERNEL, FROM SANDES ROUTINES.

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x, qx11
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy=qx*qx , qxz=qx*qzg , q2x=2*qx, qx11=qx1*qx1 )

      INTEGER*4 QN,NH,M,L0,QNNH
      REAL*8 A(QN,q2x)

      INTEGER*4 J,J0,J1,K,K0,K1,KK,L1,MOVER2,M2
      INTEGER*4 K00,K01,K10,K11
      REAL*8 ANGLE,C,IS(qx11),IU(qx11),RS(qx11),RU(qx11),S,TWOPI
      REAL*8 T00,T01,T10,T11

      DATA TWOPI /6.2831 85307 17958 64769/

      QNNH = QN/NH
      if ( QN .gt. qx11 ) stop ' R2HFTK'

      L1=L0+M
      MOVER2=(M-1)/ltwo
      M2=M*ltwo

      !print *, 'gravzz, R2HFTK'

      DO 100 K=1,q2x,M2
         K0=K+L0-lone
         K1=K+L1-lone
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(KK), SHARED(A,RS,RU),
C$OMP+            FIRSTPRIVATE(QNNH,K0,K1)
         DO KK=1,QNNH
            RS(KK)=A(KK,K0)
            RU(KK)=A(KK,K1)
            A(KK,K0)=RS(KK)+RU(KK)
            A(KK,K1)=RS(KK)-RU(KK)
         enddo
C$OMP END PARALLEL DO
  100 CONTINUE

      IF (MOVER2 .LT. lone) GO TO 400

      DO 300 J=1,MOVER2
         J0=J+lone
         J1=M-J*ltwo
         ANGLE=TWOPI*dble(J)/dble(M2)
         C=COS(ANGLE)
         S=SIN(ANGLE)

         DO 300 K0=J0,q2x,M2
            K1=K0+J1
            K10=L1+K0-lone
            K11=L1+K1-lone
            K00=L0+K0-lone
            K01=L0+K1-lone
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(KK), SHARED(A,RS,IS,RU,IU),
C$OMP+            FIRSTPRIVATE(QNNH,K00,K01,K10,K11,S,C)
            DO KK=1,QNNH
               RS(KK)=A(KK,K00)+A(KK,K01)
               IS(KK)=A(KK,K11)-A(KK,K10)
               RU(KK)=A(KK,K00)-A(KK,K01)
               IU(KK)=A(KK,K11)+A(KK,K10)
               A(KK,K00)=RS(KK)
               A(KK,K01)=IS(KK)
               A(KK,K10)=RU(KK)*C+IU(KK)*S
               A(KK,K11)=IU(KK)*C-RU(KK)*S
            enddo
C$OMP END PARALLEL DO
  300 CONTINUE

  400 CONTINUE

      IF (M.NE.(M/ltwo)*ltwo) GO TO 600
      J=M/ltwo+lone

      DO 500 K=J,q2x,M2
         K0=L0+K-lone
         K1=L1+K-lone
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(KK), SHARED(A),
C$OMP+            FIRSTPRIVATE(QNNH,K0,K1)
         DO KK=1,QNNH
            A(KK,K0)=two*A(KK,K0)
            A(KK,K1)=two*A(KK,K1)
         enddo
C$OMP END PARALLEL DO
  500 CONTINUE

  600 CONTINUE

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE PGEN(WON,ASQ)

C  FREE-BOUNDARY POTENTIAL CONVOLUTION COEFFICIENTS.
C    21 MAY 1981

      include 'qparam.cmn'
      REAL*8  WON, ASQ
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      REAL*8 RHO(qx1,qx1,q2x), HOLD(qx1,qx1,qx1), fort(qx*qx*qx,3)
      COMMON/RAYPOT/ RHO, fort, dummy
      REAL*8 dummy(qx1,qx1,qx1)

      COMMON /CPGEN/ HOLD
      REAL*8 VPOT(qx1,qx1,qx1)
      EQUIVALENCE(RHO,VPOT)

      !print *, 'gravzz, PGEN'

      COEF=WON/dble(q2x*q2x*q2x)

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(rho), PRIVATE(M,L,K,ZSQ,YSQ,XSQ),
C$OMP+            FIRSTPRIVATE(ASQ,COEF)
      DO 1 M=1,qx1
c        ZSQ=ASQ+dble((M-1)**2)
         ZSQ=    dble((M-1)**2)
      DO 1 L=1,qx1
         YSQ=ZSQ+dble((L-1)**2)
      DO 1 K=1,qx1
         XSQ=YSQ+dble((K-1)**2)
c        RHO(K,L,M) = COEF / SQRT( XSQ )
         RHO(K,L,M) = COEF / SQRT( XSQ + ASQ*exp(-xsq/asq) )
    1 CONTINUE
C$OMP END PARALLEL DO

      CALL HERMFT(RHO)

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(hold,vpot), PRIVATE(i,j,k)
      do k = 1, qx1
      do j = 1, qx1
      do i = 1, qx1
         hold(i,j,k) = vpot(i,j,k)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE HERMFT(A)

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x, qx11
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy=qx*qx, qxz=qx*qzg, q2x=2*qx, qx11=qx1*qx1 )

      REAL*8 A(qx11,q2x), TEMP(qx11)

C  LOOP FOR EACH OF THREE DIMENSIONS

      !print *, 'gravzz, HERMFT'

      KXFM=lzero
  100 CONTINUE

C  CLEAR HIGH MEMORY

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(J,KK), SHARED(A)
      DO 200 J=qx2,q2x
         DO 200 KK=1,qx11
  200       A(KK,J)=zero
C$OMP END PARALLEL DO

      M=q2x

  300 CONTINUE
      M=M/ltwo

C  CALL RADIX 2 HERMITE KERNEL

      CALL R2HFTK(qx11,lone,A,M,lone)
      IF (M.GT.lone) GO TO 300

C  BIT-REVERSED TRANSPOSITION

      j = lone
      do 30 i = 1, q2x
         if (j .gt. i) then
C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(KK), SHARED(A,TEMP),
C$OMP+            FIRSTPRIVATE(i,j)
            DO 20 KK=1,qx11
               TEMP(KK)=A(KK,j)
               A(KK,j)=A(KK,i)
               A(KK,i)=TEMP(KK)
 20         CONTINUE
C$OMP END PARALLEL DO
         endif
         m = q2x / ltwo
 21      if ( (m .ge. ltwo) .and. (j .gt. m) ) then
             j = j - m
             m = m / ltwo
         goto 21
         endif
         j = j + m
 30   continue

      CALL PARATR65

      KXFM=KXFM+lone
      IF (KXFM.LT.lthree) GO TO 100

      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE PARATR65

C PARALLEL
C  CYCLIC TRANSPOSITION OF 3-INDEX ARRAY
C    qx1,qx1,q2x VERSION TO WORK IN PGEN.  21 MAY.

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      REAL*8 A(qx1,qx1,q2x), fort(qx,qx,qx,3), temp(qx1,qx1,qx1)
      COMMON/RAYPOT/ A, fort, temp
c      write(*,*) 8*(qx1*qx1*q2x + qx*qx*qx*3 + qx1*qx1*qx1)
c      stop

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(temp,A)

      !print *, 'gravzz, PARATR65'

      DO 1 M=1,qx1
      DO 1 L=1,qx1
      DO 1 K=1,qx1
    1    temp(K,L,M) = A(K,L,M)
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(temp,A)
      DO 2 M=1,qx1
      DO 2 L=1,qx1
      DO 2 K=1,qx1
    2    A(K,L,M) = temp(L,M,K)
C$OMP END PARALLEL DO

      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE PTRANS2D

c  transposition of two '2d' arrays

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      REAL*8 A(qx,qx,qx),B(qx,qx,qx),FAKE(qx1*qx1*q2x-q2x*qx*qx),
     &           fort(qx,qx,qx,3), TEMP(qx,qx,qx)
      COMMON/RAYPOT/ A, B, FAKE, fort, TEMP, dummy
      REAL*8 dummy(qx1*qx1*qx1-qx*qx*qx)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(A,TEMP)

      !print *, 'gravzz, PTRANS2D'
      DO 1  L=1,qx
      DO 1 M=1,qzg
      DO 1 K=1,qx
    1    TEMP(K,M,L) = A(K,M,L)
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(A,TEMP)
      DO 2 K=1,qx
      DO 2 L=1,qx
      DO 2 M=1,qzg
         A(L,M,K) = TEMP(K,M,L)
    2 CONTINUE
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(B,TEMP)
      DO 3 L=1,qx
      DO 3 M=1,qzg
      DO 3 K=1,qx
    3    TEMP(K,M,L) = B(K,M,L)
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(B,TEMP)
      DO 4 K=1,qx
      DO 4 L=1,qx
      DO 4 M=1,qzg
         B(L,M,K) = TEMP(K,M,L)
    4    CONTINUE
C$OMP END PARALLEL DO

      RETURN
      END


C-----------------------------------------------------------------------

      SUBROUTINE PARATRANSP

c   PARALLEL
C   CYCLIC TRANSPOSITION OF TWO 3-INDEX ARRAYS
C     THE TWO ARRAYS ARE CONTIGUOUS IN STORAGE
C     IN COMMON/RAYPOT/

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      REAL*8 A(qx,qx,qx),B(qx,qx,qx),FAKE(qx1*qx1*q2x-q2x*qx*qx),
     &           fort(qx,qx,qx,3), temp(qx,qx,qx)
      COMMON/RAYPOT/ A, B, FAKE, fort, temp, dummy
      REAL*8 dummy(qx1*qx1*qx1-qx*qx*qx)
      integer*4 k,l,m


C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(temp,A)

      !print *, 'gravzz, PARATRANSP'

      DO 1 M=1,qx
      DO 1 L=1,qx
      DO 1 K=1,qx
         temp(K,L,M) = A(K,L,M)
    1 CONTINUE
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(temp,A)
      DO 2 M=1,qx
      DO 2 L=1,qx
      DO 2 K=1,qx
    2    A(K,L,M) = temp(L,M,K)
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(temp,B)
      DO 3 M=1,qx
      DO 3 L=1,qx
      DO 3 K=1,qx
    3    temp(K,L,M) = B(K,L,M)
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(temp,B)
      DO 4 M=1,qx
      DO 4 L=1,qx
      DO 4 K=1,qx
    4    B(K,L,M) = temp(L,M,K)
C$OMP END PARALLEL DO

      RETURN

C-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

      ENTRY PARARTRANS

C   PARALLEL
C   REVERSE TRANSFORMATION TO ENTRY TRANSP.


C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(temp,A)
      DO 5 M=1,qx
      DO 5 L=1,qx
      DO 5 K=1,qx
    5    temp(K,L,M) = A(K,L,M)
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(temp,A)
      DO 6 M=1,qx
      DO 6 L=1,qx
      DO 6 K=1,qx
    6    A(K,L,M) = temp(M,K,L)
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(temp,B)
      DO 7 M=1,qx
      DO 7 L=1,qx
      DO 7 K=1,qx
    7    temp(K,L,M) = B(K,L,M)
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(K,L,M), SHARED(temp,B)
      DO 8 M=1,qx
      DO 8 L=1,qx
      DO 8 K=1,qx
    8    B(K,L,M) = temp(M,K,L)
C$OMP END PARALLEL DO


      END

C-----------------------------------------------------------------------

      SUBROUTINE CLEAR(A)

C  CLEAR WORKING AREA

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      REAL*8 A(qxy,qx)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(KK,K), SHARED(A)

      !print *, 'gravzz, CLEAR'

      DO 10 KK=1,qx
      DO 10 K=1,qxy
   10    A(K,KK)=zero
C$OMP END PARALLEL DO

      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE COPY(A,B)

C  COPY CONTENTS OF ARRAY A INTO ARRAY B

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      REAL*8 A(qxy,qx),B(qxy,qx)

C$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(KK,K), SHARED(A,B)

    !  print *, 'gravzz, COPY'

      DO 10 KK=1,qx
      DO 10 K=1,qxy
  10     B(K,KK)=A(K,KK)
C$OMP END PARALLEL DO

      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE COPYHT(A,B)

C  COPY AND TRANSPOSE

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      real*8 a(qx,qx,qz),b(qx,qx,qx)
      integer*4 i,j,k

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(a,b), PRIVATE(i,j,k)

    !  print *, 'gravzz, COPYHT'

      do 10 i=1,qx
      do 10 j=1,qz
      do 10 k=1,qx
  10    b(k,qzg1-j,i) = a(k,i,j)
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(b), PRIVATE(i,j,k)
      do 11 i=1,qx
      do 11 j=qz+1,qzg
      do 11 k=1,qx
  11    b(k,qzg1-j,i) = zero
C$OMP END PARALLEL DO

      return

c-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

      ENTRY RCPYHT (B,A)

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(a,b), PRIVATE(i,j,k)
      do 20 i=qz+1,qzg
      do 20 j=1,qx
      do 20 k=1,qx
  20    a(k,j,qzg1-i) = -b(k,i,j)
C$OMP END PARALLEL DO

      return
      end

C-----------------------------------------------------------------------

      SUBROUTINE XSAVE(A,N)

c  array save/restore                        07 Oct 1991

      include 'qparam.cmn'
      integer*4     qxxx, qx1, q2x
      parameter ( qxxx = qx*qx*qx, qx1 = qx+1, q2x=2*qx )
      COMMON/RAYPOT/ RHO(qx1,qx1,q2x), fort(qxxx,3), dummy
      REAL*8 A(qxxx), dummy(qx1,qx1,qx1)
      INTEGER*4 N,K,L

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(a,fort), PRIVATE(i),
C$OMP+            FIRSTPRIVATE(N)

      !print *, 'gravzz, XSAVE'

      do i = 1, qxxx
         fort(i,N) = a(i)
      enddo
C$OMP END PARALLEL DO

      RETURN

C-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

      ENTRY RESTOR(A,N)

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(a,fort), PRIVATE(i),
C$OMP+            FIRSTPRIVATE(N)

      do i = 1, qxxx
         a(i) = fort(i,N)
      enddo
C$OMP END PARALLEL DO


      RETURN
      END


C-----------------------------------------------------------------------

      SUBROUTINE MIRROR(NH,A)

C  COPY 'TOP' HALF INTO 'BOTTOM' HALF

      include 'qparam.cmn'
      integer*4     qx1, qzg1, qx2, qzg2, qxy, qxz, q2x
      parameter ( qx1 = qx+1 ,  qzg1 = qzg+1 )
      parameter ( qx2 = qx+2 ,  qzg2 = qzg+2 )
      parameter ( qxy = qx*qx , qxz = qx*qzg , q2x = 2*qx )

      INTEGER*4 NH,I,J,K
      REAL*8 A(qx,qx,q2x)

C$OMP PARALLEL DO DEFAULT(NONE), SHARED(A), PRIVATE(i,j,k),
C$OMP+   FIRSTPRIVATE(NH)

      !print *, 'gravzz, MIRROR'

      DO K=1,q2x/NH
      DO J=1,qzg
      DO I=1,qx
   10    A(I,qx1-J,K) = A(I,J,K)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO

      RETURN
      END
