  PROGRAM many_resources ! population-based consumer-resource, dimension-dependent 
  IMPLICIT NONE          !  exponents in alpha constraints
  DOUBLE PRECISION, ALLOCATABLE :: X(:,:),C(:,:),X0(:),BR(:),RC(:),GA(:)
  DOUBLE PRECISION, ALLOCATABLE :: MC(:),BRD(:),S(:),PA(:),PAN(:),DX(:)
  LOGICAL, ALLOCATABLE :: LV(:) 
  INTEGER N,K,K1,K2,K3,I,IMAX,I0,IL,IT,IW,NH,NC,J1
  DOUBLE PRECISION T,DT,TMAX,TMEAS,RAN3,EPS,MUT,BRT,W,TVF,RR,R,SPEED
  DOUBLE PRECISION DTM,TV,aa,comp, DIST,XP,GW,GS,SA,POPMIN,SQ3,DCMEAS
  DOUBLE PRECISION dmin,plotcount,A1,A2,A3,AE,DR,SQ2,SQ6,SQ15,PIN,DISPCH
  DOUBLE PRECISION DTSEED, DTMERGE,TSEED,TMERGE,PT,BRW,DCC,DCMIN,MU,MK  
  DOUBLE PRECISION PTOT
  N=3 ! SPACE DIMENSIONALITY
  IMAX=500 ! MAX NUMBER OF POPULATIONS
  DR=.25D0   ! DEATH RATE
  EPS=0.005D0 ! SPREAD AROUND THE INITIAL POSITION 
  MUT=0.005D0 ! MUTATION AMPLITUDE
  MU=0.1D0    ! FRACTION OF ANCESTOR THAT GOES INTO THE OFFSPRING
  DCMIN=MUT ! MERGING DISTANCE
  DCMEAS=MUT ! MERGING DISTANCE
  I0=1 ! INITIAL NUMBER OF POPULATIONS
  PIN=1.D0 ! INITIAL POPULATION
  POPMIN=1.D-5
  DT=1.D-2 ! TIME STEP
  TMAX=2000000.D0 ! END TIME
  DTSEED=1.D0 ! TIME TO SEED A NEW CLUSTER
  DTMERGE=100.D0*DTSEED !TIME TO MERGE CLUSTERS IN SEEDING PERCLUSTER
  DTM=2000.D0 ! MEASURE TIME  SA=0.5D0 ! WIDTH OF A GAUSSIAN COMPETITION KERNEL
  ALLOCATE(X(N,IMAX),C(N,IMAX),BR(IMAX),LV(IMAX),BRD(N),S(N))
  ALLOCATE(X0(N),RC(N),MC(IMAX),PA(IMAX),DX(N),GA(N))
  GA=(/1.1D0,1.1D0,.9D0/) ! EXPONENTS IN THE CONTRAINT
  open(unit=12,file='rand.dat') ! READING THE RANDOM SEED
  read(12,*)iw
!   iw=-87615
  close(12)
  SQ2=DSQRT(2.D0)
  SQ6=DSQRT(6.D0)
  SQ15=DSQRT(1.5D0)
  SQ3=DSQRT(3.D0)
  X0=0.D0
  LV=.FALSE.
  S=1.D0 ! SETTING THE RESOURCES
  MK=0.1D0!S(1)/4 ! THE DECAY COEFFICIENT FOR RESOURCES

  X0=(/0.9D0,0.05D0,0.05D0/)

  I=I0
  DO K=1,3
     X(K,I)=X0(K)
  END DO
   LV(I)=.TRUE.
    PA(I)=3.D0/I0
  IT=I0

  DO I=1,I0 ! PUTTING EVERYBODY ON CURVED SIMPLEX
     A1=0.D0
       DO K=1,N
          A1=A1+X(K,I)**GA(K)
       END DO
        DO K=1,N
          X(K,I)=X(K,I)/A1**(1.D0/GA(K))
       END DO
  END DO

!!$  open(unit=49,file='initial_config.dat')
!!$  DO K=1,IT
!!$  IF(LV(K)) WRITE(49,*) (X(2,K)-X(1,K))/SQ2,X(3,K)*SQ15
!!$  END DO
!!$  close(49)
!!$it=4000 ! reading the previous output, it should be adjusted
!!$  open(unit=35, file='points.dat') !reading a snapshot
!!$!  read(35,*)t
!!$  DO K=1,IT
!!$  LV(K)=.TRUE.
!!$  read (35,*) (X(K1,K),K1=1,N)
!!$  END DO
!!$  close (35)

  IL=IT ! TOTAL NUMBER OF INDIVIDUALS=NUMBER OF LIVE INDIVIDUALS
  T=0.D0
  TMEAS=DTM
  TSEED=0.D0
  TMERGE=DTMERGE-DTSEED
  plotcount=0

  OPEN(UNIT=10,FILE='fvar_pop.dat')
  OPEN(UNIT=40,FILE='total_var.dat')

  DO WHILE (T.LE.TMAX.AND.IT.LT.IMAX-1) !!!! BEGINNING OF UPDATE !!!!
  BR=0.D0 ! BIRTH RATES
  BRW=0.D0
  BRD=0.D0
  DO K=1,N  ! PHI(I)
    DO I=1,IT
    BRD(K)=BRD(K)+X(K,I)*PA(I) ! PHI(I)
    END DO
    A1=S(K)+BRD(K)+MK
    IF(A1*A1-4.D0*BRD(K)*S(K).LE.0.D0)WRITE(*,*) 'ROOT'
    IF(BRD(K).EQ.0.D0)WRITE(*,*)'DENOMINATOR'
    RC(K)=0.5D0 * (A1 - DSQRT(A1*A1-4.D0*BRD(K)*S(K)))/BRD(K)
  END DO

  DO I=1,IT  ! BIRTH RATE
    DO K=1,N
       BR(I)=BR(I)+X(K,I)*RC(K)
    END DO
    BRW=BRW+BR(I)*PA(I)
 END DO
 IF (T.GE.TSEED) THEN ! SPLITTING EVENT
       TSEED=T+DTSEED
       RR=ran3(iw)
        GS=0.D0
        K=0
        DO WHILE (GS.LE.RR) ! CHOOSING THE FATHER K
          K=K+1 ! PROPORTIONAL TO ITS BIRTH RATE
          IF (LV(K)) GS=GS+BR(K)*PA(K)/BRW
       END DO
     I=1
     DO WHILE (LV(I)) ! FINDING A SPARE SPOT I
     I=I+1
     END DO
     IF(I.GT.IT) IT=I ! ADDING A NEWBORN TO THE RIGHT END OF THE LIST
     IL=IL+1 ! INCREMENTING THE NUMBER OF LIVE CLUSTERS
     LV(I)=.TRUE. ! ASSINGING A SEAT
      PA(I)=PA(K)*mu ! SPLITTING THE ANCESTOR
      PA(K)=PA(K)*(1.d0-mu)
      BR(I)=BR(K) ! TEMPORARILY THE SAME BIRTH RATE

      A2=0.D0
      DO K1=1,N ! COMPUTING DX
         DX(K1) = MUT*(2*ran3(iw)-1.D0)
         X(K1,I)=X(K1,K)+DX(K1)
         X(K1,I)=DMAX1(1.D-8,X(K1,I))
         X(K1,I)=DMIN1(1.D0,X(K1,I))
         A2=A2+X(K1,I)**GA(K1)
      END DO

      DO K2=1,N ! PUTTING X ON CONSTRAINT SURFACE
         X(K2,I)=DMIN1(1.D0,X(K2,I)/A2**(1.D0/GA(K2)))
      END DO
   END IF
      
  PTOT=0.D0
  DO I=1,IT  ! UPDATE OF POPULATION
      IF(LV(I)) THEN
         PA(I)=PA(I)+(BR(I)-DR)*PA(I)*DT ! NEW POPULATION
         PTOT=PTOT+PA(I)
         IF (PA(I).LE.POPMIN) THEN
            PA(I)=0.D0
!            WRITE(*,*)'EXTINCTION'
            LV(I)=.FALSE.
            IL=IL-1
         END IF
      END IF
   END DO
   IF(PTOT.LE.1.D-7) THEN
      WRITE(*,*)'PTOT',PTOT
      READ(*,*) A1
   END IF   
   T=T+DT

    IF(T.GE.TMERGE.OR.T.GT.TMAX) THEN ! MERGING CLOSE CLUSTERS
    TMERGE=T+DTMERGE
     DO K1=1,IT 
        IF(LV(K1)) THEN ! removing coinciding clusters
              DO K2=1,K1-1
               IF(LV(K2).AND.LV(K1)) THEN
                  DCC=0.D0
                    DO J1=1,N
                       DCC=DCC+(X(J1,K1)-X(J1,K2))**2
                    END DO
                  DCC=DSQRT(DCC)
                  IF (DCC.LE.DCMIN) THEN !LIQUIDATING CLUSTER KC
        DO J1=1,N ! PRESERVING CENTRE OF MASS POSITION
           X(J1,K2)=(X(J1,K1)*PA(K1)+X(J1,K2)*PA(K2))/(PA(K1)+PA(K2)) 
        END DO
                     LV(K1)=.FALSE.
                     PA(K2)=PA(K2)+PA(K1)
                     PA(K1)=0.D0                    
                     IL=IL-1
                  END IF
               END IF
               END DO
            END IF
        END DO
END IF ! END OF MERGING

    IF(T.GE.TMEAS) THEN ! MEASUREMENT
        TMEAS=T+DTM

     DO K1=1,IT !merging once again, now with larger cutoff  
        IF(LV(K1)) THEN ! removing coinciding clusters
              DO K2=1,K1-1
               IF(LV(K2).AND.LV(K1)) THEN
                  DCC=0.D0
                    DO J1=1,N
                       DCC=DCC+(X(J1,K1)-X(J1,K2))**2
                    END DO
                  DCC=DSQRT(DCC)
                  IF (DCC.LE.DCMEAS) THEN !LIQUIDATING CLUSTER KC
        DO J1=1,N ! PRESERVING CENTRE OF MASS POSITION
           X(J1,K2)=(X(J1,K1)*PA(K1)+X(J1,K2)*PA(K2))/(PA(K1)+PA(K2)) 
        END DO
                     LV(K1)=.FALSE.
                     PA(K2)=PA(K2)+PA(K1)
                     PA(K1)=0.D0                    
                     IL=IL-1
                  END IF
               END IF
               END DO
            END IF
        END DO



        
             WRITE(10,*) (plotcount*100,K1=1,3)
     plotcount=plotcount+1
         DO K=1,IT
            IF(LV(K)) THEN
  WRITE(10,*) (X(2,K)**GA(2)-X(1,K)**GA(1))/DSQRT(2.D0),X(3,K)**GA(3)*DSQRT(1.5D0),PA(K)
            END IF
        END DO
        WRITE(40,*) T, IL
        WRITE(*,*)  T, IL
!!$  open(unit=75, file='final_points_ran.dat') !writing last snapshot,
!!$  DO K=1,IT
!!$  IF(LV(K)) WRITE(75,*) (X(2,K)-X(1,K))/DSQRT(2.D0),X(3,K)*DSQRT(1.5D0),PA(K)
!!$  END DO
!!$  close (75) 

END IF

  END DO

   CLOSE(10)
   CLOSE(40)


  open(unit=12, file='rand.dat')
  write(12,*)nint(-ran3(iw)*1.d+6)
  close(12)         
 
  write(*,*) 'survivors, it, il', it, il
END PROGRAM many_resources

 

!==========================================================================
!  ran3 from Numerical Recipes, 2nd edition
!--------------------------------------------------------------------------
      function ran3(idum)
!==========================================================================

!         implicit real*4(m)

      integer idum
      integer mbig,seed,mz
!      real mbig,seed,mz
      double precision ran3,fac

!         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
!      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
 
      integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(55)
!      real  mj,mk,ma(55)

      save iff,inext,inextp,ma

      data iff /0/

1      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      if(ran3.le.0.or.ran3.ge.1) goto 1
      return
      end function ran3



