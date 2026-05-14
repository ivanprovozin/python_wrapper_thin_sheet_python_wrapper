c------------------------------------------------------------------------------
c---- Subprogram packages 'SUB_LUD.F' for Integral equation solution 'WSHEET.F'
c---- including Romberg 1-D numerical integral, 2-D qubic spline integral, and
c---- LU-decomposition routines.
c---- Latest revision : Sep. 24, 1997
c------------------------------------------------------------------------------
c
      SUBROUTINE QROMB(FUNC,A,B,SS)
      PARAMETER(EPS=1.e-6,JMAX=30,JMAXP=JMAX+1,K=5,KM=k-1)
      complex func, ss, s
      external func
      DIMENSION Sr(JMAXP),Hr(JMAXP),Si(JMAXP),Hi(JMAXP)
      Hr(1)=1.0
      Hi(1)=1.0
      i_r = 0
      i_i = 0
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,S,J)
        sr(j) = real(s)
        si(j) = aimag(s)
        IF (J.GE.K) THEN
          L=J-KM
          if (i_r .eq. 0) then
                CALL POLINT(Hr(L),Sr(L),K,0.0,SSr,DSSr)
                IF (ABS(DSSr).LT.EPS*ABS(SSr)) i_r = 1
          endif
          if (i_i .eq. 0) then
                CALL POLINT(Hi(L),Si(L),K,0.0,SSi,DSSi)
                IF (ABS(DSSi).LT.EPS*ABS(SSi)) i_i = 1
          endif
          IF ((i_r*i_i) .eq. 1) then
                ss = cmplx(ssr,ssi)
                RETURN
          endif
        ENDIF
        s = cmplx(sr(j),si(j))
        Hr(J+1)=0.25*Hr(J)
        Hi(J+1)=0.25*Hi(J)
11    CONTINUE
c      PAUSE 'Too many steps.'
      END
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=20) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.0)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE TRAPZD(FUNC,A,B,S,N)
      complex func, s, sum
      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=cmplx(0.0,0.0)
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
c
c-------------------------------------------------------------------------------
c
C   IMSL ROUTINE NAME   - DBCQDU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - BICUBIC SPLINE QUADRATURE
C
C   USAGE               - CALL DBCQDU (F,IFD,X,NX,Y,NY,A,B,C,D,Q,WK,IER)
C
C   ARGUMENTS    F      - NX BY NY MATRIX CONTAINING THE FUNCTION
C                           VALUES. (INPUT) F(I,J) IS THE FUNCTION VALUE
C                           AT THE POINT (X(I),Y(J)) FOR I=1,...,NX AND
C                           J=1,...,NY.
C                IFD    - ROW DIMENSION OF THE MATRIX F EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                X      - VECTOR OF LENGTH NX. (INPUT) X MUST BE
C                           ORDERED SO THAT X(I) .LT. X(I+1) FOR
C                           I=1,...,NX-1.
C                NX     - NUMBER OF ELEMENTS IN X. (INPUT) NX MUST BE
C                           .GE. 2.
C                Y      - VECTOR OF LENGTH NY. (INPUT) Y MUST BE
C                           ORDERED SO THAT Y(J) .LT. Y(J+1) FOR
C                           J=1,...,NY-1.
C                NY     - NUMBER OF ELEMENTS IN Y. (INPUT) NY MUST BE
C                           .GE. 2.
C                         NOTE - THE COORDINATE PAIRS (X(I),Y(J)),FOR
C                           I=1,...,NX AND J=1,...,NY, GIVE THE POINTS
C                           WHERE THE FUNCTION VALUES F(I,J) ARE
C                           DEFINED.
C                A,B    - X-DIRECTION LIMITS OF INTEGRATION. (INPUT)
C                C,D    - Y-DIRECTION LIMITS OF INTEGRATION. (INPUT)
C                Q      - DOUBLE INTEGRAL FROM A TO B AND C TO D.
C                           (OUTPUT)
C                WK     - WORK VECTOR OF LENGTH
C                           (NY+5)*NX+NY-1+MAX(5*NX-4,5*NY-2)
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IFD IS LESS THAN NX.
C                           IER = 130, NX IS LESS THAN 2.
C                           IER = 131, NY IS LESS THAN 2.
C                         WARNING ERROR
C                           IER = 36, A AND/OR B IS LESS THAN X(1).
C                           IER = 37, A AND/OR B IS GREATER THAN X(NX).
C                           IER = 38, C AND/OR D IS LESS THAN Y(1).
C                           IER = 39, C AND/OR D IS GREATER THAN Y(NY).
C
C   REQD. IMSL ROUTINES - ICSEVU,ICSCCU,UERSET,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      WHEN THE LIMITS OF INTEGRATION ARE OUTSIDE OF THE
C                RECTANGLE (X(1),X(NX)) X (Y(1),Y(NY)), THE
C                BICUBIC SPLINE IS EXTENDED USING THE BOUNDARY PIECES.
C                INTEGRATION IS PERFORMED OVER THE EXTENDED SPLINE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DBCQDU (F,IFD,X,NX,Y,NY,A,B,C,D,Q,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IFD,NX,NY,IER
      REAL               F(IFD,NY),X(NX),Y(NY),A,B,C,D,Q,WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IA,IAP1,IB,IC,ICP1,ID,IDM1,II,IP1,IPT,IV,
     1                   IYL,J,JER,JJER,KER,LER,MER,NCOEF,NER,NFI,NFI1,
     2                   NFL,NFTEMP,NFTMP1,NN,NWTXX,NXM1,NXX,NYL,NYLS,
     3                   NYLSP,NYM1
      REAL               AA,BB,CC,DD,DIST1,DIST2,FOUR,
     1                   FRDTRD,HALF,SUM,THIRD,V,WTYY,ZERO
      DATA               ZERO/0.0/,HALF/.5/,FOUR/4.0/
      DATA               THIRD/.3333333/,
     1                   FRDTRD/1.333333/
C                                  FIRST EXECUTABLE STATEMENT
      JER = 0
      KER = 0
      LER = 0
      MER = 0
      NER = 129
      IF (IFD.LT.NX) GO TO 140
      NER = 130
      IF (NX.LT.2) GO TO 140
      NER = 131
      IF (NY.LT.2) GO TO 140
      NER = 0
c      LEVEL = 2
c      CALL UERSET (LEVEL,LEVOLD)
      NXM1 = NX-1
      NYM1 = NY-1
C                                  FIND THE X INTERVALS
      IA = 1
      IPT = 1
      V = AMIN1(A,B)
    5 DIST1 = V-X(IA)
      DO 10 I=IA,NXM1
         IV = I
         DIST2 = V-X(I+1)
         IF (DIST2.LT.ZERO) GO TO 15
         IF (I.LT.NXM1) DIST1 = DIST2
   10 CONTINUE
      IV = NXM1
C                                  IF V .GT. X(NX) - WARNING
      IF (DIST2.GT.ZERO) KER = 37
   15 CONTINUE
C                                  CHECK FOR V .LT. X(1)
      IF (DIST1.LT.ZERO) JER = 36
      IF (IPT.EQ.2) GO TO 20
      IPT = 2
      IA = IV
      V = AMAX1(A,B)
      GO TO 5
   20 IB = IV
C                                  FIND THE Y INTERVALS
      IC = 1
      IPT = 1
      V = AMIN1(C,D)
   25 DIST1 = V-Y(IC)
      DO 30 I=IC,NYM1
         IV = I
         DIST2 = V-Y(I+1)
         IF (DIST2.LT.ZERO) GO TO 35
         IF (I.LT.NYM1) DIST1 = DIST2
   30 CONTINUE
      IV = NYM1
C                                  IF V .GT. Y(NY) - WARNING
      IF (DIST2.GT.ZERO) MER = 39
   35 CONTINUE
C                                  CHECK FOR V .LT. Y(1)
      IF (DIST1.LT.ZERO) LER = 38
      IF (IPT.EQ.2) GO TO 40
      IPT = 2
      IC = IV
      V = AMAX1(C,D)
      GO TO 25
   40 ID = IV
      AA = AMIN1(A,B)
      BB = AMAX1(A,B)
C                                  DEFINE XX(I),I=1,...,NXX
      NXX = 2*(IB-IA)+3
      WK(1) = AA
      WK(NXX) = BB
      IF (NXX.NE.3) GO TO 45
      WK(2) = (AA+BB)*HALF
      GO TO 55
   45 IAP1 = IA+1
      WK(2) = (AA+X(IAP1))*HALF
      WK(3) = X(IAP1)
      WK(NXX-1) = (X(IB)+BB)*HALF
      IF (NXX.EQ.5) GO TO 55
      NWTXX = NXX-3
      II = IAP1
      DO 50 I=4,NWTXX,2
         WK(I) = (X(II)+X(II+1))*HALF
         II = II+1
         WK(I+1) = X(II)
   50 CONTINUE
C                                  COMPUTE WTXX(I),I=1,...,NXX
   55 WK(NXX+1) = (WK(3)-WK(2))*THIRD
      WK(NXX+2) = FOUR*WK(NXX+1)
      WK(NXX+NXX) = (WK(NXX)-WK(NXX-1))*THIRD
      IF (NXX.EQ.3) GO TO 65
      NWTXX = NXX-2
      DO 60 I=3,NWTXX,2
         IP1 = I+1
         WK(NXX+I) = (WK(IP1)-WK(I-1))*THIRD
         WK(NXX+IP1) = FRDTRD*(WK(IP1)-WK(I))
   60 CONTINUE
   65 CC = AMIN1(C,D)
      DD = AMAX1(C,D)
C                                  DEFINE YL(I),I=1,...,NYL
      NYL = ID-IC+3
C                                  NYLS AND NYLSP ARE THE STARTING AND
C                                    ENDING LOCATIONS OF YL
      NYLS = NXX+NXX+1
      NYLSP = NYLS+NYL-1
      WK(NYLS) = CC
      WK(NYLSP) = DD
      IF (NYL.NE.3) GO TO 70
      WK(NYLS+1) = (CC+DD)*HALF
      GO TO 80
   70 ICP1 = IC+1
      WK(NYLS+1) = (CC+Y(ICP1))*HALF
      WK(NYLSP-1) = (Y(ID)+DD)*HALF
      IF (NYL.EQ.4) GO TO 80
      IDM1 = ID-1
      IYL = NYLS+2
      DO 75 I=ICP1,IDM1
         WK(IYL) = (Y(I)+Y(I+1))*HALF
         IYL = IYL+1
   75 CONTINUE
C                                  INTERPOLATE AT YL(I),J=1,...,NYL FOR
C                                    EACH X(I),I=1,...,NX. NFL, NFI,
C                                    NFCOEF, AND NFTEMP ARE THE STARTING
C                                    LOCATIONS OF FL, FI, COEF, AND
C                                    TEMP.
   80 NFL = NYLS+NYL
      NFI = NFL+NX*NYL
      NFI1 = NFI-1
      NCOEF = NFI+NY
      NFTEMP = NCOEF+NYM1*3
      NFTMP1 = NFTEMP+NYL-1
      DO 95 I=1,NX
C                                  MOVE F(I,J),J=1,...,NY TO FI(J)
         DO 85 J=1,NY
            WK(NFI1+J) = F(I,J)
   85    CONTINUE
C                                  INTERPOLATE
         CALL ICSCCU (Y,WK(NFI),NY,WK(NCOEF),NYM1,JJER)
C                                  EVALUATE
         CALL ICSEVU (Y,WK(NFI),NY,WK(NCOEF),NYM1,WK(NYLS),WK(NFTEMP),
     1   NYL,JJER)
C                                  MOVE FTEMP(J),J=1,...,NYL TO FL(I,J)
         NN = NYLSP+I
         DO 90 J=NFTEMP,NFTMP1
            WK(NN) = WK(J)
            NN = NN+NX
   90    CONTINUE
   95 CONTINUE
C                                  INTEGRATE
      Q = ZERO
C                                  RECOMPUTE NCOEF AND NFTEMP TO SAVE
C                                    STORAGE
      NCOEF = NFL+NX*NYL
      NFTEMP = NCOEF+NXM1*3
      NFTMP1 = NFTEMP-1-NXX
C                                  NWTXX AND NN POINT TO THE STARTING
C                                    AND ENDING LOCATIONS OF WTXX
      NWTXX = NXX+1
      NN = NXX+NXX
      NFI = NFL+NX
      NFI1 = NFL+(NYL-1)*NX
C                                  INTERPOLATE AT XX(I),I=1,NXX FOR EACH
C                                    YY WHERE YY IS THE UNION OF THE
C                                    SETS YL AND Y(J), WHERE
C                                    IC .LT. J .LT. ID
C                                    YL(1), YL(2), YL(NYL) ARE SPECIAL
C                                    CASES
C                                  COMPUTE THE WEIGHT OF YL(1)
      WTYY = (WK(NYLS+1)-WK(NYLS))*THIRD
      I = NFL
  100 SUM = ZERO
C                                  INTERPOLATE
      CALL ICSCCU (X,WK(I),NX,WK(NCOEF),NXM1,JJER)
C                                  EVALUATE
      CALL ICSEVU (X,WK(I),NX,WK(NCOEF),NXM1,WK(1),WK(NFTEMP),NXX,JJER)
C                                  COMPUTE THE SUM FOR THE INTERVAL
      DO 105 J=NWTXX,NN
         SUM = SUM+WK(J)*WK(NFTMP1+J)
  105 CONTINUE
C                                  ACCUMULATE THE INTEGRAL
      Q = Q+WTYY*SUM
      IF (I.EQ.NFI1) GO TO 115
      IF (I.EQ.NFI) GO TO 110
C                                  COMPUTE THE WEIGHT FOR YL(2)
      WTYY = FOUR*WTYY
      I = NFI
      GO TO 100
C                                  COMPUTE THE WEIGHT FOR YL(NYL)
  110 WTYY = (WK(NYLSP)-WK(NYLSP-1))*THIRD
      I = NFI1
      GO TO 100
  115 IF (NYL.EQ.3) GO TO 135
      II = ICP1
      NFL = NFI
      NFI = NYLS-1
      NFI1 = NYL-2
      DO 130 I=2,NFI1
         SUM = ZERO
C                                  INTERPOLATE
         CALL ICSCCU (X,F(1,II),NX,WK(NCOEF),NXM1,JJER)
C                                  EVALUATE
         CALL ICSEVU (X,F(1,II),NX,WK(NCOEF),NXM1,WK(1),WK(NFTEMP),NXX,
     1   JJER)
C                                  COMPUTE THE SUM FOR THE INTERVAL
         DO 120 J=NWTXX,NN
            SUM = SUM+WK(J)*WK(NFTMP1+J)
  120    CONTINUE
C                                  COMPUTE THE WEIGHT FOR Y(II)
         WTYY = (WK(NYLS+I)-WK(NFI+I))*THIRD
C                                  ACCUMULATE THE INTEGRAL
         Q = Q+WTYY*SUM
         SUM = ZERO
C                                  INTERPOLATE
         NFL = NFL+NX
         CALL ICSCCU (X,WK(NFL),NX,WK(NCOEF),NXM1,JJER)
C                                  EVALUATE
         CALL ICSEVU (X,WK(NFL),NX,WK(NCOEF),NXM1,WK(1),WK(NFTEMP),NXX,
     1   JJER)
C                                  COMPUTE THE SUM FOR THE INTERVAL
         DO 125 J=NWTXX,NN
            SUM = SUM+WK(J)*WK(NFTMP1+J)
  125    CONTINUE
C                                  COMPUTE THE WEIGHT FOR YL(I+1)
         WTYY = FRDTRD*(WK(NYLS+I)-Y(II))
C                                  ACCUMULATE THE INTEGRAL
         Q = Q+WTYY*SUM
         II = II+1
  130 CONTINUE
C                                  IF THE LIMITS OF INTEGRATION ARE
C                                    REVERSED, CHANGE THE SIGN OF THE
C                                    INTEGRAL
  135 IF (B.LT.A) Q = -Q
      IF (D.LT.C) Q = -Q
c      CALL UERSET (LEVOLD,LEVEL)
  140 IER = MAX0(JER,KER,LER,MER,NER)
 9000 CONTINUE
c      IF (NER.NE.0) CALL UERTST (NER,'DBCQDU')
c      IF (NER.NE.0) RETURN
c      IF (JER.NE.0) CALL UERTST (JER,'DBCQDU')
c      IF (KER.NE.0) CALL UERTST (KER,'DBCQDU')
c      IF (LER.NE.0) CALL UERTST (LER,'DBCQDU')
c      IF (MER.NE.0) CALL UERTST (MER,'DBCQDU')
 9005 RETURN
      END
c
C   IMSL ROUTINE NAME   - ICSEVU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EVALUATION OF A CUBIC SPLINE
C
C   USAGE               - CALL ICSEVU(X,Y,NX,C,IC,U,S,M,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSAE
C                           OF THE NX DATA POINTS (X(I),Y(I)) I=1,...,
C                           NX (INPUT). X MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                Y      - VECTOR OF LENGTH NX CONTAINING THE ORDINATES
C                           (OR FUNCTION VALUES) OF THE NX DATA POINTS
C                           (INPUT).
C                NX     - NUMBER OF ELEMENTS IN X AND Y (INPUT).
C                           NX MUST BE .GE. 2.
C                C      - SPLINE COEFFICIENTS (INPUT). C IS AN NX-1 BY
C                           3 MATRIX.
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM (INPUT).
C                           IC MUST BE .GE. NX-1.
C                U      - VECTOR OF LENGTH M CONTAINING THE ABSCISSAE
C                           OF THE M POINTS AT WHICH THE CUBIC SPLINE
C                           IS TO BE EVALUATED (INPUT).
C                S      - VECTOR OF LENGTH M (OUTPUT).
C                           THE VALUE OF THE SPLINE APPROXIMATION AT
C                           U(I) IS
C                           S(I) = ((C(J,3)*D+C(J,2))*D+C(J,1))*D+Y(J)
C                           WHERE X(J) .LE. U(I) .LT. X(J+1) AND
C                           D = U(I)-X(J).
C                M      - NUMBER OF ELEMENTS IN U AND S (INPUT).
C                IER    - ERROR PARAMETER (OUTPUT).
C                         WARNING ERROR
C                           IER = 33, U(I) IS LESS THAN X(1).
C                           IER = 34, U(I) IS GREATER THAN X(NX).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE ROUTINE ASSUMES THAT THE ABSCISSAE OF THE NX
C                DATA POINTS ARE ORDERED SUCH THAT X(I) IS LESS THAN
C                X(I+1) FOR I=1,...,NX-1. NO CHECK OF THIS CONDITION
C                IS MADE IN THE ROUTINE. UNORDERED ABSCISSAE WILL CAUSE
C                THE ALGORITHM TO PRODUCE INCORRECT RESULTS.
C            2.  THE ROUTINE GENERATES TWO WARNING ERRORS. ONE ERROR
C                OCCURS IF U(I) IS LESS THAN X(1), FOR SOME I IN THE
C                THE INTERVAL (1,M) INCLUSIVELY. THE OTHER ERROR OCCURS
C                IF U(I) IS GREATER THAN X(NX), FOR SOME I IN THE
C                INTERVAL (1,M) INCLUSIVELY.
C            3.  THE ORDINATE Y(NX) IS NOT USED BY THE ROUTINE. FOR
C                U(K) .GT. X(NX-1), THE VALUE OF THE SPLINE, S(K), IS
C                GIVEN BY
C                 S(K)=((C(NX-1,3)*D+C(NX-1,2))*D+C(NX-1,1))*D+Y(NX-1)
C                WHERE D=U(K)-X(NX-1).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSEVU  (X,Y,NX,C,IC,U,S,M,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,M,IER
      REAL               X(NX),Y(NX),C(IC,3),U(M),S(M)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JER,KER,NXM1,K
      REAL               D,DD,ZERO
      DATA               I/1/,ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      JER = 0
      KER = 0
      IF (M .LE. 0) GO TO 9005
      NXM1 = NX-1
      IF (I .GT. NXM1) I = 1
C                                  EVALUATE SPLINE AT M POINTS
      DO 40 K=1,M
C                                  FIND THE PROPER INTERVAL
         D = U(K)-X(I)
         IF (D) 5,25,15
    5    IF (I .EQ. 1) GO TO 30
         I = I-1
         D = U(K)-X(I)
         IF (D) 5,25,20
   10    I = I+1
         D = DD
   15    IF (I .GE. NX) GO TO 35
         DD = U(K)-X(I+1)
         IF (DD .GE. ZERO) GO TO 10
         IF (D .EQ. 0.0) GO TO 25
C                                  PERFORM EVALUATION
   20    S(K) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
         GO TO 40
   25    S(K) = Y(I)
         GO TO 40
C                                  WARNING - U(I) .LT. X(1)
   30    JER = 33
         GO TO 20
C                                  IF U(I) .GT. X(NX) - WARNING
   35    IF (DD .GT. ZERO) KER = 34
         D = U(K)-X(NXM1)
         I = NXM1
         GO TO 20
   40 CONTINUE
      IER = MAX0(JER,KER)
 9000 CONTINUE
c     IF (JER .GT. 0) CALL UERTST(JER,'ICSEVU')
c     IF (KER .GT. 0) CALL UERTST(KER,'ICSEVU')
 9005 RETURN
      END
c
c
C   IMSL ROUTINE NAME   - ICSCCU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CUBIC SPLINE INTERPOLATION
C                           (EASY-TO-USE VERSION)
C
C   USAGE               - CALL ICSCCU (X,Y,NX,C,IC,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSAE
C                           OF THE NX DATA POINTS (X(I),Y(I)) I=1,...,
C                           NX. (INPUT) X MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                Y      - VECTOR OF LENGTH NX CONTAINING THE ORDINATES
C                           (OR FUNCTION VALUES) OF THE NX DATA POINTS.
C                           (INPUT)
C                NX     - NUMBER OF ELEMENTS IN X AND Y. (INPUT) NX
C                           MUST BE .GE. 2.
C                C      - SPLINE COEFFICIENTS. (OUTPUT) C IS AN NX-1 BY
C                           3 MATRIX. THE VALUE OF THE SPLINE
C                           APPROXIMATION AT T IS
C                           S(T) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C                           WHERE X(I) .LE. T .LT. X(I+1) AND
C                           D = T-X(I).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. (INPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IC IS LESS THAN NX-1.
C                           IER = 130, NX IS LESS THAN 2.
C                           IER = 131, INPUT ABSCISSA ARE NOT ORDERED
C                             SO THAT X(1) .LT. X(2) ... .LT. X(NX).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSCCU (X,Y,NX,C,IC,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,IER
      REAL               X(NX),Y(NX),C(IC,3)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IM1,I,JJ,J,MM1,MP1,M,NM1,NM2
      REAL               DIVDF1,DIVDF3,DTAU,G,CNX(3)
C                                  FIRST EXECUTABLE STATEMENT
      NM1 = NX-1
      IER = 129
      IF (IC .LT. NM1) GO TO 9000
      IER = 130
      IF (NX .LT. 2) GO TO 9000
      IER = 131
      IF (NX .EQ. 2) GO TO 45
C                                  COMPUTE NOT-A-KNOT SPLINE
      DO 5 M = 2,NM1
         MM1=M-1
         C(M,2) = X(M)-X(MM1)
         IF (C(M,2).LE.0.0) GO TO 9000
         C(M,3) = (Y(M)-Y(MM1))/C(M,2)
    5 CONTINUE
      CNX(2) = X(NX)-X(NM1)
      IF (CNX(2).LE.0.0) GO TO 9000
      CNX(3) = (Y(NX)-Y(NM1))/CNX(2)
      IER = 0
      NM2 = NX-2
      IF (NX .GT. 3) GO TO 10
      C(1,3) = CNX(2)
      C(1,2) = C(2,2)+CNX(2)
      C(1,1) = ((C(2,2)+2.*C(1,2))*C(2,3)*CNX(2)+C(2,2)**2*CNX(3))
     1/C(1,2)
      GO TO 20
   10 C(1,3) = C(3,2)
      C(1,2) = C(2,2)+C(3,2)
      C(1,1) = ((C(2,2)+2.*C(1,2))*C(2,3)*C(3,2)+C(2,2)**2*C(3,3))
     1/C(1,2)
      DO 15 M=2,NM2
         MP1=M+1
         MM1=M-1
         G = -C(MP1,2)/C(MM1,3)
         C(M,1) = G*C(MM1,1)+3.*C(M,2)*C(MP1,3)+3.*C(MP1,2)*C(M,3)
         C(M,3) = G*C(MM1,2)+2.*C(M,2)+2.*C(MP1,2)
   15 CONTINUE
   20 G = -CNX(2)/C(NM2,3)
      C(NM1,1) = G*C(NM2,1)+3.*C(NM1,2)*CNX(3)+3.*CNX(2)*C(NM1,3)
      C(NM1,3) = G*C(NM2,2)+2.*C(NM1,2)+2.*CNX(2)
      IF (NX.GT.3) GO TO 25
      CNX(1)=2.*CNX(3)
      CNX(3)=1.
      G=-1./C(NM1,3)
      GO TO 30
   25 G = C(NM1,2)+CNX(2)
      CNX(1) = ((CNX(2)+2.*G)*CNX(3)*C(NM1,2)+CNX(2)**2*
     1(Y(NM1)-Y(NX-2))/C(NM1,2))/G
      G = -G/C(NM1,3)
      CNX(3) = C(NM1,2)
   30 CNX(3) = G*C(NM1,2)+CNX(3)
      CNX(1) = (G*C(NM1,1)+CNX(1))/CNX(3)
      C(NM1,1) = (C(NM1,1)-C(NM1,2)*CNX(1))/C(NM1,3)
      DO 35 JJ=1,NM2
         J = NM1-JJ
         C(J,1) = (C(J,1)-C(J,2)*C(J+1,1))/C(J,3)
   35 CONTINUE
      DO 40 I=2,NM1
         IM1 = I-1
         DTAU = C(I,2)
         DIVDF1 = (Y(I)-Y(IM1))/DTAU
         DIVDF3 = C(IM1,1)+C(I,1)-2.*DIVDF1
         C(IM1,2) = (DIVDF1-C(IM1,1)-DIVDF3)/DTAU
         C(IM1,3) = DIVDF3/DTAU**2
   40 CONTINUE
      DTAU = CNX(2)
      DIVDF1 = (Y(NX)-Y(NM1))/DTAU
      DIVDF3 = C(NM1,1)+CNX(1)-2.*DIVDF1
      C(NM1,2) = (DIVDF1-C(NM1,1)-DIVDF3)/DTAU
      C(NM1,3) = DIVDF3/DTAU**2
      GO TO 9005
   45 IF (X(1) .GE. X(2)) GO TO 9000
      IER = 0
      C(1,1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1,2) = 0.0
      C(1,3) = 0.0
      GO TO 9005
 9000 CONTINUE
c      CALL UERTST(IER,'ICSCCU')
 9005 RETURN
      END
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      complex a(np,np), b(np), sum
      dimension indx(np) 
      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II.NE.0)THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            enddo
          ELSE IF (cabs(SUM).NE.0.0) THEN
            II=I
         ENDIF
         B(I)=SUM
      enddo
      DO I=N,1,-1
         SUM=B(I)
         DO J=I+1,N
            SUM=SUM-A(I,J)*B(J)
         enddo
         B(I)=SUM/A(I,I)
      enddo
      RETURN
      END
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      PARAMETER (NMAX=2000,TINY=0.)
      complex a(np,np), sum, cdum
      dimension indx(np), vv(nmax) 
      D=1.0
      DO I=1,N
         AAMAX=0.0
         DO J=1,N
            IF (cABS(A(I,J)).GT.AAMAX) AAMAX=cABS(A(I,J))
         enddo
         IF (AAMAX.EQ.0.0) PAUSE 'Singular matrix.'
         VV(I)=1./AAMAX
      enddo
      DO J=1,N
         DO I=1,J-1
            SUM=A(I,J)
            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            enddo
            A(I,J)=SUM
         enddo
         AAMAX=0.0
         DO I=J,N
            SUM=A(I,J)
            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            enddo
            A(I,J)=SUM
            DUM=VV(I)*cABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            ENDIF
         enddo
         IF (J.NE.IMAX)THEN
            DO K=1,N
               cDUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=cDUM
            enddo
            D=-D
            VV(IMAX)=VV(J)
         ENDIF
         INDX(J)=IMAX
         IF (cabs(A(j,j)).EQ.0.0)A(N,N)=cmplx(TINY,tiny)
         IF (J.NE.N)THEN
            cDUM=1./A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*cDUM
            enddo
         endif
      enddo
      RETURN
      END
