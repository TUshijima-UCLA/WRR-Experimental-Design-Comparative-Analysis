C
C**************************  OPENIO ************************************
C
C  open the input and output files
C
C***********************************************************************
C
      SUBROUTINE OPENIO(ITERM)
C
      IMPLICIT NONE
      INTEGER  ITERM
      CHARACTER*30 FN2,FN4,FN50,FN60,FN8,FN9,FN10,FN11,FN12
      CHARACTER*30 FN13,FN14,FN15,FN20,FN21,FN35,FN36,FN37,FN61
      CHARACTER*30 FN38,FN40,FN41,FN42,FN43,FN51,FN52
C
C  call the exception handler routine for UNIX (AIX) only
C
C+      INCLUDE 'fexcp.h'
C+      CALL SIGNAL(SIGTRAP,xl_ _trce) 
C
C  open the I/O files. Define unit 4 as terminal output 
C  for the case ITERM .ne. 0 (OPNTER is for UNIX (AIX) - may need to
C  be modified for other systems (see opnter.f in 'utils' directory);
C  OPEN is for VMS).
C

      OPEN(1,FILE='sat2d.fnames',STATUS='OLD')
      READ(1,'(I6)') ITERM
      READ(1,'(A30)') FN50
      READ(1,'(A30)') FN8
      READ(1,'(A30)') FN9
      READ(1,'(A30)') FN10
      READ(1,'(A30)') FN11
      READ(1,'(A30)') FN12
      READ(1,'(A30)') FN60
C Tim edit: output file for steadystate time
      READ(1,'(A30)') FN61
C END TIM EDIT
      READ(1,'(A30)') FN15
      READ(1,'(A30)') FN21
      READ(1,'(A30)') FN36
      READ(1,'(A30)') FN40
      READ(1,'(A30)') FN41
      READ(1,'(A30)') FN42
      READ(1,'(A30)') FN51
      READ(1,'(A30)') FN52
C      print *, FN52
      IF (ITERM .EQ. 0) THEN
         READ(1,'(A30)') FN4
      END IF
      OPEN(50, FILE=FN50, STATUS='OLD')
      OPEN(8, FILE=FN8, STATUS='OLD')
      OPEN(9, FILE=FN9, STATUS='OLD')
      OPEN(10,FILE=FN10,STATUS='OLD')
      OPEN(11,FILE=FN11,STATUS='OLD')
      OPEN(12,FILE=FN12,STATUS='OLD')
      OPEN(60, FILE=FN60, STATUS='UNKNOWN')
C Tim edit: output file for steadystate time
      OPEN(61, FILE=FN61, STATUS='UNKNOWN')
C END TIM EDIT
      OPEN(15,FILE=FN15,STATUS='UNKNOWN')
      OPEN(21,FILE=FN21,STATUS='UNKNOWN')
      OPEN(36,FILE=FN36,STATUS='UNKNOWN')
      OPEN(40,FILE=FN40,STATUS='UNKNOWN')
      OPEN(41,FILE=FN41,STATUS='UNKNOWN')
      OPEN(42,FILE=FN42,STATUS='UNKNOWN')
      OPEN(51,FILE=FN51,STATUS='UNKNOWN')
      OPEN(52,FILE=FN52,STATUS='UNKNOWN')
C      OPEN(52,FILE='output/sys.sat2d',STATUS='UNKNOWN')
      IF (ITERM .EQ. 0) THEN
         OPEN(4,FILE=FN4, STATUS='UNKNOWN')
      ELSE
cc         CALL OPNTER(4)
C+         OPEN(4,FILE='SYS$OUTPUT',STATUS='NEW')
      END IF
C
      RETURN
      END

	  
C
C**************************  NODELT ************************************
C
C  average nodal input array VNOD to obtain values at each 
C  element (output array VELT)
C
C***********************************************************************
C
      SUBROUTINE NODELT(NT,TRIANG,VNOD,VELT)
C
      IMPLICIT  NONE
      INTEGER   J,IEL,INOD
      INTEGER   NT
      INTEGER   TRIANG(4,*)
      REAL*8    R3
      REAL*8    VNOD(*),VELT(*)
C
      R3=1.0D0/3.0D0
      DO IEL=1,NT
         VELT(IEL)=0.0D0
         DO J=1,3
            INOD=TRIANG(J,IEL)
            VELT(IEL)=VELT(IEL) + VNOD(INOD)
         END DO
         VELT(IEL)=VELT(IEL)*R3
      END DO
C
      RETURN
      END

C
C**************************  LOCMAS ************************************
C
C  set up LMASS, the part of the local mass matrix which is constant
C  for all elements (i.e. without the storage coefficient and without
C  the area term)
C  per il caso mass lumping si procede come segue:
C  - si sommano i coefficienti extra diagonali al coeffifiente diagonale
C  - si azzerano tutti i coefficienti extra diagonali
C
C***********************************************************************
C
      SUBROUTINE LOCMAS(LMASS,LUMP)
C
      IMPLICIT  NONE
      INTEGER   I,K
      INTEGER   LUMP
      REAL*8    R3,R6,R12
      REAL*8    LMASS(3,3)
C
      R3=1.0D0/3.0D0
      R6=1.0D0/6.0D0
      R12=1.0D0/12.0D0
      DO K=1,3
         DO I=K,3
            IF (LUMP .EQ. 0) THEN
               IF (K.EQ.I) THEN
                  LMASS(K,I)=R6
               ELSE
                  LMASS(K,I)=R12
               END IF
            ELSE
               IF (K.EQ.I) THEN
                  LMASS(K,I)=R3
               ELSE
                  LMASS(K,I)=0.0D0
               END IF
            END IF
         END DO
      END DO
      DO K=2,3
         DO I=1,K-1
            LMASS(K,I)=LMASS(I,K)
         END DO
      END DO
C
      RETURN
      END


C
C**************************  INIT0R ************************************
C
C  initialize a real array to 0.0
C
C***********************************************************************
C
      SUBROUTINE INIT0R(NUM,RVEC)
C
      IMPLICIT  NONE
      INTEGER   I
      INTEGER   NUM
      REAL*8    RVEC(*)
C
      DO I=1,NUM
         RVEC(I)=0.0D0
      END DO
C
      RETURN
      END

	  
C
C**************************  INIT0I ************************************
C
C  initialize an integer array to 0
C
C***********************************************************************
C
      SUBROUTINE INIT0I(NUM,IVEC)
C
      IMPLICIT  NONE
      INTEGER   I
      INTEGER   NUM
      INTEGER   IVEC(*)
C
      DO I=1,NUM
         IVEC(I)=0
      END DO
C
      RETURN
      END

	  
C
C*************************  DATIN  ************************************
C
C  routine per la lettura dei dati
C
C  unit  5 ---> parameters
C  unit  8 ---> grid info
C  unit  9 ---> Neumann and Dirichlet BC's
C  unit 10 ---> material properties
C  unit 11 ---> initial conditions
C
C  flag IPRT1 for output of data:
C           =0  prints parameters only (default)
C           =1  prints parameters + b.c. + geom. char.
C           =2  prints parameters + b.c. + geom. char. + grid info
C           =3  prints parameters + b.c. + geom. char. + grid info,
C               X, Y coordinate values, and then
C               terminates program execution
C
C***********************************************************************
C
      SUBROUTINE DATIN(ITERM,TRIANG,X,Y,PERMX,PERMY,ELSTOR,SPESS,
     2           CONTP,CONTQ,CONTR,
     3           TIMPRT,PTIMEP,TETAF,DELTAT,DTMAX,TMAX,
     4           DTMAGA,DTMAGM,ITMXCG,TOLCG,LUMP,IRAD,
     5           IPRT1,IPRT,INDP,
     6           NQ,NZONE,N1,NR,NPRT,N,NP,NT)
C
      IMPLICIT  NONE
      INCLUDE  'SAT2D.H'
      INTEGER   I,J,K,K1,K2
      INTEGER   ISTAP
      INTEGER   ITERM,ITMXCG,LUMP,IRAD,IPRT1,IPRT,INDP
      INTEGER   N,NT,NQ,NZONE
      INTEGER   N1,NR,NPRT,NP
      INTEGER   TRIANG(4,*),CONTP(*),CONTQ(*)
      INTEGER   CONTR(*)
      REAL*8    SUMZ
      REAL*8    TETAF,DELTAT,DTMAX,TMAX,DTMAGA,DTMAGM,TOLCG
      REAL*8    X(*),Y(*)
      REAL*8    PERMX(*),PERMY(*)
      REAL*8    ELSTOR(*),SPESS(*)
      REAL*8    TIMPRT(*),PTIMEP(*)
C
C  unit 50 input 
C
C     print *, "unit 50"
      READ(50,*) IPRT1
      READ(50,*) TETAF,LUMP,IRAD
      READ(50,*) ITMXCG,TOLCG
      READ(50,*) DELTAT,DTMAX,TMAX
      READ(50,*) DTMAGA,DTMAGM
c      READ(50,*) IPRT,NPRT,(TIMPRT(I),I=1,NPRT)
      READ(50,*) IPRT
      READ(50,*) NPRT
      READ(50,*) (TIMPRT(I),I=1,NPRT)
C      WRITE(6,*) IPRT
C      WRITE(6,*) NPRT
C      WRITE(6,*) (TIMPRT(I),I=1,NPRT)
      READ(50,*) NR
      IF (NR .NE. 0) READ(50,*)(CONTR(I),I=1,NR)
      WRITE(60,1010) IPRT1,IPRT,NPRT,NR
      WRITE(60,1000) TETAF,LUMP,IRAD
      WRITE(60,1005) ITMXCG,TOLCG
      WRITE(60,1030) DELTAT,DTMAX,TMAX
      WRITE(60,1035) DTMAGA,DTMAGM
C
C  unit 8 input
C


      READ(8,*) NZONE,N1
      READ(8,*) N,NT
      WRITE(60,1020) N,NT,NZONE,N1
C     print *, NT	  
      READ(8,*)((TRIANG(I,K),I=1,4),K=1,NT)
C     print *, "unit 8b"	  
      READ(8,*)(X(K),Y(K),K=1,N)

      IF (IPRT1 .EQ. 3) THEN
C      IF (IPRT1 .GE. 3) THEN
         DO I=1,N
            WRITE(15,1040) I,X(I),Y(I)
         END DO
         DO I=1,NT
            WRITE(15,1050) I,(TRIANG(J,I),J=1,4)
         END DO
         WRITE(60,1060) 
         WRITE(4,1060) 
         CALL CLOSIO(ITERM)
         STOP
      END IF
C
C  unit 12 input
C
C     print *, "unit 12"
      READ(12,*) INDP,ISTAP
      IF (INDP .EQ. 0) THEN
         READ(12,*) PTIMEP(1)
         IF (IPRT1 .GE. 2) WRITE(60,1070) PTIMEP(1)
         DO K=2,N
            PTIMEP(K)=PTIMEP(1)
         END DO
      ELSE
         READ(12,*)(PTIMEP(K),K=1,N)
      END IF
      IF (INDP .NE. 0 .AND. ISTAP .NE. 0 .AND. IPRT1 .GE. 1)
     1                WRITE(60,1090)(K,PTIMEP(K),K=1,N)
C
C  unit 11 input: material properties
C
C     print *, "unit 11"
      WRITE(60,1100)
      DO J=1,NZONE
         READ(11,*) PERMX(J),PERMY(J),ELSTOR(J),SPESS(J)
         IF(IPRT1.GE.2) WRITE(60,1110) J,PERMX(J), PERMY(J),
     1                         ELSTOR(J),SPESS(J)
      END DO
C
C  unit 9 input: Dirichlet B.C.'s
C
C     print *, "unit 9"
      READ(9,*) NP
      WRITE(60,1120) NP
      IF (NP .NE. 0) THEN
         READ(9,*) (CONTP(I),I=1,NP)
         IF (IPRT1 .GE. 1) THEN
            WRITE(60,1130)
            WRITE(60,1140) (CONTP(I),I=1,NP)
         END IF
      END IF
C
C  unit 10 input: Neumann B.C.'s
C
C     print *, "unit 10"
      READ(10,*) NQ
      WRITE(60,1180) NQ
      IF (NQ .NE. 0) THEN
         READ(10,*) (CONTQ(I),I=1,NQ)
         IF (IPRT1 .GE. 1) THEN
            WRITE(60,1190)
            WRITE(60,1140) (CONTQ(I),I=1,NQ)
         END IF
      END IF
C
      WRITE(60,1300) N,NT,NP
      RETURN
C
C  format statements
C
 1000 FORMAT(  5X,'TETAF  (EG: 1 BACKWARD EULER, 0.5 C-N)  = ',1PE15.5,
     1       /,5X,'LUMP   (MASS LUMPING IF NONZERO)        = ',I6,/,
     2       /,5X,'IRAD   (RADIAL FLOW IF NONZERO)         = ',I6)
 1005 FORMAT(  5X,'ITMXCG (MAX ITER FOR CG LINEAR SOLVERS) = ',I6,
     1       /,5X,'TOLCG  (TOLER. FOR CG LINEAR SOLVERS)   = ',1PE15.5)
 1010 FORMAT(  5X,'IPRT1  (FOR OUTPUT OF DATA)             = ',I6,
     1       /,5X,'  =0 prints parameters only (default)',
     2       /,5X,'  =1 prints par. + b.c. + geom. char.',
     3       /,5X,'  =2 prints par. + b.c. + geom. char.',
     4       /,5X,'     + grid info',
     5       /,5X,'  =3 prints par. + b.c. + geom. char.',
     6       /,5X,'     + grid info + X, Y, coordinate',
     7       /,5X,'     and then stops',
     8       /,5X,'IPRT   (FOR DETAILED NODAL OUTPUT)      = ',I6,
     9       /,5X,'  =0 don''t print nodal pot. head ',
     A       /,5X,'     and vel. values',
     B       /,5X,'  =1 print only nodal pot. heads',
     C       /,5X,'  =2 print nodal and el. pot. heads',
     D       /,5X,'  =3 print nodal pot. head and vel.',
     E       /,5X,'  =4 print el. pot. head and vel.',
     F       /,5X,'NPRT   (# OF TIME VALUES FOR DET OUTPUT)= ',I6,
     G       /,5X,'NR     (# OF NODES FOR PARTIAL OUTPUT)  = ',I6)
 1020 FORMAT(/,5X,'N      (NUM. NODI RETICOLO)             = ',I6,
     1       /,5X,'NT     (NUM. TRIANGOLI RETICOLO   )     = ',I6,
     2       /,5X,'NZONE  (NUMERO ZONE (MATERIAL TYPES))   = ',I6,
     3       /,5X,'N1     (NUM. MAX CONTATTI NODALI)       = ',I6)
 1030 FORMAT(/,5X,'DELTAT (INITIAL TIME STEP SIZE)         = ',1PE15.5,
     1       /,5X,'DTMAX  (MAXIMUM TIME STEP SIZE)         = ',1PE15.5,
     2       /,5X,'TMAX   (TIME AT END OF SIMULATION)      = ',1PE15.5)
 1035 FORMAT(  5X,'DTMAGA (MAG. FACTOR FOR DELTAT, ADD.)   = ',1PE15.5,
     1       /,5X,'DTMAGM (MAG. FACTOR FOR DELTAT, MULT.)  = ',1PE15.5)
 1040 FORMAT(I7,2(1PE15.6))
 1050 FORMAT(5I10)
 1060 FORMAT(//,1X,'  IPRT1=3: END OF RUN')
 1070 FORMAT(/,5X,'INITIAL POTENTIAL HEAD (UNIFORM)        = ',1PE15.5)
 1090 FORMAT(/,1X,' INITIAL POTENTIAL HEADS',/,(4(I6,2X,1PE11.3)))
 1100 FORMAT(//,3X,
     1'SATURATED HYDRAULIC CONDUCTIVITY, SPECIFIC STORAGE, AND ',
     2'THICKNESS VALUES',/,1X,
     3'  MAT.TYPE  X-PERM       Y-PERM       STORAGE      THICKNESS')
 1110 FORMAT(1X,I8,2X,4(1PE13.5))
 1120 FORMAT(/,5X,'NP (# OF DIRICHLET NODES)               = ',I6)  
 1130 FORMAT(/,5X,'DIRICHLET NODES')
 1140 FORMAT(1X,10I7)
 1180 FORMAT(/,5X,'NQ   (# OF NEUMANN NODES)               = ',I6)  
 1190 FORMAT(/,5X,'NEUMANN NODES')
 1300 FORMAT(/,5X,'N     (# OF NODES)                      = ',I6,
     1       /,5X,'NT    (# OF TRIANGLES)                  = ',I6,
     2       /,5X,'NP    (TOTAL # OF DIRICHLET NODES)      = ',I6)
      END

	  
C
C*************************   CLOSIO ************************************
C
C  chiude i files di input ed ouput
C
C***********************************************************************
C
      SUBROUTINE CLOSIO(ITERM)
C
      IMPLICIT NONE
      INTEGER  ITERM
C
      CLOSE(3)
      CLOSE(50)
      CLOSE(7)
      CLOSE(9)
      CLOSE(11)
      CLOSE(13)
      CLOSE(15)
      CLOSE(60)
      CLOSE(8)
      IF (ITERM .EQ. 0) THEN
         CLOSE(4)
      END IF
C
      RETURN
      END

	  
C
C**************************  CHKPIC ************************************
C
C  routine per controllare la correttezza della rappresentazione
C  della matrice in forma compatta nel caso simmetrico
C  calcola il numero di elementi diagonali nulli della matrice
C  ritornandolo in NDZ
C
C***********************************************************************
C
      SUBROUTINE CHKPIC(N,NTERM,TOPOL,JA,NDZ,IER)
C
      IMPLICIT NONE
      INTEGER  I,K,IS,IE,JS,JE,IEM1
      INTEGER  N,NTERM,NDZ
      INTEGER  TOPOL(*),IER(*),JA(*)
C
      IER(1)=0
C
C  controlla il numero delle righe
C
      IF (N-2) 11,1,1
C
C  controlla il numero totale di elementi della matrice
C
    1 IF(TOPOL(N+1)-1-NTERM) 12,2,12
C
C  controlla che la topologia della matrice sia corretta .
C  controlla che la lunghezza di ciascuna riga sia corretta .
C  calcola NDZ ,cioe' il numero di elementi diagonali nulli .
C
    2 NDZ=0
      DO 10 I=1,N
         IS=TOPOL(I)
         IE=TOPOL(I+1)-1
         IF(IE-IS) 3,7,5
    3    IF(IE-IS+1) 13,4,24
    4    NDZ=NDZ+1
         GO TO 10
    5    IEM1=IE-1
         DO 6 K=IS,IEM1
C
C  salta se un indice di colonna e' minore od uguale allo
C  indice che nell'ordine lo ha preceduto .
C
            IF(JA(K+1)-JA(K)) 14,14,6
    6    CONTINUE
    7    JS=JA(IS)
         JE=JA(IE)
C
C  salta se il primo indice di colonna della riga i non e' maggiore di I
C
         IF(JS-I) 15,9,8
    8    NDZ=NDZ+1
C
C  salta se l'ultimo indice di colonna della riga I e' maggiore di N
C
    9    IF(JE-N)10,10,16
   10 CONTINUE
      GO TO 20
C
C  imposta gli indicatori di errore
C
   11 IER(1)=1000
      GO TO 20
   12 IER(1)=2000
      GO TO 20
   13 IER(1)=7000
      GO TO 19
   14 IER(1)=6000
      GO TO 18
   15 IER(1)=4000
      GO TO 17
   16 IER(1)=5000
   17 IER(3)=JS
      IER(4)=JE
      GO TO 19
   18 IER(2)=K
      IER(3)=JA(K)
      IER(4)=JA(K+1)
   19 IER(5)=IS
      IER(6)=IE
      IER(7)=I
   20 IF(IER(1))22,23,22
   22 WRITE(60,100)IER(1)
   23 RETURN
C
C  questo errore non dovrebbe accadere .
C
   24 IER(1)=9061
      GO TO 19
C
  100 FORMAT(/,5X,'ROUTINE CHKPIC  CODICE DI ERRORE = ',I5)
      END

	  
C
C**************************  BASIS2 ************************************
C
C  calculate the local basis function coefficients. 
C  Basis function coefficients are divided by 2.
C
C***********************************************************************
C
      SUBROUTINE BASIS2(IP3,TRIANL,X,Z,BIL,CIL)
C
      IMPLICIT  NONE
      INTEGER   J,M,II
      INTEGER   IP3(3,3),TRIANL(4)
      REAL*8    X(*),Z(*),BIL(3),CIL(3)
C
      DO II=1,3
         J=TRIANL(IP3(II,2))
         M=TRIANL(IP3(II,3))
         BIL(II)=(Z(J) - Z(M))*0.5D0
         CIL(II)=(X(M) - X(J))*0.5D0
      END DO
C
      RETURN
      END

	  
C
C**************************  ASSF2D ************************************
C
C  assemble global stiffness and mass matrices from the local
C  contributions of each element: symmetric case, 2D
C
C***********************************************************************
C
      SUBROUTINE ASSF2D(NT,IRAD,TRIANG,TRIJA,LMASS,COEF1,COEF2,
     1                  XC,YC,PERMX,PERMZ,ELSTOR,SPESS,AREA,AREAR,BI,CI)
C
      IMPLICIT NONE
      INTEGER  K,L,IEL,KNOD,LNOD,IND,MTYPE
      INTEGER  NT,IRAD
      INTEGER  TRIANG(4,*),TRIJA(3,3,*)
      REAL*8   AKS,AE,PAK1,PAK2,PAKB,PAKC,RADCF
      REAL*8   ELSTOR(*),SPESS(*),AREA(*),AREAR(*)
      REAL*8   XC(*),YC(*),PERMX(*),PERMZ(*)
      REAL*8   COEF1(*),COEF2(*)
      REAL*8   LMASS(3,3),BI(3,*),CI(3,*)
C
      RADCF=8.0D0*DATAN(1.0D0)
      DO IEL=1,NT
         MTYPE=TRIANG(4,IEL)
         AKS=AREAR(IEL)*SPESS(MTYPE)
         AE=AREA(IEL)*ELSTOR(MTYPE)
         IF(IRAD.EQ.1) THEN
            PAK1=PERMX(MTYPE)*AKS*RADCF*XC(IEL)
            PAK2=PERMZ(MTYPE)*AKS*RADCF*XC(IEL)
            AE=AE*RADCF*XC(IEL)
         ELSE
            PAK1=PERMX(MTYPE)*AKS
            PAK2=PERMZ(MTYPE)*AKS
         END IF
         DO K=1,3
            KNOD=TRIANG(K,IEL)
            PAKB=PAK1*BI(K,IEL)
            PAKC=PAK2*CI(K,IEL)
            DO L=1,3
               LNOD=TRIANG(L,IEL)
               IF (LNOD .GE. KNOD) THEN
                  IND=TRIJA(K,L,IEL)
                  COEF1(IND)=COEF1(IND) +PAKB*BI(L,IEL) +PAKC*CI(L,IEL)
                  COEF2(IND)=COEF2(IND) + AE*LMASS(K,L)
               END IF
            END DO
         END DO
      END DO
C
C      DO K = 1,1005
C         WRITE(6,*) COEF1(K)
C      END DO
      RETURN
      END

	  
C
C**************************  AREBAS ************************************
C
C  calculate the area of each element, the area assigned to each
C  node, and the basis function coefficients. 
C  Basis function coefficients are divided by 2.
C
C***********************************************************************
C
      SUBROUTINE AREBAS(N,NT,TRIANG,IP3,BI,CI,ARENOD,AREA,AREAR,IAREA,
     1                  X,Z,ITERM)
C
      IMPLICIT  NONE
      INTEGER   J,IEL,IA1,IA2,IA3
      INTEGER   N,NT,ITERM
      INTEGER   TRIANG(4,*),IAREA(*),IP3(3,3)
      REAL*8    ARE,ARE3,R3
      REAL*8    BI(3,*),CI(3,*)
      REAL*8    ARENOD(*),AREA(*),AREAR(*),X(*),Z(*)
C
C      DATA    IP3/1,2,3,2,3,1,3,1,2/
      R3=1.0D0/3.0D0
      CALL INIT0R(N,ARENOD)
      DO IEL=1,NT
         CALL AREAS(IP3,TRIANG(1,IEL),X,Z,ARE)
         CALL BASIS2(IP3,TRIANG(1,IEL),X,Z,BI(1,IEL),CI(1,IEL))
         IF (ARE .EQ. 0.0D0) THEN
            WRITE(60,1000) IEL,(TRIANG(J,IEL),J=1,3)
            CALL CLOSIO(ITERM)
            STOP
         END IF
         ARE3=DABS(ARE)*R3
         IA1=TRIANG(1,IEL)
         IA2=TRIANG(2,IEL)
         IA3=TRIANG(3,IEL)
         ARENOD(IA1)=ARENOD(IA1) + ARE3
         ARENOD(IA2)=ARENOD(IA2) + ARE3
         ARENOD(IA3)=ARENOD(IA3) + ARE3
         IAREA(IEL)=1
         IF (ARE .LT. 0.0D0) THEN
            IAREA(IEL)=-1
            ARE=-ARE
         END IF
         AREA(IEL)=ARE
         AREAR(IEL)=1.0D0/ARE
      END DO
C
      RETURN
 1000 FORMAT(/,' ERROR IN SUBROUTINE AREBAS: ZERO AREA CALCULATED',/,
     1         '    AT ELEMENT ',I6,'  NODE NUMBERS:',3I6)
      END

	  
	  
C
C**************************  AREAS  ************************************
C
C  calculate the area of an element
C
C***********************************************************************
C
      SUBROUTINE AREAS(IP3,TRIANL,X,Z,ARE)
C
      IMPLICIT  NONE
      INTEGER   I,J,M,II
      INTEGER   IP3(3,3),TRIANL(4)
      REAL*8    A2,A3
      REAL*8    ARE
      REAL*8    X(*),Z(*)
C
      A2=0.0D0
      A3=0.0D0
      DO II=1,3
         I=TRIANL(IP3(II,1))
         J=TRIANL(IP3(II,2))
         M=TRIANL(IP3(II,3))
         A3=X(I)*Z(J) + A3
         A2=X(I)*Z(M) + A2
      END DO
      ARE=0.5D0*(A3-A2)
C
      RETURN
      END

	  
C
C**************************  TRIPIC ************************************
C
C  set up TRIJA for symmetric case
C
C***********************************************************************
C
      SUBROUTINE TRIPIC(NT,TRIANG,JA,TOPOL,TRIJA)
C
      IMPLICIT NONE
      INTEGER  I,J,II,JJ,IEL,IND
      INTEGER  NT
      INTEGER  TRIANG(4,*),TRIJA(3,3,*),JA(*),TOPOL(*)
C
      DO IEL=1,NT
         DO I=1,3
            II=TRIANG(I,IEL)
            DO J=1,3
               JJ=TRIANG(J,IEL)
               IF (JJ .GE. II) THEN
                  IND=TOPOL(II)-1
  100             IND=IND+1
                  IF (JA(IND) .NE. JJ) GO TO 100
                  TRIJA(I,J,IEL)=IND
               END IF
            END DO
         END DO
      END DO
C
      RETURN
      END

C
C**************************  STRPIC ************************************
C
C  routine per l'analisi topologica del reticolo nel caso simmetrico
C
C***********************************************************************
C
      SUBROUTINE STRPIC(N,NTERM,TRIANG,JA,TOPOL,NT,N1,IMAX,ITERM)
C
      IMPLICIT    NONE
      INTEGER     I,J,K,L,M,III,KK,KKK,MM,MCONTR
      INTEGER     N,N1,NT,NTERM,IMAX,ITERM
      INTEGER     I1(3),I2(3)
      INTEGER     TOPOL(*),JA(*),TRIANG(4,*)
C
C  fissa gli elementi diagonali a distanza costante N1
C        
      NTERM=N1*N
      CALL INIT0I(NTERM,JA)
      JA(1)=1
      DO K=1,NTERM-1
         IF (K/N1*N1 .EQ. K) JA(K+1)=K/N1+1
      END DO
C
C  analizza tutti i triangoli ( NT )
C
      DO 400 K=1,NT
         DO I=1,3
            I2(I)=TRIANG(I,K)
         END DO
C
C  ordina i nodi dell'elemento in senso crescente 
C
         DO J=1,3
            KKK=1
            KK=I2(1)
            DO I=2,3
               IF(I2(I).LT.KK) THEN
                  KK=I2(I)
                  KKK=I
               END IF
            END DO
            I1(J)=KK
            I2(KKK)=IMAX
         END DO
C  
C  genera il vettore JA
C
         DO 6 I=1,2
            J=I+1
            DO 7 L=J,3
               M=N1*(I1(I)-1)+L-J+2
               MCONTR=N1*I1(I)
    9          IF(I1(L)-JA(M))8,7,11
   11          IF(JA(M).EQ.0)GO TO 10
               M=M+1
               IF(M-MCONTR.GE.0)GO TO 99
               GO TO 9
   10          JA(M)=I1(L)
               GO TO 7
    8          MM=M
   18          MM=MM+1
               IF(MM-MCONTR.GE.0) GO TO 99
               IF(JA(MM))18,13,18
   13          JA(MM)=JA(MM-1)
               MM=MM-1
               IF(MM.GT.M)GO TO 13
               JA(M)=I1(L)
    7       CONTINUE
    6    CONTINUE
  400 CONTINUE
C
C  costruisce il vettore TOPOL
C
      TOPOL(1)=1
      M=1
      J=1
      DO 20 K=1,NTERM,N1
         DO 21 I=1,N1
            IF(JA(K+I-1).EQ.0)GO TO 21
            M=M+1
   21    CONTINUE
         J=J+1
         TOPOL(J)=M
   20 CONTINUE
C
C  compatta il vettore JA eliminando gli zeri
C
      M=0
      DO 14 K=1,NTERM
         IF(JA(K).EQ.0) GO TO 14
         M=M+1
         JA(M)=JA(K)
   14 CONTINUE
      NTERM=M
C
      RETURN
   99 III=MCONTR/N1
      WRITE(60,101)III,K
  101 FORMAT(1X,'RIGA = ',I6,' TRIANGOLO = ',I6,'  AUMENTARE N1 ')
      CALL CLOSIO(ITERM)
      STOP
      END
