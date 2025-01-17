!C
!C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!C
!C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2011, ALL RIGHTS RESERVED
!C                BY
!C         DAISUKE TAKAHASHI
!C         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!C         UNIVERSITY OF TSUKUBA
!C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!C         E-MAIL: daisuke@cs.tsukuba.ac.jp
!C
!C
!C     PARALLEL 3-D COMPLEX-TO-REAL-TO FFT ROUTINE
!C
!C     FORTRAN77 + MPI SOURCE PROGRAM
!C
!C     CALL PZDFFT3D(A,B,NX,NY,NZ,ICOMM,ME,NPU,IOPT)
!C
!C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
!C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
!C     NZ IS THE LENGTH OF THE TRANSFORMS IN THE Z-DIRECTION (INTEGER*4)
!C       ------------------------------------
!C         NX = (2**IP) * (3**IQ) * (5**IR)
!C         NY = (2**JP) * (3**JQ) * (5**JR)
!C         NZ = (2**KP) * (3**KQ) * (5**KR)
!C       ------------------------------------
!C     ICOMM IS THE COMMUNICATOR (INTEGER*4)
!C     ME IS THE RANK (INTEGER*4)
!C     NPU IS THE NUMBER OF PROCESSORS (INTEGER*4)
!C     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
!C     IOPT = +1 FOR INVERSE TRANSFORM WHERE
!C              A(NX/2+1,NY,NZ/NPU) IS COMPLEX INPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE A(*,*,BLOCK)
!C              A(NX,NY,NZ/NPU) IS REAL OUTPUT VECTOR (REAL(WP))
!C!HPF$ DISTRIBUTE A(*,*,BLOCK)
!C              B(NX/2+1,NY,NZ/NPU) IS WORK VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE B(*,*,BLOCK)
!C     IOPT = +2 FOR INVERSE TRANSFORM WHERE
!C     ME = 0   A((NX/2)/NPU+1,NY,NZ) IS COMPLEX INPUT VECTOR (COMPLEX(WP))
!C     ME > 0   A((NX/2)/NPU,NY,NZ) IS COMPLEX INPUT VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE A(BLOCK,*,*)
!C              A(NX,NY,NZ/NPU) IS REAL OUTPUT VECTOR (REAL(WP))
!C!HPF$ DISTRIBUTE A(*,*,BLOCK)
!C              B(NX/2+1,NY,NZ/NPU) IS WORK VECTOR (COMPLEX(WP))
!C!HPF$ DISTRIBUTE B(*,*,BLOCK)
!C
!C     WRITTEN BY DAISUKE TAKAHASHI
!C

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

      SUBROUTINE PZDFFT3D(A,B,NX,NY,NZ,ICOMM,ME,NPU,IOPT)
      USE CONSTANTS_MOD
      IMPLICIT REAL(WP) (A-H,O-Z)
      INCLUDE 'ffte_param.h'
      DIMENSION A(*)
      COMPLEX(WP) B(*)
      COMPLEX(WP) C((NDA3+NP)*NBLK),D(NDA3)
      COMPLEX(WP) WX(NDA3),WY(NDA3),WZ(NDA3)
      SAVE WX,WY,WZ
!C
      IF (IOPT .EQ. 0) THEN
        CALL SETTBL(WX,NX)
        CALL SETTBL(WY,NY)
        CALL SETTBL(WZ,NZ)
        RETURN
      END IF
!C
!$omp PARALLEL PRIVATE(C,D)
      CALL PZDFFT3D0(A,A,A,B,B,C,C,C,D,WX,WY,WZ,NX,NY,NZ,ICOMM,ME,NPU, &
                     IOPT)
!$omp END PARALLEL
      RETURN
      END
      SUBROUTINE PZDFFT3D0(DA,A,AXYZ,B,BXYZ,CX,CY,CZ,D,WX,WY,WZ,  &
                           NX,NY,NZ,ICOMM,ME,NPU,IOPT)
      USE CONSTANTS_MOD
      IMPLICIT REAL(WP) (A-H,O-Z)
#ifdef HAVE_MPI_GENESIS
      INCLUDE 'mpif.h'
#endif
      INCLUDE 'ffte_param.h'
      COMPLEX(WP) A(*),AXYZ(NX/2+1,NY,*)
      COMPLEX(WP) B(*),BXYZ(NX/2+1,NY,*)
      COMPLEX(WP) CX(*),CY(NY+NP,*),CZ(NZ+NP,*),D(*)
      COMPLEX(WP) WX(*),WY(*),WZ(*)
      COMPLEX(WP) TEMP
      DIMENSION DA(NX,NY,*)
      DIMENSION ISCNT(MAXNPU),ISDSP(MAXNPU),IRCNT(MAXNPU),IRDSP(MAXNPU)
      DIMENSION LNX(3),LNY(3),LNZ(3)
!C
!     DN=1.0_WP/(REAL(NX,WP)*REAL(NY,WP)*REAL(NZ,WP))
      DN=1.0_WP
!C
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
      CALL FACTOR(NZ,LNZ)
!C
      NNX=NX/NPU
      NNZ=NZ/NPU
!C
      ISCNT(1)=(NNX/2+1)*NY*NNZ
      ISDSP(1)=0
      DO 10 I=2,NPU
        ISCNT(I)=(NNX/2)*NY*NNZ
        ISDSP(I)=ISDSP(I-1)+ISCNT(I-1)
   10 CONTINUE
      IF (ME .EQ. 0) THEN
        IRCNT(1)=(NNX/2+1)*NY*NNZ
        IRDSP(1)=0
        DO 20 I=2,NPU
          IRCNT(I)=(NNX/2+1)*NY*NNZ
          IRDSP(I)=IRDSP(I-1)+IRCNT(I-1)
   20   CONTINUE
      ELSE
        IRCNT(1)=(NNX/2)*NY*NNZ
        IRDSP(1)=0
        DO 30 I=2,NPU
          IRCNT(I)=(NNX/2)*NY*NNZ
          IRDSP(I)=IRDSP(I-1)+IRCNT(I-1)
   30   CONTINUE
      END IF
!C
      IF (IOPT .EQ. 1) THEN
!$omp DO
        DO 80 K=1,NNZ
          DO 70 J=1,NY
!!DIR$ VECTOR ALIGNED
            DO 40 I=1,NNX/2+1
              B(I+(J-1)*(NNX/2+1)+(K-1)*(NNX/2+1)*NY)=AXYZ(I,J,K)
   40       CONTINUE
            DO 60 L=2,NPU
!!DIR$ VECTOR ALIGNED
              DO 50 I=1,NNX/2
                B(I+(J-1)*(NNX/2)+(K-1)*(NNX/2)*NY     &
                  +((L-2)*(NNX/2)+(NNX/2+1))*NY*NNZ)   &
               =AXYZ(I+((L-2)*(NNX/2)+(NNX/2+1)),J,K)
   50         CONTINUE
   60       CONTINUE
   70     CONTINUE
   80   CONTINUE
!$omp BARRIER
!$omp MASTER
#ifdef HAVE_MPI_GENESIS
        CALL MPI_ALLTOALLV(B,ISCNT,ISDSP,MPI_WP_COMPLEX,  &
                           A,IRCNT,IRDSP,MPI_WP_COMPLEX,  &
                           ICOMM,IERR)
#else
        DO I=1,ISCNT(1)
           A(I) = B(I)
        ENDDO
#endif
!$omp END MASTER
!$omp BARRIER
      END IF
      IF (ME .EQ. 0) THEN
!$omp DO
        DO 160 J=1,NY
          DO 150 II=1,NNX/2+1,NBLK
            DO 100 I=II,MIN0(II+NBLK-1,NNX/2+1)
!!DIR$ VECTOR ALIGNED
              DO 90 K=1,NZ
                CZ(K,I-II+1)                                     &
               =CONJG(A(I+(J-1)*(NNX/2+1)+(K-1)*(NNX/2+1)*NY))
   90         CONTINUE
  100       CONTINUE
            DO 110 I=II,MIN0(II+NBLK-1,NNX/2+1)
              CALL FFT235(CZ(1,I-II+1),D,WZ,NZ,LNZ)
  110       CONTINUE
            DO 140 L=1,NPU
              DO 130 K=1,NNZ
!!DIR$ VECTOR ALIGNED
                DO 120 I=II,MIN0(II+NBLK-1,NNX/2+1)
                  B(I+(J-1)*(NNX/2+1)+(K-1)*(NNX/2+1)*NY    &
                    +(L-1)*(NNX/2+1)*NY*NNZ)                &
                 =CZ(K+(L-1)*NNZ,I-II+1)
  120           CONTINUE
  130         CONTINUE
  140       CONTINUE
  150     CONTINUE
  160   CONTINUE
      ELSE
!$omp DO
        DO 240 J=1,NY
          DO 230 II=1,NNX/2,NBLK
            DO 180 I=II,MIN0(II+NBLK-1,NNX/2)
!!DIR$ VECTOR ALIGNED
              DO 170 K=1,NZ
                CZ(K,I-II+1)                                 &
               =CONJG(A(I+(J-1)*(NNX/2)+(K-1)*(NNX/2)*NY))
  170         CONTINUE
  180       CONTINUE
            DO 190 I=II,MIN0(II+NBLK-1,NNX/2)
              CALL FFT235(CZ(1,I-II+1),D,WZ,NZ,LNZ)
  190       CONTINUE
            DO 220 L=1,NPU
              DO 210 K=1,NNZ
!!DIR$ VECTOR ALIGNED
                DO 200 I=II,MIN0(II+NBLK-1,NNX/2)
                  B(I+(J-1)*(NNX/2)+(K-1)*(NNX/2)*NY       &
                    +(L-1)*(NNX/2)*NY*NNZ)                 &
                 =CZ(K+(L-1)*NNZ,I-II+1)
  200           CONTINUE
  210         CONTINUE
  220       CONTINUE
  230     CONTINUE
  240   CONTINUE
      END IF
!$omp BARRIER
!$omp MASTER
#ifdef HAVE_MPI_GENESIS
      CALL MPI_ALLTOALLV(B,IRCNT,IRDSP,MPI_WP_COMPLEX,   &
                         A,ISCNT,ISDSP,MPI_WP_COMPLEX,   &
                         ICOMM,IERR)
#else
      DO I=1,ISCNT(1)
         A(I)=B(I)
      ENDDO
#endif
!$omp END MASTER
!$omp BARRIER
!$omp DO
      DO 380 K=1,NNZ
        DO 300 II=1,NNX/2+1,NBLK
          DO 260 I=II,MIN0(II+NBLK-1,NNX/2+1)
!!DIR$ VECTOR ALIGNED
            DO 250 J=1,NY
              CY(J,I-II+1)=A(I+(J-1)*(NNX/2+1)+(K-1)*(NNX/2+1)*NY)
  250       CONTINUE
  260     CONTINUE
          DO 270 I=II,MIN0(II+NBLK-1,NNX/2+1)
            CALL FFT235(CY(1,I-II+1),D,WY,NY,LNY)
  270     CONTINUE
          DO 290 J=1,NY
!!DIR$ VECTOR ALIGNED
            DO 280 I=II,MIN0(II+NBLK-1,NNX/2+1)
              BXYZ(I,J,K)=CY(J,I-II+1)
  280       CONTINUE
  290     CONTINUE
  300   CONTINUE
        DO 370 L=2,NPU
          DO 360 II=1,NNX/2,NBLK
            DO 320 J=1,NY
!!DIR$ VECTOR ALIGNED
              DO 310 I=II,MIN0(II+NBLK-1,NNX/2)
                CY(J,I-II+1)                          &
               =A(I+(J-1)*(NNX/2)+(K-1)*(NNX/2)*NY    &
                  +((L-2)*(NNX/2)+(NNX/2+1))*NY*NNZ)
  310         CONTINUE
  320       CONTINUE
            DO 330 I=II,MIN0(II+NBLK-1,NNX/2)
              CALL FFT235(CY(1,I-II+1),D,WY,NY,LNY)
  330       CONTINUE
            DO 350 J=1,NY
!!DIR$ VECTOR ALIGNED
              DO 340 I=II,MIN0(II+NBLK-1,NNX/2)
                BXYZ(I+(L-2)*(NNX/2)+(NNX/2+1),J,K)=CY(J,I-II+1)
  340         CONTINUE
  350       CONTINUE
  360     CONTINUE
  370   CONTINUE
  380 CONTINUE
!$omp DO PRIVATE(TEMP)
      DO 470 K=1,NNZ
        IF (MOD(NY,2) .EQ. 0) THEN
          DO 410 J=1,NY,2
            CX(1)=DCMPLX(REAL(BXYZ(1,J,K),WP),REAL(BXYZ(1,J+1,K),WP))
!!DIR$ VECTOR ALIGNED
            DO 390 I=2,NX/2+1
              TEMP=(0.0_WP,1.0_WP)*BXYZ(I,J+1,K)
              CX(I)=BXYZ(I,J,K)+TEMP
              CX(NX-I+2)=CONJG(BXYZ(I,J,K)-TEMP)
  390       CONTINUE
            CALL FFT235(CX,D,WX,NX,LNX)
            DO 400 I=1,NX
              DA(I,J,K)=REAL(CX(I),WP)*DN
              DA(I,J+1,K)=IMAG(CX(I))*DN
  400       CONTINUE
  410     CONTINUE
        ELSE
          DO 440 J=1,NY-1,2
            CX(1)=DCMPLX(REAL(BXYZ(1,J,K),WP),REAL(BXYZ(1,J+1,K),WP))
!!DIR$ VECTOR ALIGNED
            DO 420 I=2,NX/2+1
              TEMP=(0.0_WP,1.0_WP)*BXYZ(I,J+1,K)
              CX(I)=BXYZ(I,J,K)+TEMP
              CX(NX-I+2)=CONJG(BXYZ(I,J,K)-TEMP)
  420       CONTINUE
            CALL FFT235(CX,D,WX,NX,LNX)
            DO 430 I=1,NX
              DA(I,J,K)=REAL(CX(I),WP)*DN
              DA(I,J+1,K)=IMAG(CX(I))*DN
  430       CONTINUE
  440     CONTINUE
        END IF
          CX(1)=DCMPLX(REAL(BXYZ(1,NY,K),WP),0.0_WP)
!!DIR$ VECTOR ALIGNED
          DO 450 I=2,NX/2+1
            CX(I)=BXYZ(I,NY,K)
            CX(NX-I+2)=CONJG(BXYZ(I,NY,K))
  450     CONTINUE
          CALL FFT235(CX,D,WX,NX,LNX)
          DO 460 I=1,NX
            DA(I,NY,K)=REAL(CX(I),WP)*DN
  460     CONTINUE
  470 CONTINUE
      RETURN
      END

