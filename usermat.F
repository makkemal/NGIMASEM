********************************************************************************
********************************************************************************
*ANSYS COMPATIBILITY SECTION STARTS HERE
********************************************************************************
c*******************************************************************************
c
*deck,usermat      USERDISTRIB  parallel                                gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   cutFactor, pVolDer, var2, var3, var4,
     &                   var5, var6, var7)
c*******************************************************************************
c#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ, cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 pVolDer (3),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
      EXTERNAL         usermat3d, usermatps, usermatbm, usermat1d

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
      data             var1/0.0d0/
c
c*************************************************************************
c
      IF(ncomp .GE. 4) THEN
c ***    3d, plane strain and axisymmetric example
         call usermat3d (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
      END IF
      return
      end
c***************************************************************************
*deck,usermat3d    USERDISTRIB  parallel                                gal
      subroutine usermat3d(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of 3D solid elements or plane elements
c       in plane strain or axisymmetric stress state. The plasticity
c       model is the same as TB, BISO.
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variable
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - statev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c      cutFactor(dp,sc,o)                 time step size cut-back factor 
c                                         define it if a smaller step size is wished
c                                         recommended value is 0~1
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
c#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, ONEHALF,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 ONEHALF    = 1.5d0,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 mcomp      = 6
     &                 )
c
      EXTERNAL         vzero, vmove, get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), G(mcomp),
     &                 sigDev(mcomp), JM    (mcomp,mcomp), dfds(mcomp),
     &                 sigi  (mcomp), strainEl(mcomp)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7

      DATA G/1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,0.0D0/
c
      INTEGER          i, j
      DOUBLE PRECISION pEl,   qEl,     pleq_t,  sigy_t , sigy,
     &                 dpleq, pleq, 
     &                 young, posn,    sigy0,   dsigdep, 
     &                 elast1,elast2,
     &                 twoG,  threeG,  oneOv3G, qElOv3G, threeOv2qEl, 
     &                 fratio,  con1,    con2, dperr(3)
c*************************************************************************
c
      INTEGER ndi,nshr,ntens,nstatv,nprops,noel,npt,layer,kspt,
     &  kstep,kinc
!
      REAL*8 statev(nStatev),
     &  ddsdde(ncomp,ncomp),ddsddt(ncomp),drplde(ncomp),
     &  stran(ncomp),dstran(ncomp),celent,
     &  props(nprop),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),
     &  sse,spd,scd,rpl,drpldt,predef,dpred,
     &  pnewdt
c************************************************************************
c
      keycut   = 0
      cutFactor = 0.d0
      dsigdep  = ZERO 
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young"s modulus and Poisson"s ratio, initial yield stress and others
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
      twoG     = young / (ONE+posn)
      threeG   = ONEHALF * twoG
      elast1=young*posn/((1.0D0+posn)*(1.0D0-TWO*posn))
      elast2=HALF*twoG
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = elast2
c	 
      sse = sedEl
      ndi = nDirect
      nshr = nShear
      ntens = ncomp
      nstatv = nStatev
      nprops = nProp
      time = Time
      statev = ustatev
      spd = sedPl
      stran = Strain
      dstran = dStrain
      props = prop
c
      call umat(stress,statev,ddsdde,sse,spd,
     &  ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords)
c
      dsdePL = ddsdde
      sedEl = sse
      sedPl = spd
      epseq = ustatev(1)
      ustatev = statev
      return
      end
	  
*****************
*Start of Wilke's subroutine
*****************	  
	  
	  subroutine umat(stress,statev,ddsdde,sse,spd,
     &  ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords)
!
      implicit none
!
      character*80 cmname
!
      integer ndi,nshr,ntens,nstatv,nprops,noel,npt,layer,kspt,
     &  kstep,kinc
!
      real*8 stress(ntens),statev(nstatv),
     &  ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     &  stran(ntens),dstran(ntens),time,celent,
     &  props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),
     &  sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,predef,dpred,
     &  pnewdt
!
!     IRRADIATED GRAPHITE MATERIAL
!
      integer i,j
      real*8 T_t2,T_t1,F_t2,F_t1,Eth(6),Ew(6)
      real*8 Eel(6),dEeldS(6,6),dEcpdS(6,6)
      real*8 dEcsdS(6,6),stress_t1(6),stress_t2(6),Ecp_t1(6),Ecp_t2(6)
      real*8 Ecs_t1(6),Ecs_t2(6),nrmRes,idRdS(6,6),Res(6),dResdS(6,6)
!
      call GetTempPosTime(T_t2,Coords,time+dtime)
      call GetTempPosTime(T_t1,Coords,time
      call GetFluencePosTime(F_t2,Coords,time+dtime)
      call GetFluencePosTime(F_t1,Coords,time
      call Get_Eth(T_t2,F_t2,Eth)
      call Get_Ew(T_t2,F_t2,Ew)
      Ecp_t1    = statev(1:6)
      Ecs_t1    = statev(7:12)
      stress_t1 = stress
      
      if (time.lt.dtime) then
          Ecp_t1 = 0.d0
          Ecs_t1 = 0.d0
      endif
 
********************************************************************************
********************************************************************************
*ANSYS COMPATIBILITY SECTION ENDS HERE
********************************************************************************
*******************************************************************************
      nrmRes = 1.d0
      stress_t2 = stress_t1
      do while(nrmRes.gt.1.d-12)
      
         call WarningCode(T_t2,F_t2)
         call WarningCode(T_t1,F_t1)
 
         call Get_Eel(T_t2,F_t2,stress_t2,Eel,dEeldS)
         call Get_Ecp(T_t1,F_t1,stress_t1,T_t2,F_t2,stress_t2,
     .             Ecp_t1,Ecp_t2,dEcpdS)
         call Get_Ecs(T_t1,F_t1,stress_t1,T_t2,F_t2,stress_t2,
     .             Ecs_t1,Ecs_t2,dEcsdS)

         Res = stran + dstran - props(1)*Eel - props(4)*Ew - 
     .   props(5)*Eth - props(3)*Ecs_t2 - props(2)*Ecp_t2
        
         nrmRes = 0.d0
         do i=1,6
            nrmRes = nrmRes + Res(i)**2
         enddo
         nrmRes = dsqrt(nrmRes)
         dResdS = props(1)*dEeldS + props(2)*dEcpdS + props(3)*dEcsdS
         call INV_D(dResdS, idRdS)

         do i=1,6
            do j=1,6
               stress_t2(i) = stress_t2(i) + idRdS(i,j)*Res(j)
            enddo
         enddo
      enddo

      ddsdde = idRdS
      stress = stress_t2
      statev(1:6)  = Ecp_t2
      statev(7:12) = Ecs_t2

      return
      end

!*******************************************************************
      subroutine Get_Ecp(T_t1,F_t1,stress_t1,T_t2,F_t2,stress_t2,
     .             Ecp_t1,Ecp_t2,dEcpdS)
!*******************************************************************
! This subroutine computes the primary creep at the end of the time step
! Input :  T_t1      : Temperature at start of the time step
!          T_t2      : Temperature at end of the time step
!          F_t1      : Fluence at start of the time step
!          F_t2      : Fluence at end of the time step
!          stress_t1 : Stress at start of the time step
!          stress_t2 : Stress at end of the time step
!          Ecp_t1    : Primary creep at start of the time step
! Output : Ecp_t2    : Primary creep at end of the stime step
!          dEcpdS    : Derivative of primary creep at end of time step 
!                      wrt stress at end of the time step
      implicit none
      integer i,j,counter
      real*8 T_t1,T_t2,F_t1,F_t2,stress_t1(6),stress_t2(6)
      real*8 Ecp_t1(6),Ecp_t2(6),dEcpdS(6,6)
      real*8 dF,nrmRes,Res(6),nrmRes0
      real*8 dEcp_rate_dS_t1(6,6),dEcp_rate_dS_t2(6,6)
      real*8 dEcp_rate_dE_t1,dEcp_rate_dE_t2,dResdEcp
      real*8 Ecp_rate_t1(6),Ecp_rate_t2(6)

      call GetEcp_rate(T_t1,F_t1,stress_t1,Ecp_t1,Ecp_rate_t1,
     .                 dEcp_rate_dS_t1,dEcp_rate_dE_t1)
      dF = F_t2 - F_t1
! Initial guess for primary creep at end of increment
      Ecp_t2 = Ecp_t1
      call GetEcp_rate(T_t2,F_t2,stress_t2,Ecp_t2,Ecp_rate_t2,
     .                    dEcp_rate_dS_t2,dEcp_rate_dE_t2)
      Res = Ecp_t2 - Ecp_t1 - 0.5d0*dF*(Ecp_rate_t1 + Ecp_rate_t2)
      nrmRes0 = 0.d0
      do i=1,6
          nrmRes0 = nrmRes0 + Res(i)**2
      enddo
      nrmRes0 = dsqrt(nrmRes0)
      counter = 0
      nrmRes = 1.d0
      do while((nrmRes.gt.1.d-10*nrmRes0).and.(counter.lt.20).and.
     .   (nrmRes.gt.1d-12))
         call GetEcp_rate(T_t2,F_t2,stress_t2,Ecp_t2,Ecp_rate_t2,
     .                    dEcp_rate_dS_t2,dEcp_rate_dE_t2)
         Res = Ecp_t2 - Ecp_t1 - 0.5d0*dF*(Ecp_rate_t1 + Ecp_rate_t2)
         nrmRes = 0.d0
         do i=1,6
            nrmRes = nrmRes + Res(i)**2
         enddo
         nrmRes = dsqrt(nrmRes)
         dResdEcp = 1.d0 - 0.5d0*dF*dEcp_rate_dE_t2
         do i=1,6
            Ecp_t2(i) = Ecp_t2(i) - Res(i)/dResdEcp
         enddo
         counter = counter + 1
      enddo
      dEcpdS = 0.5d0*dF*dEcp_rate_dS_t2/dResdEcp

      return
      end
!***********************************************************************
      subroutine Get_Ecs(T_t1,F_t1,stress_t1,T_t2,F_t2,stress_t2,
     .             Ecs_t1,Ecs_t2,dEcsdS)
!***********************************************************************
! This subroutine computes the secondary creep at the end of the time step
! Input :  T_t1      : temperature at start of the time step
!          T_t2      : temperature at end of the time step
!          F_t1      : fluence at start of the time step
!          F_t2      : fluence at end of the time step
!          stress_t1 : stress at start of the time step
!          stress_t2 : stress at end of the time step
!          Ecs_t1    : secondary creep at start of the time step
! Output : Ecs_t2    : secondary creep at end of the stime step
!          dEcsdS    : derivative of secondary creep at end of time step 
!                      wrt stress at end of the time step
      implicit none
      integer i,j
      real*8 T_t1,T_t2,F_t1,F_t2,stress_t1(6),stress_t2(6)
      real*8 Ecs_t1(6),Ecs_t2(6),dEcsdS(6,6),dF
      real*8 Ecs_rate_t1(6),Ecs_rate_t2(6),dEdS_t1(6,6),dEdS_t2(6,6)

      dF = F_t2 - F_t1
      call GetEcs_rate(T_t1,F_t1,stress_t1,Ecs_rate_t1,dEdS_t1)
      call GetEcs_rate(T_t2,F_t2,stress_t2,Ecs_rate_t2,dEdS_t2)
      Ecs_t2 = Ecs_t1 + 0.5d0*dF*(Ecs_rate_t1 + Ecs_rate_t2)
      dEcsdS = 0.5d0*dF*dEdS_t2
      return
      end

!************************************************************************
      subroutine Get_Eel(T,F,stress,Eel,dEeldS)
!************************************************************************
! This subroutine computes the elastic strain at the end of the time step
! Input  : T      : Temperature at end of the time step
!          F      : Fluence at the end of the time step
!          stress : Stress at the end of the time step
! Output : Eel    : Elastic strain at the end of the time step
!          dEeldS : Derivative of the elastic strain at the end of the 
!                   time step wrt the stress at the end of the time step
      implicit none
      real*8 T,F,stress(6),Eel(6),dEeldS(6,6)
      integer i,j
      call Get_invDel(T,F,dEeldS)
      Eel = 0.d0
      do i=1,6
         do j=1,6
            Eel(i) = Eel(i) + dEeldS(i,j)*stress(j)
         enddo
      enddo
      return
      end

!******************************************************************
      subroutine INV_D(A, AINV)
!******************************************************************
!  INV_D  -  Compute the inverse of a 6x6 compliance matrix.
!
!  A       = input 6x6 compliance matrix to be inverted
!  AINV    = output 6x6 inverse of compliance matrix A
!******************************************************************
      IMPLICIT NONE

      real*8, DIMENSION(6,6) :: A
      real*8, DIMENSION(6,6), INTENT(OUT) :: AINV
      real*8, DIMENSION(3,3) :: A1,A1INV,COFACTOR
      real*8 :: DET
      integer :: i,j

      do i=1,3
         do j=1,3
            A1(i,j) = A(i,j)
         enddo
      enddo

      DET =   A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)  
     &      - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1) 
     &      + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1)

      COFACTOR(1,1) = +(A1(2,2)*A1(3,3)-A1(2,3)*A1(3,2))
      COFACTOR(1,2) = -(A1(2,1)*A1(3,3)-A1(2,3)*A1(3,1))
      COFACTOR(1,3) = +(A1(2,1)*A1(3,2)-A1(2,2)*A1(3,1))
      COFACTOR(2,1) = -(A1(1,2)*A1(3,3)-A1(1,3)*A1(3,2))
      COFACTOR(2,2) = +(A1(1,1)*A1(3,3)-A1(1,3)*A1(3,1))
      COFACTOR(2,3) = -(A1(1,1)*A1(3,2)-A1(1,2)*A1(3,1))
      COFACTOR(3,1) = +(A1(1,2)*A1(2,3)-A1(1,3)*A1(2,2))
      COFACTOR(3,2) = -(A1(1,1)*A1(2,3)-A1(1,3)*A1(2,1))
      COFACTOR(3,3) = +(A1(1,1)*A1(2,2)-A1(1,2)*A1(2,1))

      A1INV = TRANSPOSE(COFACTOR) / DET

      AINV = 0.d0
      do i=1,3
         do j=1,3
            AINV(i,j) = A1INV(i,j)
         enddo
      enddo
      AINV(4,4) = 1.d0/A(4,4)
      AINV(5,5) = 1.d0/A(5,5)
      AINV(6,6) = 1.d0/A(6,6)

      RETURN
      END

!***********************************************************************
      subroutine GetEcp_rate(T,F,stress,Ecp,Ecp_rate,dEcp_rate_dS,
     .                       dEcp_rate_dE)
!***********************************************************************
! This subroutine computes the primary creep rate given temperature,
! fluence and stress
! Input  : T            : Temperature
!          F            : Fluence
!          stress       : Stress
!          Ecp          : Primary creep strain
! Output : Ecp_rate     : Primary creep strain rate (wrt fluence)
!          dEcp_rate_dS : Derivative of primary creep strain rate wrt stress
!          dEcp_rate_dE : Scalar derivative of primary creep strain rate
!                         wrt primary creep strain
      implicit none
      integer i,j
      real*8 T,F,stress(6),Ecp(6),Ecp_rate(6),alpha,G0,dEcp_rate_dE
      real*8 iDel(6,6),iDc(6,6),iDc_S(6),EoverE0_ZF,dEcp_rate_dS(6,6)

      call Get_invDel(T,0.d0,iDel)
      iDc = iDel
      stress = stress
      iDc_S = 0.d0
      do i=1,6
         do j=1,6
            iDc_S(i) = iDc_S(i) + iDc(i,j)*stress(j)
         enddo
      enddo
      Ecp_rate = -0.443997461346835d0*Ecp + 2.41778941451711d0*iDc_S
      dEcp_rate_dE = -0.443997461346835d0
      dEcp_rate_dS = 2.41778941451711d0*iDc
      return
      end
!************************************************************************
!************************************************************************
      subroutine GetEcs_rate(T,F,stress,Ecs_rate,dEcsrate_dS)
!************************************************************************
! This subroutine computes the secondary creep rate, given the current
! temperature, fluence and stress
! Input :  T           : Temperature
!          F           : Fluence
!          stress      : Stress
! Output : Ecs_rate    : Secondary creep rate
!          dEcsrate_dS : Derivative of the seondary creep rate wrt stress
      implicit none

      integer i,j
      real*8 T,F,stress(6),Ecs_rate(6)
      real*8 K,SC_a,SC_e,SC_k,beta,iDel(6,6),iDc(6,6),EoverE0_ZF
      real*8 iDc_S(6),dEcsrate_dS(6,6)

      call Get_invDel(T,0.d0,iDel)
      iDc = iDel
      stress = stress

      iDc_S = 0.d0
      do i=1,6
         do j=1,6
            iDc_S(i) = iDc_S(i) + iDc(i,j)*stress(j)
         enddo
      enddo

      Ecs_rate    = 0.81541919d0*iDc_S
      dEcsrate_dS = 0.81541919d0*iDc
      return
      end
      subroutine GetE0per(Value)
!*************************************************************************
! This subroutine constructs a subroutine that returns a constant
! Output : Warning Message Printed to the Screen 
      implicit none
      real*8 Value
      Value = 8918.93750000000d0
      return
      end
      subroutine GetE0par(Value)
!*************************************************************************
! This subroutine constructs a subroutine that returns a constant
! Output : Warning Message Printed to the Screen 
      implicit none
      real*8 Value
      Value = 9167.69230769231d0
      return
      end
      subroutine GetEoverE0per(To,Fo,Answer)
!*************************************************************************
! This subroutine computes the wigner strain in the perpendicular direction
! Input  : To      : Temperature
!          Fo      : Fluence
! Output : dL_per : Wigner strain in perpendicular direction
      implicit none
      integer i,j
      real*8 To,Fo,Answer
      real*8 T,F,Tmin,Tmax,Fmin,Fmax

      Tmin = 300.000000000000d0
      Tmax = 614.000000000000d0
      Fmin = 0.0d0
      Fmax = 47.5200000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = 10**(F*(-2.02919500156835d0*F**3*T - 3.69067427821189d0*F
     @ **3 + 8.05630863443232d0*F**2 - 1.18383567638067d0*F*T**2 -
     @ 7.24272241643419d0*F + 0.26768271083925d0*T + 3.27822398853396d0
     @ ))
      return
      end
      subroutine GetEoverE0par(To,Fo,Answer)
!*************************************************************************
! This subroutine computes the wigner strain in the perpendicular direction
! Input  : To      : Temperature
!          Fo      : Fluence
! Output : dL_per : Wigner strain in perpendicular direction
      implicit none
      integer i,j
      real*8 To,Fo,Answer
      real*8 T,F,Tmin,Tmax,Fmin,Fmax

      Tmin = 298.000000000000d0
      Tmax = 610.000000000000d0
      Fmin = 0.0d0
      Fmax = 47.7600000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = 10**(F*(-1.11508707053732d0*F**3*T - 4.08620884928165d0*F
     @ **3 + 8.9432361189593d0*F**2 - 1.24008647353985d0*F*T**2 -
     @ 7.97562684811928d0*F + 0.042604617722829d0*T + 3.46166687585701d0
     @ ))
      return
      end
      subroutine Get_Wigner_par(To,Fo,Answer)
!*************************************************************************
! This subroutine computes the wigner strain in the perpendicular direction
! Input  : To      : Temperature
!          Fo      : Fluence
! Output : dL_per : Wigner strain in perpendicular direction
      implicit none
      integer i,j
      real*8 To,Fo,Answer
      real*8 T,F,Tmin,Tmax,Fmin,Fmax

      Tmin = 298.000000000000d0
      Tmax = 610.000000000000d0
      Fmin = 0.0d0
      Fmax = 47.7600000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = 0.0153530606937785d0*F**3 + 0.349578701076742d0*F**2*T +
     @ 0.119388454793493d0*F**2 + 0.00599094058166944d0*F*T**2 -
     @ 0.153411284287507d0*F*T - 0.145297822614822d0*F +
     @ 0.0117869475771738d0*T + 0.0031840920729342d0
      return
      end
      subroutine Get_Wigner_per(To,Fo,Answer)
!*************************************************************************
! This subroutine computes the wigner strain in the perpendicular direction
! Input  : To      : Temperature
!          Fo      : Fluence
! Output : dL_per : Wigner strain in perpendicular direction
      implicit none
      integer i,j
      real*8 To,Fo,Answer
      real*8 T,F,Tmin,Tmax,Fmin,Fmax

      Tmin = 299.000000000000d0
      Tmax = 614.000000000000d0
      Fmin = 0.0d0
      Fmax = 47.5200000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = 0.0456699255847065d0*F**3 + 0.513365731196991d0*F**2*T +
     @ 0.0651624013216737d0*F**2 + 0.0208700208980301d0*F*T**2 -
     @ 0.253310245997347d0*F*T - 0.0804745837516608d0*F +
     @ 0.0200651799685113d0*T - 0.00129205950695946d0
      return
      end
      subroutine Get_CTEoCTE0par(To,Fo,Answer)
!*************************************************************************
! This subroutine computes the wigner strain in the perpendicular direction
! Input  : To      : Temperature
!          Fo      : Fluence
! Output : dL_per : Wigner strain in perpendicular direction
      implicit none
      integer i,j
      real*8 To,Fo,Answer
      real*8 T,F,Tmin,Tmax,Fmin,Fmax

      Tmin = 300.000000000000d0
      Tmax = 610.000000000000d0
      Fmin = 0.0d0
      Fmax = 38.4500000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = F*(-12.3629918266927d0*F**3 - 1.81306719440064d0*F**2*T +
     @ 28.2704331281667d0*F**2 + 0.990752925306645d0*F*T**2 +
     @ 2.87895958403812d0*F*T - 20.7106454165527d0*F -
     @ 1.81237949962376d0*T + 4.49674359309073d0)
      return
      end
      subroutine Get_CTEoCTE0per(To,Fo,Answer)
!*************************************************************************
! This subroutine computes the wigner strain in the perpendicular direction
! Input  : To      : Temperature
!          Fo      : Fluence
! Output : dL_per : Wigner strain in perpendicular direction
      implicit none
      integer i,j
      real*8 To,Fo,Answer
      real*8 T,F,Tmin,Tmax,Fmin,Fmax

      Tmin = 300.000000000000d0
      Tmax = 607.000000000000d0
      Fmin = 0.0d0
      Fmax = 38.4500000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = F*(-11.2294477015533d0*F**3 - 0.677256571597844d0*F**2*T
     @ + 25.7212675918032d0*F**2 + 0.866829645371109d0*F*T**2 +
     @ 2.02208474621037d0*F*T - 19.3955229608899d0*F -
     @ 1.82123703042884d0*T + 4.53398659717786d0)
      return
      end
      subroutine Get_CTE0par(To,Answer)
!*************************************************************************
! This subroutine computes the wigner strain in the perpendicular direction
! Input  : To      : Temperature
!          Fo      : Fluence
! Output : dL_per : Wigner strain in perpendicular direction
      implicit none
      integer i,j
      real*8 To,Answer
      real*8 T,F,Tmin,Tmax

      Tmin = 300.000000000000d0
      Tmax = 610.000000000000d0

      T = (To-Tmin)/(Tmax-Tmin)

      Answer = 10**(0.0031430377331164d0*T**2 + 0.0311030269530691d0*T -
     @ 5.369417385761d0)
      return
      end
      subroutine Get_CTE0per(To,Answer)
!*************************************************************************
! This subroutine computes the wigner strain in the perpendicular direction
! Input  : To      : Temperature
!          Fo      : Fluence
! Output : dL_per : Wigner strain in perpendicular direction
      implicit none
      integer i,j
      real*8 To,Answer
      real*8 T,F,Tmin,Tmax

      Tmin = 300.000000000000d0
      Tmax = 607.000000000000d0

      T = (To-Tmin)/(Tmax-Tmin)

      Answer = 10**(-0.0138211647738727d0*T**2 + 0.0469701716881542d0*T
     @ - 5.32664359362176d0)
      return
      end
      subroutine Get_invDel(T,F,iDel)
!*************************************************************************
! This subroutine computes the inverse elasticity tensor given the current
! temperature and fluence
! Input  : T    : Temperature
!          F    : Fluence
! Output : iDel : inverse elasticity tensor
      implicit none
      real*8 iDel(6,6)
      real*8 T,F,c0,c2,c3,T0,b1,b2,b3,b4,b5,b6,b7,b8,b9
      real*8 E10,E20,E30,nu12,nu23,nu31,G120,G230,G310
      real*8 a10,a12,a22,a23,a20,a21,k1,EoverE0per,EoverE0par
      real*8 Eper,Epar

      call GetE0par(Epar)
      call GetEoverE0par(T,F,EoverE0par)
      call GetE0per(Eper)
      call GetEoverE0per(T,F,EoverE0per)
      E10 = EoverE0par*Epar
      E20 = EoverE0par*Epar
      E30 = EoverE0per*Eper
      nu12 = 0.150000000000000d0
      nu23 = 0.150000000000000d0
      nu31 = 0.150000000000000d0
      G120 = (EoverE0par+EoverE0per)/2.d0*2000.00000000000d0
      G230 = (EoverE0par+EoverE0per)/2.d0*2000.00000000000d0
      G310 = (EoverE0par+EoverE0per)/2.d0*2000.00000000000d0

      iDel = 0.d0
      iDel(1,1) = 1.d0/E10
      iDel(2,2) = 1.d0/E20
      iDel(3,3) = 1.d0/E30
      iDel(4,4) = 1.d0/G120
      iDel(5,5) = 1.d0/G310
      iDel(6,6) = 1.d0/G230

      iDel(1,2) = -nu12/E20
      iDel(1,3) = -nu31/E30
      iDel(2,1) = -nu12/E10
      iDel(2,3) = -nu31/E30
      iDel(3,1) = -nu12/E10
      iDel(3,2) = -nu31/E20

      return
      end
      subroutine GetFluencePosTime(F,Coords,time)
!************************************************************************
! This subroutine computes the temperature as a function of time and position
! Input  : Coords   : (x,y,z) coordinate in space
!        : time     : time
! Output : F        : Fluence

      implicit none
      real*8 F,Coords(3),time,X,Y,Z

      X = Coords(1)
      Y = Coords(2)
      Z = Coords(3)
      F = 2.5d0**(-0.002d0*X)*time
      return
      end
!*************************************************************************
      subroutine GetTempPosTime(T,Coords,time)
!************************************************************************
! This subroutine computes the temperature as a function of time and position
! Input  : Coords   : (x,y,z) coordinate in space
!        : time     : time
! Output : T        : Temperature

      implicit none
      real*8 T,Coords(3),time,X,Y,Z

      X = Coords(1)
      Y = Coords(2)
      Z = Coords(3)
      T = -1.0d0/10.0d0*X + 700
      return
      end
!*************************************************************************
              subroutine Get_Ew(T,F,Ew)
!*************************************************************************
! This subroutine computes the Wigner strain at the end of the time step
! Input  : T   : Temperature
!          F   : Fluence
! Output : Eq  : Wigner strain
      implicit none
      real*8 T,F,Ew(6),dL_par,dL_per

      call Get_Wigner_par(T,F,dL_par)
      call Get_Wigner_per(T,F,dL_per)
      Ew = 0.d0
      Ew(1) = dL_par
      Ew(2) = dL_par
      Ew(3) = dL_per

      return
      end
      subroutine Get_Eth(T,F,Eth)
!*************************************************************************
! This subroutine computes the thermal strain at the end of the time step
! Input  : T   : Temperature
!          F   : Fluence
! Output : Eth : Thermal strain
      implicit none
      real*8 T,F,Eth(6),CTE0(3),Ti,b1,b2,b3,b4,b5,b6
      real*8 a0,a1,a2,CTE(3),Scale,CTE_par,CTE_per,CTE0_par,CTE0_per
      call Get_CTEoCTE0par(T,F,CTE_par)
      call Get_CTEoCTE0per(T,F,CTE_per)
      call Get_CTE0par(T,CTE0_par)
      call Get_CTE0per(T,CTE0_per)

      Eth = 0.d0
      Eth(1) = CTE0_par*CTE_par*(T-20.0000000000000d0)
      Eth(2) = CTE0_par*CTE_par*(T-20.0000000000000d0)
      Eth(3) = CTE0_per*CTE_per*(T-20.0000000000000d0)
      return
      end
      subroutine WarningCode(T,F)
!*************************************************************************
! This subroutine warns the user when extrapolation in T or F occurs
! Input  : T   : Temperature
!          F   : Fluence
! Output : Warning Message Printed to the Screen 
      implicit none
      real*8 T,F
!      if (T.lt.300.000000000000d0) then 
!          write(*,*) '*WARNING EXTRAPOLOTION* T BELOW Calibration Data'
!          write(*,*)  'Temperature=', T
!      endif
!      if (T.gt.550.000000000000d0) then 
!          write(*,*) '*WARNING EXTRAPOLOTION* T ABOVE Calibration Data'
!          write(*,*)  'Temperature=', T
!      endif
!      if (F.lt.300.000000000000d0) then 
!          write(*,*) '*WARNING EXTRAPOLOTION* F BELOW Calibration Data'
!	       write(*,*)  'Fluence=', F
!      endif
!      if (F.gt.27.9700000000000d0) then 
!          write(*,*) '*WARNING EXTRAPOLOTION* F ABOVE Calibration Data'
!	       write(*,*)  'Fluence=', F
!      endif
      return
      end
