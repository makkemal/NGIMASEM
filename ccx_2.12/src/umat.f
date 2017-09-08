!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
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
     &  stran(ntens),dstran(ntens),time(2),celent,
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
      call GetTempPosTime(T_t2,Coords,time(2)+dtime)
      call GetTempPosTime(T_t1,Coords,time(2))
      call GetFluencePosTime(F_t2,Coords,time(2)+dtime)
      call GetFluencePosTime(F_t1,Coords,time(2))
      call Get_Eth(T_t2,F_t2,Eth)
      call Get_Ew(T_t2,F_t2,Ew)

      Ecp_t1    = statev(1:6)
      Ecs_t1    = statev(7:12)
      stress_t1 = stress
      
      if (time(2).lt.dtime) then
          Ecp_t1 = 0.d0
          Ecs_t1 = 0.d0
      endif

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
      Ecp_rate = -1.1528481170908d0*Ecp + 3.09336966247907d0*iDc_S
      dEcp_rate_dE = -1.15284811709080d0
      dEcp_rate_dS = 3.09336966247907d0*iDc
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

      Ecs_rate    = 1.38616125354559d0*iDc_S
      dEcsrate_dS = 1.38616125354559d0*iDc
      return
      end
      subroutine GetE0per(Value)
!*************************************************************************
! This subroutine constructs a subroutine that returns a constant
! Output : Warning Message Printed to the Screen 
      implicit none
      real*8 Value
      Value = 8918937500.00000d0
      return
      end
      subroutine GetE0par(Value)
!*************************************************************************
! This subroutine constructs a subroutine that returns a constant
! Output : Warning Message Printed to the Screen 
      implicit none
      real*8 Value
      Value = 9167692307.69231d0
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
      Fmin = 1.51000000000000d0
      Fmax = 47.5200000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = 10**(-1.25944020758358d0*F**3 - 2.63582274787837d0*F**2*T
     @ + 0.471874638713626d0*F**2 + 0.619693846726218d0*F*T**2 +
     @ 0.12794041929499d0*F*T + 1.04578223201656d0*F +
     @ 0.382362768560357d0*T**3 - 1.17380141312693d0*T**2 +
     @ 0.799353883880987d0*T + 0.117431567083108d0)
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
      Fmin = 1.84000000000000d0
      Fmax = 47.7600000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = 10**(-0.296421474264793d0*F**3 - 0.406345680899391d0*F**2
     @ *T - 1.32647511993486d0*F**2 + 1.2041997984614d0*F*T**2 -
     @ 2.08845873871361d0*F*T + 1.99093258934403d0*F +
     @ 0.943858592229335d0*T**3 - 2.1698316048103d0*T**2 +
     @ 1.47041574884044d0*T - 0.0290942152693225d0)
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
      Fmin = 1.84000000000000d0
      Fmax = 47.7600000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = -0.0503456202866336d0*F**3 + 0.197252810406648d0*F**2*T +
     @ 0.248949747607009d0*F**2 - 0.104694581564632d0*F*T**2 +
     @ 0.0839441133248617d0*F*T - 0.224521814581155d0*F -
     @ 0.0492064563810428d0*T**3 + 0.113712475200814d0*T**2 -
     @ 0.0783113440548514d0*T + 0.0181376155571642d0
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
      Fmin = 1.51000000000000d0
      Fmax = 47.5200000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = -0.0120416989992442d0*F**3 + 0.370510740754554d0*F**2*T +
     @ 0.182381527298861d0*F**2 - 0.0863525151390481d0*F*T**2 -
     @ 0.0257149318905605d0*F*T - 0.155601368619377d0*F -
     @ 0.0351320171426914d0*T**3 + 0.0919602481363923d0*T**2 -
     @ 0.0603924305920669d0*T + 0.0140084898193773d0
      return
      end
      subroutine Get_CTE_par(To,Fo,Answer)
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
      Fmin = 2.49000000000000d0
      Fmax = 38.4500000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = 5.67761262412476d-6*F**3 + 2.63562056100359d-6*F**2*T -
     @ 3.61943697087252d-6*F**2 - 1.92366612029263d-6*F*T**2 +
     @ 1.7246170823243d-6*F*T - 5.72712295646843d-6*F -
     @ 5.67828185881912d-6*T**3 + 1.03340254086343d-5*T**2 -
     @ 6.37018065404501d-6*T + 2.52484941887536d-6
      return
      end
      subroutine Get_CTE_per(To,Fo,Answer)
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
      Fmin = 2.17000000000000d0
      Fmax = 38.4500000000000d0

      T = (To-Tmin)/(Tmax-Tmin)
      F = (Fo-Fmin)/(Fmax-Fmin)

      Answer = 8.73564829871653d-6*F**3 + 9.78096917913715d-6*F**2*T -
     @ 1.20263149989436d-5*F**2 + 2.08312596315478d-6*F*T**2 -
     @ 8.99699873250278d-6*F*T - 1.05498850844557d-9*F -
     @ 2.23612761240206d-6*T**3 + 3.38266824036031d-6*T**2 -
     @ 1.67818134455067d-6*T + 1.68730953462395d-6
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
      E10 = EoverE0per*Eper
      E20 = EoverE0per*Eper
      E30 = EoverE0par*Epar
      nu12 = 0.150000000000000d0
      nu23 = 0.150000000000000d0
      nu31 = 0.150000000000000d0
      G120 = (EoverE0par+EoverE0per)/2.d0*2000000000.00000d0
      G230 = (EoverE0par+EoverE0per)/2.d0*2000000000.00000d0
      G310 = (EoverE0par+EoverE0per)/2.d0*2000000000.00000d0

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
      F = 20.0d0*time
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
      T = 500.0d0
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
      Ew(1) = dL_per
      Ew(2) = dL_per
      Ew(3) = dL_par

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
      real*8 a0,a1,a2,CTE(3),Scale,CTE_par,CTE_per
      call Get_CTE_Par(T,F,CTE_par)
      call Get_CTE_Per(T,F,CTE_per)

      Eth = 0.d0
      Eth(1) = CTE_per*(T-20.0000000000000d0)
      Eth(2) = CTE_per*(T-20.0000000000000d0)
      Eth(3) = CTE_par*(T-20.0000000000000d0)
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
      if (T.lt.300.000000000000d0) then 
          write(*,*) '*WARNING EXTRAPOLOTION* T BELOW Calibration Data'
      endif
      if (T.gt.550.000000000000d0) then 
          write(*,*) '*WARNING EXTRAPOLOTION* T ABOVE Calibration Data'
      endif
      if (F.lt.2.49000000000000d0) then 
          write(*,*) '*WARNING EXTRAPOLOTION* F BELOW Calibration Data'
      endif
      if (F.gt.27.9700000000000d0) then 
          write(*,*) '*WARNING EXTRAPOLOTION* F ABOVE Calibration Data'
      endif
      return
      end
