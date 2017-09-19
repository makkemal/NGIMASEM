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
      stress_t1 = statev(13:18)

      nrmRes = 1.d0
      stress_t2 = stress_t1
      do while(nrmRes.gt.1.d-12)
 
         call Get_Eel(T_t2,F_t2,stress_t2,Eel,dEeldS)
         call Get_Ecp(T_t1,F_t1,stress_t1,T_t2,F_t2,stress_t2,
     .             Ecp_t1,Ecp_t2,dEcpdS)
         call Get_Ecs(T_t1,F_t1,stress_t1,T_t2,F_t2,stress_t2,
     .             Ecs_t1,Ecs_t2,dEcsdS)

         Res = stran + dstran - Eel - Ew - Eth - Ecs_t2 - Ecp_t2
         nrmRes = 0.d0
         do i=1,6
            nrmRes = nrmRes + Res(i)**2
         enddo
         nrmRes = dsqrt(nrmRes)
         dResdS = dEeldS + dEcpdS + dEcsdS
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
      statev(13:18) = stress

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
      integer i,j
      real*8 T_t1,T_t2,F_t1,F_t2,stress_t1(6),stress_t2(6)
      real*8 Ecp_t1(6),Ecp_t2(6),dEcpdS(6,6)
      real*8 dF,nrmRes,Res(6)
      real*8 dEcp_rate_dS_t1(6,6),dEcp_rate_dS_t2(6,6)
      real*8 dEcp_rate_dE_t1,dEcp_rate_dE_t2,dResdEcp
      real*8 Ecp_rate_t1(6),Ecp_rate_t2(6)

      call GetEcp_rate(T_t1,F_t1,stress_t1,Ecp_t1,Ecp_rate_t1,
     .                 dEcp_rate_dS_t1,dEcp_rate_dE_t1)
      dF = F_t2 - F_t1
! Initial guess for primary creep at end of increment
      Ecp_t2 = Ecp_t1
      nrmRes = 1.d0
      do while(nrmRes.gt.1.d-12)
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
      alpha = 1.2d0
      G0    = 0.4d0
      call Get_invDel(iDel,T,F)
      call GetEoverE0(T,0.d0,EoverE0_ZF)
      iDc = iDel*EoverE0_ZF
      iDc_S = 0.d0
      do i=1,6
         do j=1,6
            iDc_S(i) = iDc_S(i) + iDc(i,j)*stress(j)
         enddo
      enddo
      Ecp_rate = 1.d0/G0*(alpha*iDc_S - Ecp)
      dEcp_rate_dE = -1.d0/G0
      dEcp_rate_dS = 1.d0/G0*alpha*iDc
      return
      end
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

      K    = 0.23d0
      SC_a = 0.294d-5
      SC_e = 0.017d0
      SC_k = 8.62d-5
      beta = SC_a*(T-300.d0)**2 * dexp(-SC_e/(SC_k*(T-300.d0))) 
     &         + 1.d0
      
      call Get_invDel(iDel,T,F)
      call GetEoverE0(T,0.d0,EoverE0_ZF)

      iDc = iDel*EoverE0_ZF

      iDc_S = 0.d0
      do i=1,6
         do j=1,6
            iDc_S(i) = iDc_S(i) + iDc(i,j)*stress(j)
         enddo
      enddo

      Ecs_rate    = K*beta*iDc_S
      dEcsrate_dS = K*beta*iDc
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
      call Get_invDel(dEeldS,T,F)
      Eel = 0.d0
      do i=1,6
         do j=1,6
            Eel(i) = Eel(i) + dEeldS(i,j)*stress(j)
         enddo
      enddo
      return
      end
!*************************************************************************
      subroutine Get_invDel(iDel,T,F)
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
      real*8 a10,a12,a22,a23,a20,a21,k1,EoverE0

      E10 = 8.40d9
      E20 = 7.56d9
      E30 = 7.56d9
      nu12 = 0.15d0
      nu23 = 0.15d0
      nu31 = 0.15d0
      G120 = 4.69d9
      G230 = 4.45d9
      G310 = 4.69d9

      iDel = 0.d0
      iDel(1,1) = 1.d0/E10
      iDel(2,2) = 1.d0/E20
      iDel(3,3) = 1.d0/E30
      iDel(4,4) = 1.d0/G120
      iDel(5,5) = 1.d0/G310
      iDel(6,6) = 1.d0/G230
      
      iDel(1,2) = -nu12/E10
      iDel(1,3) = -nu31/E30
      iDel(2,3) = -nu23/E20
      iDel(2,1) = iDel(1,2)
      iDel(3,1) = iDel(1,3)
      iDel(3,2) = iDel(2,3)

      call GetEoverE0(T,F,EoverE0)

      iDel = EoverE0*iDel

      return
      end
!*************************************************************************
      subroutine GetEoverE0(T,F,EoverE0)
!*************************************************************************
! This subroutine computes a scalar scale factor that relates the 
! elasticiy modulus to its unirradiated value
! Input  : T       : Temperature
!          F       : Fluence
! Output : EoverE0 : Scale factor
      implicit none
      real*8 T,F,EoverE0
      real*8 c0,c2,c3,T0,b1,b2,b3,b4,b5,b6,b7,b8,b9,a10,a12,a22
      real*8 k1,a20,a21,a23

      c0 = 0.3307d0
      c2 = 7.732d-3
      c3 = -7.499d-3
      T0 = 340.d0

      a10 = c0 + ((1-c0)+c2*(T-T0))*dexp(c3*(T-T0))

      if (F.eq.0d0) then
         EoverE0 = 1.d0 + a10
         return
      endif
      b1 = 9.211d-7
      b2 = 8.941d-10
      b3 = -8.137d-5
      b4 = -3.684d-7
      b5 = 1.601d-7
      b6 = -7.877d-13
      b7 = 1.0d2
      b8 = -3.193d-2
      b9 = -5.599d-13

      a12 = b1 + b2*T + b9*T**2
      a22 = b3 + b4*T
      a23 = b5 + b6*T
      k1  = b7 + b8*T
      a20 = a10 + a22*k1**2 - 2.d0*a12*k1**3 + 2.d0*a23*k1**3
      a21 = -2.d0*a22*k1 + 3.d0*a12*k1**2 - 3.d0*a23*k1**2
      if(F.lt.k1) then
         EoverE0 = 1.d0 + a10 + a12*F**3
      else
         EoverE0 = 1.d0 + a20 + a21*F + a22*F**2+a23*F**3
      endif

      if(EoverE0.lt.1.d0) EoverE0 = 1.d0

      return
      end
!************************************************************************
      subroutine GetTempPosTime(T,Coords,time)
!************************************************************************
! This subroutime computes the temperature as a function of time and position
! Input  : Coords   : (x,y,z) coordinate in space
!        : time     : time
! Output : T        : Temperature

      implicit none
      real*8 T,Coords(3),time
      
      T  = 300.d0 + 300.d0*Coords(1)**2
!      T = 300.d0
      return
      end
!*************************************************************************
      subroutine GetFluencePosTime(F,Coords,time)
!*************************************************************************
! This subroutime computes the fluence as a function of time and position
! Input  : Coords   : (x,y,z) coordinate in space
!        : time     : time
! Output : F        : Fluence

      implicit none
      real*8 F,Coords(3),time

      F = 1.d-20*dexp(30.6263d0-9.2208*Coords(1))*
     &     3600.d0*24.d0*365.d0*time

!      F = 0.d0

      return
      end
!*************************************************************************
      subroutine Get_Eth(T,F,Eth)
!*************************************************************************
! This subroutime computes the thermal strain at the end of the time step
! Input  : T   : Temperature
!          F   : Fluence
! Output : Eth : Thermal strain
      implicit none
      real*8 T,F,Eth(6),CTE0(3),Ti,b1,b2,b3,b4,b5,b6
      real*8 a0,a1,a2,CTE(3),Scale

      CTE0(1) = 4.3d-6
      CTE0(2) = 3.87d-6
      CTE0(3) = 3.87d-6
      b1 = 7.2623d-1
      b2 = -2.8773d-4
      b3 = 1.3672d-2
      b4 = 2.5076d-5
      b5 = -7.5109d-3
      b6 = -3.2002d-5
      Ti = 20.d0
      a0 = b1 + b2*T
      a1 = b3 + b4*T
      a2 = b5 + b6*T

      Scale = (a0 + ((1.d0-a0) + a1*F))*dexp(a2*F)
      CTE(1) = CTE0(1)*Scale
      CTE(2) = CTE0(2)*Scale
      CTE(3) = CTE0(3)*Scale
      Eth = 0.d0
      Eth(1) = CTE(1)*(T-Ti)
      Eth(2) = CTE(2)*(T-Ti)
      Eth(3) = CTE(3)*(T-Ti)
      return
      end
!*************************************************************************
      subroutine Get_Ew(T,F,Ew)
!*************************************************************************
! This subroutime computes the Wigner strain at the end of the time step
! Input  : T   : Temperature
!          F   : Fluence
! Output : Eq  : Wigner strain
      implicit none
      real*8 T,F,Ew(6),dL_par,dL_per

      call Get_Wigner_par(T,F,dL_par)
      call Get_Wigner_per(T,F,dL_per)
      Ew = 0.d0
      Ew(1) = dL_par
      Ew(2) = dL_per
      Ew(3) = dL_per

      return
      end
!*************************************************************************
      subroutine Get_Wigner_par(T,F,dL_par)
!*************************************************************************
! This subroutine computes the wigner strain in the parallel direction
! Input  : T      : Temperature
!          F      : Fluence
! Output : dL_par : Wigner strain in parallel direction
      implicit none
      integer i,j
      real*8 T,F,b_par(3,3),a(3),dL_par
      b_par(1,1) = -9.3713d-2
      b_par(1,2) =  2.2227d-4
      b_par(1,3) = -1.6039d-7
      b_par(2,1) =  1.0306d-3
      b_par(2,2) = -3.87322d-6
      b_par(2,3) =  2.9755d-9
      b_par(3,1) = -3.4277d-6
      b_par(3,2) =  1.4066d-8
      b_par(3,3) = -9.6708d-12

      do i=1,3
         a(i) = 0.d0
         do j=1,3
            a(i) = a(i) + b_par(i,j)*T**(j-1)
         enddo
      enddo
      dL_par = a(1)*F + a(2)*F**2 + a(3)*F**3
      return
      end
!*************************************************************************
      subroutine Get_Wigner_per(T,F,dL_per)
!*************************************************************************
! This subroutine computes the wigner strain in the perpendicular direction
! Input  : T      : Temperature
!          F      : Fluence
! Output : dL_per : Wigner strain in perpendicular direction
      implicit none
      integer i,j
      real*8 T,F,b_per(3,3),a(3),dL_per

      b_per(1,1) = -5.6021d-2
      b_per(1,2) =  1.0798d-4
      b_per(1,3) = -9.4945d-8
      b_per(2,1) =  5.2600d-4
      b_per(2,2) = -1.7812d-6
      b_per(2,3) =  1.4581d-9
      b_per(3,1) = -1.6033d-6
      b_per(3,2) =  6.3388d-9
      b_per(3,3) = -3.7126d-12

      do i=1,3
         a(i) = 0.d0
         do j=1,3
            a(i) = a(i) + b_per(i,j)*T**(j-1)
         enddo
      enddo
      dL_per = a(1)*F + a(2)*F**2 + a(3)*F**3

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


