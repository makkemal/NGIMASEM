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

      Answer = 10**(F*(-2.02919485524931d0*F**3*T - 3.69067296269081d0*F
     @ **3 + 8.05630627594521d0*F**2 - 1.18383580569256d0*F*T**2 -
     @ 7.24272115019942d0*F + 0.267682741256171d0*T + 3.27822378551637d0
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

      Answer = 10**(F*(-1.11508717462275d0*F**3*T - 4.08620889040924d0*F
     @ **3 + 8.94323620793999d0*F**2 - 1.24008634535327d0*F*T**2 -
     @ 7.97562691470214d0*F + 0.042604583021196d0*T + 3.46166689637576d0
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

      Answer = 0.0153530586150534d0*F**3 + 0.349578699411909d0*F**2*T +
     @ 0.119388458183265d0*F**2 + 0.00599094066827001d0*F*T**2 -
     @ 0.15341128308311d0*F*T - 0.145297824185457d0*F +
     @ 0.0117869473680535d0*T + 0.00318409227515302d0
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

      Answer = 0.0456699243614627d0*F**3 + 0.513365732094782d0*F**2*T +
     @ 0.0651624030879375d0*F**2 + 0.0208700214442222d0*F*T**2 -
     @ 0.253310247001052d0*F*T - 0.0804745844722801d0*F +
     @ 0.0200651799638447d0*T - 0.00129205939267519d0
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

      Answer = 1.99112848567294d0*F**3 + 1.52818622895704d0*F**2*T -
     @ 2.22925020532312d0*F**2 + 0.547100861242998d0*F*T**2 -
     @ 1.4626664542219d0*F*T - 0.454532925706461d0*F -
     @ 0.115293825988803d0*T + 0.417211034581856d0
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

      Answer = 2.20933719030931d0*F**3 + 2.20370263505754d0*F**2*T -
     @ 3.22094006243176d0*F**2 + 0.413079561933016d0*F*T**2 -
     @ 2.06932307903771d0*F*T + 0.324400701995973d0*F -
     @ 0.0150027392087665d0*T + 0.335393164279918d0
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

      Answer = 10**(-0.0138211653718181d0*T**2 + 0.0469701719735655d0*T
     @ - 5.32664359354929d0)
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
!     if (T.lt.300.000000000000d0) then 
!          write(*,*) '*WARNING EXTRAPOLOTION* T BELOW Calibration Data'
!      endif
!      if (T.gt.550.000000000000d0) then 
!          write(*,*) '*WARNING EXTRAPOLOTION* T ABOVE Calibration Data'
!      endif
!      if (F.lt.300.000000000000d0) then 
!          write(*,*) '*WARNING EXTRAPOLOTION* F BELOW Calibration Data'
!      endif
!      if (F.gt.27.9700000000000d0) then 
!          write(*,*) '*WARNING EXTRAPOLOTION* F ABOVE Calibration Data'
!      endif
      return
      end
