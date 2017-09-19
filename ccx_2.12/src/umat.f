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
!     here, an ABAQUS umat routine can be inserted
!
!     note that reals should be double precision (REAL*8)
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
!     START EXAMPLE LINEAR ELASTIC MATERIAL
!
      integer i,j
      real*8 Tc,Tp,Fc,Fp,iDel(6,6),Eth(6),Ew(6),iDc_c(6,6),iDc_p(6,6)
      real*8 dT,dF,EoverE0_ZF,iDel_p(6,6),EoverE0_ZF_p,iDn(6,6),Dn(6,6)
      real*8 K,SC_a,SC_e,SC_k,beta_n,alpha,G0,iDc_str(6),iDcp_str(6)
      real*8 beta_p, Ecp(6),Ecs(6)
!
      call GetTempPosTime(Tc,Coords,time(2)+dtime)
      call GetTempPosTime(Tp,Coords,time(2))
      call GetFluencePosTime(Fc,Coords,time(2)+dtime)
      call GetFluencePosTime(Fp,Coords,time(2))
      dF = Fc-Fp
      dT = Tc-Tp

      K = 0.23d0
      SC_a = 0.294d-5
      SC_e = 0.017d0
      SC_k = 8.62d-5
      beta_n = SC_a*(Tc-300.d0)**2 * dexp(-SC_e/(SC_k*(Tc-300.d0))) 
     &         + 1.d0
      beta_p = SC_a*(Tp-300.d0)**2 * dexp(-SC_e/(SC_k*(Tp-300.d0))) 
     &         + 1.d0
      alpha = 1.2d0
      G0    = 0.4d0

      call Get_invDel(iDel,Tc,Fc,EoverE0_ZF)
      call Get_invDel(iDel_p,Tp,Fp,EoverE0_ZF_p)

      write(*,*) 'EoverE0 = ',EoverE0_ZF

      iDc_c = iDel*EoverE0_ZF
      iDc_p = iDel_p*EoverE0_ZF_p

      write(*,*) 'iDc_c=',iDc_c
      iDn = iDel+ (dF*alpha/(2.d0*G0+dF) + dF*K*beta_n/2.d0)*iDc_c
      call INV_D(iDn, Dn)

      call Get_Eth(Tc,Fc,Eth)
      call Get_Ew(Tc,Fc,Ew)

      Ecp = statev(1:6)
      Ecs = statev(7:12)

      iDcp_str = 0.d0
      do i=1,6
         iDcp_str(i) = 0.d0
         do j=1,6
            iDcp_str(i) = iDcp_str(i) + iDc_p(i,j)*stress(j)
         enddo
      enddo

      stran = stran+dstran-Ew-Eth-Ecs - (2.d0*G0-dF)/(2.d0*G0+dF)*Ecp
     &       -(dF*alpha/(2.d0*G0+dF)+dF*K*beta_p/2.d0)*iDcp_str

      do i=1,6
         stress(i) = 0.d0
         do j=1,6
            stress(i) = stress(i) + Dn(i,j)*(stran(j))
         enddo
      enddo

      ddsdde = Dn

      iDc_str = 0.d0
      do i=1,6
         do j=1,6
            iDc_str(i) = iDc_str(i) + iDc_c(i,j)*stress(j)
         enddo
      enddo

      Ecs = Ecs + dF*K/2.d0*(beta_p*iDcp_str+beta_n*iDc_str)
      Ecp = (Ecp + dF/(2.d0*G0)*(alpha*iDcp_str-Ecp+alpha*iDc_str))/
     &      (1.d0+dF/(2.d0*G0)) 

      statev(1:6)  = Ecp
      statev(7:12) = Ecs

      return
      end

      subroutine Get_invDel(iDel,T,F,EoverE0_ZF)
      implicit none
      real*8 iDel(6,6)
      real*8 T,F,c0,c2,c3,T0,b1,b2,b3,b4,b5,b6,b7,b8,b9
      real*8 E10,E20,E30,nu12,nu23,nu31,G120,G230,G310
      real*8 a10,a12,a22,a23,a20,a21,k1,EoverE0,EoverE0_ZF

      c0 = 0.3307d0
      c2 = 7.732d-3
      c3 = -7.499d-3
      T0 = 340.d0
      b1 = 9.211d-7
      b2 = 8.941d-10
      b3 = -8.137d-5
      b4 = -3.684d-7
      b5 = 1.601d-7
      b6 = -7.877d-13
      b7 = 1.0d2
      b8 = -3.193d-2
      b9 = -5.599d-13

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

      a10 = c0 + ((1-c0)+c2*(T-T0))*dexp(c3*(T-T0))
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

      iDel = EoverE0*iDel

      EoverE0_ZF = 1.d0 + a10

      return
      end

      subroutine GetTempPosTime(T,Coords,time)
      implicit none
      real*8 T,Coords(3),time
      
      T  = 300.d0 + 300.d0*Coords(1)**2

      return
      end

      subroutine GetFluencePosTime(F,Coords,time)
      implicit none
      real*8 F,Coords(3),time

      F = 1.d-20*dexp(30.6263d0-9.2208*Coords(1))*
     &     3600.d0*24.d0*365.d0*time

      return
      end

      subroutine Get_Eth(T,F,Eth)
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

      subroutine Get_Ew(T,F,Ew)
      implicit none
      real*8 T,F,Ew(6),b_par(3,3),b_per(3,3)
      real*8 a(3),Scale,dL_par,dL_per
      integer i,j

      b_par(1,1) = -9.3713d-2
      b_par(1,2) =  2.2227d-4
      b_par(1,3) = -1.6039d-7
      b_par(2,1) =  1.0306d-3
      b_par(2,2) = -3.87322d-6
      b_par(2,3) =  2.9755d-9
      b_par(3,1) = -3.4277d-6
      b_par(3,2) =  1.4066d-8
      b_par(3,3) = -9.6708d-12

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
            a(i) = a(i) + b_par(i,j)*T**(j-1)
         enddo
      enddo
      dL_par = a(1)*F + a(2)*F**2 + a(3)*F**3
      do i=1,3
         a(i) = 0.d0
         do j=1,3
            a(i) = a(i) + b_per(i,j)*T**(j-1)
         enddo
      enddo
      dL_per = a(1)*F + a(2)*F**2 + a(3)*F**3

      Ew = 0.d0
      Ew(1) = dL_par
      Ew(2) = dL_per
      Ew(3) = dL_per
      return
      end


      subroutine INV_D(A, AINV)
C     
C
C     !******************************************************************
C     !  INV_D  -  Compute the inverse of a 6x6 compliance matrix.
C     !
C     !  A       = input 6x6 compliance matrix to be inverted
C     !  AINV    = output 6x6 inverse of compliance matrix A
C     !  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, 
C     !            and .FALSE. if the input matrix is singular.
C     !******************************************************************

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


