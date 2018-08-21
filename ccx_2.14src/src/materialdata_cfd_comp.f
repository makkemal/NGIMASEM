!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
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
      subroutine materialdata_cfd_comp(nef,vel,shcon,nshcon,ielmatf,
     &  ntmat_,mi,cvel,vfa,cocon,ncocon,physcon,cvfa,ithermal,nface,
     &  umel,umfa,ielfa,hcfa)
!
!     calculation of material properties at elements centers and
!     face centers (compressible fluids)
!
      implicit none
!
      integer nef,i,imat,ntmat_,mi(*),ielmatf(mi(3),*),ithermal,
     &  nshcon(2,*),nface,ncocon(2,*),ielfa(4,*)
!
      real*8 t1l,vel(nef,0:7),shcon(0:3,ntmat_,*),cvel(*),vfa(0:7,*),
     &  cp,cocon(0:6,ntmat_,*),physcon(*),cvfa(*),umel(*),umfa(*),
     &  hcfa(*)
!
      intent(in) nef,shcon,nshcon,ielmatf,
     &  ntmat_,mi,cocon,ncocon,physcon,ithermal,nface,
     &  ielfa
!
      intent(inout) vel,vfa,cvel,cvfa,umel,umfa,hcfa
!     
c$omp parallel default(none)
c$omp& shared(nef,vel,ielmatf,shcon,ntmat_,nshcon,physcon,cvel,umel,
c$omp&        ithermal,nface,vfa,ielfa,cvfa,umfa,hcfa,cocon,ncocon)
c$omp& private(i,t1l,imat,cp)
!
!     element (cell) values
!
c$omp do
      do i=1,nef
         t1l=vel(i,0)
         imat=ielmatf(1,i)
!
!        density
!
c         vel(i,5)=vel(i,4)/(shcon(3,1,imat)*t1l)
!
!        heat capacity at constant volume
!
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,
     &       physcon)
!
!        cv=cp-r
!
         cvel(i)=cp-shcon(3,1,imat)
!
!        dynamic viscosity
!
         call materialdata_dvi(shcon,nshcon,imat,umel(i),t1l,ntmat_,
     &            ithermal)
      enddo
c$omp end do
!
!     facial values
!
c$omp do
      do i=1,nface
         t1l=vfa(0,i)
!
!        take the material of the first adjacent element
!
         imat=ielmatf(1,ielfa(1,i))
!
!        density
!
c         vfa(5,i)=vfa(4,i)/(shcon(3,1,imat)*t1l)
!
!        heat capacity at constant volume
!
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,
     &       physcon)
!
!        cv=cp-r
!
         cvfa(i)=cp-shcon(3,1,imat)
!
!        dynamic viscosity
!
         call materialdata_dvi(shcon,nshcon,imat,umfa(i),t1l,ntmat_,
     &            ithermal)
!
!        heat conduction
!
         call materialdata_cond(imat,ntmat_,t1l,cocon,ncocon,hcfa(i))
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
