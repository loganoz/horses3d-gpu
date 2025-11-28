#include "Includes.h"
module VariableConversion_MU
   use SMConstants
   use PhysicsStorage_MU
   use FluidData_MU
   implicit none

   private
   public   mGradientVariables
   public   GetmTwoFluidsViscosity, GetmOneFluidViscosity
   public   getVelocityGradients


   contains
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      pure subroutine mGradientVariables( nEqn, nGrad, Q, U, rho_ )
!
!        --------------------------------------------------------------
!        Returns all gradient variables EXCEPT the chemical potential,
!        to be done manually afterwards.
!        --------------------------------------------------------------
!
         !$acc routine seq
         implicit none
         integer, intent(in)                 :: nEqn, nGrad
         real(kind=RP), intent(in)           :: Q(nEqn)
         real(kind=RP), intent(out)          :: U(nGrad)
         real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local Variables
!        ---------------
!     
         real(kind=RP)  :: invSqrtRho

         invSqrtRho = 1.0_RP / sqrt(rho_)

         U = [0.0_RP, invSqrtRho, invSqrtRho, invSqrtRho, 1.0_RP] * Q

      end subroutine mGradientVariables

      pure subroutine GetmOneFluidViscosity(c, mu)
!
!        ***********************************
!           Here phi is the density, such
!           that varies linearly from the
!           density of fluid 1 to that of
!           fluid 2
!        ***********************************
!
         !$acc routine seq
         implicit none
         real(kind=RP), intent(in)   :: c
         real(kind=RP), intent(out)  :: mu
!
!        ---------------
!        Local variables
!        ---------------
!
         mu = dimensionless % mu(1)

      end subroutine GetmOneFluidViscosity

      pure subroutine GetmTwoFluidsViscosity(c, mu)
!
!        ***********************************
!           Here phi is the density, such
!           that varies linearly from the
!           density of fluid 1 to that of
!           fluid 2
!        ***********************************
!
         !$acc routine seq
         implicit none
         real(kind=RP), intent(in)   :: c
         real(kind=RP), intent(out)  :: mu
         real(kind=RP)  :: cHat

         cHat = min(max(c,0.0_RP),1.0_RP)

         mu = dimensionless % mu(1) * cHat + dimensionless % mu(2) * (1.0_RP - cHat)

      end subroutine GetmTwoFluidsViscosity
!
! /////////////////////////////////////////////////////////////////////
!
      pure subroutine getVelocityGradients(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
      !$acc routine seq
      implicit none
      !-arguments---------------------------------------------------
      real(kind=RP), intent(in)  :: Q(NCONS)
      real(kind=RP), intent(in)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD) 
      real(kind=RP), intent(out) :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
      !-local-variables---------------------------------------------
	  integer       :: i
	  real(kind=RP) :: buff
      !-------------------------------------------------------------
		 buff   = Q_x(1)
		 U_x(1) = Q_x(2)
		 U_x(2) = Q_x(3)
		 U_x(3) = Q_x(4)
		 
		 buff   = Q_y(1)
		 U_y(1) = Q_y(2)
		 U_y(2) = Q_y(3)
		 U_y(3) = Q_y(4)
		 
		 buff   = Q_z(1)
		 U_z(1) = Q_z(2)
		 U_z(2) = Q_z(3)
		 U_z(3) = Q_z(4)

   end subroutine getVelocityGradients

end module VariableConversion_MU
