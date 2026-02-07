!>
!! @file m_hb_function.f90
!! @brief Contains module m_hb_function for Herschel-Bulkley non-Newtonian viscosity
!!        with Papanastasiou regularization for improved low shear rate behavior
#:include 'macros.fpp'

!> @brief The module contains functions to compute non-Newtonian viscosity
!!        using the Herschel-Bulkley model with Papanastasiou regularization
!!
!!        Standard Herschel-Bulkley: mu = tau0/gdot + K*gdot^(nn-1)
!!        This has a singularity as gdot -> 0
!!
!!        Papanastasiou regularization: mu = (tau0/gdot)*(1 - exp(-m*gdot)) + K*gdot^(nn-1)
!!        As gdot -> 0: mu -> tau0*m (finite, no singularity)
!!        As gdot -> inf: recovers standard HB model
!!
!!        Parameter m controls the transition sharpness:
!!        - Large m (>1000): Sharp transition, closer to ideal Bingham/HB behavior
!!        - Small m (~100): Smoother transition, better numerical stability
module m_hb_function

  use m_derived_types
  use m_global_parameters

  private; public :: f_compute_hb_viscosity, &
                    f_compute_shear_rate_from_components

contains

  !> Computes Herschel-Bulkley viscosity with Papanastasiou regularization
  !!
  !! The Papanastasiou model provides a smooth, regularized approximation:
  !!   mu = (tau0/gdot)*(1 - exp(-m*gdot)) + K*gdot^(nn-1)
  !!
  !! @param tau0 Yield stress
  !! @param K_val Consistency index
  !! @param nn_val Flow behavior index
  !! @param mu_min_val Minimum viscosity limit
  !! @param mu_max_val Maximum viscosity limit
  !! @param shear_rate Shear rate magnitude
  !! @param hb_m_val Papanastasiou regularization parameter
  !! @return Viscosity
  pure function f_compute_hb_viscosity(tau0, K_val, nn_val, &
      mu_min_val, mu_max_val, shear_rate, hb_m_val) result(mu)
    $:GPU_ROUTINE(parallelism='[seq]')

    real(wp), intent(in) :: tau0, K_val, nn_val
    real(wp), intent(in) :: mu_min_val, mu_max_val
    real(wp), intent(in) :: shear_rate, hb_m_val
    real(wp) :: mu
    real(wp) :: yield_term, power_law_term, exp_term

    exp_term = exp(-hb_m_val*shear_rate)
    yield_term = tau0*(1._wp - exp_term)/shear_rate
    power_law_term = K_val*(shear_rate**(nn_val - 1._wp))

    mu = yield_term + power_law_term
    mu = min(max(mu, mu_min_val), mu_max_val)

  end function f_compute_hb_viscosity

  !> Computes shear rate from strain rate tensor components
  !! Works for 1D, 2D, and 3D cases (set D_zz, D_xz, D_yz to 0 for 2D/1D)
  !!
  !! The shear rate is: gdot = sqrt(2*D_ij*D_ij)
  !! where D_ij is the strain rate tensor
  !!
  !! @param D_xx Normal strain rate du/dx
  !! @param D_yy Normal strain rate dv/dy
  !! @param D_zz Normal strain rate dw/dz
  !! @param D_xy Shear strain rate 0.5*(du/dy + dv/dx)
  !! @param D_xz Shear strain rate 0.5*(du/dz + dw/dx)
  !! @param D_yz Shear strain rate 0.5*(dv/dz + dw/dy)
  !! @return Shear rate magnitude
  pure function f_compute_shear_rate_from_components( &
      D_xx, D_yy, D_zz, D_xy, D_xz, D_yz) result(shear_rate)
    $:GPU_ROUTINE(parallelism='[seq]')

    real(wp), intent(in) :: D_xx, D_yy, D_zz, D_xy, D_xz, D_yz
    real(wp) :: shear_rate

    ! |gdot| = sqrt(2*D_ij*D_ij) where D_ij is the symmetric strain rate tensor
    ! For symmetric tensor: 2*D_ij*D_ij = 2*(D_xx^2+D_yy^2+D_zz^2+2*(D_xy^2+D_xz^2+D_yz^2))
    shear_rate = sqrt(2._wp*(D_xx*D_xx + D_yy*D_yy + D_zz*D_zz + &
                             2._wp*(D_xy*D_xy + D_xz*D_xz + D_yz*D_yz)))

    ! Clamp for numerical safety
    shear_rate = min(max(shear_rate, 1.0e-2_wp), 1.0e5_wp)

  end function f_compute_shear_rate_from_components

end module m_hb_function
