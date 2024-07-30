!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_tis_mod
!> @brief   calculate TIS interaction energy
!! @authors Naoto Hori (NH)
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_tis_mod

  use at_pairlist_str_mod
  use at_boundary_str_mod
  use at_enefunc_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: compute_energy_tis_lstack

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_tis_lstack
  !> @brief        calculate local (=consecutive) base stacking energy
  !! @authors      NH
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] estack  : stacking energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_tis_lstack(enefunc, coord, force, virial, estack)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: estack

    ! local variables
    integer                  :: i
    integer                  :: i1, i2, i3, i4, i5, i6, i7
    integer                  :: istart, iend
    real(wp)                 :: v21(1:3), v34(1:3), v54(1:3), v56(1:3), v76(1:3)
    real(wp)                 :: f_i(1:3), f_l(1:3)
    real(wp)                 :: ediv
    real(wp)                 :: dist, ddist, dih, d
    real(wp)                 :: abs54, d5454, abs56, d5656
    real(wp)                 :: d5654
    real(wp)                 :: d7656over5656, d5456over5656
    real(wp)                 :: d3454over5454, d5654over5454
    real(wp)                 :: m(3), n(3)
    real(wp)                 :: dnn, dmm
    real(wp)                 :: sol_T, u0

    integer,         pointer :: list(:,:)
    real(wp),        pointer :: h(:), s(:), Tm(:)
    real(wp),        pointer :: Kr(:), Kphi1(:), Kphi2(:)
    real(wp),        pointer :: r0(:), phi10(:), phi20(:)
    real(wp),        pointer :: work(:,:,:)

    !
    ! Reference:
    !   Equation (3) in
    !   N.A. Denesyuk and D. Thirumalai, J Phys. Chem. B (2013) 10.1021/jp401087x
    !
    ! Potential function:
    !    U_stack = U0 / [1 + Kr(r-r0)^2 + Kphi1(phi1-phi10)^2 + Kphi2(phi2-phi20)^2]
    !
    !   P1           !    3           !
    !    \           !     \          !  Particle IDs are stored in the list
    !     S1 -- B1   !      4 -- 1    !    B1 = 1: list(1, i)
    !    /           !     /          !    B2 = 2: list(2, i)
    !   P2           !    5           !    P1 = 3: list(3, i)
    !    \           !     \          !    S1 = 4: list(4, i)
    !     S2 -- B2   !      6 -- 2    !    P2 = 5: list(5, i)
    !    /           !     /          !    S2 = 6: list(6, i)
    !   P3           !    7           !    P3 = 7: list(7, i)
    ! 
    !     r = stacking distance B1-B2 (1-2)
    !     phi1 = dihedral angle P1-S1-P2-S2 (3-4-5-6)
    !     phi2 = dihedral angle S1-P2-S2-P3 (4-5-6-7)
    !
    ! Coefficients:
    !     U0, Kr, Kphi1, Kphi2
    !
    !     U0 = -h + kB * (T - Tm) * s
    !     h, s, Tm are thermodynamic parameters of the base stacking
    !     T is the solution temperature
    !
    ! Reference values (A-type RNA):
    !     r0     = stacking distance B1-B2
    !     phi10  = dihedral angle P1-S1-P2-S2
    !     phi20  = dihedral angle S1-P2-S2-P3
    !

    call timer(TimerTISLocalStack, TimerOn)

    istart   = enefunc%istart_tis_lstack
    iend     = enefunc%iend_tis_lstack
    sol_T    = enefunc%cg_ele_sol_T

    ! use pointers
    !
    list     => enefunc%tis_lstack_list
    h     => enefunc%tis_lstack_h
    s     => enefunc%tis_lstack_s
    Tm    => enefunc%tis_lstack_Tm
    Kr    => enefunc%tis_lstack_Kr
    Kphi1 => enefunc%tis_lstack_Kphi1
    Kphi2 => enefunc%tis_lstack_Kphi2
    r0    => enefunc%tis_lstack_r0
    phi10 => enefunc%tis_lstack_phi10
    phi20 => enefunc%tis_lstack_phi20

    work  => enefunc%work_tis_stack

    work(:, :, :) = 0.0e0_wp

    ! calculation of base energy and gradient
    !
    !$omp parallel do default(none)                    &
    !$omp private(i1, i2, i3, i4, i5, i6, i7,          &
    !$omp         v21, v34, v54, v56, v76,             &
    !$omp         f_i, f_l, ediv, dist, ddist, dih, d, &
    !$omp         abs54, d5454, abs56, d5656, d5654,   &
    !$omp         d7656over5656, d5456over5656,        &
    !$omp         d3454over5454, d5654over5454,        &
    !$omp         m, n, dnn, dmm, u0)                  &
    !$omp shared(istart, iend, Kr, Kphi1, Kphi2,       &
    !$omp        r0, phi10, phi20, work, coord, list,  &
    !$omp        h, s, Tm, sol_T)                      &
    !$omp reduction(+:estack) reduction(+:virial)
    !
    do i = istart, iend

      i1 = list(1, i)
      i2 = list(2, i)
      i3 = list(3, i)
      i4 = list(4, i)
      i5 = list(5, i)
      i6 = list(6, i)
      i7 = list(7, i)

      v21(1:3) = coord(1:3, i2) - coord(1:3, i1)
      v34(1:3) = coord(1:3, i3) - coord(1:3, i4)
      v54(1:3) = coord(1:3, i5) - coord(1:3, i4)
      v56(1:3) = coord(1:3, i5) - coord(1:3, i6)
      v76(1:3) = coord(1:3, i7) - coord(1:3, i6)

      ediv = 1.0e0_wp

      !===== Distance between 1 and 2 =====
      dist = norm2(v21)
      ddist = dist - r0(i)
      ediv = ediv + Kr(i) * ddist**2

      f_i(:) = 2.0e0_wp * Kr(i) * ddist / dist * v21(:)
      work(:, 1, i) = - f_i(:)
      work(:, 2, i) = + f_i(:)

      !===== Dihedral angle 3-4-5-6 =====
      m(1) = v34(2)*v54(3) - v34(3)*v54(2)
      m(2) = v34(3)*v54(1) - v34(1)*v54(3)
      m(3) = v34(1)*v54(2) - v34(2)*v54(1)
      n(1) = v54(2)*v56(3) - v54(3)*v56(2)
      n(2) = v54(3)*v56(1) - v54(1)*v56(3)
      n(3) = v54(1)*v56(2) - v54(2)*v56(1)

      dmm = dot_product(m,m)
      dnn = dot_product(n,n)
      d5454 = dot_product(v54,v54)
      abs54 = sqrt(d5454)
      d5654 = dot_product(v56,v54)
      d3454over5454 = dot_product(v34,v54) / d5454
      d5654over5454 = d5654 / d5454

      dih = atan2(dot_product(v34,n)*abs54 , dot_product(m,n))

      d = dih - phi10(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ediv = ediv + Kphi1(i) * d**2

      f_i(:) = + 2.0e0_wp * Kphi1(i) * d * abs54 / dmm * m(:)
      f_l(:) = - 2.0e0_wp * Kphi1(i) * d * abs54 / dnn * n(:)

      work(:, 3, i) = f_i(:)
      work(:, 4, i) = (-1.0e0_wp + d3454over5454) * f_i(:) &
                     -(            d5654over5454) * f_l(:)
      work(:, 5, i) = (-1.0e0_wp + d5654over5454) * f_l(:) &
                     -(            d3454over5454) * f_i(:)
      work(:, 6, i) = f_l(:)

      !===== Dihedral angle 7-6-5-4 =====
      m(1) = v76(2)*v56(3) - v76(3)*v56(2)
      m(2) = v76(3)*v56(1) - v76(1)*v56(3)
      m(3) = v76(1)*v56(2) - v76(2)*v56(1)
      !n(1) = v56(2)*v54(3) - v56(3)*v54(2)
      !n(2) = v56(3)*v54(1) - v56(1)*v54(3)
      !n(3) = v56(1)*v54(2) - v56(2)*v54(1)
      n(:) = -n(:)

      dmm = dot_product(m,m)
      !dnn = dot_product(n,n)  !! dnn does not change.
      d5656 = dot_product(v56, v56)
      abs56 = sqrt(d5656)
      d7656over5656 = dot_product(v76,v56) / d5656
      d5456over5656 = d5654 / d5656

      dih = atan2(dot_product(v76,n)*abs56 , dot_product(m,n))

      d = dih - phi20(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ediv = ediv + Kphi2(i) * d**2

      f_i(:) = + 2.0e0_wp * Kphi2(i) * d * abs56 / dmm * m(:)
      f_l(:) = - 2.0e0_wp * Kphi2(i) * d * abs56 / dnn * n(:)

      work(:, 7, i) = work(:, 7, i) + f_i(:)
      work(:, 6, i) = work(:, 6, i) + (-1.0e0_wp + d7656over5656) * f_i(:) &
                                    - (            d5456over5656) * f_l(:)
      work(:, 5, i) = work(:, 5, i) - (            d7656over5656) * f_i(:) &
                                    + (-1.0e0_wp + d5456over5656) * f_l(:)
      work(:, 4, i) = work(:, 4, i) + f_l(:)

      !===== Total =====
      u0 = -h(i) + KBOLTZ * (sol_T - Tm(i)) * s(i)
      estack = estack + u0 / ediv

      work(:, :, i) = u0 / ediv**2 * work(:, :, i)

      !do j = 1, 3
      !  do k = j+1, 3
      !    vtmp = work(j, i) * coord(k, i_sugar)    &
      !        + work(j + 3, i) * coord(k, i_base5) &
      !        + work(j + 6, i) * coord(k, i_base3)
      !    virial(k, j) = virial(k, j) - vtmp
      !    virial(j, k) = virial(j, k) - vtmp
      !  end do
      !  vtmp =    work(j, i) * coord(j, i_sugar) &
      !      + work(j + 3, i) * coord(j, i_base5) &
      !      + work(j + 6, i) * coord(j, i_base3)
      !  virial(j,j) = virial(j,j) - vtmp
      !end do

    end do
    !$omp end parallel do

    do i = istart, iend
      force(1:3, list(1, i), 1) = force(1:3, list(1, i), 1) + work(1:3, 1, i)
      force(1:3, list(2, i), 1) = force(1:3, list(2, i), 1) + work(1:3, 2, i)
      force(1:3, list(3, i), 1) = force(1:3, list(3, i), 1) + work(1:3, 3, i)
      force(1:3, list(4, i), 1) = force(1:3, list(4, i), 1) + work(1:3, 4, i)
      force(1:3, list(5, i), 1) = force(1:3, list(5, i), 1) + work(1:3, 5, i)
      force(1:3, list(6, i), 1) = force(1:3, list(6, i), 1) + work(1:3, 6, i)
      force(1:3, list(7, i), 1) = force(1:3, list(7, i), 1) + work(1:3, 7, i)
    end do

    call timer(TimerTISLocalStack, TimerOff)

    return

  end subroutine compute_energy_tis_lstack

end module at_energy_tis_mod
