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
  use molecules_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: compute_energy_tis_lstack
  public  :: compute_energy_tis_lstack_pbc
  public  :: compute_energy_tis_mwca
  public  :: compute_energy_tis_mwca_pbc

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
  !! @note         Denesyuk and Thirumalai, J Phys. Chem. B (2013) 10.1021/jp401087x
  !
  ! Todo: virial not calculated
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_tis_lstack(enefunc, coord, force, virial, estack)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: estack

    ! local variables
    integer                  :: i, ithread 
    integer                  :: i1, i2, i3, i4, i5, i6, i7
    integer                  :: istart, iend
    real(wp)                 :: v21(1:3), v34(1:3), v54(1:3), v56(1:3), v76(1:3)
    real(wp)                 :: f_i(1:3), f_l(1:3)
    real(wp)                 :: for(1:3, 1:7)
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

    integer                  :: omp_get_thread_num

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
    list  => enefunc%tis_lstack_list
    h     => enefunc%tis_lstack_h
    s     => enefunc%tis_lstack_s
    Tm    => enefunc%tis_lstack_Tm
    Kr    => enefunc%tis_lstack_Kr
    Kphi1 => enefunc%tis_lstack_Kphi1
    Kphi2 => enefunc%tis_lstack_Kphi2
    r0    => enefunc%tis_lstack_r0
    phi10 => enefunc%tis_lstack_phi10
    phi20 => enefunc%tis_lstack_phi20

    ! calculation of local stacking energy and gradient
    !
    !$omp parallel do default(none)                    &
    !$omp private(ithread, i1, i2, i3, i4, i5, i6, i7, &
    !$omp         v21, v34, v54, v56, v76,             &
    !$omp         ediv, dist, ddist, dih, d,           &
    !$omp         f_i, f_l, for,                       &
    !$omp         abs54, d5454, abs56, d5656, d5654,   &
    !$omp         d7656over5656, d5456over5656,        &
    !$omp         d3454over5454, d5654over5454,        &
    !$omp         m, n, dnn, dmm, u0)                  &
    !$omp shared(istart, iend, Kr, Kphi1, Kphi2,       &
    !$omp        r0, phi10, phi20, coord, list,        &
    !$omp        h, s, Tm, sol_T, force)               &
    !$omp reduction(+:estack) reduction(+:virial)
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
      for(:, 1) = - f_i(:)
      for(:, 2) = + f_i(:)

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

      for(:, 3) = f_i(:)
      for(:, 4) = (-1.0e0_wp + d3454over5454) * f_i(:) &
                 -(            d5654over5454) * f_l(:)
      for(:, 5) = (-1.0e0_wp + d5654over5454) * f_l(:) &
                 -(            d3454over5454) * f_i(:)
      for(:, 6) = f_l(:)

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

      for(:, 7) =           + f_i(:)
      for(:, 6) = for(:, 6) + (-1.0e0_wp + d7656over5656) * f_i(:) &
                            - (            d5456over5656) * f_l(:)
      for(:, 5) = for(:, 5) - (            d7656over5656) * f_i(:) &
                            + (-1.0e0_wp + d5456over5656) * f_l(:)
      for(:, 4) = for(:, 4) + f_l(:)

      !===== Total =====
      u0 = -h(i) + KBOLTZ * (sol_T - Tm(i)) * s(i)
      estack = estack + u0 / ediv

      for(:, :) = u0 / ediv**2 * for(:, :)

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

#ifdef OMP
      ithread = omp_get_thread_num() + 1
#else
      ithread = 1
#endif
      force(1:3, i1, ithread) = force(1:3, i1, ithread) + for(1:3, 1)
      force(1:3, i2, ithread) = force(1:3, i2, ithread) + for(1:3, 2)
      force(1:3, i3, ithread) = force(1:3, i3, ithread) + for(1:3, 3)
      force(1:3, i4, ithread) = force(1:3, i4, ithread) + for(1:3, 4)
      force(1:3, i5, ithread) = force(1:3, i5, ithread) + for(1:3, 5)
      force(1:3, i6, ithread) = force(1:3, i6, ithread) + for(1:3, 6)
      force(1:3, i7, ithread) = force(1:3, i7, ithread) + for(1:3, 7)
    end do
    !$omp end parallel do

    call timer(TimerTISLocalStack, TimerOff)

    return

  end subroutine compute_energy_tis_lstack


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_tis_lstack_pbc
  !> @brief        calculate local (=consecutive) base stacking energy in PBC
  !! @authors      NH
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    boundary: information of boundary condition
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] estack  : stacking energy of target systems
  !! @note         Denesyuk and Thirumalai, J Phys. Chem. B (2013) 10.1021/jp401087x
  !
  ! Todo: virial not calculated
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_tis_lstack_pbc(enefunc, boundary, coord, &
                                           force, virial, estack)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_boundary),target, intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: estack

    ! local variables
    integer                  :: i, ithread
    integer                  :: i1, i2, i3, i4, i5, i6, i7
    integer                  :: istart, iend
    real(wp)                 :: v21(1:3), v34(1:3), v54(1:3), v56(1:3), v76(1:3)
    real(wp)                 :: f_i(1:3), f_l(1:3)
    real(wp)                 :: for(1:3, 1:7)
    real(wp)                 :: ediv
    real(wp)                 :: dist, ddist, dih, d
    real(wp)                 :: abs54, d5454, abs56, d5656
    real(wp)                 :: d5654
    real(wp)                 :: d7656over5656, d5456over5656
    real(wp)                 :: d3454over5454, d5654over5454
    real(wp)                 :: m(3), n(3)
    real(wp)                 :: dnn, dmm
    real(wp)                 :: sol_T, u0
    real(wp)                 :: bsize(3), half_bsize(3)

    integer,         pointer :: list(:,:)
    real(wp),        pointer :: h(:), s(:), Tm(:)
    real(wp),        pointer :: Kr(:), Kphi1(:), Kphi2(:)
    real(wp),        pointer :: r0(:), phi10(:), phi20(:)

    integer                  :: omp_get_thread_num

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
    list  => enefunc%tis_lstack_list
    h     => enefunc%tis_lstack_h
    s     => enefunc%tis_lstack_s
    Tm    => enefunc%tis_lstack_Tm
    Kr    => enefunc%tis_lstack_Kr
    Kphi1 => enefunc%tis_lstack_Kphi1
    Kphi2 => enefunc%tis_lstack_Kphi2
    r0    => enefunc%tis_lstack_r0
    phi10 => enefunc%tis_lstack_phi10
    phi20 => enefunc%tis_lstack_phi20

    bsize(1) = boundary%box_size_x
    bsize(2) = boundary%box_size_y
    bsize(3) = boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculation of local stacking energy and gradient
    !
    !$omp parallel do default(none)                    &
    !$omp private(ithread, i1, i2, i3, i4, i5, i6, i7, &
    !$omp         v21, v34, v54, v56, v76,             &
    !$omp         ediv, dist, ddist, dih, d,           &
    !$omp         f_i, f_l, for,                       &
    !$omp         abs54, d5454, abs56, d5656, d5654,   &
    !$omp         d7656over5656, d5456over5656,        &
    !$omp         d3454over5454, d5654over5454,        &
    !$omp         m, n, dnn, dmm, u0)                  &
    !$omp shared(istart, iend, Kr, Kphi1, Kphi2,       &
    !$omp        r0, phi10, phi20, coord, list,        &
    !$omp        h, s, Tm, sol_T, force)               &
    !$omp reduction(+:estack) reduction(+:virial)
    do i = istart, iend

      i1 = list(1, i)
      i2 = list(2, i)
      i3 = list(3, i)
      i4 = list(4, i)
      i5 = list(5, i)
      i6 = list(6, i)
      i7 = list(7, i)

      v21(1:3) = vpbc(coord(1:3, i2) - coord(1:3, i1))
      v34(1:3) = vpbc(coord(1:3, i3) - coord(1:3, i4))
      v54(1:3) = vpbc(coord(1:3, i5) - coord(1:3, i4))
      v56(1:3) = vpbc(coord(1:3, i5) - coord(1:3, i6))
      v76(1:3) = vpbc(coord(1:3, i7) - coord(1:3, i6))

      ediv = 1.0e0_wp

      !===== Distance between 1 and 2 =====
      dist = norm2(v21)
      ddist = dist - r0(i)
      ediv = ediv + Kr(i) * ddist**2

      f_i(:) = 2.0e0_wp * Kr(i) * ddist / dist * v21(:)
      for(:, 1) = - f_i(:)
      for(:, 2) = + f_i(:)

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

      for(:, 3) = f_i(:)
      for(:, 4) = (-1.0e0_wp + d3454over5454) * f_i(:) &
                 -(            d5654over5454) * f_l(:)
      for(:, 5) = (-1.0e0_wp + d5654over5454) * f_l(:) &
                 -(            d3454over5454) * f_i(:)
      for(:, 6) = f_l(:)

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

      for(:, 7) =           + f_i(:)
      for(:, 6) = for(:, 6) + (-1.0e0_wp + d7656over5656) * f_i(:) &
                            - (            d5456over5656) * f_l(:)
      for(:, 5) = for(:, 5) - (            d7656over5656) * f_i(:) &
                            + (-1.0e0_wp + d5456over5656) * f_l(:)
      for(:, 4) = for(:, 4) + f_l(:)

      !===== Total =====
      u0 = -h(i) + KBOLTZ * (sol_T - Tm(i)) * s(i)
      estack = estack + u0 / ediv

      for(:, :) = u0 / ediv**2 * for(:, :)

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

#ifdef OMP
      ithread = omp_get_thread_num() + 1
#else
      ithread = 1
#endif
      force(1:3, i1, ithread) = force(1:3, i1, ithread) + for(1:3, 1)
      force(1:3, i2, ithread) = force(1:3, i2, ithread) + for(1:3, 2)
      force(1:3, i3, ithread) = force(1:3, i3, ithread) + for(1:3, 3)
      force(1:3, i4, ithread) = force(1:3, i4, ithread) + for(1:3, 4)
      force(1:3, i5, ithread) = force(1:3, i5, ithread) + for(1:3, 5)
      force(1:3, i6, ithread) = force(1:3, i6, ithread) + for(1:3, 6)
      force(1:3, i7, ithread) = force(1:3, i7, ithread) + for(1:3, 7)

    end do
    !$omp end parallel do

    call timer(TimerTISLocalStack, TimerOff)

    return

    contains

    pure function vpbc(v)
      real(wp) :: vpbc(3)
      real(wp), intent(in) :: v(3)

      vpbc(1:3) = v(1:3)
      if (v(1) > half_bsize(1)) then
        vpbc(1) = v(1) - bsize(1)
      else if (v(1) < -half_bsize(1)) then
        vpbc(1) = v(1) + bsize(1)
      end if
      if (v(2) > half_bsize(2)) then
        vpbc(2) = v(2) - bsize(2)
      else if (v(2) < -half_bsize(2)) then
        vpbc(2) = v(2) + bsize(2)
      end if
      if (v(3) > half_bsize(3)) then
        vpbc(3) = v(3) - bsize(3)
      else if (v(3) < -half_bsize(3)) then
        vpbc(3) = v(3) + bsize(3)
      end if
    end function vpbc

  end subroutine compute_energy_tis_lstack_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_tis_mwca
  !> @brief        calculate mWCA energy with pairlist (NOBC)
  !! @authors      NH
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] enemwca  : mWCA energy of target systems
  !! @note         TIS modified Weeks-Chandler-Andersen potential
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_tis_mwca(enefunc, molecule, pairlist, &
                                     coord, force, virial, enemwca)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: enemwca

    ! local variables
    real(wp)                  :: a, a2, dij(3), dist
    real(wp)                  :: dr, adr2, adr4, adr8
    real(wp)                  :: grad(3), dv_dr, eps, D
    integer                   :: i, j, k, l, natom, id
    integer                   :: num_mwca, ini_mwca, fin_mwca
    integer                   :: omp_get_thread_num

    integer, pointer          :: num_mwca_calc(:,:), mwca_list(:,:)
    real(wp), pointer         :: mwca_D(:,:), mwca_eps(:,:)
    integer, pointer          :: atomtype(:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerTISmWCA, TimerOn)
    ! use pointers
    !
    natom          =  molecule%num_atoms
    num_mwca_calc  => pairlist%num_tis_mwca_calc
    mwca_list      => pairlist%tis_mwca_list
    mwca_eps       => enefunc%tis_mwca_eps
    mwca_D         => enefunc%tis_mwca_D
    atomtype      => enefunc%atom_cls
    num_mwca       = 0

    a = enefunc%tis_mwca_a
    a2 = a*a

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                    &
    !$omp firstprivate(num_mwca)                                    &
    !$omp private(ini_mwca, fin_mwca, i, k, j, l, dij, grad, dist,  &
    !$omp         D, eps, dr, adr2, adr4, adr8, dv_dr, id)          &
    !$omp shared(natom, num_mwca_calc, mwca_list, atomtype,         &
    !$omp        mwca_eps, mwca_D, coord, force, a, a2)             &
    !$omp reduction(+:virial) reduction(+:enemwca)
    !
#ifdef OMP
    id      = omp_get_thread_num() + 1
#else
    id      = 1
#endif
    do i = 1, natom-1

      ini_mwca = num_mwca + 1
      fin_mwca = num_mwca + num_mwca_calc(i, id)
      num_mwca = fin_mwca

      do k = ini_mwca, fin_mwca

        j = mwca_list(k,id)

        D = mwca_D(atomtype(i), atomtype(j))

        ! compute distance
        dij(1:3) = coord(1:3,i) - coord(1:3,j)
        dist = norm2(dij)

        ! cutoff
        !
        if (dist >= D) cycle

        eps = mwca_eps(atomtype(i), atomtype(j))

        dr = dist + a - D
        adr2 = a2 / (dr*dr)
        adr4 = adr2 * adr2
        adr8 = adr4 * adr4

        dv_dr = abs(12.0e0_wp * eps * (adr2*adr4*adr8 - adr8) * dr / a2 / dist)

        if (dv_dr > 50.0_wp) then
          dv_dr = 50.0_wp
        end if

        grad(1:3) = dv_dr * dij(1:3)
        force(1:3,i,id) = force(1:3,i,id) + grad(1:3)
        force(1:3,j,id) = force(1:3,j,id) - grad(1:3)

        !if (dist > d_inf) then
        enemwca = enemwca + eps * (adr8*adr4 - 2*adr4*adr2 + 1.0e0_wp)
        !print *, 'i, j, dist, D, a, eps, ene', i, j, dist, D, a, eps, eps * (adr8*adr4 - 2*adr4*adr2 + 1.0e0_wp)
        !print *, 'dv_dr', dv_dr
        !print *, 'coord(i)', coord(1:3, i)
        !print *, 'coord(j)', coord(1:3, j)
        !print *, 'force(i)', grad(1:3)
        !print *, 'force(j)', -grad(1:3)
        !flush(6)
        !else
        !  enemwca = enemwca + 1.0e10_wp ! High energy (to reject in Widom method)
        !endif

        ! virial
        !
        !do l = 1, 3
        !  virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
        !end do

      end do

    end do
    !$omp end parallel

    call timer(TimerTISmWCA, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_tis_mwca

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_tis_mwca_pbc
  !> @brief        calculate mWCA energy with pairlist (PBC)
  !! @authors      NH (with reference to compute_energy_general_exv_AICG2P_pbc)
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] enemwca  : mWCA energy of target systems
  !! @note         TIS modified Weeks-Chandler-Andersen potential
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

subroutine compute_energy_tis_mwca_pbc(enefunc, boundary, pairlist, &
                                     coord, force, virial, enemwca)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: enemwca

    ! local variables
    real(wp)                  :: a, a2, dij(3), dist
    real(wp)                  :: bsize(3), coord_i(3), force_tmp(3)
    real(wp)                  :: dr, adr2, adr4, adr8
    real(wp)                  :: grad(3), dv_dr, eps, D
    real(wp)                  :: ene_omp(nthread), ene_omp_tmp
    integer                   :: i, j, k, l, n_atom, id
    integer                   :: i1, i2, i3, k1
    integer                   :: num_mwca, ini_mwca, fin_mwca
    integer                   :: omp_get_thread_num

    integer, pointer          :: num_mwca_calc(:,:), mwca_list(:,:)
    real(wp), pointer         :: mwca_D(:,:), mwca_eps(:,:)
    integer, pointer          :: atomtype(:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerTISmWCA, TimerOn)

    n_atom          = enefunc%num_cg_particle_TIS_all
    num_mwca_calc  => pairlist%num_tis_mwca_calc
    mwca_list      => pairlist%tis_mwca_list
    mwca_eps       => enefunc%tis_mwca_eps
    mwca_D         => enefunc%tis_mwca_D
    atomtype      => enefunc%atom_cls

    bsize(1)       =  boundary%box_size_x
    bsize(2)       =  boundary%box_size_y
    bsize(3)       =  boundary%box_size_z

    num_mwca       = 0
    ene_omp(1:nthread) = 0.0_wp

    a = enefunc%tis_mwca_a
    a2 = a*a

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                    &
    !$omp firstprivate(num_mwca)                                    & 
    !$omp private(ini_mwca, fin_mwca, i, k, j, l, dij, grad, dist,  &
    !$omp         i1, i2, i3, k1,                                   &
    !$omp         D, eps, dr, adr2, adr4, adr8, dv_dr, id,          &
    !$omp         ene_omp_tmp, coord_i, force_tmp)                  &
    !$omp shared(n_atom, num_mwca_calc, mwca_list, bsize, atomtype,  &
    !$omp        mwca_eps, mwca_D, coord, force, a, a2, ene_omp)
    !
#ifdef OMP
    id      = omp_get_thread_num() + 1
#else
    id      = 1
#endif
    do i = 1, n_atom-1

      !proceed = .false.

      ini_mwca = num_mwca + 1
      fin_mwca = num_mwca + num_mwca_calc(i, id)
      num_mwca = fin_mwca

      !if (fin_mwca >= ini_mwca) proceed = .true.
      !if (proceed) then
      if (fin_mwca >= ini_mwca) then

        coord_i(1:3)        = coord(1:3, i)
        force_tmp(1:3)      = 0.0_wp
        !virial_tmp(1:3,1:3) = 0.0_wp
        ene_omp_tmp         = 0.0_wp

        do k = ini_mwca, fin_mwca

          k1 = mwca_list(k, id)
          j  = k1 / 27
          k1 = k1 - j*27
          i3 = k1 / 9
          k1 = k1 - i3*9
          i2 = k1 / 3
          i1 = k1 - i2*3
          i1 = i1 - 1
          i2 = i2 - 1
          i3 = i3 - 1

          D = mwca_D(atomtype(i), atomtype(j))

          ! compute distance
          dij(1)  = coord_i(1) - coord(1,j) + bsize(1)*real(i1,wp)
          dij(2)  = coord_i(2) - coord(2,j) + bsize(2)*real(i2,wp)
          dij(3)  = coord_i(3) - coord(3,j) + bsize(3)*real(i3,wp)
          !dij(1:3) = coord(1:3,i) - coord(1:3,j)
          dist = norm2(dij)

          ! cutoff
          !
          if (dist >= D) cycle

          eps = mwca_eps(atomtype(i), atomtype(j))

          dr = dist + a - D
          adr2 = a2 / (dr*dr)
          adr4 = adr2 * adr2
          adr8 = adr4 * adr4

          dv_dr = abs(12.0e0_wp * eps * (adr2*adr4*adr8 - adr8) * dr / a2 / dist)

          if (dv_dr > 50.0_wp) then
            dv_dr = 50.0_wp
          end if

          grad(1:3) = dv_dr * dij(1:3)
          !force(1:3,i,id) = force(1:3,i,id) + for(1:3)
          !force(1:3,j,id) = force(1:3,j,id) - for(1:3)

          !if (dist > d_inf) then
          ene_omp_tmp = ene_omp_tmp + eps * (adr8*adr4 - 2*adr4*adr2 + 1.0e0_wp)
          !print *, 'i, j, dist, D, a, eps, ene', i, j, dist, D, a, eps, eps * (adr8*adr4 - 2*adr4*adr2 + 1.0e0_wp)
          !print *, 'dv_dr', dv_dr
          !print *, 'coord(i)', coord_i(1:3)
          !print *, 'coord(j)', coord(1:3, j)
          !print *, 'force(i)', grad(1:3)
          !print *, 'force(j)', -grad(1:3)
          !flush(6)
          !else
          !  enemwca = enemwca + 1.0e10_wp ! High energy (to reject in Widom method)
          !endif

          ! virial
          !
          !do l = 1, 3
          !  virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
          !end do

          ! store force
          !
          force_tmp(1) = force_tmp(1) + grad(1)
          force_tmp(2) = force_tmp(2) + grad(2)
          force_tmp(3) = force_tmp(3) + grad(3)
          force(1,j,id) = force(1,j,id) - grad(1)
          force(2,j,id) = force(2,j,id) - grad(2)
          force(3,j,id) = force(3,j,id) - grad(3)

          ! virial
          !
          !virial_tmp(1,1) = virial_tmp(1,1) + dij(1)*grad(1)
          !virial_tmp(2,1) = virial_tmp(2,1) + dij(2)*grad(1)
          !virial_tmp(3,1) = virial_tmp(3,1) + dij(3)*grad(1)
          !virial_tmp(1,2) = virial_tmp(1,2) + dij(1)*grad(2)
          !virial_tmp(2,2) = virial_tmp(2,2) + dij(2)*grad(2)
          !virial_tmp(3,2) = virial_tmp(3,2) + dij(3)*grad(2)
          !virial_tmp(1,3) = virial_tmp(1,3) + dij(1)*grad(3)
          !virial_tmp(2,3) = virial_tmp(2,3) + dij(2)*grad(3)
          !virial_tmp(3,3) = virial_tmp(3,3) + dij(3)*grad(3)
        end do

        force(1:3,i,id) = force(1:3,i,id) + force_tmp(1:3)
        ene_omp(id) = ene_omp(id) + ene_omp_tmp
        !virial_omp(1:3,1:3,id) = virial_omp(1:3,1:3,id) + virial_tmp(1:3,1:3)

      end if
    end do
    !$omp end parallel

    do i = 1, nthread
      enemwca = enemwca + ene_omp(i)
      !virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,i)
    end do

    call timer(TimerTISmWCA, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_tis_mwca_pbc

end module at_energy_tis_mod
