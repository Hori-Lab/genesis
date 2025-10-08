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
  public  :: compute_energy_tis_hb
  public  :: compute_energy_tis_hb_pbc
  public  :: compute_energy_tis_harmonic_hb_pbc
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
    integer                  :: i, id 
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
    !$omp parallel default(none)                       &
    !$omp private(id, i, i1, i2, i3, i4, i5, i6, i7,   &
    !$omp         v21, v34, v54, v56, v76,             &
    !$omp         ediv, dist, ddist, dih, d,           &
    !$omp         f_i, f_l, for,                       &
    !$omp         abs54, d5454, abs56, d5656, d5654,   &
    !$omp         d7656over5656, d5456over5656,        &
    !$omp         d3454over5454, d5654over5454,        &
    !$omp         m, n, dnn, dmm, u0)                  &
    !$omp shared(istart, iend, Kr, Kphi1, Kphi2,       &
    !$omp        r0, phi10, phi20, coord, list,        &
    !$omp        h, s, Tm, sol_T, force, nthread)      &
    !$omp reduction(+:estack) reduction(+:virial)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread

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

      ediv = 1.0_wp

      !===== Distance between 1 and 2 =====
      dist = norm2(v21)
      ddist = dist - r0(i)
      ediv = ediv + Kr(i) * ddist**2

      f_i(:) = 2.0_wp * Kr(i) * ddist / dist * v21(:)
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

      f_i(:) = + 2.0_wp * Kphi1(i) * d * abs54 / dmm * m(:)
      f_l(:) = - 2.0_wp * Kphi1(i) * d * abs54 / dnn * n(:)

      for(:, 3) = f_i(:)
      for(:, 4) = (-1.0_wp + d3454over5454) * f_i(:) &
                 -(            d5654over5454) * f_l(:)
      for(:, 5) = (-1.0_wp + d5654over5454) * f_l(:) &
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

      f_i(:) = + 2.0_wp * Kphi2(i) * d * abs56 / dmm * m(:)
      f_l(:) = - 2.0_wp * Kphi2(i) * d * abs56 / dnn * n(:)

      for(:, 7) =           + f_i(:)
      for(:, 6) = for(:, 6) + (-1.0_wp + d7656over5656) * f_i(:) &
                            - (            d5456over5656) * f_l(:)
      for(:, 5) = for(:, 5) - (            d7656over5656) * f_i(:) &
                            + (-1.0_wp + d5456over5656) * f_l(:)
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

      force(1:3, i1, id+1) = force(1:3, i1, id+1) + for(1:3, 1)
      force(1:3, i2, id+1) = force(1:3, i2, id+1) + for(1:3, 2)
      force(1:3, i3, id+1) = force(1:3, i3, id+1) + for(1:3, 3)
      force(1:3, i4, id+1) = force(1:3, i4, id+1) + for(1:3, 4)
      force(1:3, i5, id+1) = force(1:3, i5, id+1) + for(1:3, 5)
      force(1:3, i6, id+1) = force(1:3, i6, id+1) + for(1:3, 6)
      force(1:3, i7, id+1) = force(1:3, i7, id+1) + for(1:3, 7)
    end do

    !$omp end parallel
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
    integer                  :: i, id
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
    real(wp)                 :: bsize(3), inv_bsize(3)

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
    inv_bsize = 1.0_wp/bsize(1:3)

    ! calculation of local stacking energy and gradient
    !
    !$omp parallel default(none)                       &
    !$omp private(id, i, i1, i2, i3, i4, i5, i6, i7,   &
    !$omp         v21, v34, v54, v56, v76,             &
    !$omp         ediv, dist, ddist, dih, d,           &
    !$omp         f_i, f_l, for,                       &
    !$omp         abs54, d5454, abs56, d5656, d5654,   &
    !$omp         d7656over5656, d5456over5656,        &
    !$omp         d3454over5454, d5654over5454,        &
    !$omp         m, n, dnn, dmm, u0)                  &
    !$omp shared(istart, iend, Kr, Kphi1, Kphi2,       &
    !$omp        r0, phi10, phi20, coord, list,        &
    !$omp        h, s, Tm, sol_T, force, inv_bsize,    &
    !$omp        nthread, bsize)                       &
    !$omp reduction(+:estack) reduction(+:virial)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread

      i1 = list(1, i)
      i2 = list(2, i)
      i3 = list(3, i)
      i4 = list(4, i)
      i5 = list(5, i)
      i6 = list(6, i)
      i7 = list(7, i)

      v21(1:3) = coord(1:3, i2) - coord(1:3, i1)
      v21(1:3) = v21(1:3) - bsize(1:3)*anint(v21(1:3)*inv_bsize(1:3))

      v34(1:3) = coord(1:3, i3) - coord(1:3, i4)
      v34(1:3) = v34(1:3) - bsize(1:3)*anint(v34(1:3)*inv_bsize(1:3))

      v54(1:3) = coord(1:3, i5) - coord(1:3, i4)
      v54(1:3) = v54(1:3) - bsize(1:3)*anint(v54(1:3)*inv_bsize(1:3))

      v56(1:3) = coord(1:3, i5) - coord(1:3, i6)
      v56(1:3) = v56(1:3) - bsize(1:3)*anint(v56(1:3)*inv_bsize(1:3))

      v76(1:3) = coord(1:3, i7) - coord(1:3, i6)
      v76(1:3) = v76(1:3) - bsize(1:3)*anint(v76(1:3)*inv_bsize(1:3))

      ediv = 1.0_wp

      !===== Distance between 1 and 2 =====
      dist = norm2(v21)
      ddist = dist - r0(i)
      ediv = ediv + Kr(i) * ddist**2

      f_i(:) = 2.0_wp * Kr(i) * ddist / dist * v21(:)
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

      f_i(:) = + 2.0_wp * Kphi1(i) * d * abs54 / dmm * m(:)
      f_l(:) = - 2.0_wp * Kphi1(i) * d * abs54 / dnn * n(:)

      for(:, 3) = f_i(:)
      for(:, 4) = (-1.0_wp + d3454over5454) * f_i(:) &
                 -(            d5654over5454) * f_l(:)
      for(:, 5) = (-1.0_wp + d5654over5454) * f_l(:) &
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

      f_i(:) = + 2.0_wp * Kphi2(i) * d * abs56 / dmm * m(:)
      f_l(:) = - 2.0_wp * Kphi2(i) * d * abs56 / dnn * n(:)

      for(:, 7) =           + f_i(:)
      for(:, 6) = for(:, 6) + (-1.0_wp + d7656over5656) * f_i(:) &
                            - (            d5456over5656) * f_l(:)
      for(:, 5) = for(:, 5) - (            d7656over5656) * f_i(:) &
                            + (-1.0_wp + d5456over5656) * f_l(:)
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

      force(1:3, i1, id+1) = force(1:3, i1, id+1) + for(1:3, 1)
      force(1:3, i2, id+1) = force(1:3, i2, id+1) + for(1:3, 2)
      force(1:3, i3, id+1) = force(1:3, i3, id+1) + for(1:3, 3)
      force(1:3, i4, id+1) = force(1:3, i4, id+1) + for(1:3, 4)
      force(1:3, i5, id+1) = force(1:3, i5, id+1) + for(1:3, 5)
      force(1:3, i6, id+1) = force(1:3, i6, id+1) + for(1:3, 6)
      force(1:3, i7, id+1) = force(1:3, i7, id+1) + for(1:3, 7)

    end do
    !$omp end parallel

    call timer(TimerTISLocalStack, TimerOff)

    return


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
    integer                   :: i, j, k, l, n_atom, id, my_id
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
    !$omp         ene_omp_tmp, coord_i, force_tmp,  my_id)          &
    !$omp shared(n_atom, num_mwca_calc, mwca_list, bsize, atomtype, &
    !$omp        mwca_eps, mwca_D, coord, force, a, a2, ene_omp,    &
    !$omp        my_city_rank, nthread)   
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1
    
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
          dij(1)  = coord_i(1) - coord(1,j) - bsize(1)*real(i1,wp)
          dij(2)  = coord_i(2) - coord(2,j) - bsize(2)*real(i2,wp)
          dij(3)  = coord_i(3) - coord(3,j) - bsize(3)*real(i3,wp)
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

          dv_dr = abs(12.0_wp * eps * (adr2*adr4*adr8 - adr8) * dr / a2 / dist)

          if (dv_dr > 50.0_wp) then
            dv_dr = 50.0_wp
          end if

          grad(1:3) = dv_dr * dij(1:3)
          !force(1:3,i,id) = force(1:3,i,id) + for(1:3)
          !force(1:3,j,id) = force(1:3,j,id) - for(1:3)

          !if (dist > d_inf) then
          ene_omp_tmp = ene_omp_tmp + eps * (adr8*adr4 - 2*adr4*adr2 + 1.0_wp)
          !print *, 'i, j, dist, D, a, eps, ene', i, j, dist, D, a, eps, eps * (adr8*adr4 - 2*adr4*adr2 + 1.0_wp)
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

  subroutine compute_energy_tis_hb(enefunc, coord, force, virial, enehb)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: enehb

    ! local variables
    integer                  :: i, id 
    integer                  :: i_atom, j_atom, i1_atom, j1_atom, i2_atom, j2_atom
    integer                  :: istart, iend
    real(wp)                 :: v12(1:3), v13(1:3), v53(1:3), v42(1:3), v46(1:3)
    real(wp)                 :: d, dih 
    real(wp)                 :: a42, a13, a12
    real(wp)                 :: cos_theta124, cos_theta312, cos_theta531, cos_theta246
    real(wp)                 :: d1212, d1313, d4242
    real(wp)                 :: d1213, d1242, d4246, d1353
    real(wp)                 :: d1213over1212, d1213over1313
    real(wp)                 :: d1242over1212, d1242over4242
    real(wp)                 :: d4246over4242, d1353over1313
    real(wp)                 :: n(3)
    real(wp)                 :: c4212(3), c1213(3), c5313(3), c4246(3)
    real(wp)                 :: c4212_abs2, c1213_abs2, c5313_abs2, c4246_abs2
    real(wp)                 :: pre, f_i(3), f_k(3), f_l(3), ex, ene, for_hb(1:3, 1:6)

    integer,         pointer :: list(:,:)
    real(wp),        pointer :: U0(:)
    real(wp),        pointer :: dist_eq(:), dist_coef(:)
    real(wp),        pointer :: ang1_eq(:), ang2_eq(:), ang_coef(:)
    real(wp),        pointer :: dih_eq(:), dih1_eq(:), dih2_eq(:), dih_coef(:)

    integer                  :: omp_get_thread_num

    call timer(TimerTISHB, TimerOn)

    istart       = enefunc%istart_tis_hb
    iend         = enefunc%iend_tis_hb

    ! use pointers
    !
    list         => enefunc%tis_hb_list
    U0           => enefunc%tis_hb_dist_U0
    dist_eq      => enefunc%tis_hb_dist_eq
    dist_coef    => enefunc%tis_hb_dist_coef
    ang1_eq      => enefunc%tis_hb_angle_ang1
    ang2_eq      => enefunc%tis_hb_angle_ang2
    ang_coef     => enefunc%tis_hb_angle_coef
    dih_eq       => enefunc%tis_hb_dihedral_dih
    dih1_eq      => enefunc%tis_hb_dihedral_dih1
    dih2_eq      => enefunc%tis_hb_dihedral_dih2
    dih_coef     => enefunc%tis_hb_dihedral_coef


    ! calculation of local stacking energy and gradient
    !
    !$omp parallel default(none)                       &
    !$omp private(id, i, i_atom, j_atom, i1_atom,      &
    !$omp         j1_atom, i2_atom, j2_atom,           &
    !$omp         d,dih,cos_theta124,cos_theta312,     &
    !$omp         cos_theta531,cos_theta246,           &
    !$omp         v12,v13,v53,v42,v46,a12,a13,a42,     &
    !$omp         d1212,d1313, for_hb, ex, ene,        &
    !$omp         pre, f_i, f_k, f_l, n,               &
    !$omp         d4242,d1213,d1242,d4246,d1353,       &
    !$omp         d1213over1212,d1213over1313,         &
    !$omp         d1242over1212,d1242over4242,         &
    !$omp         d4246over4242,d1353over1313,c4212,   &
    !$omp         c1213,c4212_abs2,c1213_abs2,         &
    !$omp         c5313, c5313_abs2, c4246, c4246_abs2)&
    !$omp shared(istart, iend, U0, dist_eq, dist_coef, &
    !$omp        coord, list, ang1_eq, ang2_eq,        &
    !$omp        ang_coef, dih_eq, dih1_eq, dih2_eq,   &
    !$omp        dih_coef, force, nthread)             &
    !$omp reduction(+:enehb) reduction(+:virial)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread
      for_hb(:,:) = 0.0_wp

      i_atom  = list(1, i) ! 1
      j_atom  = list(2, i) ! 2
      i1_atom = list(3, i) ! 3
      j1_atom = list(4, i) ! 4
      i2_atom = list(5, i) ! 5
      j2_atom = list(6, i) ! 6

      !===== Distance =====!
      v12(1:3) = coord(1:3, i_atom)  - coord(1:3, j_atom)
      d1212 = dot_product(v12,v12)
      a12 = norm2(v12)

      d = a12 - dist_eq(i)

      ! Cutoff
      if (abs(d) > 2.0_wp) then
        cycle
      end if

      ex = - dist_coef(i) * d**2
      f_i(:) = (2.0_wp * dist_coef(i) * d / a12) * v12(:)
      for_hb(:,1) = + f_i(:)
      for_hb(:,2) = - f_i(:)

      v13(1:3) = coord(1:3, i_atom)  - coord(1:3, i1_atom)
      v53(1:3) = coord(1:3, i2_atom) - coord(1:3, i1_atom)
      v42(1:3) = coord(1:3, j1_atom) - coord(1:3, j_atom)
      v46(1:3) = coord(1:3, j1_atom) - coord(1:3, j2_atom)

      d1313 = dot_product(v13,v13)
      d4242 = dot_product(v42,v42)
      a13 = sqrt(d1313)
      a42 = sqrt(d4242)
      d1213 = dot_product(v13,v12)
      d1242 = dot_product(v12,v42)
      d4246 = dot_product(v42,v46)
      d1353 = dot_product(v13,v53)
      d1213over1212 = d1213 / d1212
      d1213over1313 = d1213 / d1313
      d1242over1212 = d1242 / d1212
      d4246over4242 = d4246 / d4242
      d1242over4242 = d1242 / d4242
      d1353over1313 = d1353 / d1313

      !===== Angle of 3-1=2  =====!
      ! ENERGY
      cos_theta312 = d1213 / (a13 * a12)
      d = acos(cos_theta312) - ang1_eq(i)
      ex = ex - ang_coef(i) * d**2

      ! FORCE
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta312) > 0.99995_wp) then
          d1213 = sign(a12 * a13 * 0.99995_wp, cos_theta312)
      end if

      pre = 2.0_wp * ang_coef(i) * d / sqrt(d1313*d1212 - d1213**2)
      f_i(:) = pre * (v12(:) - (d1213over1313 * v13(:)))
      f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))
      for_hb(:,3) = for_hb(:,3) + f_i(:)
      for_hb(:,2) = for_hb(:,2) + f_k(:)
      for_hb(:,1) = for_hb(:,1) - f_i(:) - f_k(:)

      !===== Angle of 1=2-4  =====!
      ! ENERGY
      cos_theta124 = d1242 / (a12 * a42)
      d = acos(cos_theta124) - ang2_eq(i)
      ex = ex - ang_coef(i) * d**2

      ! FORCE
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta124) > 0.99995_wp) then
          d1242 = sign(a12 * a42 * 0.99995_wp, cos_theta124)
      endif

      pre = 2.0_wp * ang_coef(i) * d / sqrt(d1212*d4242 - d1242**2)
      f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1242over4242 * v42(:)))
      for_hb(:,1) = for_hb(:,1) + f_i(:)
      for_hb(:,4) = for_hb(:,4) + f_k(:)
      for_hb(:,2) = for_hb(:,2) - f_i(:) - f_k(:)

      !===== Dihedral angle among 4-2=1=3 =====!
      c4212(1) = v42(2)*v12(3) - v42(3)*v12(2)
      c4212(2) = v42(3)*v12(1) - v42(1)*v12(3)
      c4212(3) = v42(1)*v12(2) - v42(2)*v12(1)
      c1213(1) = v12(2)*v13(3) - v12(3)*v13(2)
      c1213(2) = v12(3)*v13(1) - v12(1)*v13(3)
      c1213(3) = v12(1)*v13(2) - v12(2)*v13(1)
      c4212_abs2 = dot_product(c4212,c4212)
      c1213_abs2 = dot_product(c1213,c1213)

      ! ENERGY
      dih = atan2(dot_product(v42,c1213)*a12, dot_product(c4212,c1213))
      d = dih - dih_eq(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ex = ex - dih_coef(i) * d**2

      ! FORCE
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta124) > 0.99995_wp) then
          ! |c4212|^2 = |v42|^2 * |v12|^2 * sin(theta)^2
          c4212_abs2 = d4242 * d1212 * (1.0 - 0.99995_wp**2)
      endif
      if (abs(cos_theta312) > 0.99995_wp) then
          c1213_abs2 = d1212 * d1313 * (1.0_wp - 0.99995_wp**2)
      endif

      pre = 2.0_wp * dih_coef(i) * d * a12
      f_i(:) = + pre / c4212_abs2 * c4212(:)
      f_l(:) = - pre / c1213_abs2 * c1213(:)
    
      for_hb(:,4) = for_hb(:,4) + f_i(:)
      for_hb(:,2) = for_hb(:,2) + (-1.0_wp + d1242over1212) * f_i(:) &
                                        - (  d1213over1212) * f_l(:)
      for_hb(:,1) = for_hb(:,1) + (-1.0_wp + d1213over1212) * f_l(:) &
                                        - (  d1242over1212) * f_i(:)
      for_hb(:,3) = for_hb(:,3) + f_l(:)

      !===== Dihedral angle among 5-3-1=2 =====!
      c5313(1) = v53(2) * v13(3) - v53(3) * v13(2)
      c5313(2) = v53(3) * v13(1) - v53(1) * v13(3)
      c5313(3) = v53(1) * v13(2) - v53(2) * v13(1)
      n(:) = -c1213(:)

      ! ENERGY
      dih = atan2(dot_product(v53,n)*a13 , dot_product(c5313,n))
      d = dih - dih1_eq(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ex = ex - dih_coef(i) * d**2

      ! FORCE
      c5313_abs2 = dot_product(c5313, c5313)

      cos_theta531 = d1353 / (a13 * norm2(v53))
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta531) > 0.99995_wp) then
          c5313_abs2 = dot_product(v53,v53) * d1313 * (1.0 - 0.99995_wp**2)
      endif

      pre = 2.0_wp * dih_coef(i) * d * a13
      f_i(:) = + pre / c5313_abs2 * c5313(:)
      f_l(:) = - pre / c1213_abs2 * n(:)
    
      for_hb(:,5) = for_hb(:,5) + f_i(:)
      for_hb(:,3) = for_hb(:,3) + (-1.0_wp + d1353over1313) * f_i(:) &
                                        - (  d1213over1313) * f_l(:)
      for_hb(:,1) = for_hb(:,1) + (-1.0_wp + d1213over1313) * f_l(:) &
                                        - (  d1353over1313) * f_i(:)
      for_hb(:,2) = for_hb(:,2) + f_l(:)

      !===== Dihedral angle among 1=2-4-6 =====!
      n(:) = -c4212(:)
      c4246(1) = v42(2) * v46(3) - v42(3) * v46(2)
      c4246(2) = v42(3) * v46(1) - v42(1) * v46(3)
      c4246(3) = v42(1) * v46(2) - v42(2) * v46(1)
      c4246_abs2 = dot_product(c4246, c4246)

      ! ENERGY
      dih = atan2(dot_product(v12,c4246)*a42 , dot_product(n,c4246))
      d = dih - dih2_eq(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ex = ex - dih_coef(i) * d**2

      ! FORCE
      cos_theta246 = d4246 / (a42 * norm2(v46))
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta246) > 0.99995_wp) then
          c4246_abs2 = dot_product(v46,v46) * d4242 * (1.0 - 0.99995_wp**2)
      endif

      pre = 2.0_wp * dih_coef(i) * d * a42
      f_i(:) = + pre / c4212_abs2 * n(:)
      f_l(:) = - pre / c4246_abs2 * c4246(:)
    
      for_hb(:,1) = for_hb(:,1) + f_i(:)
      for_hb(:,2) = for_hb(:,2) + (-1.0_wp + d1242over4242) * f_i(:) &
                                        - (      d4246over4242) * f_l(:)
      for_hb(:,4) = for_hb(:,4) + (-1.0_wp + d4246over4242) * f_l(:) &
                                        - (      d1242over4242) * f_i(:)
      for_hb(:,6) = for_hb(:,6) + f_l(:)

      !===== Total =====!

      ex = U0(i) * exp(ex)
      enehb = enehb + ex
      ! check this
      for_hb(:,:) = enehb * for_hb(:,:)

      force(1:3, i_atom,  id+1) = force(1:3, i_atom,  id+1) + for_hb(1:3, 1)
      force(1:3, j_atom,  id+1) = force(1:3, j_atom,  id+1) + for_hb(1:3, 2)
      force(1:3, i1_atom, id+1) = force(1:3, i1_atom, id+1) + for_hb(1:3, 3)
      force(1:3, j1_atom, id+1) = force(1:3, j1_atom, id+1) + for_hb(1:3, 4)
      force(1:3, i2_atom, id+1) = force(1:3, i2_atom, id+1) + for_hb(1:3, 5)
      force(1:3, j2_atom, id+1) = force(1:3, j2_atom, id+1) + for_hb(1:3, 6)

    end do
    !$omp end parallel

    call timer(TimerTISHB, TimerOff)

    return

  end subroutine compute_energy_tis_hb

  subroutine compute_energy_tis_hb_pbc(enefunc, boundary, coord, &
                                           force, virial, enehb)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_boundary),target, intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: enehb

    ! local variables
    integer                  :: i, id 
    integer                  :: i_atom, j_atom, i1_atom, j1_atom, i2_atom, j2_atom
    integer                  :: istart, iend
    real(wp)                 :: v12(1:3), v13(1:3), v53(1:3), v42(1:3), v46(1:3)
    real(wp)                 :: d, dih 
    real(wp)                 :: a42, a13, a12
    real(wp)                 :: cos_theta124, cos_theta312, cos_theta531, cos_theta246
    real(wp)                 :: d1212, d1313, d4242
    real(wp)                 :: d1213, d1242, d4246, d1353
    real(wp)                 :: d1213over1212, d1213over1313
    real(wp)                 :: d1242over1212, d1242over4242
    real(wp)                 :: d4246over4242, d1353over1313
    real(wp)                 :: n(3)
    real(wp)                 :: c4212(3), c1213(3), c5313(3), c4246(3)
    real(wp)                 :: c4212_abs2, c1213_abs2, c5313_abs2, c4246_abs2
    real(wp)                 :: pre, f_i(3), f_k(3), f_l(3), ex, ene, for_hb(1:3, 1:6)
    real(wp)                 :: bsize(3), inv_bsize(3)

    integer,         pointer :: list(:,:)
    real(wp),        pointer :: U0(:)
    real(wp),        pointer :: dist_eq(:), dist_coef(:)
    real(wp),        pointer :: ang1_eq(:), ang2_eq(:), ang_coef(:)
    real(wp),        pointer :: dih_eq(:), dih1_eq(:), dih2_eq(:), dih_coef(:)

    integer                  :: omp_get_thread_num

    call timer(TimerTISHB, TimerOn)

    istart       = enefunc%istart_tis_hb
    iend         = enefunc%iend_tis_hb

    ! use pointers
    !
    list         => enefunc%tis_hb_list
    U0           => enefunc%tis_hb_dist_U0
    dist_eq      => enefunc%tis_hb_dist_eq
    dist_coef    => enefunc%tis_hb_dist_coef
    ang1_eq      => enefunc%tis_hb_angle_ang1
    ang2_eq      => enefunc%tis_hb_angle_ang2
    ang_coef     => enefunc%tis_hb_angle_coef
    dih_eq       => enefunc%tis_hb_dihedral_dih
    dih1_eq      => enefunc%tis_hb_dihedral_dih1
    dih2_eq      => enefunc%tis_hb_dihedral_dih2
    dih_coef     => enefunc%tis_hb_dihedral_coef

    bsize(1) = boundary%box_size_x
    bsize(2) = boundary%box_size_y
    bsize(3) = boundary%box_size_z
    inv_bsize = 1.0_wp/bsize(1:3)

    ! calculation of local stacking energy and gradient
    !
    !$omp parallel default(none)                       &
    !$omp private(id, i, i_atom, j_atom, i1_atom,      &
    !$omp         j1_atom, i2_atom, j2_atom,           &
    !$omp         d,dih,cos_theta124,cos_theta312,     &
    !$omp         cos_theta531,cos_theta246,           &
    !$omp         v12,v13,v53,v42,v46,a12,a13,a42,     &
    !$omp         d1212,d1313, for_hb, ex, ene,        &
    !$omp         pre, f_i, f_k, f_l, n,               &
    !$omp         d4242,d1213,d1242,d4246,d1353,       &
    !$omp         d1213over1212,d1213over1313,         &
    !$omp         d1242over1212,d1242over4242,         &
    !$omp         d4246over4242,d1353over1313,c4212,   &
    !$omp         c1213,c4212_abs2,c1213_abs2,         &
    !$omp         c5313, c5313_abs2, c4246, c4246_abs2)&
    !$omp shared(istart, iend, U0, dist_eq, dist_coef, &
    !$omp        coord, list, ang1_eq, ang2_eq,        &
    !$omp        ang_coef, dih_eq, dih1_eq, dih2_eq,   &
    !$omp        dih_coef, force, nthread,             &
    !$omp        bsize, inv_bsize)                     &
    !$omp reduction(+:enehb) reduction(+:virial)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread
      for_hb(:,:) = 0.0_wp

      i_atom  = list(1, i) ! 1
      j_atom  = list(2, i) ! 2
      i1_atom = list(3, i) ! 3
      j1_atom = list(4, i) ! 4
      i2_atom = list(5, i) ! 5
      j2_atom = list(6, i) ! 6

      !===== Distance =====!
      v12(1:3) = coord(1:3, i_atom)  - coord(1:3, j_atom)
      v12(1:3) = v12(1:3) - bsize(1:3)*anint(v12(1:3)*inv_bsize(1:3))
      d1212 = dot_product(v12,v12)
      a12 = norm2(v12) 

      d = a12 - dist_eq(i)

      ! Cutoff
      if (abs(d) > 2.0_wp) then
        cycle
      end if

      ex = - dist_coef(i) * d**2
      f_i(:) = (2.0_wp * dist_coef(i) * d / a12) * v12(:)
      for_hb(:,1) = + f_i(:)
      for_hb(:,2) = - f_i(:)

      v13(1:3) = coord(1:3, i_atom)  - coord(1:3, i1_atom)
      v13(1:3) = v13(1:3) - bsize(1:3)*anint(v13(1:3)*inv_bsize(1:3))

      v53(1:3) = coord(1:3, i2_atom) - coord(1:3, i1_atom)
      v53(1:3) = v53(1:3) - bsize(1:3)*anint(v53(1:3)*inv_bsize(1:3))

      v42(1:3) = coord(1:3, j1_atom) - coord(1:3, j_atom)
      v42(1:3) = v42(1:3) - bsize(1:3)*anint(v42(1:3)*inv_bsize(1:3))

      v46(1:3) = coord(1:3, j1_atom) - coord(1:3, j2_atom)
      v46(1:3) = v46(1:3) - bsize(1:3)*anint(v46(1:3)*inv_bsize(1:3))

      d1313 = dot_product(v13,v13)
      d4242 = dot_product(v42,v42)
      a13 = sqrt(d1313)
      a42 = sqrt(d4242)
      d1213 = dot_product(v13,v12)
      d1242 = dot_product(v12,v42)
      d4246 = dot_product(v42,v46)
      d1353 = dot_product(v13,v53)
      d1213over1212 = d1213 / d1212
      d1213over1313 = d1213 / d1313
      d1242over1212 = d1242 / d1212
      d4246over4242 = d4246 / d4242
      d1242over4242 = d1242 / d4242
      d1353over1313 = d1353 / d1313

      !===== Angle of 3-1=2  =====!
      ! ENERGY
      cos_theta312 = d1213 / (a13 * a12)
      d = acos(cos_theta312) - ang1_eq(i)
      ex = ex - ang_coef(i) * d**2

      ! FORCE
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta312) > 0.99995_wp) then
          d1213 = sign(a12 * a13 * 0.99995_wp, cos_theta312)
      end if

      pre = 2.0_wp * ang_coef(i) * d / sqrt(d1313*d1212 - d1213**2)
      f_i(:) = pre * (v12(:) - (d1213over1313 * v13(:)))
      f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))
      for_hb(:,3) = for_hb(:,3) + f_i(:)
      for_hb(:,2) = for_hb(:,2) + f_k(:)
      for_hb(:,1) = for_hb(:,1) - f_i(:) - f_k(:)

      !===== Angle of 1=2-4  =====!
      ! ENERGY
      cos_theta124 = d1242 / (a12 * a42)
      d = acos(cos_theta124) - ang2_eq(i)
      ex = ex - ang_coef(i) * d**2

      ! FORCE
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta124) > 0.99995_wp) then
          d1242 = sign(a12 * a42 * 0.99995_wp, cos_theta124)
      endif

      pre = 2.0_wp * ang_coef(i) * d / sqrt(d1212*d4242 - d1242**2)
      f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1242over4242 * v42(:)))
      for_hb(:,1) = for_hb(:,1) + f_i(:)
      for_hb(:,4) = for_hb(:,4) + f_k(:)
      for_hb(:,2) = for_hb(:,2) - f_i(:) - f_k(:)

      !===== Dihedral angle among 4-2=1=3 =====!
      c4212(1) = v42(2)*v12(3) - v42(3)*v12(2)
      c4212(2) = v42(3)*v12(1) - v42(1)*v12(3)
      c4212(3) = v42(1)*v12(2) - v42(2)*v12(1)
      c1213(1) = v12(2)*v13(3) - v12(3)*v13(2)
      c1213(2) = v12(3)*v13(1) - v12(1)*v13(3)
      c1213(3) = v12(1)*v13(2) - v12(2)*v13(1)
      c4212_abs2 = dot_product(c4212,c4212)
      c1213_abs2 = dot_product(c1213,c1213)

      ! ENERGY
      dih = atan2(dot_product(v42,c1213)*a12, dot_product(c4212,c1213))
      d = dih - dih_eq(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ex = ex - dih_coef(i) * d**2

      ! FORCE
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta124) > 0.99995_wp) then
          ! |c4212|^2 = |v42|^2 * |v12|^2 * sin(theta)^2
          c4212_abs2 = d4242 * d1212 * (1.0 - 0.99995_wp**2)
      endif
      if (abs(cos_theta312) > 0.99995_wp) then
          c1213_abs2 = d1212 * d1313 * (1.0_wp - 0.99995_wp**2)
      endif

      pre = 2.0_wp * dih_coef(i) * d * a12
      f_i(:) = + pre / c4212_abs2 * c4212(:)
      f_l(:) = - pre / c1213_abs2 * c1213(:)
    
      for_hb(:,4) = for_hb(:,4) + f_i(:)
      for_hb(:,2) = for_hb(:,2) + (-1.0_wp + d1242over1212) * f_i(:) &
                                        - (  d1213over1212) * f_l(:)
      for_hb(:,1) = for_hb(:,1) + (-1.0_wp + d1213over1212) * f_l(:) &
                                        - (  d1242over1212) * f_i(:)
      for_hb(:,3) = for_hb(:,3) + f_l(:)

      !===== Dihedral angle among 5-3-1=2 =====!
      c5313(1) = v53(2) * v13(3) - v53(3) * v13(2)
      c5313(2) = v53(3) * v13(1) - v53(1) * v13(3)
      c5313(3) = v53(1) * v13(2) - v53(2) * v13(1)
      n(:) = -c1213(:)

      ! ENERGY
      dih = atan2(dot_product(v53,n)*a13 , dot_product(c5313,n))
      d = dih - dih1_eq(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ex = ex - dih_coef(i) * d**2

      ! FORCE
      c5313_abs2 = dot_product(c5313, c5313)

      cos_theta531 = d1353 / (a13 * norm2(v53))
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta531) > 0.99995_wp) then
          c5313_abs2 = dot_product(v53,v53) * d1313 * (1.0 - 0.99995_wp**2)
      endif

      pre = 2.0_wp * dih_coef(i) * d * a13
      f_i(:) = + pre / c5313_abs2 * c5313(:)
      f_l(:) = - pre / c1213_abs2 * n(:)
    
      for_hb(:,5) = for_hb(:,5) + f_i(:)
      for_hb(:,3) = for_hb(:,3) + (-1.0_wp + d1353over1313) * f_i(:) &
                                        - (  d1213over1313) * f_l(:)
      for_hb(:,1) = for_hb(:,1) + (-1.0_wp + d1213over1313) * f_l(:) &
                                        - (  d1353over1313) * f_i(:)
      for_hb(:,2) = for_hb(:,2) + f_l(:)

      !===== Dihedral angle among 1=2-4-6 =====!
      n(:) = -c4212(:)
      c4246(1) = v42(2) * v46(3) - v42(3) * v46(2)
      c4246(2) = v42(3) * v46(1) - v42(1) * v46(3)
      c4246(3) = v42(1) * v46(2) - v42(2) * v46(1)
      c4246_abs2 = dot_product(c4246, c4246)

      ! ENERGY
      dih = atan2(dot_product(v12,c4246)*a42 , dot_product(n,c4246))
      d = dih - dih2_eq(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ex = ex - dih_coef(i) * d**2

      ! FORCE
      cos_theta246 = d4246 / (a42 * norm2(v46))
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta246) > 0.99995_wp) then
          c4246_abs2 = dot_product(v46,v46) * d4242 * (1.0 - 0.99995_wp**2)
      endif

      pre = 2.0_wp * dih_coef(i) * d * a42
      f_i(:) = + pre / c4212_abs2 * n(:)
      f_l(:) = - pre / c4246_abs2 * c4246(:)
    
      for_hb(:,1) = for_hb(:,1) + f_i(:)
      for_hb(:,2) = for_hb(:,2) + (-1.0_wp + d1242over4242) * f_i(:) &
                                        - (      d4246over4242) * f_l(:)
      for_hb(:,4) = for_hb(:,4) + (-1.0_wp + d4246over4242) * f_l(:) &
                                        - (      d1242over4242) * f_i(:)
      for_hb(:,6) = for_hb(:,6) + f_l(:)

      !===== Total =====!

      ex = U0(i) * exp(ex)
      enehb = enehb + ex
      ! check this
      for_hb(:,:) = enehb * for_hb(:,:)

      force(1:3, i_atom,  id+1) = force(1:3, i_atom,  id+1) + for_hb(1:3, 1)
      force(1:3, j_atom,  id+1) = force(1:3, j_atom,  id+1) + for_hb(1:3, 2)
      force(1:3, i1_atom, id+1) = force(1:3, i1_atom, id+1) + for_hb(1:3, 3)
      force(1:3, j1_atom, id+1) = force(1:3, j1_atom, id+1) + for_hb(1:3, 4)
      force(1:3, i2_atom, id+1) = force(1:3, i2_atom, id+1) + for_hb(1:3, 5)
      force(1:3, j2_atom, id+1) = force(1:3, j2_atom, id+1) + for_hb(1:3, 6)

    end do
    !$omp end parallel

    call timer(TimerTISHB, TimerOff)

    return
  
  end subroutine compute_energy_tis_hb_pbc

  subroutine compute_energy_tis_harmonic_hb_pbc(enefunc, boundary, coord, &
                                           force, virial, enehb)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_boundary),target, intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: enehb

    ! local variables
    integer                  :: i, id 
    integer                  :: i_atom, j_atom, i1_atom, j1_atom, i2_atom, j2_atom
    integer                  :: istart, iend
    real(wp)                 :: v12(1:3), v13(1:3), v53(1:3), v42(1:3), v46(1:3)
    real(wp)                 :: d, dih 
    real(wp)                 :: a42, a13, a12
    real(wp)                 :: cos_theta124, cos_theta312, cos_theta531, cos_theta246
    real(wp)                 :: d1212, d1313, d4242
    real(wp)                 :: d1213, d1242, d4246, d1353
    real(wp)                 :: d1213over1212, d1213over1313
    real(wp)                 :: d1242over1212, d1242over4242
    real(wp)                 :: d4246over4242, d1353over1313
    real(wp)                 :: n(3)
    real(wp)                 :: c4212(3), c1213(3), c5313(3), c4246(3)
    real(wp)                 :: c4212_abs2, c1213_abs2, c5313_abs2, c4246_abs2
    real(wp)                 :: pre, f_i(3), f_k(3), f_l(3), ex, ene, for_hb(1:3, 1:6)
    real(wp)                 :: bsize(3), inv_bsize(3)
    real(wp)                 :: dist_c, ang_c, dih_c   

    integer,         pointer :: list(:,:)
    real(wp),        pointer :: U0(:)
    real(wp),        pointer :: dist_eq(:), dist_coef(:)
    real(wp),        pointer :: ang1_eq(:), ang2_eq(:), ang_coef(:)
    real(wp),        pointer :: dih_eq(:), dih1_eq(:), dih2_eq(:), dih_coef(:)

    integer                  :: omp_get_thread_num

    call timer(TimerTISHB, TimerOn)

    istart       = enefunc%istart_tis_hb
    iend         = enefunc%iend_tis_hb

    ! use pointers
    !
    list         => enefunc%tis_hb_list
    U0           => enefunc%tis_hb_dist_U0
    dist_eq      => enefunc%tis_hb_dist_eq
    dist_coef    => enefunc%tis_hb_dist_coef
    ang1_eq      => enefunc%tis_hb_angle_ang1
    ang2_eq      => enefunc%tis_hb_angle_ang2
    ang_coef     => enefunc%tis_hb_angle_coef
    dih_eq       => enefunc%tis_hb_dihedral_dih
    dih1_eq      => enefunc%tis_hb_dihedral_dih1
    dih2_eq      => enefunc%tis_hb_dihedral_dih2
    dih_coef     => enefunc%tis_hb_dihedral_coef

    bsize(1) = boundary%box_size_x
    bsize(2) = boundary%box_size_y
    bsize(3) = boundary%box_size_z
    inv_bsize = 1.0_wp/bsize(1:3)

    ! calculation of local stacking energy and gradient
    !
    !$omp parallel default(none)                       &
    !$omp private(id, i, i_atom, j_atom, i1_atom,      &
    !$omp         j1_atom, i2_atom, j2_atom,           &
    !$omp         d,dih,cos_theta124,cos_theta312,     &
    !$omp         cos_theta531,cos_theta246,           &
    !$omp         v12,v13,v53,v42,v46,a12,a13,a42,     &
    !$omp         d1212,d1313, for_hb, ex, ene,        &
    !$omp         pre, f_i, f_k, f_l, n,               &
    !$omp         d4242,d1213,d1242,d4246,d1353,       &
    !$omp         d1213over1212,d1213over1313,         &
    !$omp         d1242over1212,d1242over4242,         &
    !$omp         d4246over4242,d1353over1313,c4212,   &
    !$omp         c1213,c4212_abs2,c1213_abs2,         &
    !$omp         c5313, c5313_abs2, c4246, c4246_abs2,&
    !$omp         dist_c, ang_c, dih_c)                &
    !$omp shared(istart, iend, U0, dist_eq,            &
    !$omp        coord, list, ang1_eq, ang2_eq,        &
    !$omp        dih_eq, dih1_eq, dih2_eq,             &
    !$omp        force, nthread,                       &
    !$omp        bsize, inv_bsize,                     &
    !$omp        dist_coef, ang_coef, dih_coef)        &
    !$omp reduction(+:enehb) reduction(+:virial)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread
      for_hb(:,:) = 0.0_wp

      i_atom  = list(1, i) ! 1
      j_atom  = list(2, i) ! 2
      i1_atom = list(3, i) ! 3
      j1_atom = list(4, i) ! 4
      i2_atom = list(5, i) ! 5
      j2_atom = list(6, i) ! 6

      dist_c = -U0(i) * dist_coef(i)
      ang_c  = -U0(i) * ang_coef(i)
      dih_c  = -U0(i) * dih_coef(i)

      !===== Distance =====!
      v12(1:3) = coord(1:3, i_atom)  - coord(1:3, j_atom)
      v12(1:3) = v12(1:3) - bsize(1:3)*anint(v12(1:3)*inv_bsize(1:3))
      d1212 = dot_product(v12,v12)
      a12 = norm2(v12) 

      d = a12 - dist_eq(i)

      ex = dist_c * d**2
      f_i(:) = - (2.0_wp * dist_c * d / a12) * v12(:)
      for_hb(:,1) = + f_i(:)
      for_hb(:,2) = - f_i(:)

      v13(1:3) = coord(1:3, i_atom)  - coord(1:3, i1_atom)
      v13(1:3) = v13(1:3) - bsize(1:3)*anint(v13(1:3)*inv_bsize(1:3))

      v53(1:3) = coord(1:3, i2_atom) - coord(1:3, i1_atom)
      v53(1:3) = v53(1:3) - bsize(1:3)*anint(v53(1:3)*inv_bsize(1:3))

      v42(1:3) = coord(1:3, j1_atom) - coord(1:3, j_atom)
      v42(1:3) = v42(1:3) - bsize(1:3)*anint(v42(1:3)*inv_bsize(1:3))

      v46(1:3) = coord(1:3, j1_atom) - coord(1:3, j2_atom)
      v46(1:3) = v46(1:3) - bsize(1:3)*anint(v46(1:3)*inv_bsize(1:3))

      d1313 = dot_product(v13,v13)
      d4242 = dot_product(v42,v42)
      a13 = sqrt(d1313)
      a42 = sqrt(d4242)
      d1213 = dot_product(v13,v12)
      d1242 = dot_product(v12,v42)
      d4246 = dot_product(v42,v46)
      d1353 = dot_product(v13,v53)
      d1213over1212 = d1213 / d1212
      d1213over1313 = d1213 / d1313
      d1242over1212 = d1242 / d1212
      d4246over4242 = d4246 / d4242
      d1242over4242 = d1242 / d4242
      d1353over1313 = d1353 / d1313

      !===== Angle of 3-1=2  =====!
      ! ENERGY
      cos_theta312 = d1213 / (a13 * a12)
      d = acos(cos_theta312) - ang1_eq(i)
      ex = ex + ang_c * d**2

      ! FORCE
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta312) > 0.99995_wp) then
          d1213 = sign(a12 * a13 * 0.99995_wp, cos_theta312)
      end if

      pre = - 2.0_wp * ang_c * d / sqrt(d1313*d1212 - d1213**2)
      f_i(:) = pre * (v12(:) - (d1213over1313 * v13(:)))
      f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))
      for_hb(:,3) = for_hb(:,3) + f_i(:)
      for_hb(:,2) = for_hb(:,2) + f_k(:)
      for_hb(:,1) = for_hb(:,1) - f_i(:) - f_k(:)

      !===== Angle of 1=2-4  =====!
      ! ENERGY
      cos_theta124 = d1242 / (a12 * a42)
      d = acos(cos_theta124) - ang2_eq(i)
      ex = ex + ang_c * d**2

      ! FORCE
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta124) > 0.99995_wp) then
          d1242 = sign(a12 * a42 * 0.99995_wp, cos_theta124)
      endif

      pre = - 2.0_wp * ang_c * d / sqrt(d1212*d4242 - d1242**2)
      f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1242over4242 * v42(:)))
      for_hb(:,1) = for_hb(:,1) + f_i(:)
      for_hb(:,4) = for_hb(:,4) + f_k(:)
      for_hb(:,2) = for_hb(:,2) - f_i(:) - f_k(:)

      !===== Dihedral angle among 4-2=1=3 =====!
      c4212(1) = v42(2)*v12(3) - v42(3)*v12(2)
      c4212(2) = v42(3)*v12(1) - v42(1)*v12(3)
      c4212(3) = v42(1)*v12(2) - v42(2)*v12(1)
      c1213(1) = v12(2)*v13(3) - v12(3)*v13(2)
      c1213(2) = v12(3)*v13(1) - v12(1)*v13(3)
      c1213(3) = v12(1)*v13(2) - v12(2)*v13(1)
      c4212_abs2 = dot_product(c4212,c4212)
      c1213_abs2 = dot_product(c1213,c1213)

      ! ENERGY
      dih = atan2(dot_product(v42,c1213)*a12, dot_product(c4212,c1213))
      d = dih - dih_eq(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ex = ex + dih_c * d**2

      ! FORCE
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta124) > 0.99995_wp) then
          ! |c4212|^2 = |v42|^2 * |v12|^2 * sin(theta)^2
          c4212_abs2 = d4242 * d1212 * (1.0 - 0.99995_wp**2)
      endif
      if (abs(cos_theta312) > 0.99995_wp) then
          c1213_abs2 = d1212 * d1313 * (1.0 - 0.99995_wp**2)
      endif

      pre = - 2.0_wp * dih_c * d * a12
      f_i(:) = + pre / c4212_abs2 * c4212(:)
      f_l(:) = - pre / c1213_abs2 * c1213(:)
    
      for_hb(:,4) = for_hb(:,4) + f_i(:)
      for_hb(:,2) = for_hb(:,2) + (-1.0_wp + d1242over1212) * f_i(:) &
                                        - (  d1213over1212) * f_l(:)
      for_hb(:,1) = for_hb(:,1) + (-1.0_wp + d1213over1212) * f_l(:) &
                                        - (  d1242over1212) * f_i(:)
      for_hb(:,3) = for_hb(:,3) + f_l(:)

      !===== Dihedral angle among 5-3-1=2 =====!
      c5313(1) = v53(2) * v13(3) - v53(3) * v13(2)
      c5313(2) = v53(3) * v13(1) - v53(1) * v13(3)
      c5313(3) = v53(1) * v13(2) - v53(2) * v13(1)
      n(:) = -c1213(:)

      ! ENERGY
      dih = atan2(dot_product(v53,n)*a13 , dot_product(c5313,n))
      d = dih - dih1_eq(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ex = ex + dih_c * d**2

      ! FORCE
      c5313_abs2 = dot_product(c5313, c5313)

      cos_theta531 = d1353 / (a13 * norm2(v53))
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta531) > 0.99995_wp) then
          c5313_abs2 = dot_product(v53,v53) * d1313 * (1.0 - 0.99995_wp**2)
      endif

      pre = - 2.0_wp * dih_c * d * a13
      f_i(:) = + pre / c5313_abs2 * c5313(:)
      f_l(:) = - pre / c1213_abs2 * n(:)
    
      for_hb(:,5) = for_hb(:,5) + f_i(:)
      for_hb(:,3) = for_hb(:,3) + (-1.0_wp + d1353over1313) * f_i(:) &
                                        - (  d1213over1313) * f_l(:)
      for_hb(:,1) = for_hb(:,1) + (-1.0_wp + d1213over1313) * f_l(:) &
                                        - (  d1353over1313) * f_i(:)
      for_hb(:,2) = for_hb(:,2) + f_l(:)

      !===== Dihedral angle among 1=2-4-6 =====!
      n(:) = -c4212(:)
      c4246(1) = v42(2) * v46(3) - v42(3) * v46(2)
      c4246(2) = v42(3) * v46(1) - v42(1) * v46(3)
      c4246(3) = v42(1) * v46(2) - v42(2) * v46(1)
      c4246_abs2 = dot_product(c4246, c4246)

      ! ENERGY
      dih = atan2(dot_product(v12,c4246)*a42 , dot_product(n,c4246))
      d = dih - dih2_eq(i)
      if (d > PI) then
         d = d - 2*PI
      else if (d < -PI) then
         d = d + 2*PI
      endif
      ex = ex + dih_c * d**2

      ! FORCE
      cos_theta246 = d4246 / (a42 * norm2(v46))
      ! use maximum values at the limit angle 0.99995_wp
      if (abs(cos_theta246) > 0.99995_wp) then
          c4246_abs2 = dot_product(v46,v46) * d4242 * (1.0_wp - 0.99995_wp**2)
      endif

      pre = - 2.0_wp * dih_c * d * a42
      f_i(:) = + pre / c4212_abs2 * n(:)
      f_l(:) = - pre / c4246_abs2 * c4246(:)
    
      for_hb(:,1) = for_hb(:,1) + f_i(:)
      for_hb(:,2) = for_hb(:,2) + (-1.0_wp + d1242over4242) * f_i(:) &
                                        - (      d4246over4242) * f_l(:)
      for_hb(:,4) = for_hb(:,4) + (-1.0_wp + d4246over4242) * f_l(:) &
                                        - (      d1242over4242) * f_i(:)
      for_hb(:,6) = for_hb(:,6) + f_l(:)

      !===== Total =====!

      enehb = enehb + ex

      force(1:3, i_atom,  id+1) = force(1:3, i_atom,  id+1) + for_hb(1:3, 1)
      force(1:3, j_atom,  id+1) = force(1:3, j_atom,  id+1) + for_hb(1:3, 2)
      force(1:3, i1_atom, id+1) = force(1:3, i1_atom, id+1) + for_hb(1:3, 3)
      force(1:3, j1_atom, id+1) = force(1:3, j1_atom, id+1) + for_hb(1:3, 4)
      force(1:3, i2_atom, id+1) = force(1:3, i2_atom, id+1) + for_hb(1:3, 5)
      force(1:3, j2_atom, id+1) = force(1:3, j2_atom, id+1) + for_hb(1:3, 6)

    end do
    !$omp end parallel

    call timer(TimerTISHB, TimerOff)

    return
  
  end subroutine compute_energy_tis_harmonic_hb_pbc

end module at_energy_tis_mod
