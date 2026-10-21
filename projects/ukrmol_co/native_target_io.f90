! Read-only diagnostic driver over installed UKRmol+/GBTOlib entry points.
! No integral evaluation, Hamiltonian construction or diagonalization is called.
program native_target_io
    use precisn_gbl, only: wp
    use const_gbl, only: line_len, one_elham, overlap_ints, two_el_ints, set_verbosity_level
    use mpi_gbl, only: mpi_mod_start, mpi_mod_finalize
    use atomic_basis_gbl, only: atomic_orbital_basis_obj
    use molecular_basis_gbl, only: molecular_orbital_basis_obj
    use parallel_arrays_gbl, only: p2d_array_obj
    use integral_storage_gbl, only: integral_options_obj
    use scatci_routines, only: cirmat
    implicit none

    character(len=line_len) :: mode, input, csfdir, output, path
    character(len=3) :: unit_text
    integer :: i, j, k, l, g, ordinal, n, err, u, v, column_h, column_s, idx(1)
    integer :: dims(8) = [15548,24366,15196,24540,15196,24540,15044,24564]
    integer :: units(8) = [710,730,711,731,712,732,713,733]
    integer :: roots(8) = [13,13,8,8,8,8,8,8], target(4) = [7,3,3,0]
    integer :: counts(8), mapping(13), two(2,1), four(4,1), nnuc
    integer, allocatable :: phases(:), sets(:), states(:)
    real(wp), allocatable :: vectors(:), energies(:), geometry(:,:)
    real(wp) :: overlap(13,13), h(13,13), eri(13,13,13,13), enuc, distance
    logical, allocatable :: continuum(:)
    type(atomic_orbital_basis_obj), target :: ao
    type(molecular_orbital_basis_obj) :: mo
    type(p2d_array_obj), target :: one_electron, two_electron
    type(integral_options_obj) :: options

    call get_command_argument(1, mode)
    call get_command_argument(2, input)
    call get_command_argument(3, csfdir)
    call get_command_argument(4, output)
    if (len_trim(output) == 0) error stop 'four arguments required'
    call mpi_mod_start
    call set_verbosity_level(1)
    open(newunit=u, file=trim(output), status='new', form='unformatted', action='write')

    select case (trim(mode))
    case ('target')
        open(unit=26, file=trim(input), status='old', form='unformatted', action='read')
        do g=1,8
            n = dims(g)
            allocate(phases(n), sets(roots(g)), states(roots(g)))
            allocate(vectors(n*roots(g)), energies(roots(g)))
            write(unit_text,'(i3)') units(g)
            path = trim(csfdir)//'/fort.'//unit_text
            open(newunit=v, file=trim(path), status='old', form='unformatted', action='read')
            read(v)
            read(v) phases
            close(v)
            sets = g
            states = [(i, i=1,roots(g))]
            call cirmat(1, 26, sets, vectors, energies, states, [n], 6, [roots(g)], phases, n, 2)
            write(u) g, n, roots(g)
            write(u) energies
            write(u) vectors
            deallocate(phases, sets, states, vectors, energies)
        end do
        close(26)
        close(u)
        write(*,'(a)') 'CO_NATIVE_TARGET_READ_COMPLETE'
    case ('integrals')
        err = mo%read_basis_and_integrals(input, ao, one_electron, two_electron, options)
        if (err /= 0) error stop 'saved integral read failed'
        counts = 0
        mapping = 0
        n = 0
        if (mo%no_irr > 8) error stop 'unexpected symmetry inventory'
        do g=1,mo%no_irr
            counts(g) = mo%get_number_of_orbitals(g)
            call mo%get_continuum_flags(g, continuum)
            if (g <= 4) then
                if (counts(g) < target(g)) error stop 'missing target orbitals'
                if (any(continuum(:target(g)))) error stop 'target flagged as continuum'
                if (.not. all(continuum(target(g)+1:))) error stop 'unexpected extra target orbital'
                do ordinal=1,target(g)
                    n = n+1
                    mapping(n) = mo%get_absolute_index(ordinal,g)
                end do
            else
                if (counts(g) /= 0) error stop 'unexpected occupied symmetry'
            end if
        end do
        if (n /= 13) error stop 'target orbital count mismatch'
        column_h = one_electron%find_column_matching_name(one_elham)
        column_s = one_electron%find_column_matching_name(overlap_ints)
        if (min(column_h,column_s) <= 0) error stop 'missing stored h/overlap column'
        do j=1,13
            do i=1,13
                two(:,1) = [mapping(i),mapping(j)]
                idx = mo%integral_index(one_elham,two,options%two_p_continuum)
                if (idx(1) <= 0) error stop 'invalid stored one-electron index'
                h(i,j) = one_electron%a(idx(1),column_h)
                overlap(i,j) = one_electron%a(idx(1),column_s)
                do l=1,13
                    do k=1,13
                        four(:,1) = [mapping(i),mapping(j),mapping(k),mapping(l)]
                        idx = mo%integral_index(two_el_ints,four,options%two_p_continuum)
                        eri(i,j,k,l) = 0.0_wp
                        if (idx(1) > 0) eri(i,j,k,l) = two_electron%a(idx(1),1)
                    end do
                end do
            end do
        end do
        nnuc = mo%symmetry_data%no_nuc
        allocate(geometry(4,nnuc))
        do i=1,nnuc
            geometry(:3,i) = mo%symmetry_data%nucleus(i)%center
            geometry(4,i) = mo%symmetry_data%nucleus(i)%charge
            if (mo%symmetry_data%nucleus(i)%ecp%n_core /= 0) error stop 'unexpected ECP'
        end do
        enuc = 0.0_wp
        do j=2,nnuc
            do i=1,j-1
                if (abs(geometry(4,i)*geometry(4,j)) < tiny(enuc)) cycle
                distance = sqrt(sum((geometry(:3,i)-geometry(:3,j))**2))
                if (distance <= 0.0_wp) error stop 'coincident charged nuclei'
                enuc = enuc + geometry(4,i)*geometry(4,j)/distance
            end do
        end do
        write(u) 13, mo%number_of_functions, nnuc
        write(u) counts, mapping
        write(u) geometry
        write(u) overlap
        write(u) h
        write(u) eri
        write(u) enuc
        close(u)
        err = mo%final()
        if (err /= 0) error stop 'MO finalization failed'
        err = one_electron%final()
        if (err /= 0) error stop 'one-electron finalization failed'
        err = two_electron%final()
        if (err /= 0) error stop 'two-electron finalization failed'
        write(*,'(a)') 'CO_NATIVE_INTEGRAL_READ_COMPLETE'
    case default
        error stop 'unknown read-only mode'
    end select
    call mpi_mod_finalize
end program native_target_io
