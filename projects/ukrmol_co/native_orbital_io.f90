! Read saved GBTO AO/MO coefficients and stored overlaps; never evaluate integrals.
program native_orbital_io
    use ieee_arithmetic, only: ieee_is_finite
    use precisn_gbl, only: wp, cfp
    use const_gbl, only: line_len, overlap_ints, set_verbosity_level
    use mpi_gbl, only: mpi_mod_start, mpi_mod_finalize
    use atomic_basis_gbl, only: atomic_orbital_basis_obj, CGTO_shell_data_obj
    use molecular_basis_gbl, only: molecular_orbital_basis_obj
    use parallel_arrays_gbl, only: p2d_array_obj
    use integral_storage_gbl, only: integral_options_obj
    implicit none

    character(len=line_len) :: input, output
    character(len=32) :: magic = 'CO_ORBITAL_EXPORT_V1'
    type(atomic_orbital_basis_obj), target :: ao
    type(molecular_orbital_basis_obj) :: mo
    type(p2d_array_obj), target :: one, two
    type(integral_options_obj) :: options
    type(CGTO_shell_data_obj), allocatable :: shells(:)
    real(cfp), allocatable :: coefficients(:,:), overlaps(:,:), geometry(:,:)
    integer, allocatable :: counts(:), mapping(:,:), flags(:), ao_map(:,:), descriptors(:,:)
    logical, allocatable :: continuum(:)
    integer :: mode, nao, nmo, nnuc, nsym, nshell, i, j, g, k, u, err, column_s, pair(2,1), idx(1)
    real(cfp) :: radius

    call get_command_argument(1, input)
    call get_command_argument(2, output)
    if (len_trim(output) == 0) error stop 'input and fresh output required'
    if (wp /= 8 .or. cfp /= 8) error stop 'matching double GBTO modules required'
    call mpi_mod_start
    call set_verbosity_level(1)
    mode = 1
    if (trim(input) == '--self-test') then
        mode = 0
        nao=3; nmo=2; nnuc=0; nsym=1; nshell=1; radius=18
        allocate(coefficients(nao,nmo), overlaps(nmo,nmo), geometry(4,0))
        allocate(counts(1), mapping(3,nmo), flags(nmo), ao_map(2,nao), descriptors(6,1))
        coefficients=0; coefficients(1,1)=1; coefficients(2,2)=1
        overlaps=0; overlaps(1,1)=1; overlaps(2,2)=1
        counts=[2]; mapping(:,1)=[1,1,1]; mapping(:,2)=[2,1,2]; flags=[0,1]
        ao_map(1,:)=1; ao_map(2,:)=[1,2,3]; descriptors(:,1)=[1,1,0,1,3,1]
    else
        err = mo%read_basis_and_integrals(input, ao, one, two, options)
        if (err /= 0) error stop 'saved basis/integral read failed'
        if (ao%contains_btos()) error stop 'this export version supports CGTO-only AO bases'
        call ao%get_all_CGTO_shells(shells, nshell)
        call mo%get_orbital_coefficient_matrix(coefficients)
        nao=size(coefficients,1); nmo=size(coefficients,2)
        nnuc=mo%symmetry_data%no_nuc; nsym=mo%no_irr; radius=options%a
        if (nao > 10000 .or. nmo > 2000) error stop 'export dimension cap'
        allocate(overlaps(nmo,nmo), geometry(4,nnuc), counts(nsym), mapping(3,nmo), flags(nmo))
        ao_map=ao%indices_to_shells
        descriptors=ao%shell_descriptor
        if (size(descriptors,2) /= nshell .or. any(descriptors(1,:) /= 1)) error stop 'CGTO inventory'
        k=0
        do g=1,nsym
            counts(g)=mo%get_number_of_orbitals(g)
            call mo%get_continuum_flags(g, continuum)
            do j=1,counts(g)
                k=k+1
                mapping(:,k)=[mo%get_absolute_index(j,g),g,j]
                flags(k)=merge(1,0,continuum(j))
            end do
        end do
        if (k /= nmo) error stop 'MO inventory'
        column_s=one%find_column_matching_name(overlap_ints)
        if (column_s <= 0) error stop 'missing stored overlap column'
        do j=1,nmo
            do i=1,nmo
                pair(:,1)=[mapping(1,i),mapping(1,j)]
                idx=mo%integral_index(overlap_ints,pair,options%two_p_continuum)
                if (idx(1) <= 0) error stop 'invalid stored overlap index'
                overlaps(i,j)=one%a(idx(1),column_s)
            end do
        end do
        do i=1,nnuc
            geometry(:3,i)=mo%symmetry_data%nucleus(i)%center
            geometry(4,i)=mo%symmetry_data%nucleus(i)%charge
            if (mo%symmetry_data%nucleus(i)%ecp%n_core /= 0) error stop 'unexpected ECP'
        end do
    end if
    if (.not. all(ieee_is_finite(coefficients)) .or. .not. all(ieee_is_finite(overlaps))) &
        error stop 'nonfinite saved coefficient/overlap'
    open(newunit=u,file=trim(output),status='new',form='unformatted',action='write')
    write(u) magic
    write(u) mode,nao,nmo,nnuc,nsym,nshell
    write(u) radius
    write(u) counts
    write(u) mapping
    write(u) flags
    write(u) coefficients
    write(u) overlaps
    write(u) ao_map
    write(u) descriptors
    write(u) geometry
    if (mode == 0) then
        write(u) 1,1,0
        write(u) 0.0_cfp,0.0_cfp,0.0_cfp,1.0_cfp
        write(u) 1.0_cfp
        write(u) 1.0_cfp
        write(u) 1.0_cfp
    else
        do i=1,nshell
            write(u) shells(i)%l,shells(i)%number_of_primitives,merge(1,0,shells(i)%is_continuum())
            write(u) shells(i)%center,shells(i)%norm
            write(u) shells(i)%exponents
            write(u) shells(i)%contractions
            write(u) shells(i)%norms
        end do
    end if
    close(u)
    if (mode == 1) then
        err=mo%final()
        if (err /= 0) error stop 'MO finalization'
        err=one%final()
        if (err /= 0) error stop 'one-electron storage finalization'
        err=two%final()
        if (err /= 0) error stop 'two-electron storage finalization'
    end if
    write(*,'(a)') 'CO_SAVED_ORBITAL_EXPORT_COMPLETE'
    call mpi_mod_finalize
end program native_orbital_io
