! Test the audit predicate independently of SCATCI and its data.
module association_probe
    implicit none
    type configuration
        type(determinant), allocatable :: orbitals(:)
        integer :: id
    end type
    type determinant
        type(configuration), pointer :: parent
        integer :: id
    end type
contains
    subroutine audit_polymorphic(values)
        class(configuration), target, intent(in) :: values(:)
        integer :: i
        do i = 1, size(values)
            write(*,*) 'CLASS element', i, 'associated=', &
                associated(values(i)%orbitals(1)%parent, values(i)), &
                'actual=', loc(values(i)%orbitals(1)%parent), 'expected=', loc(values(i))
            if (loc(values(i)%orbitals(1)%parent) /= loc(values(i))) error stop 1
        end do
    end subroutine
    subroutine audit_nonpolymorphic(values)
        type(configuration), target, intent(in) :: values(:)
        integer :: i
        do i = 1, size(values)
            write(*,*) 'TYPE element', i, 'associated=', &
                associated(values(i)%orbitals(1)%parent, values(i)), &
                'actual=', loc(values(i)%orbitals(1)%parent), 'expected=', loc(values(i))
            if (.not. associated(values(i)%orbitals(1)%parent, values(i))) error stop 2
        end do
    end subroutine
    subroutine verify_address_guard(values)
        type(configuration), target, intent(inout) :: values(:)
        integer :: i
        do i = 1, size(values)
            if (loc(values(i)%orbitals(1)%parent) /= loc(values(i))) error stop 3
        end do
        values(2)%orbitals(1)%parent => values(1)
        if (loc(values(2)%orbitals(1)%parent) == loc(values(2))) error stop 4
        write(*,*) 'ADDRESS guard detects deliberately redirected parent: PASS'
        values(2)%orbitals(1)%parent => values(2)
        nullify(values(2)%orbitals(1)%parent)
        if (loc(values(2)%orbitals(1)%parent) == loc(values(2))) error stop 5
        write(*,*) 'ADDRESS guard detects null parent: PASS'
        values(2)%orbitals(1)%parent => values(2)
    end subroutine
end module

program probe
    use association_probe
    implicit none
    type(configuration), allocatable, target :: values(:)
    integer :: i
    allocate(values(4))
    do i = 1, size(values)
        values(i)%id = i
        allocate(values(i)%orbitals(1))
        values(i)%orbitals(1)%parent => values(i)
        values(i)%orbitals(1)%id = i
    end do
    call audit_nonpolymorphic(values)
    call audit_polymorphic(values)
    call verify_address_guard(values)
end program
