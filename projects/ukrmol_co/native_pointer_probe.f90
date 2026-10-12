! Reduced probe of the retained TARGET-dummy pointers used by MPI-SCATCI.
! Compile with the campaign's pinned compiler; this contains no scattering physics.
module pointer_probe
    implicit none
    type manager
        integer :: marker = 17
    end type
    type csf
        type(manager), pointer :: owner
        type(csf), pointer :: parent
        real(8), allocatable :: coefficients(:)
    end type
    type target_ci
        real(8), allocatable :: coefficients(:)
    end type
    type holder
        class(csf), pointer :: configurations(:)
        class(manager), pointer :: factory
        class(target_ci), pointer :: targets(:)
    end type
contains
    subroutine create(factory, configurations)
        class(manager), target, intent(inout) :: factory
        type(csf), allocatable, target, intent(out) :: configurations(:)
        integer :: i
        allocate(configurations(4))
        do i = 1, 4
            configurations(i)%owner => factory
            configurations(i)%parent => configurations(i)
            allocate(configurations(i)%coefficients(2))
            configurations(i)%coefficients = i
        end do
    end subroutine
    subroutine retain(h, configurations, factory, targets)
        type(holder), intent(out) :: h
        class(csf), target, intent(in) :: configurations(:)
        class(manager), target, intent(in) :: factory
        class(target_ci), target, intent(in) :: targets(:)
        h%configurations => configurations
        h%factory => factory
        h%targets => targets
    end subroutine
    subroutine original()
        type(manager) :: factory
        type(csf), allocatable :: configurations(:)
        type(target_ci), allocatable :: targets(:)
        type(holder) :: h
        call create(factory, configurations)
        allocate(targets(2))
        call retain(h, configurations, factory, targets)
        write(*,*) 'ORIGINAL configuration address match:', &
            loc(h%configurations) == loc(configurations)
        write(*,*) 'ORIGINAL factory address match:', loc(h%factory) == loc(factory)
        write(*,*) 'ORIGINAL target address match:', loc(h%targets) == loc(targets)
        write(*,*) 'ORIGINAL parent address match:', &
            loc(configurations(4)%parent) == loc(configurations(4))
        write(*,*) 'ORIGINAL owner address match:', &
            loc(configurations(4)%owner) == loc(factory)
    end subroutine
    subroutine explicit_targets()
        type(manager), target :: factory
        type(csf), allocatable, target :: configurations(:)
        type(target_ci), allocatable, target :: targets(:)
        type(holder) :: h
        call create(factory, configurations)
        allocate(targets(2))
        call retain(h, configurations, factory, targets)
        write(*,*) 'TARGET configuration address match:', &
            loc(h%configurations) == loc(configurations)
        write(*,*) 'TARGET factory address match:', loc(h%factory) == loc(factory)
        write(*,*) 'TARGET target address match:', loc(h%targets) == loc(targets)
        write(*,*) 'TARGET parent address match:', &
            loc(configurations(4)%parent) == loc(configurations(4))
        write(*,*) 'TARGET owner address match:', &
            loc(configurations(4)%owner) == loc(factory)
        if (loc(h%configurations) /= loc(configurations)) error stop 1
        if (loc(h%factory) /= loc(factory)) error stop 2
        if (loc(h%targets) /= loc(targets)) error stop 3
        if (loc(configurations(4)%parent) /= loc(configurations(4))) error stop 4
        if (loc(configurations(4)%owner) /= loc(factory)) error stop 5
    end subroutine
end module
program probe
    use pointer_probe
    call original()
    call explicit_targets()
end program
