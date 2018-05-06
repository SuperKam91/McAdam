!> This allows for a simple C interface... 


module interfaces_module
    use utils_module, only: dp

    implicit none
    interface run_polychord
        module procedure run_polychord_full
    end interface run_polychord

contains


    subroutine run_polychord_full(loglikelihood, settings_in)
        use settings_module,          only: program_settings,initialise_settings
        use random_module,            only: initialise_random
        use nested_sampling_module,   only: NestedSampling
#ifdef MPI
        use mpi_module,               only: initialise_mpi, finalise_mpi
        use mpi,                      only: MPI_COMM_WORLD
#endif
        implicit none

        interface
            function loglikelihood(hyp,theta,phi)
                import :: dp
                real(dp), intent(in), dimension(:) :: hyp
                real(dp), intent(out), dimension(:) :: theta
                real(dp), intent(out), dimension(:) :: phi
                real(dp) :: loglikelihood
            end function loglikelihood
        end interface

        type(program_settings),intent(in)    :: settings_in
        type(program_settings)               :: settings 

        real(dp), dimension(4) :: output_info

#ifdef MPI
        call initialise_mpi
#endif
        call initialise_random()
        settings = settings_in
        call initialise_settings(settings)   
#ifdef MPI
        output_info = NestedSampling(loglikelihood,settings,MPI_COMM_WORLD) 
        call finalise_mpi
#else
        output_info = NestedSampling(loglikelihood,settings,0) 
#endif

    end subroutine run_polychord_full
end module interfaces_module
