module kind_def

  integer, parameter :: dp = kind(1.0D0)
  integer, parameter :: rk = SELECTED_REAL_KIND(12,200)
  !integer :: num_params  = 2

end module kind_def
