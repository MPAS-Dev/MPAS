program mpas

   use mpas_subdriver

   implicit none

   call mpas_init()

   call mpas_run() 

   call mpas_finalize()

   stop

end program mpas
