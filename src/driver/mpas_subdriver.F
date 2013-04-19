module mpas_subdriver

   use mpas_framework
   use mpas_core
   use mpas_kind_types

   type (dm_info), pointer :: dminfo
   type (domain_type), pointer :: domain
   type (io_output_object), save :: output_obj
   integer :: output_frame


   contains


   subroutine mpas_init()
 
      implicit none

      real (kind=RKIND) :: dt
      character(len=StrKIND) :: timeStamp

      !
      ! Initialize infrastructure
      !
      call mpas_framework_init(dminfo, domain)

      call mpas_timer_start("total time")
      call mpas_timer_start("initialize")

      call mpas_input_state_for_domain(domain)


      !
      ! Initialize core
      !
      call mpas_core_init(domain, timeStamp)

      call mpas_timer_stop("initialize")


      !
      ! Set up output streams to be written to by the MPAS core
      !
      output_frame = 1

      if(config_frames_per_outfile > 0) then
         call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp))
      else
         call mpas_output_state_init(output_obj, domain, "OUTPUT")         
      end if


   end subroutine mpas_init


   subroutine mpas_run()

      implicit none

      call mpas_core_run(domain, output_obj, output_frame)

   end subroutine mpas_run


   subroutine mpas_finalize()
   
      implicit none

      !
      ! Finalize output streams
      !
      call mpas_output_state_finalize(output_obj, domain % dminfo)


      !
      ! Finalize core
      !
      call mpas_core_finalize(domain)

      call mpas_timer_stop("total time")
      call mpas_timer_write()


      !
      ! Finalize infrastructure
      !
      call mpas_framework_finalize(dminfo, domain)

   end subroutine mpas_finalize

end module mpas_subdriver
