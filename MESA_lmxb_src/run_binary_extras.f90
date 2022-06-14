! ***********************************************************************
!
!    Copyright (C) 2012  Bill Paxton
!
!    this file is part of mesa.
!
!    mesa is free software; you can redistribute it and/or modify
!    it under the terms of the gnu general library public license as published
!    by the free software foundation;either version 2 of the license, or
!    (at your option) any later version.
!
!    mesa is distributed in the hope that it will be useful, 
!    but without any warranty; without even the implied warranty of
!    merchantability or fitness for a particular purpose.  see the
!    gnu library general public license for more details.
!
!    you should have received a copy of the gnu library general public license
!    along with this software; if not, write to the free software
!    foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! *********************************************************************** 
        module run_binary_extras 

        use star_lib
        use star_def
        use const_def
        use const_def
        use chem_def
        use num_lib
        use binary_def
        use math_lib
        
        implicit none
        
        contains
        
        subroutine extras_binary_controls(binary_id, ierr)
            integer :: binary_id
            integer, intent(out) :: ierr
            type (binary_info), pointer :: b
            ierr = 0
            call binary_ptr(binary_id, b, ierr)
            if (ierr /= 0) then
                write(*,*) 'failed in binary_ptr'
                return
            end if
            ! write(*,*) 'hello from extra_binary_controls'
            b% other_mdot_edd => mdot_edd_routine
            b% other_jdot_mb => jdot_mb_routine

            b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
            b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

            b% extras_binary_startup => extras_binary_startup
            b% extras_binary_check_model => extras_binary_check_model
            b% extras_binary_finish_step => extras_binary_finish_step
            b% extras_binary_after_evolve => extras_binary_after_evolve
            b% warn_binary_extra =.false.
        end subroutine extras_binary_controls

        subroutine mdot_edd_routine(binary_id, mdot_edd, mdot_edd_eta, ierr)
            use const_def, only: dp
            integer, intent(in) :: binary_id
            real(dp), intent(out) :: mdot_edd
            real(dp), intent(out) :: mdot_edd_eta
            integer, intent(out) :: ierr
            type (binary_info), pointer :: b
            ierr = 0
            call binary_ptr(binary_id, b, ierr)
            if (ierr /= 0) then
                write(*,*) 'failed in binary_ptr'
                return
            end if

            ! changing mdot_edd from default MESA from BH to NS

            ! mdot_edd = 4*pi*clight*b% s1% x_ctrl(5)/(0.2*(1+b% s1% surface_h1))
            ! hard coding in the radius of 11.5km into the equation results in the next line.
            mdot_edd_eta = 0.1
            mdot_edd = 2.1666d18 / ((1.d0 + b% s1% surface_h1)) 
            ! write (*,*) "Modified Mdot_edd = ", mdot_edd
            ! write (*,*) "mdot_system_transfer = ", b% mdot_system_transfer(b% a_i)
            ! write (*,*) " "

        end subroutine mdot_edd_routine

        subroutine check_jdot_mb_conditions(b, s, apply_jdot_mb, qconv_env)
            type (binary_info), pointer :: b
            type (star_info), pointer :: s
            logical, intent(out) :: apply_jdot_mb
            real(dp), intent(out) :: qconv_env

            real(dp) :: qrad_core
            integer :: i, k, id

            include 'formats.inc'

            ! calculate how much of inner region is convective
            qrad_core = 0d0
            do k = s% nz, 1, -1
                if (s% q(k) > b% jdot_mb_qlim_for_check_rad_core .and. & 
                   (qrad_core == 0d0 .or. s% mixing_type(k) /= convective_mixing)) exit
                if (s% mixing_type(k) == convective_mixing) &
                   qrad_core = qrad_core + s% dq(k)
            end do

            ! calculate how much of the envelope
            qconv_env = 0d0
            do k = 1, s% nz
                if (s% q(k) < b% jdot_mb_qlim_for_check_conv_env .and. &
                   (qconv_env == 0d0 .or. s% mixing_type(k) /= convective_mixing)) exit
                if (s% mixing_type(k) == convective_mixing) &
                   qconv_env = qconv_env + s% dq(k)
            end do

            apply_jdot_mb = .true.
            if (qconv_env < b% jdot_mb_min_qconv_env) then
                apply_jdot_mb = .false.
            return
            end if

            if (qconv_env > b% jdot_mb_max_qconv_env) then
                apply_jdot_mb = .false.
            return
            end if

            if (qrad_core > b% jdot_mb_max_qrad_core) then
                apply_jdot_mb = .false.
            return
        end if
         
        end subroutine check_jdot_mb_conditions

        subroutine jdot_mb_routine(binary_id, ierr)
            integer, intent(in) :: binary_id
            integer, intent(out) :: ierr
            integer :: k, nz
            type (binary_info), pointer :: b
            logical :: apply_jdot_mb
            type (star_info), pointer :: s
            real(dp) :: turnover_time, tt_temp, tt_temp_scaled, tt_old, tt_diff
            real(dp) :: vel, vel_ratio, vel_diff, upper_lim, lower_lim, scaled_vel
            real(dp) :: eps_nuc_lim, eps_nuc
            real(dp) :: dr, tau_lim, delta_mag_chk
            real(dp) :: rsun4, qconv_env, two_pi_div_p3, two_pi_div_p2, K2
            real(dp) :: tt_ratio, tt4
            real(dp) :: rot_ratio, rot4
            real(dp) :: rad4
            real(dp) :: v_esc2, v_mod2
            real(dp) :: alfven_no_R, R_alfven
            real(dp) :: jdot_mb_old, jdot_mb_new
            real(dp) :: conv_env_r, conv_env_m, sonic_cross_time, mag_field, f2
            common/ old_var/ tt_old
            logical :: conv_env_found
            ierr = 0
            call binary_ptr(binary_id, b, ierr)
            if (ierr /= 0) then
                write(*,*) 'failed in binary_ptr'
                return
            end if

            ! write (*,*) " "
            ! write (*,*) "====================================================="
            ! write (*,*) 'doing jdot'
            ! write (*,*) "====================================================="
            ! write (*,*) " "

        ! INITIALIZE THE VARIABLES

            s => b% s_donor
            nz = s% nz
            vel_ratio = s% x_ctrl(1)
            tau_lim = s% x_ctrl(2)

            conv_env_found = .false.

            turnover_time = 0.0
            tt_temp = 0.0
            tt_temp_scaled = 0.0

            sonic_cross_time = 0.0
            conv_env_m = 0.0
            conv_env_r = 0.0
            mag_field = 0.0 
            f2 = 0.0

            eps_nuc_lim = 1.0d-2
            vel_diff = 0.0
            scaled_vel = 0.0

        ! INITIAL TURNOVER TIME CALCULATION

            do k = nz, 1, -1 ! beginning of do loop to calculate convective turnover time

                eps_nuc = s% eps_nuc(k)
                ! check if the cell we are looping through satisfies our convection criteria
                if ((s% gradr(k) .gt. s% grada(k)) .and. (eps_nuc .lt. eps_nuc_lim)) then
                    ! toggle the boolean to begin integration
                    conv_env_found = .true.
                end if

                ! only enter this portion if the convective boolean is true
                ! this loop will go from the innermost cell that is convective to 
                ! the surface. This is to try and smooth through any numeric issues
                ! with convective zones appearing and disappearing in MESA.
                if (conv_env_found) then

                    ! loop to calculate the size of the cell, the innermost cell
                    ! needs special consideration as it is above the core
                    if (k .lt. s% nz) then
                        dr = (s% r(k) - s% r(k + 1))
                    else
                        dr = (s% r(k) - s% R_center)
                    end if
                    
                    ! determine the convective velocity inside each given cell
                    if (s% mixing_type(k) == convective_mixing) then

                        ! need to ensure that the convective velocity is within
                        ! our defined limits, if they are outside of these limits
                        ! set them to be the max/min value allowed.
                        vel = s% conv_vel(k)
                        lower_lim = vel_ratio * s% csound(k)
                        upper_lim = 1.0 * s% csound(k)

                        if (vel .lt. lower_lim) then
                            vel = lower_lim
                        else if (vel .gt. upper_lim) then
                            vel = upper_lim
                        end if
                    
                    ! if the cell isnt defined by MESA to be convective take the
                    ! convective velocity to be equal to sound speed
                    else
                        vel = s% csound(k)
                    end if

                    ! Final check involving the opacity of the given cell. If the 
                    ! cell isn't near the surface (low tau) then include it in our integration
                    if (s% tau(k) .gt. tau_lim) then
                        sonic_cross_time = sonic_cross_time + (dr / s% csound(k))
                        conv_env_r = conv_env_r + dr
                        conv_env_m = conv_env_m + s% dm(k)
                        tt_temp = tt_temp + (dr / vel)
                    end if
                end if

            end do ! end of do loop to calculate convective turnover time

            ! reset the boolean just in case
            conv_env_found = .false.

        ! TURNOVER TIME CHECK, THIS IS TO TRY AND AVOID LARGE CHANGES

            ! simply set the turnover time to the internal variable calculated above
            turnover_time = tt_temp

            if (s% model_number .gt. 1) then
                ! calculate the variables used to check if our system is rapidly evolving
                tt_diff = abs(tt_old - tt_temp) / tt_old
                delta_mag_chk = s% dt / tt_old

                ! write (*,*) "tt_diff = ", tt_diff
                ! write (*,*) "delta_mag = ", delta_mag_chk
                ! write (*,*) "turnover_time = ", turnover_time
                ! write (*,*) "tt_old = ", tt_old

                ! check if timesteps are very small or if the relative change is very large
                if (tt_diff .gt. delta_mag_chk) then 
                    write (*,*) "large change, adjusting accordingly"
                    turnover_time = tt_old + (tt_temp - tt_old) * min((s% dt / tt_old), 0.5)
                    mag_field = (turnover_time / 2.8d6) * (2073600. / b% period)
 
                end if ! end of timestep/relative change check
            end if

            ! remember the current values to be used as comparison in the next step

            tt_old = turnover_time

        ! MAGNETIC BRAKING CALCULATION

            b% jdot_mb = 0
            rsun4 = pow4(rsun)

            call check_jdot_mb_conditions(b, b% s_donor, apply_jdot_mb, qconv_env)

            two_pi_div_p3 = (2.0*pi/b% period)*(2.0*pi/b% period)*(2.0*pi/b% period)
            two_pi_div_p2 = (2.0*pi/b% period)*(2.0*pi/b% period)

            K2 = 0.07 * 0.07

            ! use the formula from rappaport, verbunt, and joss.  apj, 275, 713-731. 1983.
            if (apply_jdot_mb .or. b% keep_mb_on) then

                jdot_mb_old = -3.8d-30*b% m(b% d_i)*rsun4* &            
                    pow(min(b% r(b% d_i),b% rl(b% d_i))/rsun,b% magnetic_braking_gamma)* &
                    two_pi_div_p3

                tt_ratio = turnover_time / 2.8d6
                tt4 = pow4(tt_ratio)
                ! write (*,*) "tt4 = ", tt4

                rot_ratio = (2073600. / b% period )
                rot4 = pow4(rot_ratio)
                ! write (*,*) "rot4 = ", rot4

                rad4 = pow4(b% r(b% d_i))
                ! write (*,*) "rad4 = ", rad4

                v_esc2 = 2.0 * standard_cgrav * b% m(b% d_i) / b% r(b% d_i)
                v_mod2 = v_esc2 + 2.0 * two_pi_div_p2 * b% r(b% d_i) * b% r(b% d_i) / K2 
                ! write (*,*) "G = ", standard_cgrav
                ! write (*,*) "M = ", b% m(b% d_i)
                ! write (*,*) "R = ", b% r(b% d_i)
                ! write (*,*) "v_esc2 = ", v_esc2
                ! write (*,*) "v_mod2 = ", v_mod2

                alfven_no_R = rad4 * rot4 * tt4 / (b% mdot_system_wind(b% d_i) * b% mdot_system_wind(b% d_i)) * (1.0 / v_mod2)

                R_alfven = b% r(b% d_i) * alfven_no_R**(1.d0/3.d0)

                ! HALVED VALUE FOR REF CHECK
                jdot_mb_new = (1.0/3.0) * (2.0*pi/b% period) * b% mdot_system_wind(b% d_i) * R_alfven * R_alfven 
                ! write (*,*) "jdot_mb_new = ", jdot_mb_new
                ! write (*,*) "jdot_mb_old = ", jdot_mb_old

                b% jdot_mb = 2*jdot_mb_new
                
            end if

            if (b% point_mass_i == 0 .and. b% include_accretor_mb .and. &
                (apply_jdot_mb .or. b% keep_mb_on)) then
                f2 = 1.0
                b% jdot_mb = b% jdot_mb - &
                                    3.8d-30*b% m(b% a_i)*rsun4* &
                                    pow(min(b% r(b% a_i),b% rl(b% a_i))/rsun,b% magnetic_braking_gamma)* &
                                    two_pi_div_p3
            end if

            s% xtra(1) = turnover_time
            s% xtra(2) = mag_field
            s% xtra(3) = conv_env_r
            s% xtra(4) = conv_env_m
            s% xtra(5) = sonic_cross_time
            s% xtra(6) = R_alfven                    

        end subroutine jdot_mb_routine

        integer function how_many_extra_binary_history_columns(binary_id)
            use binary_def, only: binary_info
            integer, intent(in) :: binary_id
            how_many_extra_binary_history_columns = 6
        end function how_many_extra_binary_history_columns
        
        subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
            use const_def, only: dp
            type (binary_info), pointer :: b
            type (star_info), pointer :: s
            integer, intent(in) :: binary_id
            integer, intent(in) :: n
            character (len=maxlen_binary_history_column_name) :: names(n)
            real(dp) :: vals(n)
            integer, intent(out) :: ierr
            real(dp) :: beta
            ierr = 0
            call binary_ptr(binary_id, b, ierr)
            if (ierr /= 0) then
                write(*,*) 'failed in binary_ptr'
                return
            end if

            s => b% s_donor
            
            names(1) = "turnover_time"
            vals(1) = s% xtra(1)

            names(2) = "mag_field"
            vals(2) = s% xtra(2)

            names(3) = "conv_env_r"
            vals(3) = s% xtra(3)

            names(4) = "conv_env_m"
            vals(4) = s% xtra(4)

            names(5) = "sonic_cross_time"
            vals(5) = s% xtra(5)

            names(6) = "R_alfven"
            vals(6) = s% xtra(6)

        end subroutine data_for_extra_binary_history_columns
        
        
        integer function extras_binary_startup(binary_id,restart,ierr)
            type (binary_info), pointer :: b
            integer, intent(in) :: binary_id
            integer, intent(out) :: ierr
            logical, intent(in) :: restart
            call binary_ptr(binary_id, b, ierr)
            if (ierr /= 0) then ! failure in  binary_ptr
                return
            end if
            
            ! b% s1% job% warn_run_star_extras = .false.
            ! extras_binary_startup = keep_going

            write(*,*) "starting modified MB"
            if (b% angular_momentum_j <= 0) then
                write(*,*) "angular_momentum_j <= 0", b% angular_momentum_j
                b% s1% dt_next = min(b% s1% dt * 0.50, b% s1% dt_next)
                extras_binary_startup = retry
            end if
        end function  extras_binary_startup
        
        !Return either retry,backup,keep_going or terminate
        integer function extras_binary_check_model(binary_id)
            type (binary_info), pointer :: b
            integer, intent(in) :: binary_id
            integer :: ierr
            real(dp) :: j_check1
            call binary_ptr(binary_id, b, ierr)
            if (ierr /= 0) then ! failure in  binary_ptr
                return
            end if

            ! write (*,*) " "
            ! write (*,*) "====================================================="
            ! write (*,*) "checking model"
            ! write (*,*) "====================================================="
            ! write (*,*) " "

            ! write (*,*) "angular momentum before jdot = ", b% angular_momentum_j

            if (b% angular_momentum_j <= 0) then
                write(*,*) "bad angular momentum"
                b% s1% dt_next = min(b% s1% dt * 0.50, b% s1% dt_next)
                extras_binary_check_model = retry
            end if
          
            extras_binary_check_model = keep_going
        end function extras_binary_check_model
        
        
        ! returns either keep_going or terminate.
        ! note: cannot request retry or backup; extras_check_model can do that.
        integer function extras_binary_finish_step(binary_id)
            type (binary_info), pointer :: b
            integer, intent(in) :: binary_id
            integer :: ierr
            call binary_ptr(binary_id, b, ierr)
            if (ierr /= 0) then ! failure in binary_ptr
                return
            end if  
            extras_binary_finish_step = keep_going
            
        end function extras_binary_finish_step
        
        subroutine extras_binary_after_evolve(binary_id, ierr)
            type (binary_info), pointer :: b
            integer, intent(in) :: binary_id
            integer, intent(out) :: ierr
            call binary_ptr(binary_id, b, ierr)
            if (ierr /= 0) then ! failure in binary_ptr
                return
            end if
 
        end subroutine extras_binary_after_evolve      
        
        end module run_binary_extras
