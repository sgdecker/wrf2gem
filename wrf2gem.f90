! The wrf2gem program reads a WRF history file in netCDF format and outputs
! specified variables (from convert.nl) into a new GEMPAK file.


!   Copyright 2015 Steven G. Decker

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.


! Identifiers within this program use the following conventions:
! o Variables are typed like numOuts
! o Constants are typed like GasConstant
! o Functions and subroutines generally have underscores.

program wrf2gem
  use wrf2gem_subs, only: get_out_fields, get_num_files, open_wrf_file,  &
                          get_grid_info, update_grid_info, mem_error,  &
                          get_times, close_wrf_file
  use gempak,       only: init_gem, create_gemfile, close_and_exit_gem
  use registry,     only: init_reg, output_var
  implicit none

  !
  ! Initialization
  !

  real,    dimension(3)          :: ang
  real,    dimension(2)          :: lat, lon
  integer, dimension(:), pointer :: outList
  integer :: numOuts, pb, pt, dp, numFiles, iCount, ncid, nx, ny, nz,  &
             numTimes, proj, timesPer12hr, i, igdfln, status
  character(len=80), dimension(:), pointer     :: outFiles
  character(len=15), dimension(:), allocatable :: times
  logical :: ok

  !
  ! Execution
  !

  nullify(outList, outFiles)

  ! Find out from namelist which fields we should output to GEMPAK
  call get_out_fields(numOuts, outList, pb, pt, dp)
  call get_num_files(numFiles, outFiles)
  if (numOuts < 1 .or. numFiles < 1) stop "Nothing to output!"

  ! Do this for each WRF history file
  do iCount = 1, numFiles
     write (*,"(3a)") "Processing ", trim(outFiles(iCount)), "..."

     ! Open WRF history file
     call open_wrf_file(outFiles(iCount), ncid)

     if (iCount == 1) then
        ! Determine grid size, map projection, and number of output times
        ! Note nx, ny are based on number of grid boxes, i.e., cross points
        call get_grid_info(ncid, nx, ny, nz, numTimes, proj, lat, lon, ang)
     else
        ! Verify that grid projection matches the first file, and get new
        ! number of times in file.
        call update_grid_info(ncid, nx, ny, nz, proj, lat, lon, ang,  &
             numTimes, ok)
        if (.not. ok) then
           write (*,"(3a)") "Skipping ", trim(outFiles(iCount)), " because of&
                & grid navigation mismatch with previous netCDF files."
           call close_wrf_file(ncid)
           cycle
        end if
     end if
     allocate (times(numTimes), stat=status)
     call mem_error(status, 1, "wrf2gem.f90")
     call get_times(iCount, outFiles(1:numFiles), ncid, times, timesPer12hr)
     
     ! (Re)initialize registry, since timesPer12hr might have changed
     call init_reg(timesPer12hr, pb, pt, dp)
     
     ! Initialize GEMPAK and create or add to GEMPAK file
     call init_gem
     call create_gemfile(proj, nx, ny, lat, lon, ang, igdfln, ok)
     if (.not. ok) then
        write (*,"(3a)") "Skipping ", trim(outFiles(iCount)), " because of&
             & grid navigation mismatch with GEMPAK file."
     else
        ! Main loop over output fields
        do i = 1, numOuts
           call output_var(ncid, igdfln, nx, ny, nz, outList(i), iCount,  &
                outFiles(1:numFiles), times)
        end do
     end if

     ! Clean up
     call close_wrf_file(ncid)
     call close_and_exit_gem(igdfln)
     deallocate (times, stat=status)
     call mem_error(status, 2, "wrf2gem.f90")
  end do

  print *, "Successful completion!"
end program wrf2gem
