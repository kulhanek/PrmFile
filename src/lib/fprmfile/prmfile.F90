!===============================================================================
! PrmFileLib - Parametr file parsing and manipulation library
!-------------------------------------------------------------------------------
! (c) 2007 Petr Kulhanek (kulhanek@enzim.hu)
!          Institute of Enzymology, Karolina ut 29, Budapest H-1113, Hungary
!-------------------------------------------------------------------------------
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, 
!    Boston, MA  02110-1301  USA
!===============================================================================

module prmfile

use prmfile_dat
use prmfile_core

implicit none
contains

!===============================================================================
! subroutine prmfile_init(prmfile)
! ------------------------------------------------------------------------------
! initialize object
!===============================================================================

subroutine prmfile_init(prmfile)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 !------------------------------------------------------------------------------

 prmfile%FirstGroup      => NULL()
 prmfile%CurrentGroup    => NULL()
 prmfile%CurrentSection  => NULL()
 prmfile%CurrentLine     => NULL()
 prmfile%FieldPosition   = 0

end subroutine prmfile_init

!===============================================================================
! logical function prmfile_read(prmfile,filename)
! ------------------------------------------------------------------------------
! this function load file and decompose it to groups, sections and lines, which
! will be stored to memory accessbile via prmfile structure
!===============================================================================

logical function prmfile_read(prmfile,filename)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: filename
 ! -----------------------------------------------
 integer                            :: stat
 !------------------------------------------------------------------------------

 ! clear previous data -----------------
 call prmfile_clear(prmfile)

#if defined(__IBM__) || defined(__IBMC__)
 ! this is necessary on systems with xlf compiler
 ! control file is usually opened more times
 ! this on xlf systems leads to an error with
 ! iostat = 23
 ! 'Attempt to connect a file that is already connected to another unit.'
 call setrteopts("multconn=yes")
#endif

 ! open file ---------------------------
 open(unit=PRMFILE_UNIT, file = filename, status = 'old', form='formatted', &
      action='read', iostat = stat, access='sequential')
 if(stat .ne. 0) then
    prmfile_read = .false.
    return
 end if
 prmfile%Name = filename
 ! process file ------------------------
 if( prmfile_parse(prmfile) ) then
    prmfile_read = .true.
 else
    prmfile_read = .false.
    call prmfile_clear(prmfile)
 end if

 ! close file --------------------------
 close(PRMFILE_UNIT)

 return

end function prmfile_read

!===============================================================================
! subroutine prmfile_clear(prmfile)
! ------------------------------------------------------------------------------
! this subroutine remove all data associated with prmfile
!===============================================================================

subroutine prmfile_clear(prmfile)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 ! -----------------------------------------------
 type(GROUP_TYPE), pointer          :: pgroup,cgroup
 type(SECTION_TYPE), pointer        :: psection,csection
 type(LINE_TYPE), pointer           :: pline,cline
 !------------------------------------------------------------------------------

 !loop over groups
 pgroup => prmfile%FirstGroup

 do while(associated(pgroup))
    cgroup => pgroup

    ! loop over sections
    psection => cgroup%FirstSection

    do while(associated(psection))
        csection => psection

        ! loop over lines
        pline => csection%FirstLine

        do while(associated(pline))
            cline => pline
            pline => cline%Next
            deallocate(cline)
        end do

        psection => csection%Next
        deallocate(csection)

    end do

    pgroup => cgroup%Next
    deallocate(cgroup)

 end do

 prmfile%Name = ''
 nullify(prmfile%FirstGroup)
 nullify(prmfile%CurrentGroup)
 nullify(prmfile%CurrentSection)
 nullify(prmfile%CurrentLine)
 prmfile%FieldPosition = 0

 return

end subroutine prmfile_clear

!===============================================================================
! subroutine prmfile_dump(prmfile,ounit,unprocessed)
! ------------------------------------------------------------------------------
! it prints parameter file prmfile to ounit stream
! optional parameter unproccessed determines if only onprocessed items
! are printed
! this subroutine does not influence the iterator states
!===============================================================================

subroutine prmfile_dump(prmfile,ounit,unprocessed)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 integer                            :: ounit
 logical, intent(in), optional      :: unprocessed
 ! -----------------------------------------------
 type(GROUP_TYPE), pointer          :: pgroup
 type(SECTION_TYPE), pointer        :: psection
 type(LINE_TYPE), pointer           :: pline
 logical                            :: only_unprocessed
 !------------------------------------------------------------------------------

 only_unprocessed = .false.
 if(present(unprocessed)) only_unprocessed = unprocessed

 !loop over groups
 pgroup => prmfile%FirstGroup

 do while(associated(pgroup))
    if( .not. only_unprocessed ) then
        write(ounit,100) trim(pgroup%name)
    else
        if( prmfile_count_ulines_in_grp(pgroup) .gt. 0 ) then
            write(ounit,100) trim(pgroup%name)
        end if
    end if

    ! loop over sections
    psection => pgroup%FirstSection

    do while(associated(psection))
        if( .not. only_unprocessed ) then
            write(ounit,200) trim(psection%name)
        else
            if( prmfile_count_ulines_in_sec(psection) .gt. 0 ) then
                write(ounit,200) trim(psection%name)
            end if
        end if

        ! loop over lines
        pline => psection%FirstLine

        do while(associated(pline))
            if( .not. only_unprocessed ) then
                write(ounit,300) trim(pline%text)
            else
                if( .not. pline%Processed ) then
                    write(ounit,300) trim(pline%text)
                end if
            end if
            pline => pline%Next
         end do

        psection => psection%Next
    end do

    pgroup => pgroup%Next
 end do

 return

100 format('{',A,'}')
200 format('[',A,']')
300 format(A)

end subroutine prmfile_dump

!===============================================================================
! subroutine prmfile_dump_group(prmfile,ounit,gname,unprocessed)
! ------------------------------------------------------------------------------
! it prints gname group to ounit stream
! optional parameter unproccessed determines if only onprocessed items
! are printed
! this subroutine does not influence the iterator states
!===============================================================================

subroutine prmfile_dump_group(prmfile,ounit,gname,unprocessed)

 implicit none
 type(PRMFILE_TYPE)                     :: prmfile
 integer                                :: ounit
 character(*),optional                  :: gname
 logical, intent(in), optional          :: unprocessed
 ! -----------------------------------------------
 type(GROUP_TYPE), pointer              :: pgroup
 type(SECTION_TYPE), pointer            :: psection
 type(LINE_TYPE), pointer               :: pline
 logical                                :: only_unprocessed
 character(len=PRMFILE_MAX_GROUP_NAME)  :: lgname
 character(len=PRMFILE_MAX_GROUP_NAME)  :: lname
 !------------------------------------------------------------------------------

 only_unprocessed = .false.
 if(present(unprocessed)) only_unprocessed = unprocessed

 !loop over groups
 pgroup => prmfile%FirstGroup

 lgname = gname
 call prmfile_upcase(lgname)

 do while(associated(pgroup))
    lname = pgroup%Name
    call prmfile_upcase(lname)
    if( trim(lname) .eq. trim(lgname) ) exit
    pgroup => pgroup%Next
 end do

 if(associated(pgroup)) then
    if( .not. only_unprocessed ) then
        write(ounit,100) trim(pgroup%name)
    else
        if( prmfile_count_ulines_in_grp(pgroup) .gt. 0 ) then
            write(ounit,100) trim(pgroup%name)
        end if
    end if

    ! loop over sections
    psection => pgroup%FirstSection

    do while(associated(psection))
        if( .not. only_unprocessed ) then
            write(ounit,200) trim(psection%name)
        else
            if( prmfile_count_ulines_in_sec(psection) .gt. 0 ) then
                write(ounit,200) trim(psection%name)
            end if
        end if

        ! loop over lines
        pline => psection%FirstLine

        do while(associated(pline))
            if( .not. only_unprocessed ) then
                write(ounit,300) trim(pline%text)
            else
                if( .not. pline%Processed ) then
                    write(ounit,300) trim(pline%text)
                end if
            end if
            pline => pline%Next
         end do

        psection => psection%Next
    end do
 end if

 return

100 format('{',A,'}')
200 format('[',A,']')
300 format(A)

end subroutine prmfile_dump_group

!===============================================================================
! integer function prmfile_count_ulines(prmfile)
! ------------------------------------------------------------------------------
! it returns total number of unprocessed lines
!===============================================================================

integer function prmfile_count_ulines(prmfile,gname)

 implicit none
 type(PRMFILE_TYPE)                     :: prmfile
 character(*),optional                  :: gname
 ! -----------------------------------------------
 type(GROUP_TYPE), pointer              :: pgroup
 type(SECTION_TYPE), pointer            :: psection
 type(LINE_TYPE), pointer               :: pline
 character(len=PRMFILE_MAX_GROUP_NAME)  :: lgname
 character(len=PRMFILE_MAX_GROUP_NAME)  :: lname
 !------------------------------------------------------------------------------

 prmfile_count_ulines = 0

 pgroup => prmfile%FirstGroup

 if( present(gname) ) then
    lgname = gname
    call prmfile_upcase(lgname)
    do while(associated(pgroup))
        lname = pgroup%Name
        call prmfile_upcase(lname)
        if( trim(lname) .eq. trim(lgname) ) exit
        pgroup => pgroup%Next
    end do
    if( associated(pgroup) ) then
        prmfile_count_ulines = prmfile_count_ulines_in_grp(pgroup)
    end if
    return
 end if

 !loop over groups
 do while(associated(pgroup))
    ! loop over sections
    psection => pgroup%FirstSection
    do while(associated(psection))
        ! loop over lines
        pline => psection%FirstLine
        do while(associated(pline))
            if( .not. pline%Processed ) then
                prmfile_count_ulines = prmfile_count_ulines + 1
            end if
            pline => pline%Next
        end do
        psection => psection%Next
    end do
    pgroup => pgroup%Next
 end do

 return

end function prmfile_count_ulines

!===============================================================================
! subroutine prmfile_first_group(prmfile)
! ------------------------------------------------------------------------------
! it rewinds file to the first group
!===============================================================================

logical function prmfile_first_group(prmfile)

 implicit none
 type(PRMFILE_TYPE)                         :: prmfile
 !------------------------------------------------------------------------------

 nullify(prmfile%CurrentGroup)
 nullify(prmfile%CurrentSection)
 nullify(prmfile%CurrentLine)
 prmfile%FieldPosition = 0

 prmfile%CurrentGroup => prmfile%FirstGroup
 if( associated(prmfile%CurrentGroup) ) then
    prmfile%CurrentSection => prmfile%CurrentGroup%FirstSection
    if( associated(prmfile%CurrentSection) )then
        prmfile%CurrentLine => prmfile%CurrentSection%FirstLine
    end if
 end if

 prmfile_first_group = associated(prmfile%CurrentGroup)
 if( prmfile_first_group ) then
    prmfile%CurrentGroup%Processed = .true.
 end if

 return

end function prmfile_first_group

!===============================================================================
! logical function prmfile_open_group(prmfile,name)
! ------------------------------------------------------------------------------
! it opens FIRST group with name in the file
!===============================================================================

logical function prmfile_open_group(prmfile,name)

 implicit none
 type(PRMFILE_TYPE)                         :: prmfile
 character(*)                               :: name
 !------------------------------------------------
 character(len=PRMFILE_MAX_GROUP_NAME)      :: uname
 character(len=PRMFILE_MAX_GROUP_NAME)      :: lname
 !------------------------------------------------------------------------------

 prmfile_open_group = .false.

 if( .not. prmfile_first_group(prmfile) ) return

 uname = name
 call prmfile_upcase(uname)

 do while(associated(prmfile%CurrentGroup))
    lname = prmfile%CurrentGroup%Name
    call prmfile_upcase(lname)
    if( trim(lname) .eq. trim(uname) ) then
        prmfile%CurrentSection => prmfile%CurrentGroup%FirstSection
        if( associated(prmfile%CurrentGroup%FirstSection) ) then
            prmfile%CurrentLine => prmfile%CurrentGroup%FirstSection%FirstLine
        end if
        prmfile%CurrentGroup%Processed = .true.
        prmfile_open_group = .true.
        return
    end if
    prmfile%CurrentGroup => prmfile%CurrentGroup%Next
 end do

 return

end function prmfile_open_group

!===============================================================================
! logical function prmfile_next_group(prmfile)
! ------------------------------------------------------------------------------
! it moves to next group
!===============================================================================

logical function prmfile_next_group(prmfile)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 !------------------------------------------------------------------------------

 if( associated(prmfile%CurrentGroup) ) then
    prmfile%CurrentGroup => prmfile%CurrentGroup%Next
 end if

 nullify(prmfile%CurrentSection)
 nullify(prmfile%CurrentLine)
 prmfile%FieldPosition = 0

 if( associated(prmfile%CurrentGroup) ) then
    prmfile%CurrentSection => prmfile%CurrentGroup%FirstSection
    if( associated(prmfile%CurrentSection)  ) then
        prmfile%CurrentLine => prmfile%CurrentSection%FirstLine
    end if
 end if

 prmfile_next_group = associated(prmfile%CurrentGroup)

 if( prmfile_next_group ) then
    prmfile%CurrentGroup%Processed = .true.
 end if

 return

end function prmfile_next_group

!===============================================================================
! integer function prmfile_count_group(prmfile)
! ------------------------------------------------------------------------------
! it counts number of sections in current group
!===============================================================================

integer function prmfile_count_group(prmfile)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 !------------------------------------------------
 type(SECTION_TYPE), pointer        :: psection
 !------------------------------------------------------------------------------

 prmfile_count_group = 0

 nullify(psection)
 if( associated(prmfile%CurrentGroup) ) then
    psection => prmfile%CurrentGroup%FirstSection
 end if

 do while(associated(psection))
    prmfile_count_group = prmfile_count_group + 1
    psection => psection%Next
 end do

 return

end function prmfile_count_group

!===============================================================================
! logical function prmfile_get_group_name(prmfile,name)
! ------------------------------------------------------------------------------
! it returns name of current group
!===============================================================================

logical function prmfile_get_group_name(prmfile,name)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: name
 !------------------------------------------------------------------------------

 if( associated(prmfile%CurrentGroup) ) then
    name = prmfile%CurrentGroup%Name
 end if

 prmfile_get_group_name = associated(prmfile%CurrentGroup)

 return

end function prmfile_get_group_name

!===============================================================================
! function prmfile_first_section(prmfile)
! ------------------------------------------------------------------------------
! it rewinds file to the first section of current group
!===============================================================================

logical function prmfile_first_section(prmfile)

 implicit none
 type(PRMFILE_TYPE)                         :: prmfile
 !------------------------------------------------------------------------------

 nullify(prmfile%CurrentSection)
 nullify(prmfile%CurrentLine)
 prmfile%FieldPosition = 0

 if(associated(prmfile%CurrentGroup)) then
    prmfile%CurrentSection => prmfile%CurrentGroup%FirstSection
    if(associated(prmfile%CurrentSection)) then
        prmfile%CurrentLine => prmfile%CurrentSection%FirstLine
    end if
 end if

 prmfile_first_section = associated(prmfile%CurrentSection)

 if( prmfile_first_section ) then
    prmfile%CurrentSection%Processed = .true.
 end if

 return

end function prmfile_first_section

!===============================================================================
! logical function prmfile_open_section(prmfile,name)
! ------------------------------------------------------------------------------
! it opens FIRST section with name in  if( prmfile_first_section ) then
!===============================================================================

logical function prmfile_open_section(prmfile,name)

 implicit none
 type(PRMFILE_TYPE)                         :: prmfile
 character(*)                               :: name
 !------------------------------------------------
 character(len=PRMFILE_MAX_SECTION_NAME)    :: uname
 character(len=PRMFILE_MAX_SECTION_NAME)    :: lname
 !------------------------------------------------------------------------------

 prmfile_open_section = .false.

 if( .not. prmfile_first_section(prmfile) ) return

 uname = name
 call prmfile_upcase(uname)

 do while(associated(prmfile%CurrentSection))
    lname = prmfile%CurrentSection%Name
    call prmfile_upcase(lname)
    if( trim(lname) .eq. trim(uname) ) then
        if(associated(prmfile%CurrentSection)) then
            prmfile%CurrentLine => prmfile%CurrentSection%FirstLine
        end if
        prmfile%CurrentSection%Processed = .true.
        prmfile_open_section = .true.
        return
    end if
    prmfile%CurrentSection => prmfile%CurrentSection%Next
 end do

 return

end function prmfile_open_section

!===============================================================================
! integer function prmfile_open_section_and_count(prmfile,name)
! ------------------------------------------------------------------------------
! it opens section and return number of lines in it
!===============================================================================

integer function prmfile_open_section_and_count(prmfile,name)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: name
 !------------------------------------------------------------------------------

 if(prmfile_open_section(prmfile,name)) then
    prmfile_open_section_and_count = prmfile_count_section(prmfile)
 else
    prmfile_open_section_and_count = 0
 end if

 return

end function prmfile_open_section_and_count

!===============================================================================
! logical function prmfile_next_section(prmfile)
! ------------------------------------------------------------------------------
! it moves to next section
!===============================================================================

logical function prmfile_next_section(prmfile)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 !------------------------------------------------------------------------------

 if( associated(prmfile%CurrentSection) ) then
    prmfile%CurrentSection => prmfile%CurrentSection%Next
 end if

 nullify(prmfile%CurrentLine)
 prmfile%FieldPosition = 0

 if( associated(prmfile%CurrentSection) ) then
    prmfile%CurrentLine => prmfile%CurrentSection%FirstLine
 end if 

 prmfile_next_section = associated(prmfile%CurrentSection)

 if( prmfile_next_section ) then
    prmfile%CurrentSection%Processed = .true.
 end if

 return

end function prmfile_next_section

!===============================================================================
! integer function prmfile_count_section(prmfile)
! ------------------------------------------------------------------------------
! it counts number of line in current section
!===============================================================================

integer function prmfile_count_section(prmfile)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 !------------------------------------------------
 type(LINE_TYPE), pointer           :: pline
 !------------------------------------------------------------------------------

 prmfile_count_section = 0

 nullify(pline)
 if( associated(prmfile%CurrentSection) ) then
    pline => prmfile%CurrentSection%FirstLine
 end if

 do while(associated(pline))
    prmfile_count_section = prmfile_count_section + 1
    pline => pline%Next
 end do

 return

end function prmfile_count_section

!===============================================================================
! logical function prmfile_get_section_name(prmfile,name)
! ------------------------------------------------------------------------------
! it returns name of current section
!===============================================================================

logical function prmfile_get_section_name(prmfile,name)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: name
 !------------------------------------------------------------------------------

 if( associated(prmfile%CurrentSection) ) then
    name = prmfile%CurrentSection%Name
 end if

 prmfile_get_section_name = associated(prmfile%CurrentSection)

 return

end function prmfile_get_section_name

!===============================================================================
! logical function prmfile_first_line(prmfile)
! ------------------------------------------------------------------------------
! it rewinds file to the first line of current section of current group
!===============================================================================

logical function prmfile_first_line(prmfile)

 implicit none
 type(PRMFILE_TYPE)                         :: prmfile
 !------------------------------------------------------------------------------

 nullify(prmfile%CurrentLine)
 prmfile%FieldPosition = 0

 if(associated(prmfile%CurrentSection)) then
    prmfile%CurrentLine => prmfile%CurrentSection%FirstLine
 end if

 prmfile_first_line = associated(prmfile%CurrentSection)

 return

end function prmfile_first_line

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_next_line(prmfile)

 implicit none
 type(PRMFILE_TYPE)                         :: prmfile
 !------------------------------------------------------------------------------

 prmfile%FieldPosition = 0

 if(associated(prmfile%CurrentSection)) then
    if( associated(prmfile%CurrentLine) ) then
        prmfile%CurrentLine => prmfile%CurrentLine%Next
    end if
 end if

 prmfile_next_line = associated(prmfile%CurrentLine)

 return

end function prmfile_next_line

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_integer_by_key(prmfile,key,value)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: key
 integer                            :: value
 !------------------------------------------------
 character(len=PRMFILE_MAX_VALUE)   :: locvalue
 integer                            :: stat
 ! -----------------------------------------------------------------------------

 ! find key --------------------------------------
 if( prmfile_find_key(prmfile,key,locvalue) ) then
    ! try to decode data
    read(locvalue, fmt=*, iostat=stat) value
    if (stat == 0 ) then
        prmfile_get_integer_by_key = .true.
        call prmfile_set_kline_as_processed(prmfile)
        return
    end if
 end if

 ! do not change value if default argument is not provided
 prmfile_get_integer_by_key = .false.

 return

end function prmfile_get_integer_by_key

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_real_by_key(prmfile,key,value)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: key
 real                               :: value
 !------------------------------------------------
 character(len=PRMFILE_MAX_VALUE)   :: locvalue
 integer                            :: stat
 ! -----------------------------------------------------------------------------

 ! find key --------------------------------------
 if( prmfile_find_key(prmfile,key,locvalue) ) then
    ! try to decode data
    read(locvalue, fmt=*, iostat=stat) value
    if (stat == 0 ) then
        prmfile_get_real_by_key = .true.
        call prmfile_set_kline_as_processed(prmfile)
        return
    end if
 end if

 ! do not change value if default argument is not provided
 prmfile_get_real_by_key = .false.

 return

end function prmfile_get_real_by_key

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_real8_by_key(prmfile,key,value)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: key
 real(8)                            :: value
 !------------------------------------------------
 character(len=PRMFILE_MAX_VALUE)   :: locvalue
 integer                            :: stat
 ! -----------------------------------------------------------------------------

 ! find key --------------------------------------
 if( prmfile_find_key(prmfile,key,locvalue) ) then
    ! try to decode data
    read(locvalue, fmt=*, iostat=stat) value
    if (stat == 0 ) then
        prmfile_get_real8_by_key = .true.
        call prmfile_set_kline_as_processed(prmfile)
        return
    end if
 end if

 ! do not change value if default argument is not provided
 prmfile_get_real8_by_key = .false.

 return

end function prmfile_get_real8_by_key

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_string_by_key(prmfile,key,value)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: key
 character(*)                       :: value
 !------------------------------------------------
 character(len=PRMFILE_MAX_VALUE)   :: locvalue
 integer                            :: sbeg,send
 ! -----------------------------------------------------------------------------

 ! do not change value if default argument is not provided
 prmfile_get_string_by_key = .false.

 ! find key --------------------------------------
 if( .not. prmfile_find_key(prmfile,key,locvalue) ) then
    return
 end if

 sbeg = 1
 send = len(trim(locvalue))

 ! handle quotations
 if( send .ge. 2 ) then
    if( (locvalue(1:1) .eq. '''') .and. (locvalue(send:send) .eq. '''')  ) then
        sbeg = sbeg + 1
        send = send - 1
    else if( (locvalue(1:1) .eq. '"') .and. (locvalue(send:send) .eq. '"') ) then
        sbeg = sbeg + 1
        send = send - 1
    end if
 end if

 if( sbeg .le. send ) then
    value = locvalue(sbeg:send)
 else
    value = ''
 end if

 call prmfile_set_kline_as_processed(prmfile)
 prmfile_get_string_by_key = .true.
 return

end function prmfile_get_string_by_key

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_logical_by_key(prmfile,key,value)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: key
 logical                            :: value
 !------------------------------------------------
 character(len=PRMFILE_MAX_VALUE)   :: locvalue
 ! -----------------------------------------------------------------------------

 prmfile_get_logical_by_key = .false.

 if( prmfile_find_key(prmfile,key,locvalue) ) then
    if(     trim(locvalue) == 'on' &
       .or. trim(locvalue) == 'ON' &
       .or. trim(locvalue) == 'true' &
       .or. trim(locvalue) == 'TRUE' &
       .or. trim(locvalue) == 'T' &
       .or. trim(locvalue) =='1') then
        value = .true.
        prmfile_get_logical_by_key = .true.
        call prmfile_set_kline_as_processed(prmfile)
    else if(     trim(locvalue) == 'off' &
            .or. trim(locvalue) == 'OFF' &
            .or. trim(locvalue) == 'false' &
            .or. trim(locvalue) == 'FALSE' &
            .or. trim(locvalue) == 'F' &
            .or. trim(locvalue) =='0') then
        value = .false.
        prmfile_get_logical_by_key = .true.
        call prmfile_set_kline_as_processed(prmfile)
    end if
 end if

 return

end function prmfile_get_logical_by_key

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_line(prmfile,line)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: line
 !------------------------------------------------------------------------------

 prmfile%FieldPosition = 0

 !if there's a line
 if(associated(prmfile%CurrentLine)) then
    line = prmfile%CurrentLine%Text
    prmfile%CurrentLine%Processed = .true.
    prmfile%CurrentLine => prmfile%CurrentLine%Next
    prmfile_get_line = .true.
 else
    line = ''
    prmfile_get_line = .false.
 end if
 return

end function prmfile_get_line

!-------------------------------------------------------------------------------

logical function prmfile_get_current_line(prmfile,line)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: line
 !------------------------------------------------------------------------------

 prmfile%FieldPosition = 0

 !if there's a line
 if(associated(prmfile%CurrentLine)) then
    line = prmfile%CurrentLine%Text
    prmfile%CurrentLine%Processed = .true.
    prmfile_get_current_line = .true.
 else
    line = ''
    prmfile_get_current_line = .false.
 end if
 return

end function prmfile_get_current_line

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

integer function prmfile_max_enum(prmfile,section,count_out)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: section
 integer, optional                  :: count_out
 !------------------------------------------------
 integer                            :: max_enum, enum, count
 logical                            :: dummy
 !------------------------------------------------------------------------------

 max_enum = 0
 count = 0

 if( prmfile_open_section(prmfile,section) ) then
    do while(prmfile_get_int(prmfile,enum))
        if(enum > max_enum) max_enum = enum
    end do
    count = prmfile_count_section(prmfile)
    dummy = prmfile_first_line(prmfile)                 ! and rewind
 end if

 prmfile_max_enum = max_enum
 if(present(count_out)) count_out = count

 return

end function prmfile_max_enum

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_int_int(prmfile,key,value)

 implicit none
 type(PRMFILE_TYPE)                 ::  prmfile
 integer,intent(out)                ::  key
 integer                            ::  value
 !------------------------------------------------
 character(len=PRMFILE_MAX_LINE)    ::  line
 integer                            ::  stat
 !-----------------------------------------------------------------------------

 prmfile_get_int_int = .false.

 if( .not. prmfile_get_line(prmfile,line) ) return

 read(line, fmt=*, iostat=stat) key, value
 prmfile_get_int_int = stat == 0

end function prmfile_get_int_int

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_int(prmfile,value)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 integer                            :: value
 !------------------------------------------------
 character(len=PRMFILE_MAX_LINE)    :: line
 integer                            :: stat
 ! -----------------------------------------------------------------------------

 prmfile_get_int = .false.

 if( .not. prmfile_get_line(prmfile,line) ) return

 read(line, fmt=*, iostat=stat) value

 if(stat .eq. 0) then
        prmfile_get_int = .true.
 end if

 return

end function prmfile_get_int

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_init_field_by_key(prmfile,key)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: key
 !------------------------------------------------
 character(len=2)                   :: ws
 integer                            :: llen,ekey,svalue
 character(len=PRMFILE_MAX_LINE)    :: value
 ! -----------------------------------------------------------------------------

 prmfile_init_field_by_key = .false.

 ! find key
 if( .not. prmfile_find_key(prmfile,key,value) ) return
 if( .not. associated(prmfile%CurrentLine) ) return

 ! init white space characters
 ws(1:1) = char(9)    ! tabulator
 ws(2:2) = ' '        ! space

 ! we need whole line
 value = prmfile%CurrentLine%Text

 llen = len(value)

 ! find end of key
 ekey = scan(value,ws)
 if( ekey .le. 0 ) return       ! line is only key, this is not permitted

 ! find beggining of value
 svalue = verify(value(ekey:llen),ws)

 prmfile%FieldPosition = ekey+svalue-1

 prmfile_init_field_by_key = .true.

 return

end function prmfile_init_field_by_key

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_field_by_key(prmfile,field)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: field
 !------------------------------------------------
 character(len=2)                   :: ws
 integer                            :: slen,istart,iend
 ! -----------------------------------------------------------------------------

 ! init white space characters
 ws(1:1) = char(9)    ! tabulator
 ws(2:2) = ' '        ! space

 prmfile_get_field_by_key = .false.

 if( .not. associated(prmfile%CurrentLine) ) return      ! no data to process
 if( prmfile%FieldPosition .le. 0 ) return ! key was not found

 slen = len_trim(prmfile%CurrentLine%Text)

 istart = verify(prmfile%CurrentLine%Text(prmfile%FieldPosition:slen),ws)
 if( istart .le. 0 ) then
    ! this cannot happen
    stop 'Fatal error in prmfile_get_field'
 end if
 iend  = scan(prmfile%CurrentLine%Text(prmfile%FieldPosition+istart-1:slen),ws)
 if( iend .le. 0 ) then
    iend = slen - prmfile%FieldPosition+istart + 1
 end if
 field = prmfile%CurrentLine%Text(prmfile%FieldPosition+istart-1:prmfile%FieldPosition+istart+iend-2)

 ! correct position
 prmfile%FieldPosition = prmfile%FieldPosition+istart+iend-1

 if( prmfile%FieldPosition .gt. slen ) then
    ! set line as procesessed
    prmfile%CurrentLine%Processed = .true.
    ! now new field
    prmfile%FieldPosition = 0
 end if

 prmfile_get_field_by_key = .true.
 return

end function prmfile_get_field_by_key

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_field(prmfile,field)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: field
 !------------------------------------------------
 character(len=2)                   :: ws
 integer                            :: slen,istart,iend
 ! -----------------------------------------------------------------------------

 ! init white space characters
 ws(1:1) = char(9)    ! tabulator
 ws(2:2) = ' '        ! space

 prmfile_get_field = .false.

 if( .not. associated(prmfile%CurrentLine) ) then
    return      ! no data to process
 end if

 if( prmfile%FieldPosition .eq. 0 ) then
    prmfile%FieldPosition = 1       ! first record
 end if

 slen = len_trim(prmfile%CurrentLine%Text)

 istart = verify(prmfile%CurrentLine%Text(prmfile%FieldPosition:slen),ws)
 if( istart .le. 0 ) then
    ! this cannot happen
    stop 'Fatal error in prmfile_get_field'
 end if
 iend  = scan(prmfile%CurrentLine%Text(prmfile%FieldPosition+istart-1:slen),ws)
 if( iend .le. 0 ) then
    iend = slen - prmfile%FieldPosition+istart + 1
 end if
 field = prmfile%CurrentLine%Text(prmfile%FieldPosition+istart-1:prmfile%FieldPosition+istart+iend-2)

 ! correct position
 prmfile%FieldPosition = prmfile%FieldPosition+istart+iend-1

 if( prmfile%FieldPosition .gt. slen ) then
    ! set line as procesessed
    prmfile%CurrentLine%Processed = .true.
    ! and start new one
    prmfile%FieldPosition = 1
    prmfile%CurrentLine => prmfile%CurrentLine%Next
 end if

 prmfile_get_field = .true.
 return

end function prmfile_get_field

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_get_field_on_line(prmfile,field)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 character(*)                       :: field
 !------------------------------------------------
 character(len=2)                   :: ws
 integer                            :: slen,istart,iend
 ! -----------------------------------------------------------------------------

 ! init white space characters
 ws(1:1) = char(9)    ! tabulator
 ws(2:2) = ' '        ! space

 prmfile_get_field_on_line = .false.

 if( .not. associated(prmfile%CurrentLine) ) then
    return      ! no data to process
 end if

 if( prmfile%FieldPosition .eq. 0 ) then
    prmfile%FieldPosition = 1       ! first record
 end if

 slen = len_trim(prmfile%CurrentLine%Text)

 if( prmfile%FieldPosition .gt. slen ) then
    return
 end if

 istart = verify(prmfile%CurrentLine%Text(prmfile%FieldPosition:slen),ws)
 if( istart .le. 0 ) then
    ! this cannot happen
    stop 'Fatal error in prmfile_get_field'
 end if
 iend  = scan(prmfile%CurrentLine%Text(prmfile%FieldPosition+istart-1:slen),ws)
 if( iend .le. 0 ) then
    iend = slen - prmfile%FieldPosition+istart + 1
 end if
 field = prmfile%CurrentLine%Text(prmfile%FieldPosition+istart-1:prmfile%FieldPosition+istart+iend-2)

 ! correct position
 prmfile%FieldPosition = prmfile%FieldPosition+istart+iend-1

 if( prmfile%FieldPosition .gt. slen ) then
    ! set line as procesessed
    prmfile%CurrentLine%Processed = .true.
 end if

 prmfile_get_field_on_line = .true.

 return

end function prmfile_get_field_on_line

!===============================================================================
! subroutine prmfile_set_sec_as_processed(prmfile)
! ------------------------------------------------------------------------------
! set all lines in current section as processed
!===============================================================================

subroutine prmfile_set_sec_as_processed(prmfile)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 ! --------------------------------------------
 type(LINE_TYPE), pointer           :: pline
 !-----------------------------------------------------------------------------

 nullify(pline)

 if( associated(prmfile%CurrentSection) ) then
    pline =>  prmfile%CurrentSection%FirstLine
 end if

 do while(associated(pline))
    pline%Processed = .true.
    pline => pline%Next
 end do

 return

end subroutine prmfile_set_sec_as_processed

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

character(len=3) function prmfile_onoff(value)

 implicit none
 logical      :: value
 ! ----------------------------------------------------------------------------

 if(value) then 
    prmfile_onoff = ' on'
 else
    prmfile_onoff = 'off'
 endif

end function prmfile_onoff

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module prmfile
