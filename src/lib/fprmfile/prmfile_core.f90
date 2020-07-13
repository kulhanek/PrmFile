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

module prmfile_core

use prmfile_dat

implicit none
contains

!===============================================================================
! logical function prmfile_parse(prmfile)
! ------------------------------------------------------------------------------
! it parses prmfile to individual groups, sections and lines
!===============================================================================

logical function prmfile_parse(prmfile)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 ! -----------------------------------------------
 character(len=PRMFILE_MAX_LINE)    :: line
 integer                            :: lineid, stat
 !------------------------------------------------------------------------------

 lineid = 1
 do while(.true.)
    line = ''

    ! read line --------------------
    read(unit=PRMFILE_UNIT,fmt='(a)',iostat=stat) line
    if( stat .gt. 0 ) then
        prmfile_parse = .false.
        return
    end if

    ! parse line -------------------
    if( .not. prmfile_parse_line(prmfile,lineid,line) ) then
        prmfile_parse = .false.
        return
    end if

    if( stat .lt. 0 ) then          ! EOF
        prmfile_parse = .true.
        return
    end if

    lineid = lineid + 1
 end do

 return

end function prmfile_parse

!===============================================================================
! logical function prmfile_parse(prmfile,line)
! ------------------------------------------------------------------------------
! it parses line to individual groups, sections and lines
!===============================================================================

logical function prmfile_parse_line(prmfile,lineid,line)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 integer                            :: lineid
 character(len=PRMFILE_MAX_LINE)    :: line
 ! -----------------------------------------------
 character(len=2)                   :: ws
 integer                            :: tstart,tend
 !------------------------------------------------------------------------------

 ! init white space characters
 ws(1:1) = char(9)    ! tabulator
 ws(2:2) = ' '        ! space

 tstart = verify(line,ws)

 ! skip if empty line
 if( tstart .eq. 0 ) then
    prmfile_parse_line = .true.
    return
 end if

 ! determine type of line
 select case(line(tstart:tstart))
    case('!','#','*')
        ! it is comment skip it
        prmfile_parse_line = .true.
        return
    case('{')
        ! it is group
        tend = index(line,'}')
        if(tend .le. tstart+1) then
            call add_error(prmfile,lineid,'Group name is not terminated or it is empty string!')
            prmfile_parse_line = .false.
            return
        end if
        prmfile_parse_line = prmfile_add_group(prmfile,lineid,line(tstart+1:tend-1))
        return
    case('[')
        ! it is section
        tend = index(line,']')
        if(tend .le. tstart+1) then
            call add_error(prmfile,lineid,'Section name is not terminated or it is empty string!')
            prmfile_parse_line = .false.
            return
        end if
        prmfile_parse_line = prmfile_add_section(prmfile,lineid,line(tstart+1:tend-1))
        return
    case default
        ! it is line
        prmfile_parse_line = prmfile_add_line(prmfile,lineid,line)
        return
 end select

end function prmfile_parse_line

!===============================================================================
! logical function prmfile_add_group(prmfile,lineid,name)
! ------------------------------------------------------------------------------
! it adds group with name to prmfile
!===============================================================================

logical function prmfile_add_group(prmfile,lineid,name)

 implicit none
 type(PRMFILE_TYPE)                :: prmfile
 integer                            :: lineid
 character(*)                       :: name
 ! -----------------------------------------------
 type(GROUP_TYPE), pointer          :: pgroup
 integer                            :: stat
 !------------------------------------------------------------------------------

 ! check group name length
 if( len(name) .gt. PRMFILE_MAX_GROUP_NAME ) then
    call add_error(prmfile,lineid,'Group name is too long (Longer then PRMFILE_MAX_GROUP_NAME)!')
    prmfile_add_group = .false.
    return
 end if

 ! allocate new group
 allocate(pgroup,stat=stat)
 if( stat .ne. 0 ) then
    call add_error(prmfile,lineid,'Unable to allocate memory for group!')
    prmfile_add_group = .false.
    return
 end if

 ! init it
 pgroup%Name = name
 pgroup%Processed = .false.
 nullify(pgroup%Next)
 nullify(pgroup%FirstSection)

 ! and connect it to prmfile
 if( .not. associated(prmfile%FirstGroup) ) then
    prmfile%FirstGroup  => pgroup
    prmfile%CurrentGroup => pgroup
 else
    prmfile%CurrentGroup%Next => pgroup
    prmfile%CurrentGroup => pgroup
 end if

 nullify(prmfile%CurrentSection)
 nullify(prmfile%CurrentLine)

 prmfile_add_group = .true.
 return

end function prmfile_add_group

!===============================================================================
! logical function prmfile_add_section(prmfile,lineid,name)
! ------------------------------------------------------------------------------
! it adds section with name to prmfile and currect group
! if none group exists then generic {MAIN} is created automatically
!===============================================================================

logical function prmfile_add_section(prmfile,lineid,name)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 integer                            :: lineid
 character(*)                       :: name
 ! -----------------------------------------------
 type(SECTION_TYPE), pointer        :: psection
 integer                            :: stat
 !------------------------------------------------------------------------------

 ! check group name length
 if( len(name) .gt. PRMFILE_MAX_SECTION_NAME ) then
    call add_error(prmfile,lineid,'Section name is too long (Longer then PRMFILE_MAX_SECTION_NAME)!')
    prmfile_add_section = .false.
    return
 end if

 ! add first group if it does not exist
 if( .not. associated(prmfile%FirstGroup) ) then
    if( .not. prmfile_add_group(prmfile,lineid,PRMFILE_MAIN_GROUP_NAME) ) then
        prmfile_add_section = .false.
        return
    end if
 end if

 ! allocate new section
 allocate(psection,stat=stat)
 if( stat .ne. 0 ) then
    call add_error(prmfile,lineid,'Unable to allocate memory for section!')
    prmfile_add_section = .false.
    return
 end if

 ! init section
 psection%Name = name
 psection%Processed = .false.
 nullify(psection%Next)
 nullify(psection%FirstLine)

 nullify(prmfile%CurrentLine)

 ! and connect it to current group
 if( associated(prmfile%CurrentSection) ) then
    prmfile%CurrentSection%Next => psection
    prmfile%CurrentSection => psection
 else
    prmfile%CurrentGroup%FirstSection => psection
    prmfile%CurrentSection => psection
 end if

 prmfile_add_section = .true.
 return

end function prmfile_add_section

!===============================================================================
! logical function prmfile_add_line(prmfile,lineid,name)
! ------------------------------------------------------------------------------
! it adds line to prmfile and currect section
! if none section exists then en error occures
!===============================================================================

logical function prmfile_add_line(prmfile,lineid,line)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 integer                            :: lineid
 character(*)                       :: line
 ! -----------------------------------------------
 type(LINE_TYPE), pointer           :: pline
 character(len=2)                   :: ws
 character(len=3)                   :: comm
 character(len=2)                   :: quot
 character                          :: fquot
 integer                            :: tstart,tend,pos,stat,qpos
 !------------------------------------------------------------------------------

 ! find the first and last character of line

 ! init white space characters
 ws(1:1)   = char(9)    ! tabulator
 ws(2:2)   = ' '        ! space
 comm(1:1) = '!'
 comm(2:2) = '*'
 comm(3:3) = '#'
 quot(1:1) = '"'
 quot(2:2) = ''''

 ! beggining
 tstart = verify(line,ws)

 ! try to detect beginning of comment
 pos = tstart
 tend = tstart
 do while(pos .le. len(line))
    qpos = scan(line(pos:),quot)
    if( qpos .ne. 0 ) then
        fquot = line(qpos:qpos)
        ! quotation found move after it
        pos = pos + qpos
        if( pos .gt. len(line) ) then
            tend = len(line)
            exit
        end if
        tend = pos
        ! find end
        qpos = scan(line(pos:),fquot)
        if( qpos .ne. 0 ) then
            pos = pos + qpos
            tend = pos
        else
            ! not found take everything to the end of line
            pos = len(line)
            tend = pos
            exit
        end if
    end if
    pos = scan(line(pos:),comm)
    if( pos .eq. 0 ) then
        ! it is not found
        tend = len(line)
        exit
    end if
    ! tend - points to the leading character of comment
    tend = tend + pos - 1
    ! check if it is real comment
    ! it is the comment only if the preceeding character is a white space
    if( scan(line(tend-1:tend-1),ws) .ne. 0 ) then
        tend = tend - 1
        exit
    end if
    ! correct tend
    tend = tend + 1
    pos = tend
 end do

 ! now remove white spaces between line end and the remaining part
 tend = verify(line(1:tend),ws,back=.true.)

 ! check if there is active section
 if( .not. associated(prmfile%CurrentSection) ) then
    call add_error(prmfile,lineid,'None section is active to add this line!')
    prmfile_add_line = .false.
    return
 end if

 ! allocate new line
 allocate(pline,stat=stat)
 if( stat .ne. 0 ) then
    call add_error(prmfile,lineid,'Unable to allocate memory for line!')
    prmfile_add_line = .false.
    return
 end if

 ! init section
 pline%Text = line(tstart:tend)
 pline%Processed = .false.
 nullify(pline%Next)

 ! and connect it to current section
 if( associated(prmfile%CurrentLine) ) then
    prmfile%CurrentLine%Next => pline
    prmfile%CurrentLine => pline
 else
    prmfile%CurrentSection%FirstLine => pline
    prmfile%CurrentLine => pline
 end if

 prmfile_add_line = .true.
 return

end function prmfile_add_line

!===============================================================================
! subroutine add_error(prmfile,lineid,error)
! ------------------------------------------------------------------------------
! add error message to prmfile error buffer
!===============================================================================

subroutine add_error(prmfile,lineid,error)

 implicit none
 type(PRMFILE_TYPE)                 :: prmfile
 integer                            :: lineid
 character(*)                       :: error
 !------------------------------------------------------------------------------

 write(*,'(A,A)')  '>>> ERROR: ',trim(error)
 write(*,'(A,A)')  '            Name : ',trim(prmfile%name)
 write(*,'(A,I6)') '            Line : ',lineid
 stop
 return

end subroutine add_error

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_find_key(prmfile,key,value)

 implicit none
 type(PRMFILE_TYPE)                     :: prmfile
 character(*), intent(in)               :: key
 character(*), intent(out)              :: value
 !------------------------------------------------
 character(len=PRMFILE_MAX_KEY_NAME)    :: ikey
 character(len=PRMFILE_MAX_KEY_NAME)    :: lkey
 character(len=PRMFILE_MAX_VALUE)       :: lvalue
 ! -----------------------------------------------------------------------------

 ! go to beginning of section --------------------
 nullify(prmfile%CurrentLine)
 prmfile%FieldPosition = 0

 if(associated(prmfile%CurrentSection)) then
    prmfile%CurrentLine => prmfile%CurrentSection%FirstLine
 end if

 ikey = key
 call prmfile_upcase(ikey)

 do while(associated(prmfile%CurrentLine))

    if( .not. prmfile_split_line(prmfile%CurrentLine%Text,lkey,lvalue) ) then
        prmfile_find_key = .false.
        return
    end if

    call prmfile_upcase(lkey)

    ! and compare
    if( trim(ikey) .eq. trim(lkey) ) then
        prmfile_find_key = .true.
        value = lvalue
        return
    end if

    prmfile%CurrentLine => prmfile%CurrentLine%Next
 end do
 
 prmfile_find_key = .false.
 return

end function prmfile_find_key

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine prmfile_set_kline_as_processed(prmfile)

 implicit none
 type(PRMFILE_TYPE)                     :: prmfile
 ! -----------------------------------------------------------------------------

 if(associated(prmfile%CurrentLine)) then
    prmfile%CurrentLine%Processed = .true.
 end if

 return

end subroutine prmfile_set_kline_as_processed

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine prmfile_upcase(string)

 implicit none
 character(*), intent(inout)    :: string
 !------------------------------------------------
 integer                        :: i,c
 ! ---------------------------------------------------------------------------

 do i=1, len_trim(string)
    c = ichar(string(i:i))
    if(c .ge. ichar('a') .and. c .le. ichar('z') ) c = c - ichar('a') + ichar('A')
    string(i:i) = char(c)
 end do

end subroutine prmfile_upcase

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function prmfile_split_line(line,key,value)

 implicit none
 character(*), intent(in)   :: line
 character(*), intent(out)  :: key
 character(*), intent(out)  :: value
 !------------------------------------------------
 character(len=2)           :: ws
 integer                    :: skey,ekey,svalue,llen
 ! -----------------------------------------------------------------------------

 ! init white space characters
 ws(1:1) = char(9)    ! tabulator
 ws(2:2) = ' '        ! space

 prmfile_split_line = .false.
 key = ''
 value = ''

 llen = len(line)

 ! key starts always on 1 (due to parsing of file)
 skey = 1
 ! find end of key
 ekey = scan(line,ws)
 if( ekey .le. 0 ) return       ! line is only key, this is not permitted

 ! extract key
 key = line(skey:ekey-1)

 ! find beggining of value
 svalue = verify(line(ekey:llen),ws)

 ! value cannot be only empty string because of parsing
 value = line(ekey+svalue-1:llen)

 prmfile_split_line = .true.

 return

end function prmfile_split_line

!===============================================================================
! integer function prmfile_count_ulines_in_grp(pgroup)
! ------------------------------------------------------------------------------
! it returns total number of unprocessed lines in given group
!===============================================================================

integer function prmfile_count_ulines_in_grp(igroup)

 implicit none
 type(GROUP_TYPE), pointer              :: igroup
 ! -----------------------------------------------
 type(GROUP_TYPE), pointer              :: pgroup
 type(SECTION_TYPE), pointer            :: psection
 type(LINE_TYPE), pointer               :: pline
 !------------------------------------------------------------------------------

 prmfile_count_ulines_in_grp = 0

 pgroup => igroup

 if(associated(pgroup)) then
    ! loop over sections
    psection => pgroup%FirstSection
    do while(associated(psection))
        ! loop over lines
        pline => psection%FirstLine
        do while(associated(pline))
            if( .not. pline%Processed ) then
                prmfile_count_ulines_in_grp = prmfile_count_ulines_in_grp + 1
            end if
            pline => pline%Next
        end do
        psection => psection%Next
    end do
 end if

 return

end function prmfile_count_ulines_in_grp

!===============================================================================
! integer function prmfile_count_ulines_in_sec(pgroup)
! ------------------------------------------------------------------------------
! it returns total number of unprocessed lines in given section
!===============================================================================

integer function prmfile_count_ulines_in_sec(isection)

 implicit none
 type(SECTION_TYPE), pointer            :: isection
 ! -----------------------------------------------
 type(SECTION_TYPE), pointer            :: psection
 type(LINE_TYPE), pointer               :: pline
 !------------------------------------------------------------------------------

 prmfile_count_ulines_in_sec = 0

 psection => isection

 if(associated(psection)) then
    ! loop over lines
    pline => psection%FirstLine
    do while(associated(pline))
        if( .not. pline%Processed ) then
            prmfile_count_ulines_in_sec = prmfile_count_ulines_in_sec + 1
        end if
        pline => pline%Next
    end do
 end if

 return

end function prmfile_count_ulines_in_sec

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module prmfile_core
