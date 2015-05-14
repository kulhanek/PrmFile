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

module prmfile_dat

implicit none

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

! these constants determine maximal limits

integer,parameter       :: PRMFILE_UNIT                = 3245
integer,parameter       :: PRMFILE_MAX_LINE            = 65534
integer,parameter       :: PRMFILE_MAX_PATH            = 255
integer,parameter       :: PRMFILE_MAX_GROUP_NAME      = 80
integer,parameter       :: PRMFILE_MAX_SECTION_NAME    = 80
integer,parameter       :: PRMFILE_MAX_KEY_NAME        = 80
integer,parameter       :: PRMFILE_MAX_VALUE           = PRMFILE_MAX_LINE
character(*),parameter  :: PRMFILE_MAIN_GROUP_NAME     = 'MAIN'

!-------------------------------------------------------------------------------
! NULL() is necessary for -auto option of compiler

! prmfile main class used by every prmfile subroutine

type PRMFILE_TYPE
    character(len=PRMFILE_MAX_PATH)         :: Name                         ! name of file
    type(GROUP_TYPE), pointer               :: FirstGroup                   ! first group in file
    type(GROUP_TYPE), pointer               :: CurrentGroup                 ! current group
    type(SECTION_TYPE), pointer             :: CurrentSection               ! current section
    type(LINE_TYPE), pointer                :: CurrentLine                  ! current line
    integer                                 :: FieldPosition                ! used by prmfile_get_field
end type PRMFILE_TYPE

!-------------------------------------------------------------------------------

type GROUP_TYPE
    character(len=PRMFILE_MAX_GROUP_NAME)   :: Name                         ! name of group
    type(GROUP_TYPE), pointer               :: Next                         ! next group
    type(SECTION_TYPE), pointer             :: FirstSection                 ! first section
    logical                                 :: Processed
end type GROUP_TYPE

!-------------------------------------------------------------------------------

type SECTION_TYPE
    character(len=PRMFILE_MAX_GROUP_NAME)   :: Name                         ! name of section
    type(SECTION_TYPE), pointer             :: Next                         ! next section
    type(LINE_TYPE), pointer                :: FirstLine                    ! first line
    logical                                 :: Processed
end type SECTION_TYPE

!-------------------------------------------------------------------------------

type LINE_TYPE
    character(len=PRMFILE_MAX_LINE)         :: Text                         ! text of line
    logical                                 :: Processed                    ! processed status
    type(LINE_TYPE), pointer                :: Next                         ! next line in section
end type LINE_TYPE

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module prmfile_dat
