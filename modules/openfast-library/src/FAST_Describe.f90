!**********************************************************************************************************************************
! FAST_Solver.f90, FAST_Subs.f90, FAST_Lin.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
! FAST_Prog.f90, FAST_Library.f90, FAST_Prog.c are different drivers for this code.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of FAST.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************
MODULE FAST_Describe
    
    USE NWTC_Library, ONLY: ReKi, R8Ki, SiKi, ProgDesc, BITS_IN_ADDR, num2lstr, getnvd
    Use VersionInfo
     
   IMPLICIT NONE

CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------
!> This function returns a string describing the glue code and some of the compilation options we're using.
FUNCTION GetVersion(ThisProgVer)

   ! Passed Variables:

   TYPE(ProgDesc), INTENT( IN    ) :: ThisProgVer     !< program name/date/version description
   CHARACTER(1024)                 :: GetVersion      !< String containing a description of the compiled precision.
   
   CHARACTER(200)                  :: git_commit

#ifdef COMPILE_SIMULINK
   LOGICAL, PARAMETER :: Cmpl4SFun = .TRUE.                            ! Is the module being compiled as an S-Function for Simulink?
#else
   LOGICAL, PARAMETER :: Cmpl4SFun = .FALSE.                           ! Is the module being compiled as an S-Function for Simulink?
#endif

#ifdef COMPILE_LABVIEW
   LOGICAL, PARAMETER :: Cmpl4LV = .TRUE.                            ! Is the module being compiled for Labview?
#else
   LOGICAL, PARAMETER :: Cmpl4LV = .FALSE.                           ! Is the module being compiled for Labview?
#endif

   GetVersion = TRIM(GetNVD(ThisProgVer))//', compiled'

   IF ( Cmpl4SFun )  THEN     ! FAST has been compiled as an S-Function for Simulink
      GetVersion = TRIM(GetVersion)//' as a DLL S-Function for Simulink'
   ELSEIF ( Cmpl4LV )  THEN     ! FAST has been compiled as a DLL for Labview
      GetVersion = TRIM(GetVersion)//' as a DLL for LabVIEW'
   ENDIF   
   
   GetVersion = TRIM(GetVersion)//' as a '//TRIM(Num2LStr(BITS_IN_ADDR))//'-bit application using'
   
   ! determine precision

      IF ( ReKi == SiKi )  THEN     ! Single precision
         GetVersion = TRIM(GetVersion)//' single'
      ELSEIF ( ReKi == R8Ki )  THEN ! Double precision
         GetVersion = TRIM(GetVersion)// ' double'
      ELSE                          ! Unknown precision
         GetVersion = TRIM(GetVersion)//' unknown'
      ENDIF
      

!   GetVersion = TRIM(GetVersion)//' precision with '//OS_Desc
   GetVersion = TRIM(GetVersion)//' precision'

   ! add git info
   git_commit = QueryGitVersion()
   GetVersion = TRIM(GetVersion)//' at commit '//git_commit

   RETURN
END FUNCTION GetVersion

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine parses and compiles the relevant version and compile data for a givne program
subroutine GetProgramMetadata(ThisProgVer, name, version, git_commit, architecture, precision)

   TYPE(ProgDesc), INTENT(IN ) :: ThisProgVer     !< program name/date/version description
   character(200), intent(out) :: name, version
   character(200), intent(out) :: git_commit, architecture, precision
   
   name = trim(ThisProgVer%Name)
   version = trim(ThisProgVer%Ver)
   
   git_commit = QueryGitVersion()

   architecture = TRIM(Num2LStr(BITS_IN_ADDR))//' bit'
   
   if (ReKi == SiKi) then
     precision = 'single'
   else if (ReKi == R8Ki) then
     precision = 'double'
   else
     precision = 'unknown'
   end if
   
end subroutine GetProgramMetadata

end module