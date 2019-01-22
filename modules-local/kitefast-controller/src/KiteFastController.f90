!**********************************************************************************************************************************
!> ## SC
!! The KiteFastController  module implements a controller for the KiteFAST code. 
!! KiteFastController_Types will be auto-generated by the FAST registry program, based on the variables specified in the
!! KiteFastController_Registry.txt file.
!!
! ..................................................................................................................................
!! ## LICENSING 
!! Copyright (C) 2018  National Renewable Energy Laboratory
!!
!!    This file is part of KiteFAST.
!!
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!! You may obtain a copy of the License at
!!
!!     http://www.apache.org/licenses/LICENSE-2.0
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!**********************************************************************************************************************************
module KiteFastController

   use KiteFastController_Types
   use NWTC_Library
   use, intrinsic :: ISO_C_Binding
   
   implicit none
   private
   
   integer,        parameter    :: IntfStrLen  = 1025       ! length of strings through the C interface
   type(ProgDesc), parameter    :: KFC_Ver = ProgDesc( 'KiteFastController', '', '' )

      !> Definition of the DLL Interface for the SuperController
      !! 
   abstract interface
      subroutine KFC_DLL_Init_PROC ( errStat, errMsg )  BIND(C)
         use, intrinsic :: ISO_C_Binding
         integer(C_INT),         intent(  out) :: errStat             !< error status code (uses NWTC_Library error codes)
         character(kind=C_CHAR), intent(inout) :: errMsg          (*) !< Error Message from DLL to simulation code        
      end subroutine KFC_DLL_Init_PROC   
   end interface   

   abstract interface
      subroutine KFC_DLL_Step_PROC ( dcm_g2b_c, pqr_c, acc_norm_c, Xg_c, Vg_c, Vb_c, Ag_c, Ab_c, rho_c, apparent_wind_c, &
         tether_forceb_c, wind_g_c, kFlapA_c, Motor_c, errStat, errMsg )  BIND(C)
         use, intrinsic :: ISO_C_Binding
         real(C_DOUBLE),         intent(in   ) :: dcm_g2b_c(9)      
         real(C_DOUBLE),         intent(in   ) :: pqr_c(3)          
         real(C_DOUBLE),         intent(in   ) :: acc_norm_c    
         real(C_DOUBLE),         intent(in   ) :: Xg_c(3)           
         real(C_DOUBLE),         intent(in   ) :: Vg_c(3)           
         real(C_DOUBLE),         intent(in   ) :: Vb_c(3)           
         real(C_DOUBLE),         intent(in   ) :: Ag_c(3)           
         real(C_DOUBLE),         intent(in   ) :: Ab_c(3)           
         real(C_DOUBLE),         intent(in   ) :: rho_c          
         real(C_DOUBLE),         intent(in   ) :: apparent_wind_c(3)
         real(C_DOUBLE),         intent(in   ) :: tether_forceb_c(3) 
         real(C_DOUBLE),         intent(in   ) :: wind_g_c(3) 
         real(C_DOUBLE),         intent(  out) :: kFlapA_c(10)                      
         real(C_DOUBLE),         intent(  out) :: Motor_c(8)
         integer(C_INT),         intent(inout) :: errStat           !< error status code (uses NWTC_Library error codes)
         character(kind=C_CHAR), intent(inout) :: errMsg(1025)      !< Error Message from DLL to simulation code        
      end subroutine KFC_DLL_Step_PROC   
   end interface   
 
   abstract interface
      subroutine KFC_DLL_END_PROC ( errStat, errMsg )  BIND(C)
         use, intrinsic :: ISO_C_Binding
         integer(C_INT),         intent(  out) :: errStat             !< error status code (uses NWTC_Library error codes)
         character(kind=C_CHAR), intent(inout) :: errMsg          (*) !< Error Message from DLL to simulation code        
      end subroutine KFC_DLL_END_PROC   
   end interface   

   public :: KFC_Init                     ! Initialization routine
   public :: KFC_End                      ! Ending routine (includes clean up)
   public :: KFC_Step                     ! Routine for computing outputs and internally updating states
  
   contains   
   
   
   subroutine KFC_End(p, errStat, errMsg)

      type(KFC_ParameterType),        intent(inout)  :: p               !< Parameters
      integer(IntKi),                 intent(  out)  :: errStat         !< Error status of the operation
      character(*),                   intent(  out)  :: errMsg          !< Error message if ErrStat /= ErrID_None

         ! local variables
      character(*), parameter                        :: routineName = 'KFC_End'
      integer(IntKi)                                 :: errStat2       ! The error status code
      character(ErrMsgLen)                           :: errMsg2        ! The error message, if an error occurred
      procedure(KFC_DLL_END_PROC),pointer            :: DLL_KFC_End_Subroutine       ! The address of the controller cc_end procedure in the DLL
      character(kind=C_CHAR)                         :: errMsg_c(IntfStrLen)
      errStat = ErrID_None
      errMsg= ''
      
      if (.not. p%useDummy) then
            ! Call the DLL's end subroutine:
         call C_F_PROCPOINTER( p%DLL_Trgt%ProcAddr(3), DLL_KFC_End_Subroutine) 
         call DLL_KFC_End_Subroutine ( errStat, errMsg_c ) 
         call c_to_fortran_string(errMsg_c, errMsg)
      
         print *, " KFC_End errStat - ", errStat, " errMsg - ", trim(errMsg)

         ! TODO: Check errors
      
            ! Free the library
         call FreeDynamicLib( p%DLL_Trgt, errStat2, errMsg2 )  
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )
      end if
      
   end subroutine KFC_End

   subroutine KFC_Init(InitInp, u, p, y, interval, InitOut, errStat, errMsg )

      type(KFC_InitInputType),      intent(in   )  :: InitInp     !< Input data for initialization routine
      type(KFC_InputType),          intent(inout)  :: u           !< An initial guess for the input
      type(KFC_ParameterType),      intent(  out)  :: p           !< Parameters
      type(KFC_OutputType),         intent(  out)  :: y           !< Initial system outputs 
      real(DbKi),                   intent(inout)  :: interval    !< Coupling interval in seconds: 
                                                                  !<   Input is the timestep size requested by caller, returned is the Controller's required timestep
      type(KFC_InitOutputType),     intent(  out)  :: InitOut     !< Output for initialization routine
      integer(IntKi),               intent(  out)  :: errStat     !< Error status of the operation
      character(*),                 intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
 
         ! local variables
      character(*), parameter                 :: routineName = 'KFC_Init'
      integer(IntKi)                          :: errStat2                     ! The error status code
      character(ErrMsgLen)                    :: errMsg2                      ! The error message, if an error occurred
      procedure(KFC_DLL_Init_PROC),pointer    :: DLL_KFC_Init_Subroutine       ! The address of the controller cc_init procedure in the DLL
      
      integer(IntKi)                          :: nParams
      character(kind=C_CHAR)                  :: errMsg_c(IntfStrLen)
      
      errStat2 = ErrID_None
      errMsg2  = ''
   
      call DispNVD( KFC_Ver )  ! Display the version of this interface
      
         ! Check that key Kite model components match the requirements of this controller interface.
      !=============================================================================================
      ! NOTE: GJH: Perhaps a better design is to let the actual controller (shared object) determine if numFlaps and numPylons and interval are acceptable
      !
      !if (InitInp%numFlaps /= 3) call SetErrStat( ErrID_Fatal, 'The current KiteFAST controller interface requires numFlaps = 3', errStat, errMsg, routineName )
      !if (InitInp%numPylons /= 2) call SetErrStat( ErrID_Fatal, 'The current KiteFAST controller interface requires numPylons = 2', errStat, errMsg, routineName )
      !if (.not. EqualRealNos(interval, 0.01_DbKi)) call SetErrStat( ErrID_Fatal, 'The current KiteFAST controller interface requires DT = 0.01 seconds', errStat, errMsg, routineName )
      !   if (errStat >= AbortErrLev ) return
      !=============================================================================================  
      
               ! Set the module's parameters
      p%numFlaps  = InitInp%numFlaps
      p%numPylons = InitInp%numPylons
      p%DT        = interval
      p%useDummy  = InitInp%useDummy
      
         ! allocate the inputs and outputs
      call AllocAry( u%SPyAeroTorque, 2, p%numPylons, 'u%SPyAeroTorque', errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )
      call AllocAry( u%PPyAeroTorque, 2, p%numPylons, 'u%PPyAeroTorque', errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )
      call AllocAry( y%SPyGenTorque,  2, p%numPylons, 'y%SPyGenTorque',  errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )
      call AllocAry( y%PPyGenTorque,  2, p%numPylons, 'y%PPyGenTorque',  errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )   
      call AllocAry( y%SPyRtrSpd,     2, p%numPylons, 'y%SPyRtrSpd',     errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )
      call AllocAry( y%PPyRtrSpd,     2, p%numPylons, 'y%PPyRtrSpd',     errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )  
       call AllocAry( y%SPyRtrAcc,     2, p%numPylons, 'y%SPyRtrAcc',     errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )
      call AllocAry( y%PPyRtrAcc,     2, p%numPylons, 'y%PPyRtrAcc',     errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )  
      call AllocAry( y%SPyBldPitch,   2, p%numPylons, 'y%SPyBldPitch',   errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )
      call AllocAry( y%PPyBldPitch,   2, p%numPylons, 'y%PPyBldPitch',   errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )   
      call AllocAry( y%SFlp, p%numFlaps, 'y%SFlp', errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )
      call AllocAry( y%PFlp, p%numFlaps, 'y%PFlp', errStat2, errMsg2 )   
         
      if (errStat >= AbortErrLev ) return

      
      ! Are we simply using a dummy controller?  If so, we will skip trying to call into a DLL/SO      
      
      if (.not. p%useDummy) then
     
            ! Define and load the DLL:
         p%DLL_Trgt%FileName = InitInp%DLL_FileName

         p%DLL_Trgt%ProcName = "" ! initialize all procedures to empty so we try to load only one
         p%DLL_Trgt%ProcName(1) = 'kfc_dll_init'
         p%DLL_Trgt%ProcName(2) = 'kfc_dll_step'
         p%DLL_Trgt%ProcName(3) = 'kfc_dll_end'
      
         call LoadDynamicLib ( p%DLL_Trgt, errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, routineName )
         if (errStat >= AbortErrLev ) return

         ! Now that the library is loaded, call the controller's kfc_dll_init routine

            ! Call the DLL (first associate the address from the procedure in the DLL with the subroutine):
         call C_F_PROCPOINTER( p%DLL_Trgt%ProcAddr(1), DLL_KFC_Init_Subroutine) 
      
   ! TODO: jjonkman's plan doc assumes that the initial outputs are returned by KFC_Init(), but we aren't doing that here.  GJH 12/19/18
            ! Can we modify the following to send the controller numFlaps, numPylons, and interval and let the controller throw an error and/or change interval as needed? GJH 12/19/18
         ! also add Irot for each rotor.
         call DLL_KFC_Init_Subroutine ( errStat, errMsg_c ) 
      
         call c_to_fortran_string(errMsg_c, errMsg)
         print *, " KFC_Init errStat - ", errStat, " errMsg - ", trim(errMsg)
         ! TODO: Check errors
         print *, " debug marker - pre errStat >= Abort"
         if (errStat >= AbortErrLev ) return
         print *, " debug marker - post errStat >= Abort"
         
         ! TODO: obtain initial outputs from the DLL and set them
           ! Set outputs to zero for now
         y%SPyGenTorque = 0.0_ReKi
         y%PPyGenTorque = 0.0_ReKi
         y%SPyRtrSpd    = 0.0_ReKi
         y%PPyRtrSpd    = 0.0_ReKi
         y%SPyRtrAcc    =   0.0_ReKi  ! rad/s^2
         y%PPyRtrAcc    =   0.0_ReKi  ! rad/s^2
         y%SPyBldPitch  = 0.0_ReKi
         y%PPyBldPitch  = 0.0_ReKi
         y%SFlp         = 0.0_ReKi
         y%PFlp         = 0.0_ReKi
         y%Rudr         = 0.0_ReKi
         y%SElv         = 0.0_ReKi
         y%PElv         = 0.0_ReKi
        
      else
           ! Set outputs to zero except for RtrSpd which is set to be constant for the dummy controller
         y%SPyGenTorque = 0.0_ReKi
         y%PPyGenTorque = 0.0_ReKi
            ! TODO: Determine what would be a realistic dummy set of speed and the correct signs for each rotor
         y%SPyRtrSpd    = 180.0_ReKi  ! rad/s
         y%PPyRtrSpd    = 180.0_ReKi  ! rad/s
         y%SPyRtrAcc    =   0.0_ReKi  ! rad/s^2
         y%PPyRtrAcc    =   0.0_ReKi  ! rad/s^2
         y%SPyBldPitch  = 0.0_ReKi
         y%PPyBldPitch  = 0.0_ReKi
         y%SFlp         = 0.0_ReKi
         y%PFlp         = 0.0_ReKi
         y%Rudr         = 0.0_ReKi
         y%SElv         = 0.0_ReKi
         y%PElv         = 0.0_ReKi
         
      end if
      
      
      
      
   end subroutine KFC_Init

   subroutine KFC_Step(t, u, p, y, errStat, errMsg )
      real(DbKi),                    intent(in   )  :: t           !< Current simulation time in seconds
      type(KFC_InputType),           intent(in   )  :: u           !< Inputs at Time t
      type(KFC_ParameterType),       intent(in   )  :: p           !< Parameters
      type(KFC_OutputType),          intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                   !!   nectivity information does not have to be recalculated)
      integer(IntKi),                intent(  out)  :: errStat     !< Error status of the operation
      character(*),                  intent(  out)  :: errMsg      !< Error message if ErrStat /= ErrID_None
   
      character(*), parameter                       :: routineName = 'KFC_Step'
      integer(IntKi)                                :: errStat2       ! The error status code
      character(ErrMsgLen)                          :: errMsg2        ! The error message, if an error occurred     
      procedure(KFC_DLL_Step_PROC),pointer          :: DLL_KFC_Step_Subroutine              ! The address of the supercontroller sc_calcoutputs procedure in the DLL
      real(C_DOUBLE)                                :: dcm_g2b_c(9)      
      real(C_DOUBLE)                                :: pqr_c(3)          
      real(C_DOUBLE)                                :: acc_norm_c    
      real(C_DOUBLE)                                :: Xg_c(3)           
      real(C_DOUBLE)                                :: Vg_c(3)           
      real(C_DOUBLE)                                :: Vb_c(3)           
      real(C_DOUBLE)                                :: Ag_c(3)           
      real(C_DOUBLE)                                :: Ab_c(3)           
      real(C_DOUBLE)                                :: rho_c          
      real(C_DOUBLE)                                :: apparent_wind_c(3)
      real(C_DOUBLE)                                :: tether_forceb_c(3) 
      real(C_DOUBLE)                                :: wind_g_c(3)       
      character(kind=C_CHAR)                        :: errMsg_c(IntfStrLen)
      real(C_DOUBLE)                                :: kFlapA_c(10)                   
      real(C_DOUBLE)                                :: Motor_c(8)

      errStat2 = ErrID_None
      errMsg2  = ''
      
         ! Cast and massage inputs to match DLL datatypes
      dcm_g2b_c       = reshape(u%dcm_g2b,(/9/))
      pqr_c           = u%pqr
      acc_norm_c      = u%acc_norm
      Xg_c            = u%Xg
      Vg_c            = u%Vg
      Vb_c            = u%Vb
      Ag_c            = u%Ag
      Ab_c            = u%Ab
      rho_c           = u%rho
      apparent_wind_c = u%apparent_wind
      tether_forceb_c = u%tether_forceb
      wind_g_c        = u%wind_g
      
      
      if (.not. p%useDummy) then
            ! Call the DLL (first associate the address from the procedure in the DLL with the subroutine):
         call C_F_PROCPOINTER( p%DLL_Trgt%ProcAddr(2), DLL_KFC_Step_Subroutine) 
         call DLL_KFC_Step_Subroutine ( dcm_g2b_c, pqr_c, acc_norm_c, Xg_c, Vg_c, Vb_c, Ag_c, Ab_c, rho_c, apparent_wind_c, tether_forceb_c, wind_g_c, kFlapA_c, Motor_c, errStat, errMsg_c ) 
         call c_to_fortran_string(errMsg_c, errMsg)

         print *, " KFC_Step errStat - ", errStat, " errMsg - ", trim(errMsg)

            ! Convert the controller outputs into the KiteFAST Fortran-style controller outputs
         y%SFlp(1) = kFlapA_c(5)
         y%SFlp(2) = kFlapA_c(7)
         y%SFlp(3) = kFlapA_c(8)
         y%PFlp(1) = kFlapA_c(4)
         y%PFlp(2) = kFlapA_c(2)
         y%PFlp(3) = kFlapA_c(1)
         y%Rudr(:) = kFlapA_c(10)
         y%SElv(:) = kFlapA_c(9)
         y%PElv(:) = kFlapA_c(9)   
      
         y%SPyGenTorque(1,1) = Motor_c(7)  ! starboard top rotor, pylon 1 (inboard)
         y%SPyGenTorque(2,1) = Motor_c(2)  ! starboard bottom rotor, pylon 1 (inboard)
         y%SPyGenTorque(1,2) = Motor_c(8)  ! starboard top rotor, pylon 2 (outboard)
         y%SPyGenTorque(2,2) = Motor_c(1)  ! starboard bottom rotor, pylon 2 (outboard)
         y%PPyGenTorque(1,1) = Motor_c(6)  ! port top rotor, pylon 1 (inboard)
         y%PPyGenTorque(2,1) = Motor_c(3)  ! port bottom rotor, pylon 1 (inboard)
         y%PPyGenTorque(1,2) = Motor_c(5)  ! port top rotor, pylon 2 (outboard)
         y%PPyGenTorque(2,2) = Motor_c(4)  ! port bottom rotor, pylon 2 (outboard)
 
   !TODO: How to we obtain rotor speeds?
         y%SPyRtrSpd(1,1) = 0.0  !RtrSpd_c(7)  ! starboard top rotor, pylon 1 (inboard)
         y%SPyRtrSpd(2,1) = 0.0  !RtrSpd_c(2)  ! starboard bottom rotor, pylon 1 (inboard)
         y%SPyRtrSpd(1,2) = 0.0  !RtrSpd_c(8)  ! starboard top rotor, pylon 2 (outboard)
         y%SPyRtrSpd(2,2) = 0.0  !RtrSpd_c(1)  ! starboard bottom rotor, pylon 2 (outboard)
         y%PPyRtrSpd(1,1) = 0.0  !RtrSpd_c(6)  ! port top rotor, pylon 1 (inboard)
         y%PPyRtrSpd(2,1) = 0.0  !RtrSpd_c(3)  ! port bottom rotor, pylon 1 (inboard)
         y%PPyRtrSpd(1,2) = 0.0  !RtrSpd_c(5)  ! port top rotor, pylon 2 (outboard)
         y%PPyRtrSpd(2,2) = 0.0  !RtrSpd_c(4)  ! port bottom rotor, pylon 2 (outboard)

   ! TODO: Are we still receiving rotor accelerations from controller?
         y%SPyRtrAcc(1,1) = 0.0_ReKi ! starboard top rotor, pylon 1 (inboard)
         y%SPyRtrAcc(2,1) = 0.0_ReKi ! starboard bottom rotor, pylon 1 (inboard)
         y%SPyRtrAcc(1,2) = 0.0_ReKi ! starboard top rotor, pylon 2 (outboard)
         y%SPyRtrAcc(2,2) = 0.0_ReKi ! starboard bottom rotor, pylon 2 (outboard)
         y%PPyRtrAcc(1,1) = 0.0_ReKi ! port top rotor, pylon 1 (inboard)
         y%PPyRtrAcc(2,1) = 0.0_ReKi ! port bottom rotor, pylon 1 (inboard)
         y%PPyRtrAcc(1,2) = 0.0_ReKi ! port top rotor, pylon 2 (outboard)
         y%PPyRtrAcc(2,2) = 0.0_ReKi ! port bottom rotor, pylon 2 (outboard)

         ! Currently blade pitch is not being set by controller and was initialized to 0.0
         y%SPyBldPitch  = 0.0_ReKi
         y%PPyBldPitch  = 0.0_ReKi
            ! TODO Error checking
      else
         
         ! TODO: Determine what would be a realistic dummy set of speed and the correct signs for each rotor
            ! NOTE: Speed should match the settings used in the Init routine.
         y%SPyRtrSpd(1,:) = 180.0  ! starboard top rotor, all pylons 
         y%SPyRtrSpd(2,:) = 180.0  ! starboard bottom rotor, all pylons
         y%PPyRtrSpd(1,:) = 180.0  ! port top rotor, all pylons
         y%PPyRtrSpd(2,:) = 180.0  ! port bottom rotor, all pylons
         
            ! Zero rotor acceleration
         y%SPyRtrAcc(1,:) = 0.0_ReKi ! starboard top rotor, all pylons 
         y%SPyRtrAcc(2,:) = 0.0_ReKi ! starboard bottom rotor, all pylons
         y%PPyRtrAcc(1,:) = 0.0_ReKi ! port top rotor, all pylons
         y%PPyRtrAcc(2,:) = 0.0_ReKi ! port bottom rotor, all pylons
         

            ! Currently blade pitch is not being set by controller and was initialized to 0.0
         y%SPyBldPitch  = 0.0_ReKi
         y%PPyBldPitch  = 0.0_ReKi
         
            ! Set GenTorque = -AeroTorque
         y%SPyGenTorque = -u%SPyAeroTorque  ! starboard rotors
         y%PPyGenTorque = -u%PPyAeroTorque  ! port rotors
        
            ! All flag commands are constant for the dummy controller
         y%SFlp = 0.0_ReKi
         y%PFlp = 0.0_ReKi
         y%Rudr = 0.0_ReKi
         y%SElv = 0.0_ReKi
         y%PElv = 0.0_ReKi
         
      end if
      
   end subroutine KFC_Step
   
   subroutine c_to_fortran_string(input, output)
      character(kind=C_CHAR), intent(in) :: input(IntfStrLen)
      character(*), intent(out) :: output
      character(1024) :: temp_string
      temp_string = transfer(input(1:1024), output)
      call RemoveNullChar(temp_string)
      output = trim(temp_string)
   end subroutine

end module KiteFastController
