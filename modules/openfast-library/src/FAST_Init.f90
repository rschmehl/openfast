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
MODULE FAST_Initialization
   
   use FAST_Describe
   USE FAST_ModTypes
   USE FAST_IO
   USE FAST_Linear, ONLY: Init_Lin
   USE FAST_Solver,       ONLY: initmodulemappings
   USE AeroDyn,           ONLY: AD_Init
   USE AeroDyn14,         ONLY: AD14_Init
   USE InflowWind,        ONLY: InflowWind_Init
   USE ElastoDyn,         ONLY: ED_Init
   USE BeamDyn,           ONLY: BD_Init
   USE FEAMooring,        ONLY: FEAM_Init
   USE MoorDyn,           ONLY: MD_Init
   USE MAP,               ONLY: MAP_Init
   USE OrcaFlexInterface, ONLY: Orca_Init
   USE HydroDyn,          ONLY: HydroDyn_Init
   USE IceDyn,            ONLY: IceD_Init
   USE IceFloe,           ONLY: IceFloe_Init
   USE ServoDyn,          ONLY: SrvD_Init, Cmpl4SFun, Cmpl4LV
   USE SubDyn,            ONLY: SD_Init
   USE OpenFOAM,          ONLY: Init_OpFM
   USE SuperController,   ONLY: Init_SC
   Use ExtPtfm_MCKF,      ONLY: ExtPtfm_Init

   USE BeamDyn_Types, ONLY: BD_STATIC_ANALYSIS, BD_MESH_QP
   use OpenFOAM_Types
   USE FAST_Types, ONLY: FAST_OutputFileType, FAST_ModuleMapType, FAST_ParameterType, FAST_MiscVarType
   USE FAST_Types, ONLY: FAST_TurbineType, FAST_ExternInitType
   USE FAST_Types, ONLY: NumModules, sensortype_none, IceD_MaxLegs
   USE FAST_Types, ONLY: Module_None, Module_IfW, Module_MAP, Module_IceF, Module_OpFM, Module_SD, Module_MD, Module_Orca
   USE FAST_Types, ONLY: Module_IceD, Module_SrvD, Module_ED, Module_AD14, Module_FEAM, Module_ExtPtfm, Module_HD
   USE FAST_Types, ONLY: Module_AD, Module_BD
   USE FAST_Types, ONLY: ElastoDyn_Data, BeamDyn_Data, ServoDyn_Data, AeroDyn14_Data, AeroDyn_Data, InflowWind_Data, OpenFOAM_Data
   USE FAST_Types, ONLY: SuperController_Data, HydroDyn_Data, SubDyn_Data, ExtPtfm_Data, MAP_Data, FEAMooring_Data, MoorDyn_Data
   USE FAST_Types, ONLY: OrcaFlex_Data, IceFloe_Data, IceDyn_Data
   USE FAST_Types, ONLY: ED_InitInputType , ED_InitOutputType
   USE FAST_Types, ONLY: BD_InitInputType , BD_InitOutputType
   USE FAST_Types, ONLY: SrvD_InitInputType , SrvD_InitOutputType
   USE FAST_Types, ONLY: AD14_InitInputType , AD14_InitOutputType
   USE FAST_Types, ONLY: AD_InitInputType , AD_InitOutputType
   USE FAST_Types, ONLY: InflowWind_InitInputType, InflowWind_InitOutputType
   USE FAST_Types, ONLY: OpFM_InitInputType, OpFM_InitOutputType
   USE FAST_Types, ONLY: SC_InitInputType, SC_InitOutputType
   USE FAST_Types, ONLY: HydroDyn_InitInputType, HydroDyn_InitOutputType
   USE FAST_Types, ONLY: SD_InitInputType, SD_InitOutputType
   USE FAST_Types, ONLY: ExtPtfm_InitInputType, ExtPtfm_InitOutputType
   USE FAST_Types, ONLY: MAP_InitInputType, MAP_InitOutputType
   USE FAST_Types, ONLY: FEAM_InitInputType, FEAM_InitOutputType
   USE FAST_Types, ONLY: MD_InitInputType, MD_InitOutputType
   USE FAST_Types, ONLY: Orca_InitInputType, Orca_InitOutputType
   USE FAST_Types, ONLY: IceFloe_InitInputType, IceFloe_InitOutputType
   USE FAST_Types, ONLY: IceD_InitInputType, IceD_InitOutputType

   IMPLICIT NONE

CONTAINS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INITIALIZATION ROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> a wrapper routine to call FAST_Initialize a the full-turbine simulation level (makes easier to write top-level driver)
SUBROUTINE FAST_InitializeAll_T( t_initial, TurbID, Turbine, ErrStat, ErrMsg, InFile, ExternInitData )

   REAL(DbKi),                        INTENT(IN   ) :: t_initial      !< initial time
   INTEGER(IntKi),                    INTENT(IN   ) :: TurbID         !< turbine Identifier (1-NumTurbines)
   TYPE(FAST_TurbineType),            INTENT(INOUT) :: Turbine        !< all data for one instance of a turbine
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat        !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   CHARACTER(*),             OPTIONAL,INTENT(IN   ) :: InFile         !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)   
   TYPE(FAST_ExternInitType),OPTIONAL,INTENT(IN   ) :: ExternInitData !< Initialization input data from an external source (Simulink)
   
   Turbine%TurbID = TurbID  
   
   
   IF (PRESENT(InFile)) THEN
      IF (PRESENT(ExternInitData)) THEN
         CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC,&
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg, InFile, ExternInitData )
      ELSE         
         CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC, &
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg, InFile  )
      END IF
   ELSE
      CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC, &
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg )
   END IF
   
         
END SUBROUTINE FAST_InitializeAll_T
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to call Init routine for each module. This routine sets all of the init input data for each module.
SUBROUTINE FAST_InitializeAll( t_initial, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, SC, HD, SD, ExtPtfm, &
                               MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg, InFile, ExternInitData )

   use ElastoDyn_Parameters, only: Method_RK4

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(SuperController_Data), INTENT(INOUT) :: SC                !< SuperController data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   CHARACTER(*), OPTIONAL,   INTENT(IN   ) :: InFile              !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)
   
   TYPE(FAST_ExternInitType), OPTIONAL, INTENT(IN) :: ExternInitData !< Initialization input data from an external source (Simulink)
   
   ! local variables      
   CHARACTER(1024)                         :: InputFile           !< A CHARACTER string containing the name of the primary FAST input file

   TYPE(ED_InitInputType)                  :: InitInData_ED       ! Initialization input data
   TYPE(ED_InitOutputType)                 :: InitOutData_ED      ! Initialization output data
   
   TYPE(BD_InitInputType)                  :: InitInData_BD       ! Initialization input data
   TYPE(BD_InitOutputType), ALLOCATABLE    :: InitOutData_BD(:)   ! Initialization output data
                                           
   TYPE(SrvD_InitInputType)                :: InitInData_SrvD     ! Initialization input data
   TYPE(SrvD_InitOutputType)               :: InitOutData_SrvD    ! Initialization output data
                                           
   TYPE(AD14_InitInputType)                :: InitInData_AD14     ! Initialization input data
   TYPE(AD14_InitOutputType)               :: InitOutData_AD14    ! Initialization output data
                                           
   TYPE(AD_InitInputType)                  :: InitInData_AD       ! Initialization input data
   TYPE(AD_InitOutputType)                 :: InitOutData_AD      ! Initialization output data
      
   TYPE(InflowWind_InitInputType)          :: InitInData_IfW      ! Initialization input data
   TYPE(InflowWind_InitOutputType)         :: InitOutData_IfW     ! Initialization output data
   
   TYPE(OpFM_InitInputType)                :: InitInData_OpFM     ! Initialization input data
   TYPE(OpFM_InitOutputType)               :: InitOutData_OpFM    ! Initialization output data
      
   TYPE(SC_InitInputType)                  :: InitInData_SC       ! Initialization input data
   TYPE(SC_InitOutputType)                 :: InitOutData_SC      ! Initialization output data

   TYPE(HydroDyn_InitInputType)            :: InitInData_HD       ! Initialization input data
   TYPE(HydroDyn_InitOutputType)           :: InitOutData_HD      ! Initialization output data
                                           
   TYPE(SD_InitInputType)                  :: InitInData_SD       ! Initialization input data
   TYPE(SD_InitOutputType)                 :: InitOutData_SD      ! Initialization output data
                                           
   TYPE(ExtPtfm_InitInputType)             :: InitInData_ExtPtfm  ! Initialization input data
   TYPE(ExtPtfm_InitOutputType)            :: InitOutData_ExtPtfm ! Initialization output data
                                           
   TYPE(MAP_InitInputType)                 :: InitInData_MAP      ! Initialization input data
   TYPE(MAP_InitOutputType)                :: InitOutData_MAP     ! Initialization output data
                                           
   TYPE(FEAM_InitInputType)                :: InitInData_FEAM     ! Initialization input data
   TYPE(FEAM_InitOutputType)               :: InitOutData_FEAM    ! Initialization output data
                                           
   TYPE(MD_InitInputType)                  :: InitInData_MD       ! Initialization input data
   TYPE(MD_InitOutputType)                 :: InitOutData_MD      ! Initialization output data
             
   TYPE(Orca_InitInputType)                :: InitInData_Orca     ! Initialization input data
   TYPE(Orca_InitOutputType)               :: InitOutData_Orca    ! Initialization output data
   
   TYPE(IceFloe_InitInputType)             :: InitInData_IceF     ! Initialization input data
   TYPE(IceFloe_InitOutputType)            :: InitOutData_IceF    ! Initialization output data
                                           
   TYPE(IceD_InitInputType)                :: InitInData_IceD     ! Initialization input data
   TYPE(IceD_InitOutputType)               :: InitOutData_IceD    ! Initialization output data (each instance will have the same output channels)
       
   REAL(ReKi)                              :: AirDens             ! air density for initialization/normalization of OpenFOAM data
   REAL(DbKi)                              :: dt_IceD             ! tmp dt variable to ensure IceDyn doesn't specify different dt values for different legs (IceDyn instances)
   REAL(DbKi)                              :: dt_BD               ! tmp dt variable to ensure BeamDyn doesn't specify different dt values for different instances
   INTEGER(IntKi)                          :: ErrStat2
   INTEGER(IntKi)                          :: IceDim              ! dimension we're pre-allocating for number of IceDyn legs/instances
   INTEGER(IntKi)                          :: I                   ! generic loop counter
   INTEGER(IntKi)                          :: k                   ! blade loop counter
   logical                                 :: CallStart
   
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
                                           
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_InitializeAll'       
   
   
   !..........
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   y_FAST%UnSum = -1                                                    ! set the summary file unit to -1 to indicate it's not open
   y_FAST%UnOu  = -1                                                    ! set the text output file unit to -1 to indicate it's not open
   y_FAST%UnGra = -1                                                    ! set the binary graphics output file unit to -1 to indicate it's not open
   
   p_FAST%WrVTK = VTK_Unknown                                           ! set this so that we can potentially output VTK information on initialization error
   y_FAST%VTK_LastWaveIndx = 1                                          ! Start looking for wave data at the first index
   y_FAST%VTK_count = 0                                                 ! first VTK file has 0 as output      
   y_FAST%n_Out = 0                                                     ! set the number of ouptut channels to 0 to indicate there's nothing to write to the binary file
   p_FAST%ModuleInitialized = .FALSE.                                   ! (array initialization) no modules are initialized 
   
      ! Get the current time
   CALL DATE_AND_TIME ( Values=m_FAST%StrtTime )                        ! Let's time the whole simulation
   CALL CPU_TIME ( m_FAST%UsrTime1 )                                    ! Initial time (this zeros the start time when used as a MATLAB function)
   m_FAST%UsrTime1 = MAX( 0.0_ReKi, m_FAST%UsrTime1 )                   ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned
   

   m_FAST%t_global        = t_initial - 20.                             ! initialize this to a number < t_initial for error message in ProgAbort
   m_FAST%calcJacobian    = .TRUE.                                      ! we need to calculate the Jacobian
   m_FAST%NextJacCalcTime = m_FAST%t_global                             ! We want to calculate the Jacobian on the first step
   p_FAST%TDesc           = ''

   if (present(ExternInitData)) then
      CallStart = .not. ExternInitData%FarmIntegration ! .and. ExternInitData%TurbineID == 1
      if (ExternInitData%TurbineID > 0) p_FAST%TDesc = 'T'//trim(num2lstr(ExternInitData%TurbineID)) 
   else
      CallStart = .true.
   end if
           
   
      ! Init NWTC_Library, display copyright and version information:
   if (CallStart) then
      AbortErrLev = ErrID_Fatal                                 ! Until we read otherwise from the FAST input file, we abort only on FATAL errors
      CALL FAST_ProgStart( FAST_Ver )
      p_FAST%WrSttsTime = .TRUE.
   else
      ! if we don't call the start data (e.g., from FAST.Farm), we won't override AbortErrLev either 
      CALL DispNVD( FAST_Ver )
      p_FAST%WrSttsTime = .FALSE.
   end if
   
   IF (PRESENT(InFile)) THEN
      p_FAST%UseDWM = .FALSE.
      InputFile = InFile
   ELSE
      CALL GetInputFileName(InputFile,p_FAST%UseDWM,ErrStat2,ErrMsg2)            
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
   END IF
   
   ! ... Open and read input files ...
   ! also, set turbine reference position for graphics output
   if (PRESENT(ExternInitData)) then
      p_FAST%TurbinePos = ExternInitData%TurbinePos
      
      if (ExternInitData%FarmIntegration) then ! we're integrating with FAST.Farm
         CALL FAST_Init( p_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2, ExternInitData%TMax, OverrideAbortLev=.false., RootName=ExternInitData%RootName )         
      else
         CALL FAST_Init( p_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2, ExternInitData%TMax, ExternInitData%TurbineID )  ! We have the name of the input file and the simulation length from somewhere else (e.g. Simulink)         
      end if
      
   else
      p_FAST%TurbinePos = 0.0_ReKi
      CALL FAST_Init( p_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2 )                       ! We have the name of the input file from somewhere else (e.g. Simulink)
   end if
         
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
      
      
   !...............................................................................................................................  
      
   p_FAST%dt_module = p_FAST%dt ! initialize time steps for each module   

   ! ........................
   ! initialize ElastoDyn (must be done first)
   ! ........................
   
   ALLOCATE( ED%Input( p_FAST%InterpOrder+1 ), ED%InputTimes( p_FAST%InterpOrder+1 ), ED%Output( p_FAST%InterpOrder+1 ),STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating ED%Input, ED%Output, and ED%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   InitInData_ED%Linearize = p_FAST%Linearize
   InitInData_ED%InputFile = p_FAST%EDFile
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      InitInData_ED%ADInputFile = p_FAST%AeroFile
   ELSE
      InitInData_ED%ADInputFile = ""
   END IF
   
   InitInData_ED%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_ED))
   InitInData_ED%CompElast     = p_FAST%CompElast == Module_ED

   CALL ED_Init( InitInData_ED, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                  ED%Output(1), ED%m, p_FAST%dt_module( MODULE_ED ), InitOutData_ED, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   p_FAST%ModuleInitialized(Module_ED) = .TRUE.
   CALL SetModuleSubstepTime(Module_ED, p_FAST, y_FAST, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      ! bjj: added this check per jmj; perhaps it would be better in ElastoDyn, but I'll leave it here for now:
   IF ( p_FAST%TurbineType == Type_Offshore_Floating ) THEN
      IF ( ED%p%TowerBsHt < 0.0_ReKi .AND. .NOT. EqualRealNos( ED%p%TowerBsHt, 0.0_ReKi ) ) THEN
         CALL SetErrStat(ErrID_Fatal,"ElastoDyn TowerBsHt must not be negative for floating offshore systems.",ErrStat,ErrMsg,RoutineName)
      END IF      
   END IF   

   allocate( y_FAST%Lin%Modules(MODULE_ED)%Instance(1), stat=ErrStat2)
   if (ErrStat2 /= 0 ) then
      call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(ED).", ErrStat, ErrMsg, RoutineName )
   else
   
      if (allocated(InitOutData_ED%LinNames_y)) call move_alloc(InitOutData_ED%LinNames_y,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%Names_y)
      if (allocated(InitOutData_ED%LinNames_x)) call move_alloc(InitOutData_ED%LinNames_x,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%Names_x)
      if (allocated(InitOutData_ED%LinNames_u)) call move_alloc(InitOutData_ED%LinNames_u,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%Names_u)
      if (allocated(InitOutData_ED%RotFrame_y)) call move_alloc(InitOutData_ED%RotFrame_y,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%RotFrame_y)
      if (allocated(InitOutData_ED%RotFrame_x)) call move_alloc(InitOutData_ED%RotFrame_x,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%RotFrame_x)
      if (allocated(InitOutData_ED%RotFrame_u)) call move_alloc(InitOutData_ED%RotFrame_u,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%RotFrame_u)
      if (allocated(InitOutData_ED%IsLoad_u  )) call move_alloc(InitOutData_ED%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%IsLoad_u  )
         
      if (allocated(InitOutData_ED%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%NumOutputs = size(InitOutData_ED%WriteOutputHdr)
   end if

   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
   
   ! ........................
   ! initialize BeamDyn 
   ! ........................
   IF ( p_FAST%CompElast == Module_BD ) THEN      
      p_FAST%nBeams = InitOutData_ED%NumBl          ! initialize number of BeamDyn instances = number of blades      
   ELSE
      p_FAST%nBeams = 0
   END IF

   ALLOCATE( BD%Input( p_FAST%InterpOrder+1, p_FAST%nBeams ), BD%InputTimes( p_FAST%InterpOrder+1, p_FAST%nBeams ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating BD%Input and BD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF  
                        
   ALLOCATE( BD%x(           p_FAST%nBeams,2), &
             BD%xd(          p_FAST%nBeams,2), &
             BD%z(           p_FAST%nBeams,2), &
             BD%OtherSt(     p_FAST%nBeams,2), &
             BD%p(           p_FAST%nBeams  ), &
             BD%u(           p_FAST%nBeams  ), &
             BD%y(           p_FAST%nBeams  ), &
             BD%m(           p_FAST%nBeams  ), &
             InitOutData_BD( p_FAST%nBeams  ), &
                                             STAT = ErrStat2 )                                                  
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating BeamDyn state, input, and output data.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF        
   
   IF (p_FAST%CompElast == Module_BD) THEN

      InitInData_BD%DynamicSolve = .TRUE.       ! FAST can only couple to BeamDyn when dynamic solve is used.

      InitInData_BD%Linearize = p_FAST%Linearize
      InitInData_BD%gravity      = (/ 0.0_ReKi, 0.0_ReKi, -InitOutData_ED%Gravity /)       ! "Gravitational acceleration" m/s^2
      
         ! now initialize BeamDyn for all beams
      dt_BD = p_FAST%dt_module( MODULE_BD )
                        
      InitInData_BD%HubPos = ED%Output(1)%HubPtMotion%Position(:,1)
      InitInData_BD%HubRot = ED%Output(1)%HubPtMotion%RefOrientation(:,:,1)
            
      p_FAST%BD_OutputSibling = .true.
      
      allocate( y_FAST%Lin%Modules(MODULE_BD)%Instance(p_FAST%nBeams), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(BD).", ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      end if

      DO k=1,p_FAST%nBeams
         InitInData_BD%RootName     = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_BD))//TRIM( Num2LStr(k) )
         
         
         InitInData_BD%InputFile    = p_FAST%BDBldFile(k)
         
         InitInData_BD%GlbPos       = ED%Output(1)%BladeRootMotion(k)%Position(:,1)          ! {:}    - - "Initial Position Vector of the local blade coordinate system"
         InitInData_BD%GlbRot       = ED%Output(1)%BladeRootMotion(k)%RefOrientation(:,:,1)  ! {:}{:} - - "Initial direction cosine matrix of the local blade coordinate system"
         
         InitInData_BD%RootDisp     = ED%Output(1)%BladeRootMotion(k)%TranslationDisp(:,1)   ! {:}    - - "Initial root displacement"
         InitInData_BD%RootOri      = ED%Output(1)%BladeRootMotion(k)%Orientation(:,:,1)     ! {:}{:} - - "Initial root orientation"
         InitInData_BD%RootVel(1:3) = ED%Output(1)%BladeRootMotion(k)%TranslationVel(:,1)    ! {:}    - - "Initial root velocities and angular veolcities"                  
         InitInData_BD%RootVel(4:6) = ED%Output(1)%BladeRootMotion(k)%RotationVel(:,1)       ! {:}    - - "Initial root velocities and angular veolcities"                  
                           
         CALL BD_Init( InitInData_BD, BD%Input(1,k), BD%p(k),  BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), &
                           BD%OtherSt(k,STATE_CURR), BD%y(k),  BD%m(k), dt_BD, InitOutData_BD(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)  
            
         !bjj: we're going to force this to have the same timestep because I don't want to have to deal with n BD modules with n timesteps.
         IF ( k == 1 ) THEN
            p_FAST%dt_module( MODULE_BD ) = dt_BD
            
            p_FAST%ModuleInitialized(Module_BD) = .TRUE. ! this really should be once per BD instance, but BD doesn't care so I won't go through the effort to track this
            CALL SetModuleSubstepTime(Module_BD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)            
         ELSEIF ( .NOT. EqualRealNos( p_FAST%dt_module( MODULE_BD ),dt_BD )) THEN
            CALL SetErrStat(ErrID_Fatal,"All instances of BeamDyn (one per blade) must have the same time step.",ErrStat,ErrMsg,RoutineName)
         END IF
                       
         ! BeamDyn shouldn't be run in static mode when coupled with FAST
         if (BD%p(k)%analysis_type == BD_STATIC_ANALYSIS) then ! static
            CALL SetErrStat(ErrID_Fatal,"BeamDyn cannot perform static analysis when coupled with FAST.",ErrStat,ErrMsg,RoutineName)
         end if

            ! We're going to do fewer computations if the BD input and output meshes that couple to AD are siblings:
         if (BD%p(k)%BldMotionNodeLoc /= BD_MESH_QP) p_FAST%BD_OutputSibling = .false.
      
         if (ErrStat>=AbortErrLev) exit !exit this loop so we don't get p_FAST%nBeams of the same errors
         
         if (allocated(InitOutData_BD(k)%LinNames_y)) call move_alloc(InitOutData_BD(k)%LinNames_y, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%Names_y )
         if (allocated(InitOutData_BD(k)%LinNames_x)) call move_alloc(InitOutData_BD(k)%LinNames_x, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%Names_x )
         if (allocated(InitOutData_BD(k)%LinNames_u)) call move_alloc(InitOutData_BD(k)%LinNames_u, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%Names_u )
         if (allocated(InitOutData_BD(k)%RotFrame_y)) call move_alloc(InitOutData_BD(k)%RotFrame_y, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%RotFrame_y )
         if (allocated(InitOutData_BD(k)%RotFrame_x)) call move_alloc(InitOutData_BD(k)%RotFrame_x, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%RotFrame_x )
         if (allocated(InitOutData_BD(k)%RotFrame_u)) call move_alloc(InitOutData_BD(k)%RotFrame_u, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%RotFrame_u )
         if (allocated(InitOutData_BD(k)%IsLoad_u  )) call move_alloc(InitOutData_BD(k)%IsLoad_u  , y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%IsLoad_u   )
         
         if (allocated(InitOutData_BD(k)%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%NumOutputs = size(InitOutData_BD(k)%WriteOutputHdr)
         
      END DO
               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF           
      
   END IF         
      

   ! ........................
   ! initialize AeroDyn 
   ! ........................
   ALLOCATE( AD14%Input( p_FAST%InterpOrder+1 ), AD14%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating AD14%Input and AD14%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
     
   ALLOCATE( AD%Input( p_FAST%InterpOrder+1 ), AD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating AD%Input and AD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
      
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
               
      CALL AD_SetInitInput(InitInData_AD14, InitOutData_ED, ED%Output(1), p_FAST, ErrStat2, ErrMsg2)            ! set the values in InitInData_AD14
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                                       
      CALL AD14_Init( InitInData_AD14, AD14%Input(1), AD14%p, AD14%x(STATE_CURR), AD14%xd(STATE_CURR), AD14%z(STATE_CURR), &
                     AD14%OtherSt(STATE_CURR), AD14%y, AD14%m, p_FAST%dt_module( MODULE_AD14 ), InitOutData_AD14, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_AD14) = .TRUE.            
      CALL SetModuleSubstepTime(Module_AD14, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
         ! bjj: this really shouldn't be in the FAST glue code, but I'm going to put this check here so people don't use an invalid model 
         !    and send me emails to debug numerical issues in their results.
      IF ( AD14%p%TwrProps%PJM_Version .AND. p_FAST%TurbineType == Type_Offshore_Floating ) THEN
         CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 tower influence model "NEWTOWER" is invalid for models of floating offshore turbines.',ErrStat,ErrMsg,RoutineName)
      END IF         
            
      AirDens = InitOutData_AD14%AirDens
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
      
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      
      
         ! set initialization data for AD
      CALL AllocAry( InitInData_AD%BladeRootPosition,      3, InitOutData_ED%NumBl, 'InitInData_AD%BladeRootPosition', errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry( InitInData_AD%BladeRootOrientation,3, 3, InitOutData_ED%NumBl, 'InitInData_AD%BladeRootOrientation', errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
      InitInData_AD%Gravity            = InitOutData_ED%Gravity      
      InitInData_AD%Linearize          = p_FAST%Linearize
      InitInData_AD%InputFile          = p_FAST%AeroFile
      InitInData_AD%NumBlades          = InitOutData_ED%NumBl
      InitInData_AD%RootName           = p_FAST%OutFileRoot
      InitInData_AD%HubPosition        = ED%Output(1)%HubPtMotion%Position(:,1)
      InitInData_AD%HubOrientation     = ED%Output(1)%HubPtMotion%RefOrientation(:,:,1)
      
      do k=1,InitOutData_ED%NumBl
         InitInData_AD%BladeRootPosition(:,k)      = ED%Output(1)%BladeRootMotion(k)%Position(:,1)
         InitInData_AD%BladeRootOrientation(:,:,k) = ED%Output(1)%BladeRootMotion(k)%RefOrientation(:,:,1)
      end do
      
            
      CALL AD_Init( InitInData_AD, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                    AD%OtherSt(STATE_CURR), AD%y, AD%m, p_FAST%dt_module( MODULE_AD ), InitOutData_AD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_AD) = .TRUE.            
      CALL SetModuleSubstepTime(Module_AD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                               
      allocate( y_FAST%Lin%Modules(MODULE_AD)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(AD).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(InitOutData_AD%LinNames_u)) call move_alloc(InitOutData_AD%LinNames_u,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%Names_u )
         if (allocated(InitOutData_AD%LinNames_y)) call move_alloc(InitOutData_AD%LinNames_y,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%Names_y )
         if (allocated(InitOutData_AD%LinNames_z)) call move_alloc(InitOutData_AD%LinNames_z,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%Names_z )
         if (allocated(InitOutData_AD%RotFrame_u)) call move_alloc(InitOutData_AD%RotFrame_u,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%RotFrame_u )
         if (allocated(InitOutData_AD%RotFrame_y)) call move_alloc(InitOutData_AD%RotFrame_y,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%RotFrame_y )
         if (allocated(InitOutData_AD%RotFrame_z)) call move_alloc(InitOutData_AD%RotFrame_z,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%RotFrame_z )
         if (allocated(InitOutData_AD%IsLoad_u  )) call move_alloc(InitOutData_AD%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%IsLoad_u   )
         
         if (allocated(InitOutData_AD%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%NumOutputs = size(InitOutData_AD%WriteOutputHdr)
      end if
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
      
      AirDens = InitOutData_AD%AirDens
      
   ELSE
      AirDens = 0.0_ReKi
   END IF ! CompAero
   
               
   ! ........................
   ! initialize InflowWind
   ! ........................   
   ALLOCATE( IfW%Input( p_FAST%InterpOrder+1 ), IfW%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IfW%Input and IfW%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
              
   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      
      InitInData_IfW%Linearize        = p_FAST%Linearize
      InitInData_IfW%InputFileName    = p_FAST%InflowFile
      InitInData_IfW%RootName         = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IfW))
      InitInData_IfW%UseInputFile     = .TRUE.
   
      InitInData_IfW%NumWindPoints = 0      
      IF ( p_FAST%CompServo == Module_SrvD ) InitInData_IfW%NumWindPoints = InitInData_IfW%NumWindPoints + 1
      IF ( p_FAST%CompAero  == Module_AD14 ) THEN
         InitInData_IfW%NumWindPoints = InitInData_IfW%NumWindPoints + InitOutData_ED%NumBl * AD14%Input(1)%InputMarkers(1)%NNodes + AD14%Input(1)%Twr_InputMarkers%NNodes
      ELSEIF ( p_FAST%CompAero  == Module_AD ) THEN
         InitInData_IfW%NumWindPoints = InitInData_IfW%NumWindPoints + AD%Input(1)%TowerMotion%NNodes
         DO k=1,InitOutData_ED%NumBl
            InitInData_IfW%NumWindPoints = InitInData_IfW%NumWindPoints + AD%Input(1)%BladeMotion(k)%NNodes
         END DO
      END IF
      
      ! lidar        
      InitInData_IfW%lidar%Tmax                   = p_FAST%TMax
      InitInData_IfW%lidar%HubPosition            = ED%Output(1)%HubPtMotion%Position(:,1) 
      IF ( PRESENT(ExternInitData) ) THEN
         InitInData_IfW%Use4Dext = ExternInitData%FarmIntegration

         if (InitInData_IfW%Use4Dext) then
            InitInData_IfW%FDext%n      = ExternInitData%windGrid_n
            InitInData_IfW%FDext%delta  = ExternInitData%windGrid_delta
            InitInData_IfW%FDext%pZero  = ExternInitData%windGrid_pZero
         end if
         
         ! bjj: these lidar inputs should come from an InflowWind input file; I'm hard coding them here for now
         InitInData_IfW%lidar%SensorType          = ExternInitData%SensorType   
         InitInData_IfW%lidar%LidRadialVel        = ExternInitData%LidRadialVel   
         InitInData_IfW%lidar%RotorApexOffsetPos  = 0.0         
         InitInData_IfW%lidar%NumPulseGate        = 0
      ELSE
         InitInData_IfW%lidar%SensorType          = SensorType_None
         InitInData_IfW%Use4Dext                  = .false.
      END IF
                                     
      CALL InflowWind_Init( InitInData_IfW, IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR),  &
                     IfW%OtherSt(STATE_CURR), IfW%y, IfW%m, p_FAST%dt_module( MODULE_IfW ), InitOutData_IfW, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_IfW) = .TRUE.            
      CALL SetModuleSubstepTime(Module_IfW, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      allocate( y_FAST%Lin%Modules(MODULE_IfW)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(IfW).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(InitOutData_IfW%LinNames_y)) call move_alloc(InitOutData_IfW%LinNames_y,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%Names_y )
         if (allocated(InitOutData_IfW%LinNames_u)) call move_alloc(InitOutData_IfW%LinNames_u,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%Names_u )
         if (allocated(InitOutData_IfW%RotFrame_y)) call move_alloc(InitOutData_IfW%RotFrame_y,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%RotFrame_y )
         if (allocated(InitOutData_IfW%RotFrame_u)) call move_alloc(InitOutData_IfW%RotFrame_u,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%RotFrame_u )
         if (allocated(InitOutData_IfW%IsLoad_u  )) call move_alloc(InitOutData_IfW%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%IsLoad_u   )

         if (allocated(InitOutData_IfW%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%NumOutputs = size(InitOutData_IfW%WriteOutputHdr)
      end if
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
      
   ELSEIF ( p_FAST%CompInflow == Module_OpFM ) THEN
      
      IF ( PRESENT(ExternInitData) ) THEN
         InitInData_OpFM%NumSC2Ctrl = ExternInitData%NumSC2Ctrl
         InitInData_OpFM%NumCtrl2SC = ExternInitData%NumCtrl2SC  
         InitInData_OpFM%NumActForcePtsBlade = ExternInitData%NumActForcePtsBlade
         InitInData_OpFM%NumActForcePtsTower = ExternInitData%NumActForcePtsTower 
      ELSE
         CALL SetErrStat( ErrID_Fatal, 'OpenFOAM integration can be used only with external input data (not the stand-alone executable).', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN         
      END IF
      InitInData_OpFM%BladeLength = InitOutData_ED%BladeLength
      InitInData_OpFM%TowerHeight = InitOutData_ED%TowerHeight
      InitInData_OpFM%TowerBaseHeight = InitOutData_ED%TowerBaseHeight
      ALLOCATE(InitInData_OpFM%StructBldRNodes( SIZE(InitOutData_ED%BldRNodes)),  STAT=ErrStat2)
      InitInData_OpFM%StructBldRNodes(:) = InitOutData_ED%BldRNodes(:)
      ALLOCATE(InitInData_OpFM%StructTwrHNodes( SIZE(InitOutData_ED%TwrHNodes)),  STAT=ErrStat2)
      InitInData_OpFM%StructTwrHNodes(:) = InitOutData_ED%TwrHNodes(:)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating OpFM%InitInput.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
         ! set up the data structures for integration with OpenFOAM
      CALL Init_OpFM( InitInData_OpFM, p_FAST, AirDens, AD14%Input(1), AD%Input(1), InitOutData_AD, AD%y, ED%Output(1), OpFM, InitOutData_OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
                  
      !bjj: fix me!!! to do
      InitOutData_IfW%WindFileInfo%MWS = 0.0_ReKi
      
   ELSE
      InitOutData_IfW%WindFileInfo%MWS = 0.0_ReKi
   END IF   ! CompInflow
   
   ! ........................
   ! initialize SuperController
   ! ........................   
      IF ( PRESENT(ExternInitData) ) THEN
         InitInData_SC%NumSC2Ctrl = ExternInitData%NumSC2Ctrl
         InitInData_SC%NumCtrl2SC = ExternInitData%NumCtrl2SC  
      ELSE
         InitInData_SC%NumSC2Ctrl = 0
         InitInData_SC%NumCtrl2SC = 0
      END IF
      
         ! set up the data structures for integration with supercontroller
      CALL Init_SC( InitInData_SC, SC, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       

   ! ........................
   ! some checks for AeroDyn14's Dynamic Inflow with Mean Wind Speed from InflowWind:
   ! (DO NOT COPY THIS CODE!)
   ! bjj: AeroDyn14 should not need this rule of thumb; it should check the instantaneous values when the code runs
   ! ........................   
   
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      IF (AD14%p%DynInfl) THEN               
         IF ( InitOutData_IfW%WindFileInfo%MWS  < 8.0 ) THEN
            CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 "DYNINFL" InfModel is invalid for models with wind speeds less than 8 m/s.',ErrStat,ErrMsg,RoutineName)
            !CALL SetErrStat(ErrID_Info,'Estimated average inflow wind speed is less than 8 m/s. Dynamic Inflow will be turned off.',ErrStat,ErrMess,RoutineName )
         END IF
      END IF      
   END IF
   
   
   ! ........................
   ! initialize ServoDyn 
   ! ........................
   ALLOCATE( SrvD%Input( p_FAST%InterpOrder+1 ), SrvD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating SrvD%Input and SrvD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      InitInData_SrvD%InputFile     = p_FAST%ServoFile
      InitInData_SrvD%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_SrvD))
      InitInData_SrvD%NumBl         = InitOutData_ED%NumBl
      InitInData_SrvD%gravity       = InitOutData_ED%gravity
      InitInData_SrvD%r_N_O_G       = ED%Input(1)%NacelleLoads%Position(:,1)
      InitInData_SrvD%r_TwrBase     = InitOutData_ED%TwrBasePos
      InitInData_SrvD%TMax          = p_FAST%TMax
      InitInData_SrvD%AirDens       = AirDens
      InitInData_SrvD%AvgWindSpeed  = InitOutData_IfW%WindFileInfo%MWS
      InitInData_SrvD%Linearize     = p_FAST%Linearize
      
      IF ( PRESENT(ExternInitData) ) THEN
         InitInData_SrvD%NumSC2Ctrl = ExternInitData%NumSC2Ctrl
         InitInData_SrvD%NumCtrl2SC = ExternInitData%NumCtrl2SC  
      ELSE
         InitInData_SrvD%NumSC2Ctrl = 0
         InitInData_SrvD%NumCtrl2SC = 0
      END IF      
            
      CALL AllocAry(InitInData_SrvD%BlPitchInit, InitOutData_ED%NumBl, 'BlPitchInit', ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      InitInData_SrvD%BlPitchInit   = InitOutData_ED%BlPitch
      CALL SrvD_Init( InitInData_SrvD, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                      SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, p_FAST%dt_module( MODULE_SrvD ), InitOutData_SrvD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      p_FAST%ModuleInitialized(Module_SrvD) = .TRUE.

      !IF ( InitOutData_SrvD%CouplingScheme == ExplicitLoose ) THEN ...  bjj: abort if we're doing anything else!

      CALL SetModuleSubstepTime(Module_SrvD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      !! initialize SrvD%y%ElecPwr and SrvD%y%GenTq because they are one timestep different (used as input for the next step)?
                  
      allocate( y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(SrvD).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(InitOutData_SrvD%LinNames_y)) call move_alloc(InitOutData_SrvD%LinNames_y,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%Names_y )
         if (allocated(InitOutData_SrvD%LinNames_u)) call move_alloc(InitOutData_SrvD%LinNames_u,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%Names_u )
         if (allocated(InitOutData_SrvD%RotFrame_y)) call move_alloc(InitOutData_SrvD%RotFrame_y,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%RotFrame_y )
         if (allocated(InitOutData_SrvD%RotFrame_u)) call move_alloc(InitOutData_SrvD%RotFrame_u,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%RotFrame_u )
         if (allocated(InitOutData_SrvD%IsLoad_u  )) call move_alloc(InitOutData_SrvD%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%IsLoad_u   )

         if (allocated(InitOutData_SrvD%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%NumOutputs = size(InitOutData_SrvD%WriteOutputHdr)
      end if
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
   ! ........................
   ! some checks for AeroDyn and ElastoDyn inputs with the high-speed shaft brake hack in ElastoDyn:
   ! (DO NOT COPY THIS CODE!)
   ! ........................   
         ! bjj: this is a hack to get high-speed shaft braking in FAST v8
      
      IF ( InitOutData_SrvD%UseHSSBrake ) THEN
         IF ( p_FAST%CompAero == Module_AD14 ) THEN
            IF ( AD14%p%DYNINFL ) THEN
               CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 "DYNINFL" InfModel is invalid for models with high-speed shaft braking.',ErrStat,ErrMsg,RoutineName)
            END IF
         END IF
         

         IF ( ED%p%method == Method_RK4 ) THEN ! bjj: should be using ElastoDyn's Method_ABM4 Method_AB4 parameters
            CALL SetErrStat(ErrID_Fatal,'ElastoDyn must use the AB4 or ABM4 integration method to implement high-speed shaft braking.',ErrStat,ErrMsg,RoutineName)
         ENDIF
      END IF ! InitOutData_SrvD%UseHSSBrake
      
      
   END IF

   ! ........................
   ! set some VTK parameters required before HydroDyn init (so we can get wave elevations for visualization)
   ! ........................
   
      ! get wave elevation data for visualization
   if ( p_FAST%WrVTK > VTK_None ) then   
      call SetVTKParameters_B4HD(p_FAST, InitOutData_ED, InitInData_HD, BD, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF       
   end if
   
   
   ! ........................
   ! initialize HydroDyn 
   ! ........................
   ALLOCATE( HD%Input( p_FAST%InterpOrder+1 ), HD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating HD%Input and HD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   IF ( p_FAST%CompHydro == Module_HD ) THEN

      InitInData_HD%Gravity       = InitOutData_ED%Gravity
      InitInData_HD%UseInputFile  = .TRUE.
      InitInData_HD%InputFile     = p_FAST%HydroFile
      InitInData_HD%OutRootName   = p_FAST%OutFileRoot
      InitInData_HD%TMax          = p_FAST%TMax
      InitInData_HD%hasIce        = p_FAST%CompIce /= Module_None
            
      
         ! if wave field needs an offset, modify these values (added at request of SOWFA developers):
      InitInData_HD%PtfmLocationX = p_FAST%TurbinePos(1) 
      InitInData_HD%PtfmLocationY = p_FAST%TurbinePos(2)
      
      CALL HydroDyn_Init( InitInData_HD, HD%Input(1), HD%p,  HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), &
                          HD%OtherSt(STATE_CURR), HD%y, HD%m, p_FAST%dt_module( MODULE_HD ), InitOutData_HD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_HD) = .TRUE.
      CALL SetModuleSubstepTime(Module_HD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
   END IF   ! CompHydro

   ! ........................
   ! initialize SubDyn or ExtPtfm_MCKF
   ! ........................
   ALLOCATE( SD%Input( p_FAST%InterpOrder+1 ), SD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating SD%Input and SD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   ALLOCATE( ExtPtfm%Input( p_FAST%InterpOrder+1 ), ExtPtfm%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating ExtPtfm%Input and ExtPtfm%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   IF ( p_FAST%CompSub == Module_SD ) THEN
          
      IF ( p_FAST%CompHydro == Module_HD ) THEN
         InitInData_SD%WtrDpth = InitOutData_HD%WtrDpth
      ELSE
         InitInData_SD%WtrDpth = 0.0_ReKi
      END IF
            
      InitInData_SD%g             = InitOutData_ED%Gravity     
      !InitInData_SD%UseInputFile = .TRUE. 
      InitInData_SD%SDInputFile   = p_FAST%SubFile
      InitInData_SD%RootName      = p_FAST%OutFileRoot
      InitInData_SD%TP_RefPoint   = ED%Output(1)%PlatformPtMesh%Position(:,1)  ! bjj: not sure what this is supposed to be 
      InitInData_SD%SubRotateZ    = 0.0                                        ! bjj: not sure what this is supposed to be 
      
            
      CALL SD_Init( InitInData_SD, SD%Input(1), SD%p,  SD%x(STATE_CURR), SD%xd(STATE_CURR), SD%z(STATE_CURR),  &
                    SD%OtherSt(STATE_CURR), SD%y, SD%m, p_FAST%dt_module( MODULE_SD ), InitOutData_SD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_SD) = .TRUE.
      CALL SetModuleSubstepTime(Module_SD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF   
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN

      InitInData_ExtPtfm%InputFile = p_FAST%SubFile
!      InitInData_ExtPtfm%RootName  = trim(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_ExtPtfm))
      InitInData_ExtPtfm%Linearize = p_FAST%Linearize
      
      
      CALL ExtPtfm_Init( InitInData_ExtPtfm, ExtPtfm%Input(1), ExtPtfm%p,  &
                         ExtPtfm%x(STATE_CURR), ExtPtfm%xd(STATE_CURR), ExtPtfm%z(STATE_CURR),  ExtPtfm%OtherSt(STATE_CURR), &
                         ExtPtfm%y, ExtPtfm%m, p_FAST%dt_module( MODULE_ExtPtfm ), InitOutData_ExtPtfm, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(MODULE_ExtPtfm) = .TRUE.
      CALL SetModuleSubstepTime(MODULE_ExtPtfm, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF   
      
   END IF

   ! ------------------------------
   ! initialize CompMooring modules 
   ! ------------------------------
   ALLOCATE( MAPp%Input( p_FAST%InterpOrder+1 ), MAPp%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating MAPp%Input and MAPp%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
   ALLOCATE( MD%Input( p_FAST%InterpOrder+1 ), MD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating MD%Input and MD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF   
   ALLOCATE( FEAM%Input( p_FAST%InterpOrder+1 ), FEAM%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating FEAM%Input and FEAM%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF   
   ALLOCATE( Orca%Input( p_FAST%InterpOrder+1 ), Orca%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating Orca%Input and Orca%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF   
      
   ! ........................
   ! initialize MAP 
   ! ........................
   IF (p_FAST%CompMooring == Module_MAP) THEN
      !bjj: until we modify this, MAP requires HydroDyn to be used. (perhaps we could send air density from AeroDyn or something...)
      
      CALL WrScr(NewLine) !bjj: I'm printing two blank lines here because MAP seems to be writing over the last line on the screen.
      

!      InitInData_MAP%rootname          =  p_FAST%OutFileRoot        ! Output file name 
      InitInData_MAP%gravity           =  InitOutData_ED%Gravity    ! This need to be according to g used in ElastoDyn
      InitInData_MAP%sea_density       =  InitOutData_HD%WtrDens    ! This needs to be set according to seawater density in HydroDyn
      InitInData_MAP%depth             =  InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
                  
   ! differences for MAP++
      InitInData_MAP%file_name         =  p_FAST%MooringFile        ! This needs to be set according to what is in the FAST input file. 
      InitInData_MAP%summary_file_name =  TRIM(p_FAST%OutFileRoot)//'.MAP.sum'        ! Output file name 
      InitInData_MAP%depth             = -InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            

      
      CALL MAP_Init( InitInData_MAP, MAPp%Input(1), MAPp%p,  MAPp%x(STATE_CURR), MAPp%xd(STATE_CURR), MAPp%z(STATE_CURR), MAPp%OtherSt, &
                      MAPp%y, p_FAST%dt_module( MODULE_MAP ), InitOutData_MAP, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_MAP) = .TRUE.
      CALL SetModuleSubstepTime(Module_MAP, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF              
   ! ........................
   ! initialize MoorDyn 
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
                        
      InitInData_MD%FileName  = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file. 
      InitInData_MD%RootName  = p_FAST%OutFileRoot
      
      InitInData_MD%PtfmInit  = InitOutData_ED%PlatformPos !ED%x(STATE_CURR)%QT(1:6)   ! initial position of the platform !bjj: this should come from InitOutData_ED, not x_ED
      InitInData_MD%g         = InitOutData_ED%Gravity     ! This need to be according to g used in ElastoDyn 
      InitInData_MD%rhoW      = InitOutData_HD%WtrDens     ! This needs to be set according to seawater density in HydroDyn      
      InitInData_MD%WtrDepth  = InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            
      CALL MD_Init( InitInData_MD, MD%Input(1), MD%p, MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), &
                    MD%OtherSt(STATE_CURR), MD%y, MD%m, p_FAST%dt_module( MODULE_MD ), InitOutData_MD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_MD) = .TRUE.
      CALL SetModuleSubstepTime(Module_MD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   ! ........................
   ! initialize FEAM 
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
            
      InitInData_FEAM%InputFile   = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file. 
      InitInData_FEAM%RootName    = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_FEAM))
      
      InitInData_FEAM%PtfmInit    = InitOutData_ED%PlatformPos !ED%x(STATE_CURR)%QT(1:6)   ! initial position of the platform !bjj: this should come from InitOutData_ED, not x_ED
      InitInData_FEAM%NStepWave   = 1                          ! an arbitrary number > 0 (to set the size of the wave data, which currently contains all zero values)     
      InitInData_FEAM%gravity     = InitOutData_ED%Gravity     ! This need to be according to g used in ElastoDyn 
      InitInData_FEAM%WtrDens     = InitOutData_HD%WtrDens     ! This needs to be set according to seawater density in HydroDyn      
!      InitInData_FEAM%depth       =  InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            
      CALL FEAM_Init( InitInData_FEAM, FEAM%Input(1), FEAM%p,  FEAM%x(STATE_CURR), FEAM%xd(STATE_CURR), FEAM%z(STATE_CURR), &
                      FEAM%OtherSt(STATE_CURR), FEAM%y, FEAM%m, p_FAST%dt_module( MODULE_FEAM ), InitOutData_FEAM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_FEAM) = .TRUE.
      CALL SetModuleSubstepTime(Module_FEAM, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF              
   ! ........................
   ! initialize OrcaFlex Interface 
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_Orca) THEN
            
      InitInData_Orca%InputFile = p_FAST%MooringFile
      InitInData_Orca%RootName  = p_FAST%OutFileRoot
      InitInData_Orca%TMax      = p_FAST%TMax 
                  
      CALL Orca_Init( InitInData_Orca, Orca%Input(1), Orca%p,  Orca%x(STATE_CURR), Orca%xd(STATE_CURR), Orca%z(STATE_CURR), Orca%OtherSt(STATE_CURR), &
                      Orca%y, Orca%m, p_FAST%dt_module( MODULE_Orca ), InitOutData_Orca, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(MODULE_Orca) = .TRUE.
      CALL SetModuleSubstepTime(MODULE_Orca, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   END IF

   ! ------------------------------
   ! initialize CompIce modules 
   ! ------------------------------
   ALLOCATE( IceF%Input( p_FAST%InterpOrder+1 ), IceF%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceF%Input and IceF%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF  
      
      ! We need this to be allocated (else we have issues passing nonallocated arrays and using the first index of Input(),
      !   but we don't need the space of IceD_MaxLegs if we're not using it. 
   IF ( p_FAST%CompIce /= Module_IceD ) THEN   
      IceDim = 1
   ELSE
      IceDim = IceD_MaxLegs
   END IF
      
      ! because there may be multiple instances of IceDyn, we'll allocate arrays for that here
      ! we could allocate these after 
   ALLOCATE( IceD%Input( p_FAST%InterpOrder+1, IceDim ), IceD%InputTimes( p_FAST%InterpOrder+1, IceDim ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceD%Input and IceD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF  
      
     ALLOCATE( IceD%x(           IceDim,2), &
               IceD%xd(          IceDim,2), &
               IceD%z(           IceDim,2), &
               IceD%OtherSt(     IceDim,2), &
               IceD%p(           IceDim  ), &
               IceD%u(           IceDim  ), &
               IceD%y(           IceDim  ), &
               IceD%m(           IceDim  ), &
                                             STAT = ErrStat2 )                                                  
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceD state, input, and output data.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF      
         
         
   ! ........................
   ! initialize IceFloe 
   ! ........................
   IF ( p_FAST%CompIce == Module_IceF ) THEN
                      
      InitInData_IceF%InputFile     = p_FAST%IceFile
      InitInData_IceF%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceF))     
      InitInData_IceF%simLength     = p_FAST%TMax  !bjj: IceFloe stores this as single-precision (ReKi) TMax is DbKi
      InitInData_IceF%MSL2SWL       = InitOutData_HD%MSL2SWL
      InitInData_IceF%gravity       = InitOutData_ED%Gravity
      
      CALL IceFloe_Init( InitInData_IceF, IceF%Input(1), IceF%p,  IceF%x(STATE_CURR), IceF%xd(STATE_CURR), IceF%z(STATE_CURR), &
                         IceF%OtherSt(STATE_CURR), IceF%y, IceF%m, p_FAST%dt_module( MODULE_IceF ), InitOutData_IceF, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_IceF) = .TRUE.
      CALL SetModuleSubstepTime(Module_IceF, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
              
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF              
   ! ........................
   ! initialize IceDyn 
   ! ........................
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN  
      
      InitInData_IceD%InputFile     = p_FAST%IceFile
      InitInData_IceD%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceD))//'1'     
      InitInData_IceD%MSL2SWL       = InitOutData_HD%MSL2SWL      
      InitInData_IceD%WtrDens       = InitOutData_HD%WtrDens    
      InitInData_IceD%gravity       = InitOutData_ED%Gravity
      InitInData_IceD%TMax          = p_FAST%TMax
      InitInData_IceD%LegNum        = 1
      
      CALL IceD_Init( InitInData_IceD, IceD%Input(1,1), IceD%p(1),  IceD%x(1,STATE_CURR), IceD%xd(1,STATE_CURR), IceD%z(1,STATE_CURR), &
                      IceD%OtherSt(1,STATE_CURR), IceD%y(1), IceD%m(1), p_FAST%dt_module( MODULE_IceD ), InitOutData_IceD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      p_FAST%ModuleInitialized(Module_IceD) = .TRUE.
      CALL SetModuleSubstepTime(Module_IceD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
         ! now initialize IceD for additional legs (if necessary)
      dt_IceD           = p_FAST%dt_module( MODULE_IceD )
      p_FAST%numIceLegs = InitOutData_IceD%numLegs     
      
      IF (p_FAST%numIceLegs > IceD_MaxLegs) THEN
         CALL SetErrStat(ErrID_Fatal,'IceDyn-FAST coupling is supported for up to '//TRIM(Num2LStr(IceD_MaxLegs))//' legs, but ' &
                           //TRIM(Num2LStr(p_FAST%numIceLegs))//' legs were specified.',ErrStat,ErrMsg,RoutineName)
      END IF
                  

      DO i=2,p_FAST%numIceLegs  ! basically, we just need IceDyn to set up its meshes for inputs/outputs and possibly initial values for states
         InitInData_IceD%LegNum = i
         InitInData_IceD%RootName = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceD))//TRIM(Num2LStr(i))     
         
         CALL IceD_Init( InitInData_IceD, IceD%Input(1,i), IceD%p(i),  IceD%x(i,STATE_CURR), IceD%xd(i,STATE_CURR), IceD%z(i,STATE_CURR), &
                            IceD%OtherSt(i,STATE_CURR), IceD%y(i), IceD%m(i), dt_IceD, InitOutData_IceD, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
         !bjj: we're going to force this to have the same timestep because I don't want to have to deal with n IceD modules with n timesteps.
         IF (.NOT. EqualRealNos( p_FAST%dt_module( MODULE_IceD ),dt_IceD )) THEN
            CALL SetErrStat(ErrID_Fatal,"All instances of IceDyn (one per support-structure leg) must be the same",ErrStat,ErrMsg,RoutineName)
         END IF
      END DO
            
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF           
      
   END IF   
   

   ! ........................
   ! Set up output for glue code (must be done after all modules are initialized so we have their WriteOutput information)
   ! ........................

   CALL FAST_InitOutput( p_FAST, y_FAST, InitOutData_ED, InitOutData_BD, InitOutData_SrvD, InitOutData_AD14, InitOutData_AD, &
                         InitOutData_IfW, InitOutData_OpFM, InitOutData_HD, InitOutData_SD, InitOutData_ExtPtfm, InitOutData_MAP, &
                         InitOutData_FEAM, InitOutData_MD, InitOutData_Orca, InitOutData_IceF, InitOutData_IceD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data
   ! -------------------------------------------------------------------------

   CALL InitModuleMappings(p_FAST, ED, BD, AD14, AD, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF      
      
   ! -------------------------------------------------------------------------
   ! Initialize for linearization:
   ! -------------------------------------------------------------------------
   if ( p_FAST%Linearize ) then      
      call Init_Lin(p_FAST, y_FAST, m_FAST, AD, InitOutData_ED%NumBl, ErrStat2, ErrMsg2)      
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if      
   end if
   
      
   ! -------------------------------------------------------------------------
   ! Initialize data for VTK output
   ! -------------------------------------------------------------------------
   if ( p_FAST%WrVTK > VTK_None ) then
      call SetVTKParameters(p_FAST, InitOutData_ED, InitOutData_AD, InitInData_HD, InitOutData_HD, ED, BD, AD, HD, ErrStat2, ErrMsg2)      
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if
   
   ! -------------------------------------------------------------------------
   ! Write initialization data to FAST summary file:
   ! -------------------------------------------------------------------------
   
   CALL FAST_WrSum( p_FAST, y_FAST, MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
   
   ! -------------------------------------------------------------------------
   ! other misc variables initialized here:
   ! -------------------------------------------------------------------------
      
   m_FAST%t_global        = t_initial
   m_FAST%NextLinTimeIndx = 1 
         
   ! Initialize external inputs for first step  
   if ( p_FAST%CompServo == MODULE_SrvD ) then      
      m_FAST%ExternInput%GenTrq     = SrvD%Input(1)%ExternalGenTrq !0.0_ReKi
      m_FAST%ExternInput%ElecPwr    = SrvD%Input(1)%ExternalElecPwr
      m_FAST%ExternInput%YawPosCom  = SrvD%Input(1)%ExternalYawPosCom
      m_FAST%ExternInput%YawRateCom = SrvD%Input(1)%ExternalYawRateCom
      m_FAST%ExternInput%HSSBrFrac  = SrvD%Input(1)%ExternalHSSBrFrac
      
      do i=1,SIZE(SrvD%Input(1)%ExternalBlPitchCom)
         m_FAST%ExternInput%BlPitchCom(i) = SrvD%Input(1)%ExternalBlPitchCom(i)
      end do   
   end if
   
   m_FAST%ExternInput%LidarFocus = 1.0_ReKi  ! make this non-zero (until we add the initial position in the InflowWind input file)
         
   
   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................      
   CALL Cleanup()
   
CONTAINS
   SUBROUTINE Cleanup()
   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................
   
      CALL ED_DestroyInitInput(  InitInData_ED,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL ED_DestroyInitOutput( InitOutData_ED, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      CALL BD_DestroyInitInput(  InitInData_BD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)         
      IF (ALLOCATED(InitOutData_BD)) THEN
         DO i=1,p_FAST%nBeams
            CALL BD_DestroyInitOutput( InitOutData_BD(i), ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         END DO
         DEALLOCATE(InitOutData_BD)
      END IF
                     
      CALL AD14_DestroyInitInput(  InitInData_AD14,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AD14_DestroyInitOutput( InitOutData_AD14, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL AD_DestroyInitInput(  InitInData_AD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AD_DestroyInitOutput( InitOutData_AD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      CALL InflowWind_DestroyInitInput(  InitInData_IfW,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL InflowWind_DestroyInitOutput( InitOutData_IfW, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      CALL OpFM_DestroyInitInput(  InitInData_OpFM,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL OpFM_DestroyInitOutput( InitOutData_OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)         
         
      CALL SrvD_DestroyInitInput(  InitInData_SrvD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL SrvD_DestroyInitOutput( InitOutData_SrvD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL HydroDyn_DestroyInitInput(  InitInData_HD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL HydroDyn_DestroyInitOutput( InitOutData_HD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL SD_DestroyInitInput(  InitInData_SD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL SD_DestroyInitOutput( InitOutData_SD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      CALL ExtPtfm_DestroyInitInput(  InitInData_ExtPtfm,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL ExtPtfm_DestroyInitOutput( InitOutData_ExtPtfm, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      CALL MAP_DestroyInitInput(  InitInData_MAP,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL MAP_DestroyInitOutput( InitOutData_MAP, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      CALL FEAM_DestroyInitInput(  InitInData_FEAM,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL FEAM_DestroyInitOutput( InitOutData_FEAM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL MD_DestroyInitInput(  InitInData_MD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL MD_DestroyInitOutput( InitOutData_MD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  
      CALL Orca_DestroyInitInput(  InitInData_Orca,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL Orca_DestroyInitOutput( InitOutData_Orca, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  
      CALL IceFloe_DestroyInitInput(  InitInData_IceF,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL IceFloe_DestroyInitOutput( InitOutData_IceF, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      CALL IceD_DestroyInitInput(  InitInData_IceD,  ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL IceD_DestroyInitOutput( InitOutData_IceD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
   
   END SUBROUTINE Cleanup

END SUBROUTINE FAST_InitializeAll
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine checks for command-line arguments, gets the root name of the input files
!! (including full path name), and creates the names of the output files.
SUBROUTINE FAST_Init( p, y_FAST, t_initial, InputFile, ErrStat, ErrMsg, TMax, TurbID, OverrideAbortLev, RootName )

      IMPLICIT                        NONE

   ! Passed variables

   TYPE(FAST_ParameterType), INTENT(INOUT)         :: p                 !< The parameter data for the FAST (glue-code) simulation
   TYPE(FAST_OutputFileType),INTENT(INOUT)         :: y_FAST            !< The output data for the FAST (glue-code) simulation
   REAL(DbKi),               INTENT(IN)            :: t_initial         !< the beginning time of the simulation
   INTEGER(IntKi),           INTENT(OUT)           :: ErrStat           !< Error status
   CHARACTER(*),             INTENT(OUT)           :: ErrMsg            !< Error message
   CHARACTER(*),             INTENT(IN)            :: InputFile         !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)
   REAL(DbKi),               INTENT(IN), OPTIONAL  :: TMax              !< the length of the simulation (from Simulink or FAST.Farm)
   INTEGER(IntKi),           INTENT(IN), OPTIONAL  :: TurbID            !< an ID for naming the tubine output file
   LOGICAL,                  INTENT(IN), OPTIONAL  :: OverrideAbortLev  !< whether or not we should override the abort error level (e.g., FAST.Farm)
   CHARACTER(*),             INTENT(IN), OPTIONAL  :: RootName          !< A CHARACTER string containing the root name of FAST output files, overriding normal naming convention
      ! Local variables

   INTEGER                      :: i                                    ! loop counter
   !CHARACTER(1024)              :: DirName                              ! A CHARACTER string containing the path of the current working directory


   LOGICAL                      :: OverrideAbortErrLev  
   CHARACTER(*), PARAMETER      :: RoutineName = "FAST_Init"
   
   INTEGER(IntKi)               :: ErrStat2
   CHARACTER(1024)              :: ErrMsg2
   
      ! Initialize some variables
   ErrStat = ErrID_None
   ErrMsg = ''
   
   IF (PRESENT(OverrideAbortLev)) THEN
      OverrideAbortErrLev = OverrideAbortLev
   ELSE
      OverrideAbortErrLev = .true.
   END IF
   

   
   !...............................................................................................................................
   ! Set the root name of the output files based on the input file name
   !...............................................................................................................................
   
   if (present(RootName)) then
      p%OutFileRoot = RootName
   else         
         ! Determine the root name of the primary file (will be used for output files)
      CALL GetRoot( InputFile, p%OutFileRoot )
      IF ( Cmpl4SFun )  p%OutFileRoot = TRIM( p%OutFileRoot )//'.SFunc'
      IF ( PRESENT(TurbID) ) THEN
         IF ( TurbID > 0 ) THEN
            p%OutFileRoot = TRIM( p%OutFileRoot )//'.T'//TRIM(Num2LStr(TurbID))
         END IF
      END IF
   
   end if
   
   
   !...............................................................................................................................
   ! Initialize the module name/date/version info:
   !...............................................................................................................................

   DO i=1,NumModules
      y_FAST%Module_Ver(i)%Date = 'unknown date'
      y_FAST%Module_Ver(i)%Ver  = 'unknown version'
   END DO       
   y_FAST%Module_Ver( Module_IfW    )%Name = 'InflowWind'
   y_FAST%Module_Ver( Module_OpFM   )%Name = 'OpenFOAM integration'
   y_FAST%Module_Ver( Module_ED     )%Name = 'ElastoDyn'
   y_FAST%Module_Ver( Module_BD     )%Name = 'BeamDyn'
   y_FAST%Module_Ver( Module_AD14   )%Name = 'AeroDyn14'
   y_FAST%Module_Ver( Module_AD     )%Name = 'AeroDyn'
   y_FAST%Module_Ver( Module_SrvD   )%Name = 'ServoDyn'
   y_FAST%Module_Ver( Module_HD     )%Name = 'HydroDyn'
   y_FAST%Module_Ver( Module_SD     )%Name = 'SubDyn'
   y_FAST%Module_Ver( Module_ExtPtfm)%Name = 'ExtPtfm_MCKF'
   y_FAST%Module_Ver( Module_MAP    )%Name = 'MAP'
   y_FAST%Module_Ver( Module_FEAM   )%Name = 'FEAMooring'
   y_FAST%Module_Ver( Module_MD     )%Name = 'MoorDyn'
   y_FAST%Module_Ver( Module_Orca   )%Name = 'OrcaFlexInterface'
   y_FAST%Module_Ver( Module_IceF   )%Name = 'IceFloe'
   y_FAST%Module_Ver( Module_IceD   )%Name = 'IceDyn'
         
   y_FAST%Module_Abrev( Module_IfW    ) = 'IfW'
   y_FAST%Module_Abrev( Module_OpFM   ) = 'OpFM'
   y_FAST%Module_Abrev( Module_ED     ) = 'ED'
   y_FAST%Module_Abrev( Module_BD     ) = 'BD'
   y_FAST%Module_Abrev( Module_AD14   ) = 'AD'
   y_FAST%Module_Abrev( Module_AD     ) = 'AD'
   y_FAST%Module_Abrev( Module_SrvD   ) = 'SrvD'
   y_FAST%Module_Abrev( Module_HD     ) = 'HD'
   y_FAST%Module_Abrev( Module_SD     ) = 'SD'
   y_FAST%Module_Abrev( Module_ExtPtfm) = 'ExtPtfm'
   y_FAST%Module_Abrev( Module_MAP    ) = 'MAP'
   y_FAST%Module_Abrev( Module_FEAM   ) = 'FEAM'
   y_FAST%Module_Abrev( Module_MD     ) = 'MD'
   y_FAST%Module_Abrev( Module_Orca   ) = 'Orca'
   y_FAST%Module_Abrev( Module_IceF   ) = 'IceF'
   y_FAST%Module_Abrev( Module_IceD   ) = 'IceD'   
   
   p%n_substeps = 1                                                ! number of substeps for between modules and global/FAST time
   p%BD_OutputSibling = .false.
   
   !...............................................................................................................................
   ! Read the primary file for the glue code:
   !...............................................................................................................................
   CALL FAST_ReadPrimaryFile( InputFile, p, OverrideAbortErrLev, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      
      ! overwrite TMax if necessary)
   IF (PRESENT(TMax)) THEN
      p%TMax = TMax
      !p%TMax = MAX( TMax, p%TMax )
   END IF
   
   IF ( ErrStat >= AbortErrLev ) RETURN


   p%KMax = 1                 ! after more checking, we may put this in the input file...
   !IF (p%CompIce == Module_IceF) p%KMax = 2
   p%SizeJac_Opt1 = 0  ! initialize this vector to zero; after we figure out what size the ED/SD/HD/BD meshes are, we'll fill this
   
   p%numIceLegs = 0           ! initialize number of support-structure legs in contact with ice (IceDyn will set this later)
   
   p%nBeams = 0               ! initialize number of BeamDyn instances (will be set later)
   
      ! determine what kind of turbine we're modeling:
   IF ( p%CompHydro == Module_HD ) THEN
      IF ( p%CompSub == Module_SD ) THEN
         p%TurbineType = Type_Offshore_Fixed
      ELSE
         p%TurbineType = Type_Offshore_Floating
      END IF
   ELSEIF ( p%CompMooring == Module_Orca ) THEN
      p%TurbineType = Type_Offshore_Floating
   !bjj: what about ExtPtfm_MCKF ???
   ELSE      
      p%TurbineType = Type_LandBased
   END IF   
         
    
   p%n_TMax_m1  = CEILING( ( (p%TMax - t_initial) / p%DT ) ) - 1 ! We're going to go from step 0 to n_TMax (thus the -1 here)

   if (p%TMax < 1.0_DbKi) then ! log10(0) gives floating point divide-by-zero error
      p%TChanLen = 10
   else
      p%TChanLen = max( 10, int(log10(p%TMax))+7 )
   end if
   p%OutFmt_t = 'F'//trim(num2lstr( p%TChanLen ))//'.4' ! 'F10.4'    
    
   !...............................................................................................................................
   ! Do some error checking on the inputs (validation):
   !...............................................................................................................................   
   call ValidateInputData(p, ErrStat2, ErrMsg2)    
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
    

   
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   
   RETURN
END SUBROUTINE FAST_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the output for the glue code, including writing the header for the primary output file.
SUBROUTINE FAST_InitOutput( p_FAST, y_FAST, InitOutData_ED, InitOutData_BD, InitOutData_SrvD, InitOutData_AD14, InitOutData_AD, &
                            InitOutData_IfW, InitOutData_OpFM, InitOutData_HD, InitOutData_SD, InitOutData_ExtPtfm, InitOutData_MAP, &
                            InitOutData_FEAM, InitOutData_MD, InitOutData_Orca, InitOutData_IceF, InitOutData_IceD, ErrStat, ErrMsg )

   IMPLICIT NONE

      ! Passed variables
   TYPE(FAST_ParameterType),       INTENT(IN)           :: p_FAST                                !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(INOUT)        :: y_FAST                                !< Glue-code simulation outputs

   TYPE(ED_InitOutputType),        INTENT(IN)           :: InitOutData_ED                        !< Initialization output for ElastoDyn
   TYPE(BD_InitOutputType),        INTENT(IN)           :: InitOutData_BD(:)                     !< Initialization output for BeamDyn (each instance)
   TYPE(SrvD_InitOutputType),      INTENT(IN)           :: InitOutData_SrvD                      !< Initialization output for ServoDyn
   TYPE(AD14_InitOutputType),      INTENT(IN)           :: InitOutData_AD14                      !< Initialization output for AeroDyn14
   TYPE(AD_InitOutputType),        INTENT(IN)           :: InitOutData_AD                        !< Initialization output for AeroDyn
   TYPE(InflowWind_InitOutputType),INTENT(IN)           :: InitOutData_IfW                       !< Initialization output for InflowWind
   TYPE(OpFM_InitOutputType),      INTENT(IN)           :: InitOutData_OpFM                      !< Initialization output for OpenFOAM
   TYPE(HydroDyn_InitOutputType),  INTENT(IN)           :: InitOutData_HD                        !< Initialization output for HydroDyn
   TYPE(SD_InitOutputType),        INTENT(IN)           :: InitOutData_SD                        !< Initialization output for SubDyn
   TYPE(ExtPtfm_InitOutputType),   INTENT(IN)           :: InitOutData_ExtPtfm                   !< Initialization output for ExtPtfm_MCKF
   TYPE(MAP_InitOutputType),       INTENT(IN)           :: InitOutData_MAP                       !< Initialization output for MAP
   TYPE(Orca_InitOutputType),      INTENT(IN)           :: InitOutData_Orca                      !< Initialization output for OrcaFlex interface
   TYPE(FEAM_InitOutputType),      INTENT(IN)           :: InitOutData_FEAM                      !< Initialization output for FEAMooring
   TYPE(MD_InitOutputType),        INTENT(IN)           :: InitOutData_MD                        !< Initialization output for MoorDyn
   TYPE(IceFloe_InitOutputType),   INTENT(IN)           :: InitOutData_IceF                      !< Initialization output for IceFloe
   TYPE(IceD_InitOutputType),      INTENT(IN)           :: InitOutData_IceD                      !< Initialization output for IceDyn

   INTEGER(IntKi),                 INTENT(OUT)          :: ErrStat                               !< Error status
   CHARACTER(*),                   INTENT(OUT)          :: ErrMsg                                !< Error message corresponding to ErrStat


      ! Local variables.

   INTEGER(IntKi)                   :: I, J                                            ! Generic index for DO loops.
   INTEGER(IntKi)                   :: indxLast                                        ! The index of the last value to be written to an array
   INTEGER(IntKi)                   :: indxNext                                        ! The index of the next value to be written to an array
   INTEGER(IntKi)                   :: NumOuts                                         ! number of channels to be written to the output file(s)



   !......................................................
   ! Set the description lines to be printed in the output file
   !......................................................
   y_FAST%FileDescLines(1)  = 'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//TRIM(GetVersion(FAST_Ver))
   y_FAST%FileDescLines(2)  = 'linked with ' //' '//TRIM(GetNVD(NWTC_Ver            ))  ! we'll get the rest of the linked modules in the section below
   y_FAST%FileDescLines(3)  = 'Description from the FAST input file: '//TRIM(p_FAST%FTitle)
   
   !......................................................
   ! We'll fill out the rest of FileDescLines(2), 
   ! and save the module version info for later use, too:
   !......................................................

   y_FAST%Module_Ver( Module_ED )   = InitOutData_ED%Ver
   y_FAST%FileDescLines(2)          = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_ED )  ))

   IF ( p_FAST%CompElast == Module_BD )  THEN
      y_FAST%Module_Ver( Module_BD ) = InitOutData_BD(1)%Ver ! call copy routine for this type if it every uses dynamic memory     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_BD ))) 
   END IF   
   
   
   IF ( p_FAST%CompInflow == Module_IfW )  THEN
      y_FAST%Module_Ver( Module_IfW ) = InitOutData_IfW%Ver ! call copy routine for this type if it every uses dynamic memory     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IfW ))) 
   ELSEIF ( p_FAST%CompInflow == Module_OpFM )  THEN
      y_FAST%Module_Ver( Module_OpFM ) = InitOutData_OpFM%Ver ! call copy routine for this type if it every uses dynamic memory     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_OpFM ))) 
   END IF   
   
   IF ( p_FAST%CompAero == Module_AD14 )  THEN
      y_FAST%Module_Ver( Module_AD14  ) = InitOutData_AD14%Ver     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_AD14  ) ))                  
   ELSEIF ( p_FAST%CompAero == Module_AD )  THEN
      y_FAST%Module_Ver( Module_AD  ) = InitOutData_AD%Ver     
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_AD  ) ))                  
   END IF

   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      y_FAST%Module_Ver( Module_SrvD ) = InitOutData_SrvD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SrvD )))
   END IF
         
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      y_FAST%Module_Ver( Module_HD )   = InitOutData_HD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_HD )))
   END IF

   IF ( p_FAST%CompSub == Module_SD ) THEN
      y_FAST%Module_Ver( Module_SD )   = InitOutData_SD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SD )))
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      y_FAST%Module_Ver( Module_ExtPtfm )   = InitOutData_ExtPtfm%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_ExtPtfm )))
   END IF

   IF ( p_FAST%CompMooring == Module_MAP ) THEN
      y_FAST%Module_Ver( Module_MAP )   = InitOutData_MAP%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_MAP )))
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      y_FAST%Module_Ver( Module_MD )   = InitOutData_MD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_MD )))
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      y_FAST%Module_Ver( Module_FEAM )   = InitOutData_FEAM%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_FEAM )))
   ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
      y_FAST%Module_Ver( Module_Orca )   = InitOutData_Orca%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_Orca)))
   END IF   
   
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      y_FAST%Module_Ver( Module_IceF )   = InitOutData_IceF%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IceF )))
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      y_FAST%Module_Ver( Module_IceD )   = InitOutData_IceD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IceD )))   
   END IF      
   
   !......................................................
   ! Set the number of output columns from each module
   !......................................................
   y_FAST%numOuts = 0    ! Inintialize entire array
   
   
   
   !y_FAST%numOuts(Module_InfW)  = 3  !hack for now: always output 3 wind speeds at hub-height
   IF ( ALLOCATED( InitOutData_IfW%WriteOutputHdr  ) ) y_FAST%numOuts(Module_IfW)  = SIZE(InitOutData_IfW%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_OpFM%WriteOutputHdr ) ) y_FAST%numOuts(Module_OpFM) = SIZE(InitOutData_OpFM%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_ED%WriteOutputHdr   ) ) y_FAST%numOuts(Module_ED)   = SIZE(InitOutData_ED%WriteOutputHdr)
do i=1,p_FAST%nBeams
   IF ( ALLOCATED( InitOutData_BD(i)%WriteOutputHdr) ) y_FAST%numOuts(Module_BD)   = y_FAST%numOuts(Module_BD) + SIZE(InitOutData_BD(i)%WriteOutputHdr)
end do   
                                                       y_FAST%numOuts(Module_AD14) = 0
                                                       
   IF ( ALLOCATED( InitOutData_AD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_AD)     = SIZE(InitOutData_AD%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_SrvD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_SrvD)   = SIZE(InitOutData_SrvD%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_HD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_HD)     = SIZE(InitOutData_HD%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_SD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_SD)     = SIZE(InitOutData_SD%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_ExtPtfm%WriteOutputHdr) ) y_FAST%numOuts(Module_ExtPtfm)= SIZE(InitOutData_ExtPtfm%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_MAP%WriteOutputHdr    ) ) y_FAST%numOuts(Module_MAP)    = SIZE(InitOutData_MAP%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_FEAM%WriteOutputHdr   ) ) y_FAST%numOuts(Module_FEAM)   = SIZE(InitOutData_FEAM%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_MD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_MD)     = SIZE(InitOutData_MD%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_Orca%WriteOutputHdr   ) ) y_FAST%numOuts(Module_Orca)   = SIZE(InitOutData_Orca%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_IceF%WriteOutputHdr   ) ) y_FAST%numOuts(Module_IceF)   = SIZE(InitOutData_IceF%WriteOutputHdr)
   IF ( ALLOCATED( InitOutData_IceD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_IceD)   = SIZE(InitOutData_IceD%WriteOutputHdr)*p_FAST%numIceLegs         
   
   !......................................................
   ! Initialize the output channel names and units
   !......................................................
   NumOuts   = 1 + SUM( y_FAST%numOuts )

   CALL AllocAry( y_FAST%ChannelNames,NumOuts, 'ChannelNames', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( y_FAST%ChannelUnits,NumOuts, 'ChannelUnits', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN

   y_FAST%ChannelNames(1) = 'Time'
   y_FAST%ChannelUnits(1) = '(s)'

   indxLast = 1
   indxNext = 2

   IF ( y_FAST%numOuts(Module_IfW) > 0_IntKi ) THEN  
      indxLast = indxNext + y_FAST%numOuts(Module_IfW) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_IfW%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_IfW%WriteOutputUnt      
      indxNext = indxLast + 1
   ELSEIF ( y_FAST%numOuts(Module_OpFM) > 0_IntKi ) THEN  
      indxLast = indxNext + y_FAST%numOuts(Module_OpFM) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_OpFM%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_OpFM%WriteOutputUnt      
      indxNext = indxLast + 1
   END IF


   IF ( y_FAST%numOuts(Module_ED) > 0_IntKi ) THEN !ElastoDyn
      indxLast = indxNext + y_FAST%numOuts(Module_ED) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_ED%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_ED%WriteOutputUnt
      indxNext = indxLast + 1
   END IF

   IF ( y_FAST%numOuts(Module_BD) > 0_IntKi ) THEN !BeamDyn
      do i=1,p_FAST%nBeams
         if ( allocated(InitOutData_BD(i)%WriteOutputHdr) ) then            
            do j=1,size(InitOutData_BD(i)%WriteOutputHdr) 
               y_FAST%ChannelNames(indxNext) = 'B'//TRIM(Num2Lstr(i))//trim(InitOutData_BD(i)%WriteOutputHdr(j))
               y_FAST%ChannelUnits(indxNext) = InitOutData_BD(i)%WriteOutputUnt(j)
               indxNext = indxNext + 1
            end do ! j            
         end if         
      end do                 
   END IF
   
   
      ! none for AeroDyn14 
   
   IF ( y_FAST%numOuts(Module_AD) > 0_IntKi ) THEN !AeroDyn
      indxLast = indxNext + y_FAST%numOuts(Module_AD) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_AD%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_AD%WriteOutputUnt
      indxNext = indxLast + 1
   END IF

   
   IF ( y_FAST%numOuts(Module_SrvD) > 0_IntKi ) THEN !ServoDyn
      indxLast = indxNext + y_FAST%numOuts(Module_SrvD) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_SrvD%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_SrvD%WriteOutputUnt
      indxNext = indxLast + 1
   END IF

   IF ( y_FAST%numOuts(Module_HD) > 0_IntKi ) THEN !HydroDyn
      indxLast = indxNext + y_FAST%numOuts(Module_HD) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_HD%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_HD%WriteOutputUnt
      indxNext = indxLast + 1
   END IF

   
   IF ( y_FAST%numOuts(Module_SD) > 0_IntKi ) THEN !SubDyn
      indxLast = indxNext + y_FAST%numOuts(Module_SD) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_SD%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_SD%WriteOutputUnt
      indxNext = indxLast + 1
   ELSEIF ( y_FAST%numOuts(Module_ExtPtfm) > 0_IntKi ) THEN !ExtPtfm_MCKF
      indxLast = indxNext + y_FAST%numOuts(Module_ExtPtfm) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_ExtPtfm%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_ExtPtfm%WriteOutputUnt
      indxNext = indxLast + 1
   END IF

   
   IF ( y_FAST%numOuts(Module_MAP) > 0_IntKi ) THEN !MAP
      indxLast = indxNext + y_FAST%numOuts(Module_MAP) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_MAP%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_MAP%WriteOutputUnt
      indxNext = indxLast + 1
   ELSEIF ( y_FAST%numOuts(Module_MD) > 0_IntKi ) THEN !MoorDyn
      indxLast = indxNext + y_FAST%numOuts(Module_MD) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_MD%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_MD%WriteOutputUnt
      indxNext = indxLast + 1
   ELSEIF ( y_FAST%numOuts(Module_FEAM) > 0_IntKi ) THEN !FEAMooring
      indxLast = indxNext + y_FAST%numOuts(Module_FEAM) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_FEAM%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_FEAM%WriteOutputUnt
      indxNext = indxLast + 1
   ELSEIF ( y_FAST%numOuts(Module_Orca) > 0_IntKi ) THEN !OrcaFlex
      indxLast = indxNext + y_FAST%numOuts(Module_Orca) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_Orca%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_Orca%WriteOutputUnt
      indxNext = indxLast + 1
   END IF
   
   
   IF ( y_FAST%numOuts(Module_IceF) > 0_IntKi ) THEN !IceFloe
      indxLast = indxNext + y_FAST%numOuts(Module_IceF) - 1
      y_FAST%ChannelNames(indxNext:indxLast) = InitOutData_IceF%WriteOutputHdr
      y_FAST%ChannelUnits(indxNext:indxLast) = InitOutData_IceF%WriteOutputUnt
      indxNext = indxLast + 1
   ELSEIF ( y_FAST%numOuts(Module_IceD) > 0_IntKi ) THEN !IceDyn
      DO I=1,p_FAST%numIceLegs         
         DO J=1,SIZE(InitOutData_IceD%WriteOutputHdr) 
            y_FAST%ChannelNames(indxNext) =TRIM(InitOutData_IceD%WriteOutputHdr(J))//'L'//TRIM(Num2Lstr(I))  !bjj: do we want this "Lx" at the end?
            y_FAST%ChannelUnits(indxNext) = InitOutData_IceD%WriteOutputUnt(J)
            indxNext = indxNext + 1
         END DO ! J
      END DO ! I
   END IF   
      
   
   !......................................................
   ! Open the text output file and print the headers
   !......................................................

   IF (p_FAST%WrTxtOutFile) THEN

      CALL GetNewUnit( y_FAST%UnOu, ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL OpenFOutFile ( y_FAST%UnOu, TRIM(p_FAST%OutFileRoot)//'.out', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

         ! Add some file information:

      WRITE (y_FAST%UnOu,'(/,A)')  TRIM( y_FAST%FileDescLines(1) )
      WRITE (y_FAST%UnOu,'(1X,A)') TRIM( y_FAST%FileDescLines(2) )
      WRITE (y_FAST%UnOu,'()' )    !print a blank line
      WRITE (y_FAST%UnOu,'(A)'   ) TRIM( y_FAST%FileDescLines(3) )
      WRITE (y_FAST%UnOu,'()' )    !print a blank line


         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................

      CALL WrFileNR ( y_FAST%UnOu, y_FAST%ChannelNames(1) )

      DO I=2,NumOuts
         CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//y_FAST%ChannelNames(I) )
      ENDDO ! I

      WRITE (y_FAST%UnOu,'()')

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................

      CALL WrFileNR ( y_FAST%UnOu, y_FAST%ChannelUnits(1) )

      DO I=2,NumOuts
         CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//y_FAST%ChannelUnits(I) )
      ENDDO ! I

      WRITE (y_FAST%UnOu,'()')

   END IF

   !......................................................
   ! Allocate data for binary output file
   !......................................................
   IF (p_FAST%WrBinOutFile) THEN

         ! calculate the size of the array of outputs we need to store
      y_FAST%NOutSteps = CEILING ( (p_FAST%TMax - p_FAST%TStart) / p_FAST%DT_OUT ) + 1

      CALL AllocAry( y_FAST%AllOutData, NumOuts-1, y_FAST%NOutSteps, 'AllOutData', ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

      IF ( p_FAST%WrBinMod == FileFmtID_WithTime ) THEN   ! we store the entire time array
         CALL AllocAry( y_FAST%TimeData, y_FAST%NOutSteps, 'TimeData', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN
      ELSE  
         CALL AllocAry( y_FAST%TimeData, 2_IntKi, 'TimeData', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

         y_FAST%TimeData(1) = 0.0_DbKi           ! This is the first output time, which we will set later
         y_FAST%TimeData(2) = p_FAST%DT_out      ! This is the (constant) time between subsequent writes to the output file
      END IF

      y_FAST%n_Out = 0  !number of steps actually written to the file

   END IF

   y_FAST%VTK_count = 0  ! first VTK file has 0 as output

RETURN
END SUBROUTINE FAST_InitOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets up the information needed to initialize AeroDyn, then initializes AeroDyn
SUBROUTINE AD_SetInitInput(InitInData_AD14, InitOutData_ED, y_ED, p_FAST, ErrStat, ErrMsg)

   ! Passed variables:
   TYPE(AD14_InitInputType),INTENT(INOUT) :: InitInData_AD14  !< The initialization input to AeroDyn14
   TYPE(ED_InitOutputType), INTENT(IN)    :: InitOutData_ED   !< The initialization output from structural dynamics module
   TYPE(ED_OutputType),     INTENT(IN)    :: y_ED             !< The outputs of the structural dynamics module (meshes with position/RefOrientation set)
   TYPE(FAST_ParameterType),INTENT(IN)    :: p_FAST           !< The parameters of the glue code
   INTEGER(IntKi)                         :: ErrStat          !< Error status of the operation
   CHARACTER(*)                           :: ErrMsg           !< Error message if ErrStat /= ErrID_None

      ! Local variables

   !TYPE(AD_InitOptions)       :: ADOptions                  ! Options for AeroDyn

   INTEGER                    :: K


   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! Set up the AeroDyn parameters
   InitInData_AD14%ADFileName   = p_FAST%AeroFile
   InitInData_AD14%OutRootName  = p_FAST%OutFileRoot
   InitInData_AD14%WrSumFile    = p_FAST%SumPrint      
   InitInData_AD14%NumBl        = InitOutData_ED%NumBl
   InitInData_AD14%UseDWM       = p_FAST%UseDWM
   
   InitInData_AD14%DWM%IfW%InputFileName   = p_FAST%InflowFile
   
      ! Hub position and orientation (relative here, but does not need to be)

   InitInData_AD14%TurbineComponents%Hub%Position(:)      = y_ED%HubPtMotion14%Position(:,1) - y_ED%HubPtMotion14%Position(:,1)  ! bjj: was 0; mesh was changed by adding p_ED%HubHt to 3rd component
   InitInData_AD14%TurbineComponents%Hub%Orientation(:,:) = y_ED%HubPtMotion14%RefOrientation(:,:,1)
   InitInData_AD14%TurbineComponents%Hub%TranslationVel   = 0.0_ReKi ! bjj: we don't need this field
   InitInData_AD14%TurbineComponents%Hub%RotationVel      = 0.0_ReKi ! bjj: we don't need this field

      ! Blade root position and orientation (relative here, but does not need to be)

   IF (.NOT. ALLOCATED( InitInData_AD14%TurbineComponents%Blade ) ) THEN
      ALLOCATE( InitInData_AD14%TurbineComponents%Blade( InitInData_AD14%NumBl ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_AD%TurbineComponents%Blade.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   DO K=1, InitInData_AD14%NumBl
      InitInData_AD14%TurbineComponents%Blade(K)%Position        = y_ED%BladeRootMotion14%Position(:,K)
      InitInData_AD14%TurbineComponents%Blade(K)%Orientation     = y_ED%BladeRootMotion14%RefOrientation(:,:,K)
      InitInData_AD14%TurbineComponents%Blade(K)%TranslationVel  = 0.0_ReKi ! bjj: we don't need this field
      InitInData_AD14%TurbineComponents%Blade(K)%RotationVel     = 0.0_ReKi ! bjj: we don't need this field      
   END DO
  

      ! Blade length
   IF (p_FAST%CompElast == Module_ED) THEN  ! note, we can't get here if we're using BeamDyn....
      InitInData_AD14%TurbineComponents%BladeLength = InitOutData_ED%BladeLength
   END IF
   
   
      ! Tower mesh ( here only because we currently need line2 meshes to contain the same nodes/elements )
      
   InitInData_AD14%NumTwrNodes = y_ED%TowerLn2Mesh%NNodes - 2
   IF (.NOT. ALLOCATED( InitInData_AD14%TwrNodeLocs ) ) THEN
      ALLOCATE( InitInData_AD14%TwrNodeLocs( 3, InitInData_AD14%NumTwrNodes ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_AD%TwrNodeLocs.'
         RETURN
      ELSE
         ErrStat = ErrID_None
      END IF
   END IF   
   
   IF ( InitInData_AD14%NumTwrNodes > 0 ) THEN
      InitInData_AD14%TwrNodeLocs = y_ED%TowerLn2Mesh%Position(:,1:InitInData_AD14%NumTwrNodes)  ! ED has extra nodes at beginning and top and bottom of tower
   END IF
   
      ! hub height         
   InitInData_AD14%HubHt = InitOutData_ED%HubHt
             

   RETURN
END SUBROUTINE AD_SetInitInput

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine is called at the start (or restart) of a FAST program (or FAST.Farm). It initializes the NWTC subroutine library,
!! displays the copyright notice, and displays some version information (including addressing scheme and precision).
SUBROUTINE FAST_ProgStart(ThisProgVer)
   TYPE(ProgDesc), INTENT(IN) :: ThisProgVer     !< program name/date/version description
   character(200) :: name, version
   character(200) :: git_commit, architecture, precision
   character(200) :: execution_date, execution_time, execution_zone
   
   ! ... Initialize NWTC Library (open console, set pi constants) ...
   ! sets the pi constants, open console for output, etc...
   CALL NWTC_Init( ProgNameIN=ThisProgVer%Name, EchoLibVer=.FALSE. )
   
   ! Display the copyright notice
   CALL DispCopyrightLicense( ThisProgVer )
   
   ! Display the program metadata
   call GetProgramMetadata(ThisProgVer, name, version, git_commit, architecture, precision)
   
   call wrscr(trim(name)//'-'//trim(git_commit))
   call wrscr('Compile Info:')
   call wrscr(' - Architecture: '//trim(architecture))
   call wrscr(' - Precision: '//trim(precision))
   call wrscr(' - Date: '//__DATE__)
   call wrscr(' - Time: '//__TIME__)
   ! use iso_fortran_env for compiler_version() and compiler_options()
   ! call wrscr(' - Compiler: '//trim(compiler_version()))
   ! call wrscr(' - Options: '//trim(compiler_options()))

   call date_and_time(execution_date, execution_time, execution_zone) 
   
   call wrscr('Execution Info:')
   call wrscr(' - Date: '//trim(execution_date(5:6)//'/'//execution_date(7:8)//'/'//execution_date(1:4)))
   call wrscr(' - Time: '//trim(execution_time(1:2)//':'//execution_time(3:4)//':'//execution_time(5:6))//trim(execution_zone))
   
   call wrscr('')
   
END SUBROUTINE FAST_ProgStart

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the number of subcycles (substeps) for modules at initialization, checking to make sure that their requested 
!! time step is valid.
SUBROUTINE SetModuleSubstepTime(ModuleID, p_FAST, y_FAST, ErrStat, ErrMsg)
   INTEGER(IntKi),           INTENT(IN   ) :: ModuleID            !< ID of the module to check time step and set
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

      
   ErrStat = ErrID_None
   ErrMsg  = "" 
   
   IF ( EqualRealNos( p_FAST%dt_module( ModuleID ), p_FAST%dt ) ) THEN
      p_FAST%n_substeps(ModuleID) = 1
   ELSE
      IF ( p_FAST%dt_module( ModuleID ) > p_FAST%dt ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = "The "//TRIM(y_FAST%Module_Ver(ModuleID)%Name)//" module time step ("//&
                          TRIM(Num2LStr(p_FAST%dt_module( ModuleID )))// &
                    " s) cannot be larger than FAST time step ("//TRIM(Num2LStr(p_FAST%dt))//" s)."
      ELSE
            ! calculate the number of subcycles:
         p_FAST%n_substeps(ModuleID) = NINT( p_FAST%dt / p_FAST%dt_module( ModuleID ) )
            
            ! let's make sure THE module DT is an exact integer divisor of the global (FAST) time step:
         IF ( .NOT. EqualRealNos( p_FAST%dt, p_FAST%dt_module( ModuleID ) * p_FAST%n_substeps(ModuleID) )  ) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = "The "//TRIM(y_FAST%Module_Ver(ModuleID)%Name)//" module time step ("//&
                              TRIM(Num2LStr(p_FAST%dt_module( ModuleID )))// &
                              " s) must be an integer divisor of the FAST time step ("//TRIM(Num2LStr(p_FAST%dt))//" s)."
         END IF
            
      END IF
   END IF      
                 
   RETURN
      
END SUBROUTINE SetModuleSubstepTime


END MODULE