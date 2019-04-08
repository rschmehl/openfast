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
! Unless required by applicable law or agreed to in viting, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************
MODULE FAST_VTK

   USE FAST_ModTypes

IMPLICIT NONE

CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------
!> This function builds the path for the vtk directory based on the output file root
FUNCTION get_vtkdir_path( out_file_root )
   CHARACTER(1024) :: get_vtkdir_path
   CHARACTER(*), INTENT(IN) :: out_file_root
   INTEGER(IntKi) :: last_separator_index
   
   ! get the directory of the primary input file (i.e. the case directory); Windows can have either forward or backward slashes (compare with GetPath())
   
   last_separator_index =      index(out_file_root, '/', back=.true.)
   last_separator_index = max( index(out_file_root, '\', back=.true.), last_separator_index )
   
   if (last_separator_index==0) then
      get_vtkdir_path = '.'//PathSep//'vtk'
   else
      get_vtkdir_path = trim(out_file_root(1 : last_separator_index) // 'vtk')
   end if
END FUNCTION
!----------------------------------------------------------------------------------------------------------------------------------
!> This function builds the path for the vtk root file name based on the output file root
FUNCTION get_vtkroot_path( out_file_root )
   CHARACTER(1024) :: get_vtkroot_path
   CHARACTER(*), INTENT(IN) :: out_file_root
   INTEGER(IntKi) :: last_separator_index
   INTEGER(IntKi) :: path_length

   last_separator_index =      index(out_file_root, '/', back=.true.)
   last_separator_index = max( index(out_file_root, '\', back=.true.), last_separator_index )

   get_vtkroot_path = trim( get_vtkdir_path(out_file_root) ) // PathSep &
                      // out_file_root( last_separator_index + 1 :)
END FUNCTION
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets up some of the information needed for plotting VTK surfaces. It initializes only the data needed before 
!! HD initialization. (HD needs some of this data so it can return the wave elevation data we want.)
SUBROUTINE SetVTKParameters_B4HD(p_FAST, InitOutData_ED, InitInData_HD, BD, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType),     INTENT(INOUT) :: p_FAST           !< The parameters of the glue code
   TYPE(ED_InitOutputType),      INTENT(IN   ) :: InitOutData_ED   !< The initialization output from structural dynamics module
   TYPE(HydroDyn_InitInputType), INTENT(INOUT) :: InitInData_HD    !< The initialization input to HydroDyn
   TYPE(BeamDyn_Data),           INTENT(IN   ) :: BD               !< BeamDyn data
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None

      
   REAL(SiKi)                              :: BladeLength, Width, WidthBy2
   REAL(SiKi)                              :: dx, dy                
   INTEGER(IntKi)                          :: i, j, n
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SetVTKParameters_B4HD'
   
         
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! Get radius for ground (blade length + hub radius):
   if ( p_FAST%CompElast == Module_BD ) then  
      BladeLength = TwoNorm(BD%y(1)%BldMotion%Position(:,1) - BD%y(1)%BldMotion%Position(:,BD%y(1)%BldMotion%Nnodes))
   else
      BladeLength = InitOutData_ED%BladeLength 
   end if
   p_FAST%VTK_Surface%GroundRad =  BladeLength + InitOutData_ED%HubRad 

   !........................................................................................................
   ! We don't use the rest of this routine for stick-figure output
   if (p_FAST%VTK_Type /= VTK_Surf) return  
   !........................................................................................................
      
      ! initialize wave elevation data:
   if ( p_FAST%CompHydro == Module_HD ) then
      
      p_FAST%VTK_surface%NWaveElevPts(1) = 25
      p_FAST%VTK_surface%NWaveElevPts(2) = 25
            
      call allocAry( InitInData_HD%WaveElevXY, 2, p_FAST%VTK_surface%NWaveElevPts(1)*p_FAST%VTK_surface%NWaveElevPts(2), 'WaveElevXY', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return

      Width = p_FAST%VTK_Surface%GroundRad * VTK_GroundFactor
      dx = Width / (p_FAST%VTK_surface%NWaveElevPts(1) - 1)
      dy = Width / (p_FAST%VTK_surface%NWaveElevPts(2) - 1)
            
      WidthBy2 = Width / 2.0_SiKi
      n = 1
      do i=1,p_FAST%VTK_surface%NWaveElevPts(1)
         do j=1,p_FAST%VTK_surface%NWaveElevPts(2)
            InitInData_HD%WaveElevXY(1,n) = dx*(i-1) - WidthBy2 !+ p_FAST%TurbinePos(1) ! HD takes p_FAST%TurbinePos into account already
            InitInData_HD%WaveElevXY(2,n) = dy*(j-1) - WidthBy2 !+ p_FAST%TurbinePos(2)
            n = n+1
         end do
      end do
      
   end if
         
      
END SUBROUTINE SetVTKParameters_B4HD
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets up the information needed for plotting VTK surfaces.
SUBROUTINE SetVTKParameters(p_FAST, InitOutData_ED, InitOutData_AD, InitInData_HD, InitOutData_HD, ED, BD, AD, HD, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType),     INTENT(INOUT) :: p_FAST           !< The parameters of the glue code
   TYPE(ED_InitOutputType),      INTENT(IN   ) :: InitOutData_ED   !< The initialization output from structural dynamics module
   TYPE(AD_InitOutputType),      INTENT(INOUT) :: InitOutData_AD   !< The initialization output from AeroDyn
   TYPE(HydroDyn_InitInputType), INTENT(INOUT) :: InitInData_HD    !< The initialization input to HydroDyn
   TYPE(HydroDyn_InitOutputType),INTENT(INOUT) :: InitOutData_HD   !< The initialization output from HydroDyn
   TYPE(ElastoDyn_Data),         INTENT(IN   ) :: ED               !< ElastoDyn data
   TYPE(BeamDyn_Data),           INTENT(IN   ) :: BD               !< BeamDyn data
   TYPE(AeroDyn_Data),           INTENT(IN   ) :: AD               !< AeroDyn data
   TYPE(HydroDyn_Data),          INTENT(IN   ) :: HD               !< HydroDyn data
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None

   REAL(SiKi)                              :: RefPoint(3), RefLengths(2)               
   REAL(SiKi)                              :: x, y                
   REAL(SiKi)                              :: TwrDiam_top, TwrDiam_base, TwrRatio, TwrLength
   INTEGER(IntKi)                          :: topNode, baseNode
   INTEGER(IntKi)                          :: tipNode, rootNode, cylNode
   INTEGER(IntKi)                          :: NumBl, k
   CHARACTER(1024)                         :: VTK_path
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SetVTKParameters'

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! create the VTK directory if it does not exist
   call MKDIR( get_vtkdir_path(p_FAST%OutFileRoot) )

   ! initialize the vtk data
   p_FAST%VTK_Surface%NumSectors = 18
   p_FAST%VTK_Surface%HubRad     = InitOutData_ED%HubRad
   ! NOTE: we set p_FAST%VTK_Surface%GroundRad in SetVTKParameters_B4HD

   ! write the ground or seabed reference polygon:
   VTK_path = get_vtkroot_path( p_FAST%OutFileRoot )
   RefPoint = p_FAST%TurbinePos
   if (p_FAST%CompHydro == MODULE_HD) then
      RefLengths = p_FAST%VTK_Surface%GroundRad*VTK_GroundFactor/2.0_SiKi
      
      ! note that p_FAST%TurbinePos(3) must be 0 for offshore turbines
      RefPoint(3) = p_FAST%TurbinePos(3) - InitOutData_HD%WtrDpth      
      call WrVTK_Ground ( RefPoint, RefLengths, trim(VTK_path) // '.SeabedSurface', ErrStat2, ErrMsg2 )   
      
      RefPoint(3) = p_FAST%TurbinePos(3) - InitOutData_HD%MSL2SWL    
      call WrVTK_Ground ( RefPoint, RefLengths, trim(VTK_path) // '.StillWaterSurface', ErrStat2, ErrMsg2 )       
   else
      RefLengths = p_FAST%VTK_Surface%GroundRad !array = scalar
      call WrVTK_Ground ( RefPoint, RefLengths, trim(VTK_path) // '.GroundSurface', ErrStat2, ErrMsg2 )         
   end if
   
   
   !........................................................................................................
   ! We don't use the rest of this routine for stick-figure output
   if (p_FAST%VTK_Type /= VTK_Surf) return  
   !........................................................................................................
            
      ! we're going to create a box using these dimensions
   y  =          ED%Output(1)%HubPtMotion%Position(3,  1) - ED%Output(1)%NacelleMotion%Position(3,  1)
   x  = TwoNorm( ED%Output(1)%HubPtMotion%Position(1:2,1) - ED%Output(1)%NacelleMotion%Position(1:2,1) ) - InitOutData_ED%HubRad
   
   p_FAST%VTK_Surface%NacelleBox(:,1) = (/ -x,  y, 0.0_SiKi /)
   p_FAST%VTK_Surface%NacelleBox(:,2) = (/  x,  y, 0.0_SiKi /) 
   p_FAST%VTK_Surface%NacelleBox(:,3) = (/  x, -y, 0.0_SiKi /)
   p_FAST%VTK_Surface%NacelleBox(:,4) = (/ -x, -y, 0.0_SiKi /) 
   p_FAST%VTK_Surface%NacelleBox(:,5) = (/ -x, -y, 2*y      /)
   p_FAST%VTK_Surface%NacelleBox(:,6) = (/  x, -y, 2*y      /) 
   p_FAST%VTK_Surface%NacelleBox(:,7) = (/  x,  y, 2*y      /)
   p_FAST%VTK_Surface%NacelleBox(:,8) = (/ -x,  y, 2*y      /) 
   
   !.......................
   ! tapered tower
   !.......................
      
   CALL AllocAry(p_FAST%VTK_Surface%TowerRad,ED%Output(1)%TowerLn2Mesh%NNodes,'VTK_Surface%TowerRad',ErrStat2,ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN
   
   topNode   = ED%Output(1)%TowerLn2Mesh%NNodes - 1
   baseNode  = ED%Output(1)%TowerLn2Mesh%NNodes  
   TwrLength = TwoNorm( ED%Output(1)%TowerLn2Mesh%position(:,topNode) - ED%Output(1)%TowerLn2Mesh%position(:,baseNode) ) ! this is the assumed length of the tower
   TwrRatio  = TwrLength / 87.6_SiKi  ! use ratio of the tower length to the length of the 5MW tower
   TwrDiam_top  = 3.87*TwrRatio
   TwrDiam_base = 6.0*TwrRatio
   
   TwrRatio = 0.5 * (TwrDiam_top - TwrDiam_base) / TwrLength
   do k=1,ED%Output(1)%TowerLn2Mesh%NNodes
      TwrLength = TwoNorm( ED%Output(1)%TowerLn2Mesh%position(:,k) - ED%Output(1)%TowerLn2Mesh%position(:,baseNode) ) 
      p_FAST%VTK_Surface%TowerRad(k) = 0.5*TwrDiam_Base + TwrRatio*TwrLength
   end do
   
   !.......................
   ! blade surfaces
   !.......................
   NumBl = SIZE(ED%Output(1)%BladeRootMotion,1)
   allocate(p_FAST%VTK_Surface%BladeShape(NumBl),stat=ErrStat2)
   if (errStat2/=0) then
      call setErrStat(ErrID_Fatal,'Error allocating VTK_Surface%BladeShape.',ErrStat,ErrMsg,RoutineName)
      return
   end if
            
   IF ( p_FAST%CompAero == Module_AD ) THEN  ! These meshes may have airfoil data associated with nodes...

      IF (ALLOCATED(InitOutData_AD%BladeShape)) THEN
         do k=1,NumBl   
            call move_alloc( InitOutData_AD%BladeShape(k)%AirfoilCoords, p_FAST%VTK_Surface%BladeShape(k)%AirfoilCoords )
         end do
      ELSE
#ifndef USE_DEFAULT_BLADE_SURFACE
         call setErrStat(ErrID_Fatal,'Cannot do surface visualization without airfoil coordinates defined in AeroDyn.',ErrStat,ErrMsg,RoutineName)
         return
      END IF
   ELSE
      call setErrStat(ErrID_Fatal,'Cannot do surface visualization without using AeroDyn.',ErrStat,ErrMsg,RoutineName)
      return
   END IF      
#else
      ! AD used without airfoil coordinates specified

         rootNode = 1
      
         DO K=1,NumBl   
            tipNode  = AD%Input(1)%BladeMotion(K)%NNodes
            cylNode  = min(3,AD%Input(1)%BladeMotion(K)%Nnodes)
         
            call SetVTKDefaultBladeParams(AD%Input(1)%BladeMotion(K), p_FAST%VTK_Surface%BladeShape(K), tipNode, rootNode, cylNode, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF (ErrStat >= AbortErrLev) RETURN
         END DO                           
      END IF
      
   ELSE IF ( p_FAST%CompElast == Module_BD ) THEN
      rootNode = 1      
      DO K=1,NumBl   
         tipNode  = BD%y(k)%BldMotion%NNodes
         cylNode  = min(3,BD%y(k)%BldMotion%NNodes)
         
         call SetVTKDefaultBladeParams(BD%y(k)%BldMotion, p_FAST%VTK_Surface%BladeShape(K), tipNode, rootNode, cylNode, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      END DO      
   ELSE
      DO K=1,NumBl   
         rootNode = ED%Output(1)%BladeLn2Mesh(K)%NNodes     
         tipNode  = ED%Output(1)%BladeLn2Mesh(K)%NNodes-1
         cylNode  = min(2,ED%Output(1)%BladeLn2Mesh(K)%NNodes)
         
         call SetVTKDefaultBladeParams(ED%Output(1)%BladeLn2Mesh(K), p_FAST%VTK_Surface%BladeShape(K), tipNode, rootNode, cylNode, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      END DO  
   END IF   
#endif 
   
   
   !.......................
   ! wave elevation 
   !.......................

   !bjj: interpolate here instead of each time step?
   if ( allocated(InitOutData_HD%WaveElevSeries) ) then
      call move_alloc( InitInData_HD%WaveElevXY, p_FAST%VTK_Surface%WaveElevXY )
      call move_alloc( InitOutData_HD%WaveElevSeries, p_FAST%VTK_Surface%WaveElev )
      
         ! put the following lines in loops to avoid stack-size issues:
      do k=1,size(p_FAST%VTK_Surface%WaveElevXY,2)
         p_FAST%VTK_Surface%WaveElevXY(:,k) = p_FAST%VTK_Surface%WaveElevXY(:,k) + p_FAST%TurbinePos(1:2)
      end do
         
      !do k=1,size(p_FAST%VTK_Surface%WaveElev,2)
      !   p_FAST%VTK_Surface%WaveElev(:,k) = p_FAST%VTK_Surface%WaveElev(:,k) + p_FAST%TurbinePos(3)  ! not sure this is really accurate if p_FAST%TurbinePos(3) is non-zero
      !end do
      
   end if
   
   !.......................
   ! morison surfaces
   !.......................
   
   IF ( HD%Input(1)%Morison%DistribMesh%Committed ) THEN      
      
      call move_alloc(InitOutData_HD%Morison%Morison_Rad, p_FAST%VTK_Surface%MorisonRad)
      
   END IF
   
END SUBROUTINE SetVTKParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine comes up with some default airfoils for blade surfaces for a given blade mesh, M.
SUBROUTINE SetVTKDefaultBladeParams(M, BladeShape, tipNode, rootNode, cylNode, ErrStat, ErrMsg)

   TYPE(MeshType),               INTENT(IN   ) :: M                !< The Mesh the defaults should be calculated for
   TYPE(FAST_VTK_BLSurfaceType), INTENT(INOUT) :: BladeShape       !< BladeShape to set to default values
   INTEGER(IntKi),               INTENT(IN   ) :: rootNode         !< Index of root node (innermost node) for this mesh
   INTEGER(IntKi),               INTENT(IN   ) :: tipNode          !< Index of tip node (outermost node) for this mesh
   INTEGER(IntKi),               INTENT(IN   ) :: cylNode          !< Index of last node to have a cylinder shape
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None

      
   REAL(SiKi)                                  :: bladeLength, chord, pitchAxis
   REAL(SiKi)                                  :: bladeLengthFract, bladeLengthFract2, ratio, posLength ! temporary quantities               
   REAL(SiKi)                                  :: cylinderLength, x, y, angle               
   INTEGER(IntKi)                              :: i, j
   INTEGER(IntKi)                              :: ErrStat2
   CHARACTER(ErrMsgLen)                        :: ErrMsg2
   CHARACTER(*), PARAMETER                     :: RoutineName = 'SetVTKDefaultBladeParams'
   
   !Note: jmj does not like this default option

   integer, parameter :: N = 66
   
   ! default airfoil shape coordinates; uses S809 values from http://wind.nrel.gov/airfoils/Shapes/S809_Shape.html:   
   real, parameter, dimension(N) :: xc=(/ 1.0,0.996203,0.98519,0.967844,0.945073,0.917488,0.885293,0.848455,0.80747,0.763042,0.715952,0.667064,0.617331,0.56783,0.519832,0.474243,0.428461,0.382612,0.33726,0.29297,0.250247,0.209576,0.171409,0.136174,0.104263,0.076035,0.051823,0.03191,0.01659,0.006026,0.000658,0.000204,0.0,0.000213,0.001045,0.001208,0.002398,0.009313,0.02323,0.04232,0.065877,0.093426,0.124111,0.157653,0.193738,0.231914,0.271438,0.311968,0.35337,0.395329,0.438273,0.48192,0.527928,0.576211,0.626092,0.676744,0.727211,0.776432,0.823285,0.86663,0.905365,0.938474,0.965086,0.984478,0.996141,1.0 /)
   real, parameter, dimension(N) :: yc=(/ 0.0,0.000487,0.002373,0.00596,0.011024,0.017033,0.023458,0.03028,0.037766,0.045974,0.054872,0.064353,0.074214,0.084095,0.093268,0.099392,0.10176,0.10184,0.10007,0.096703,0.091908,0.085851,0.078687,0.07058,0.061697,0.052224,0.042352,0.032299,0.02229,0.012615,0.003723,0.001942,-0.00002,-0.001794,-0.003477,-0.003724,-0.005266,-0.011499,-0.020399,-0.030269,-0.040821,-0.051923,-0.063082,-0.07373,-0.083567,-0.092442,-0.099905,-0.105281,-0.108181,-0.108011,-0.104552,-0.097347,-0.086571,-0.073979,-0.060644,-0.047441,-0.0351,-0.024204,-0.015163,-0.008204,-0.003363,-0.000487,0.000743,0.000775,0.00029,0.0 /)

   call AllocAry(BladeShape%AirfoilCoords, 2, N, M%NNodes, 'BladeShape%AirfoilCoords', ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN
         
   ! Chord length and pitch axis location are given by scaling law
   bladeLength       = TwoNorm( M%position(:,tipNode) - M%Position(:,rootNode) )
   cylinderLength    = TwoNorm( M%Position(:,cylNode) - M%Position(:,rootNode) )
   bladeLengthFract  = 0.22*bladeLength
   bladeLengthFract2 = bladeLength-bladeLengthFract != 0.78*bladeLength
   
   DO i=1,M%Nnodes
      posLength = TwoNorm( M%Position(:,i) - M%Position(:,rootNode) )
         
      IF (posLength .LE. bladeLengthFract) THEN
         ratio     = posLength/bladeLengthFract
         chord     =  (0.06 + 0.02*ratio)*bladeLength
         pitchAxis =   0.25 + 0.125*ratio
      ELSE
         chord     = (0.08 - 0.06*(posLength-bladeLengthFract)/bladeLengthFract2)*bladeLength
         pitchAxis = 0.375
      END IF
         
      IF (posLength .LE. cylinderLength) THEN 
         ! create a cylinder for this node
         
         chord = chord/2.0_SiKi
         
         DO j=1,N
            ! normalized x,y coordinates for airfoil
            x = yc(j)
            y = xc(j) - 0.5
                     
            angle = ATAN2( y, x)
         
               ! x,y coordinates for cylinder
            BladeShape%AirfoilCoords(1,j,i) = chord*COS(angle) ! x (note that "chord" is really representing chord/2 here)
            BladeShape%AirfoilCoords(2,j,i) = chord*SIN(angle) ! y (note that "chord" is really representing chord/2 here)
         END DO                                                     
         
      ELSE
         ! create an airfoil for this node
            
         DO j=1,N                  
            ! normalized x,y coordinates for airfoil, assuming an upwind turbine
            x = yc(j)
            y = xc(j) - pitchAxis
                  
               ! x,y coordinates for airfoil
            BladeShape%AirfoilCoords(1,j,i) =  chord*x
            BladeShape%AirfoilCoords(2,j,i) =  chord*y                        
         END DO
         
      END IF
      
   END DO ! nodes on mesh
         
END SUBROUTINE SetVTKDefaultBladeParams
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes the ground or seabed reference surface information in VTK format.
!! see VTK file information format for XML, here: http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
SUBROUTINE WrVTK_Ground ( RefPoint, HalfLengths, FileRootName, ErrStat, ErrMsg )
      
   REAL(SiKi),      INTENT(IN)           :: RefPoint(3)     !< reference point (plane will be created around it)
   REAL(SiKi),      INTENT(IN)           :: HalfLengths(2)  !< half of the X-Y lengths of plane surrounding RefPoint
   CHARACTER(*),    INTENT(IN)           :: FileRootName    !< Name of the file to write the output in (excluding extension)
   
   INTEGER(IntKi),  INTENT(OUT)          :: ErrStat         !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),    INTENT(OUT)          :: ErrMsg          !< Error message associated with the ErrStat


   ! local variables
   INTEGER(IntKi)                        :: Un            ! fortran unit number
   INTEGER(IntKi)                        :: ix            ! loop counters
   CHARACTER(1024)                       :: FileName
   INTEGER(IntKi), parameter             :: NumberOfPoints = 4
   INTEGER(IntKi), parameter             :: NumberOfLines = 0
   INTEGER(IntKi), parameter             :: NumberOfPolys = 1
        
   INTEGER(IntKi)                        :: ErrStat2 
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*),PARAMETER                :: RoutineName = 'WrVTK_Ground'
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !.................................................................
   ! write the data that potentially changes each time step:
   !.................................................................
      
   ! PolyData (.vtp) - Serial vtkPolyData (unstructured) file
   FileName = TRIM(FileRootName)//'.vtp'
      
   call WrVTK_header( FileName, NumberOfPoints, NumberOfLines, NumberOfPolys, Un, ErrStat2, ErrMsg2 )    
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
         
! points (nodes, augmented with NumSegments):   
      WRITE(Un,'(A)')         '      <Points>'
      WRITE(Un,'(A)')         '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
               
      WRITE(Un,VTK_AryFmt) RefPoint(1) + HalfLengths(1) , RefPoint(2) + HalfLengths(2), RefPoint(3)
      WRITE(Un,VTK_AryFmt) RefPoint(1) + HalfLengths(1) , RefPoint(2) - HalfLengths(2), RefPoint(3)
      WRITE(Un,VTK_AryFmt) RefPoint(1) - HalfLengths(1) , RefPoint(2) - HalfLengths(2), RefPoint(3)
      WRITE(Un,VTK_AryFmt) RefPoint(1) - HalfLengths(1) , RefPoint(2) + HalfLengths(2), RefPoint(3)
            
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Points>'
  
                  
      WRITE(Un,'(A)')         '      <Polys>'      
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'         
      WRITE(Un,'('//trim(num2lstr(NumberOfPoints))//'(i7))') (ix, ix=0,NumberOfPoints-1)                   
      WRITE(Un,'(A)')         '        </DataArray>'      
      
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="offsets" format="ascii">'            
      WRITE(Un,'(i7)') NumberOfPoints
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Polys>'      
            
      call WrVTK_footer( Un )       
                     
END SUBROUTINE WrVTK_Ground
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes all the committed meshes to VTK-formatted files. It doesn't bother with returning an error code.
SUBROUTINE WrVTK_AllMeshes(p_FAST, y_FAST, MeshMapData, ED, BD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_ModuleMapType), INTENT(IN   ) :: MeshMapData         !< Data for mapping between modules

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(IN   ) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(IN   ) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(IN   ) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop


   logical                                 :: outputFields        ! flag to determine if we want to output the HD mesh fields
   INTEGER(IntKi)                          :: NumBl, k, Twidth
   CHARACTER(1024)                         :: VTK_path
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_AllMeshes'

   ! Calculate the number of digits for the maximum number of output steps to be written.
   ! This will be used to pad the write-out step in the VTK filename with zeros in calls to MeshWrVTK()
   if ( (p_FAST%n_VTKTime>0) .and. (p_FAST%n_TMax_m1+1>0) ) then
      Twidth = CEILING( log10( real(p_FAST%n_TMax_m1+1, ReKi) / p_FAST%n_VTKTime ) ) + 1
   else
      Twidth = 1
   endif
   
   NumBl = 0
   if (allocated(ED%Output)) then
      if (allocated(ED%Output(1)%BladeRootMotion)) then
         NumBl = SIZE(ED%Output(1)%BladeRootMotion)      
      end if
   end if

   VTK_path = get_vtkroot_path( p_FAST%OutFileRoot )
   
   
! I'm first going to just put all of the meshes that get mapped together, then decide if we're going to print/plot them all
         
!  ElastoDyn
   if (allocated(ED%Output) .and. allocated(ED%Input)) then
   
         !  ElastoDyn outputs (motions)
      DO K=1,NumBl        
         !%BladeLn2Mesh(K) used only when not BD (see below)
         call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%BladeRootMotion(K), trim(VTK_path)//'.ED_BladeRootMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth )
      END DO
      
      call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%TowerLn2Mesh, trim(VTK_path)//'.ED_TowerLn2Mesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth )     

! these will get output with their sibling input meshes
      !call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%HubPtMotion, trim(p_FAST%OutFileRoot)//'.ED_HubPtMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
      !call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%NacelleMotion, trim(p_FAST%OutFileRoot)//'.ED_NacelleMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
      !call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%PlatformPtMesh, trim(p_FAST%OutFileRoot)//'.ED_PlatformPtMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
      
         !  ElastoDyn inputs (loads)
      ! %BladePtLoads used only when not BD (see below)
      call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%TowerPtLoads, trim(VTK_path)//'.ED_TowerPtLoads', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, ED%Output(1)%TowerLn2Mesh )     
      call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%HubPtLoad, trim(VTK_path)//'.ED_Hub', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, ED%Output(1)%HubPtMotion )
      call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%NacelleLoads, trim(VTK_path)//'.ED_Nacelle' ,y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, ED%Output(1)%NacelleMotion )     
      call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%PlatformPtMesh, trim(VTK_path)//'.ED_PlatformPtMesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, ED%Output(1)%PlatformPtMesh )     
   end if
   
   
!  BeamDyn
   IF ( p_FAST%CompElast == Module_BD .and. allocated(BD%Input) .and. allocated(BD%y)) THEN
            
      do K=1,NumBl        
            ! BeamDyn inputs
         !call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%RootMotion, trim(p_FAST%OutFileRoot)//'.BD_RootMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )
         call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%HubMotion, trim(VTK_path)//'.BD_HubMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth )    
      end do
      if (allocated(MeshMapData%y_BD_BldMotion_4Loads)) then
         do K=1,NumBl 
            call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%DistrLoad, trim(VTK_path)//'.BD_DistrLoad'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, MeshMapData%y_BD_BldMotion_4Loads(k) )
            ! skipping PointLoad
         end do
      elseif (p_FAST%BD_OutputSibling) then
         do K=1,NumBl
            call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%DistrLoad, trim(p_FAST%OutFileRoot)//'.BD_Blade'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, BD%y(k)%BldMotion )
            ! skipping PointLoad
         end do
      end if
      
      do K=1,NumBl
            ! BeamDyn outputs
         call MeshWrVTK(p_FAST%TurbinePos, BD%y(k)%ReactionForce, trim(p_FAST%OutFileRoot)//'.BD_ReactionForce_RootMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, BD%Input(1,k)%RootMotion )
      end do  
      
      if (.not. p_FAST%BD_OutputSibling) then !otherwise this mesh has been put with the DistrLoad mesh
         do K=1,NumBl
               ! BeamDyn outputs
         call MeshWrVTK(p_FAST%TurbinePos, BD%y(k)%BldMotion, trim(p_FAST%OutFileRoot)//'.BD_BldMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth )
      end do  
      end if
      
      
   ELSE if (allocated(ED%Input) .and. allocated(ED%Output)) then
      ! ElastoDyn
      DO K=1,NumBl        
         call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%BladeLn2Mesh(K), trim(VTK_path)//'.ED_BladeLn2Mesh_motion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth )
         call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%BladePtLoads(K), trim(VTK_path)//'.ED_BladePtLoads'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, ED%Output(1)%BladeLn2Mesh(K) )
      END DO      
   END IF
            
!  ServoDyn
   if (allocated(SrvD%Input)) then
      IF ( SrvD%Input(1)%NTMD%Mesh%Committed ) THEN         
         !call MeshWrVTK(p_FAST%TurbinePos, SrvD%Input(1)%NTMD%Mesh, trim(p_FAST%OutFileRoot)//'.SrvD_NTMD_Motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
         call MeshWrVTK(p_FAST%TurbinePos, SrvD%y%NTMD%Mesh, trim(VTK_path)//'.SrvD_NTMD', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, SrvD%Input(1)%TTMD%Mesh )                
      END IF      
      IF ( SrvD%Input(1)%TTMD%Mesh%Committed ) THEN 
         !call MeshWrVTK(p_FAST%TurbinePos, SrvD%Input(1)%TTMD%Mesh, trim(p_FAST%OutFileRoot)//'.SrvD_TTMD_Motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
         call MeshWrVTK(p_FAST%TurbinePos, SrvD%y%TTMD%Mesh, trim(VTK_path)//'.SrvD_TTMD', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, SrvD%Input(1)%TTMD%Mesh )         
      END IF   
   end if
   
      
!  AeroDyn   
   IF ( p_FAST%CompAero == Module_AD .and. allocated(AD%Input)) THEN 
               
      if (allocated(AD%Input(1)%BladeRootMotion)) then      
      
         DO K=1,NumBl   
            call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%BladeRootMotion(K), trim(VTK_path)//'.AD_BladeRootMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth )     
            !call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%BladeMotion(K), trim(p_FAST%OutFileRoot)//'.AD_BladeMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
         END DO            
         call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%HubMotion, trim(VTK_path)//'.AD_HubMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth )     
         !call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%TowerMotion, trim(p_FAST%OutFileRoot)//'.AD_TowerMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
               
         DO K=1,NumBl   
            call MeshWrVTK(p_FAST%TurbinePos, AD%y%BladeLoad(K), trim(VTK_path)//'.AD_Blade'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, AD%Input(1)%BladeMotion(k) )     
         END DO            
         call MeshWrVTK(p_FAST%TurbinePos, AD%y%TowerLoad, trim(VTK_path)//'.AD_Tower', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, AD%Input(1)%TowerMotion )     
         
      end if
      
   END IF
   
! HydroDyn            
   IF ( p_FAST%CompHydro == Module_HD .and. allocated(HD%Input)) THEN       
      !call MeshWrVTK(p_FAST%TurbinePos, HD%Input(1)%Mesh, trim(p_FAST%OutFileRoot)//'.HD_Mesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
      !call MeshWrVTK(p_FAST%TurbinePos, HD%Input(1)%Morison%LumpedMesh, trim(p_FAST%OutFileRoot)//'.HD_MorisonLumped_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
      !call MeshWrVTK(p_FAST%TurbinePos, HD%Input(1)%Morison%DistribMesh, trim(p_FAST%OutFileRoot)//'.HD_MorisonDistrib_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
      
      if (p_FAST%CompSub == Module_NONE) then
         call MeshWrVTK(p_FAST%TurbinePos, HD%y%AllHdroOrigin, trim(VTK_path)//'.HD_AllHdroOrigin', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, HD%Input(1)%Mesh )
         outputFields = .false.
      else         
         call MeshWrVTK(p_FAST%TurbinePos, HD%y%Mesh, trim(VTK_path)//'.HD_Mesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, HD%Input(1)%Mesh )
         outputFields = p_FAST%VTK_fields
      end if
      call MeshWrVTK(p_FAST%TurbinePos, HD%y%Morison%LumpedMesh, trim(VTK_path)//'.HD_MorisonLumped', y_FAST%VTK_count, outputFields, ErrStat2, ErrMsg2, Twidth, HD%Input(1)%Morison%LumpedMesh )     
      call MeshWrVTK(p_FAST%TurbinePos, HD%y%Morison%DistribMesh, trim(VTK_path)//'.HD_MorisonDistrib', y_FAST%VTK_count, outputFields, ErrStat2, ErrMsg2, Twidth, HD%Input(1)%Morison%DistribMesh )     
      
                  
   END IF
   
! SubDyn   
   IF ( p_FAST%CompSub == Module_SD .and. allocated(SD%Input)) THEN
      !call MeshWrVTK(p_FAST%TurbinePos, SD%Input(1)%TPMesh, trim(p_FAST%OutFileRoot)//'.SD_TPMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
      call MeshWrVTK(p_FAST%TurbinePos, SD%Input(1)%LMesh, trim(VTK_path)//'.SD_LMesh_y2Mesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, SD%y%y2Mesh )     
      
      call MeshWrVTK(p_FAST%TurbinePos, SD%y%y1Mesh, trim(VTK_path)//'.SD_y1Mesh_TPMesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, SD%Input(1)%TPMesh )     
      !call MeshWrVTK(p_FAST%TurbinePos, SD%y%y2Mesh, trim(p_FAST%OutFileRoot)//'.SD_y2Mesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )        
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm .and. allocated(ExtPtfm%Input)) THEN
      call MeshWrVTK(p_FAST%TurbinePos, ExtPtfm%y%PtfmMesh, trim(VTK_path)//'.ExtPtfm', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, ExtPtfm%Input(1)%PtfmMesh )     
   END IF     
       
! MAP
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
      if (allocated(MAPp%Input)) then
         call MeshWrVTK(p_FAST%TurbinePos, MAPp%y%PtFairleadLoad, trim(VTK_path)//'.MAP_PtFairlead', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, MAPp%Input(1)%PtFairDisplacement )     
         !call MeshWrVTK(p_FAST%TurbinePos, MAPp%Input(1)%PtFairDisplacement, trim(p_FAST%OutFileRoot)//'.MAP_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )        
      end if
      
! MoorDyn      
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      if (allocated(MD%Input)) then
         call MeshWrVTK(p_FAST%TurbinePos, MD%y%PtFairleadLoad, trim(VTK_path)//'.MD_PtFairlead', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, MD%Input(1)%PtFairleadDisplacement )     
         !call MeshWrVTK(p_FAST%TurbinePos, MD%Input(1)%PtFairleadDisplacement, trim(p_FAST%OutFileRoot)//'.MD_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )        
      end if
      
! FEAMooring                   
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      if (allocated(FEAM%Input)) then
         call MeshWrVTK(p_FAST%TurbinePos, FEAM%y%PtFairleadLoad, trim(VTK_path)//'.FEAM_PtFairlead', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, FEAM%Input(1)%PtFairleadDisplacement )     
         !call MeshWrVTK(p_FAST%TurbinePos, FEAM%Input(1)%PtFairleadDisplacement, trim(p_FAST%OutFileRoot)//'.FEAM_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )        
      end if
      
! Orca      
   ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
      if (allocated(Orca%Input)) then
         call MeshWrVTK(p_FAST%TurbinePos, Orca%y%PtfmMesh, trim(VTK_path)//'.Orca_PtfmMesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, Orca%Input(1)%PtfmMesh )     
         !call MeshWrVTK(p_FAST%TurbinePos, Orca%Input(1)%PtfmMesh, trim(p_FAST%OutFileRoot)//'.Orca_PtfmMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )
      end if
   END IF
            
         
! IceFloe      
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      if (allocated(IceF%Input)) then
         call MeshWrVTK(p_FAST%TurbinePos, IceF%y%iceMesh, trim(VTK_path)//'.IceF_iceMesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, IceF%Input(1)%iceMesh )     
         !call MeshWrVTK(p_FAST%TurbinePos, IceF%Input(1)%iceMesh, trim(p_FAST%OutFileRoot)//'.IceF_iceMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )
      end if
      
! IceDyn
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      if (allocated(IceD%Input)) then
            
         DO k = 1,p_FAST%numIceLegs
            call MeshWrVTK(p_FAST%TurbinePos, IceD%y(k)%PointMesh, trim(VTK_path)//'.IceD_PointMesh'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, IceD%Input(1,k)%PointMesh )     
            !call MeshWrVTK(p_FAST%TurbinePos, IceD%Input(1,k)%PointMesh, trim(p_FAST%OutFileRoot)//'.IceD_PointMesh_motion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )
         END DO
      end if
      
   END IF
   
   
END SUBROUTINE WrVTK_AllMeshes 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes a minimal subset of meshes (enough to visualize the turbine) to VTK-formatted files. It doesn't bother with 
!! returning an error code.
SUBROUTINE WrVTK_BasicMeshes(p_FAST, y_FAST, MeshMapData, ED, BD, AD14, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_ModuleMapType), INTENT(IN   ) :: MeshMapData         !< Data for mapping between modules

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(IN   ) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(IN   ) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop

   logical                                 :: OutputFields
   INTEGER(IntKi)                          :: NumBl, k, Twidth
   CHARACTER(1024)                         :: VTK_path
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_BasicMeshes'

   ! Calculate the number of digits for the maximum number of output steps to be written.
   ! This will be used to pad the write-out step in the VTK filename with zeros in calls to MeshWrVTK()
   if ( (p_FAST%n_VTKTime>0) .and. (p_FAST%n_TMax_m1+1>0) ) then
      Twidth = CEILING( log10( real(p_FAST%n_TMax_m1+1, ReKi) / p_FAST%n_VTKTime ) ) + 1
   else
      Twidth = 1
   endif
   
   
   NumBl = 0
   if (allocated(ED%Output(1)%BladeRootMotion)) then
      NumBl = SIZE(ED%Output(1)%BladeRootMotion)    
   end if

   VTK_path = get_vtkroot_path( p_FAST%OutFileRoot )

! Nacelle
   call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%NacelleMotion, trim(VTK_path)//'.ED_Nacelle', y_FAST%VTK_count, &
                  p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, Sib=ED%Input(1)%NacelleLoads )     
               
! Hub
   call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%HubPtMotion, trim(VTK_path)//'.ED_Hub', y_FAST%VTK_count, &
                  p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, Sib=ED%Input(1)%HubPtLoad )     
   
! Blades
   IF ( p_FAST%CompAero == Module_AD ) THEN  ! These meshes may have airfoil data associated with nodes...
      DO K=1,NumBl   
         call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%BladeMotion(K), trim(VTK_path)//'.AD_Blade'//trim(num2lstr(k)), &
                        y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, Sib=AD%y%BladeLoad(K) )     
      END DO                  
   ELSE IF ( p_FAST%CompElast == Module_BD ) THEN
      DO K=1,NumBl                 
         call MeshWrVTK(p_FAST%TurbinePos, BD%y(k)%BldMotion, trim(VTK_path)//'.BD_BldMotion'//trim(num2lstr(k)), &
                        y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth )         
      END DO  
   ELSE
      DO K=1,NumBl        
         call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%BladeLn2Mesh(K), trim(VTK_path)//'.ED_BladeLn2Mesh_motion'//trim(num2lstr(k)), &
                        y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth )
      END DO  
   END IF   
         
! Tower motions
   call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%TowerLn2Mesh, trim(VTK_path)//'.ED_TowerLn2Mesh_motion', &
                  y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth )     
   
   
! Substructure   
!   call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%PlatformPtMesh, trim(p_FAST%OutFileRoot)//'.ED_PlatformPtMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
!   IF ( p_FAST%CompSub == Module_SD ) THEN
!     call MeshWrVTK(p_FAST%TurbinePos, SD%Input(1)%TPMesh, trim(p_FAST%OutFileRoot)//'.SD_TPMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )     
!      call MeshWrVTK(p_FAST%TurbinePos, SD%y%y2Mesh, trim(p_FAST%OutFileRoot)//'.SD_y2Mesh_motion', y_FAST%VTK_count, ErrStat2, ErrMsg2 )        
!   END IF     
      
   IF ( p_FAST%CompHydro == Module_HD ) THEN 
      
      if (p_FAST%CompSub == Module_NONE) then
         call MeshWrVTK(p_FAST%TurbinePos, HD%y%AllHdroOrigin, trim(VTK_path)//'.HD_AllHdroOrigin', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, Twidth, HD%Input(1)%Mesh )
         outputFields = .false.
      else         
         OutputFields = p_FAST%VTK_fields
      end if
      
      call MeshWrVTK(p_FAST%TurbinePos, HD%Input(1)%Morison%DistribMesh, trim(VTK_path)//'.HD_MorisonDistrib', &
                     y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, Twidth, Sib=HD%y%Morison%DistribMesh )           
   END IF
   
   
! Mooring Lines?            
!   IF ( p_FAST%CompMooring == Module_MAP ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, MAPp%Input(1)%PtFairDisplacement, trim(p_FAST%OutFileRoot)//'.MAP_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )        
!   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, MD%Input(1)%PtFairleadDisplacement, trim(p_FAST%OutFileRoot)//'.MD_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )        
!   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, FEAM%Input(1)%PtFairleadDisplacement, trim(p_FAST%OutFileRoot)//'FEAM_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )        
!   END IF
         
   
END SUBROUTINE WrVTK_BasicMeshes 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes a minimal subset of meshes with surfaces to VTK-formatted files. It doesn't bother with 
!! returning an error code.
SUBROUTINE WrVTK_Surfaces(t_global, p_FAST, y_FAST, MeshMapData, ED, BD, AD14, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code (only because we're updating VTK_LastWaveIndx)
   TYPE(FAST_ModuleMapType), INTENT(IN   ) :: MeshMapData         !< Data for mapping between modules

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(IN   ) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(IN   ) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop


   logical, parameter                      :: OutputFields = .FALSE. ! due to confusion about what fields mean on a surface, we are going to just output the basic meshes if people ask for fields
   INTEGER(IntKi)                          :: NumBl, k, Twidth
   CHARACTER(1024)                         :: VTK_path
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_Surfaces'

   ! Calculate the number of digits for the maximum number of output steps to be written.
   ! This will be used to pad the write-out step in the VTK filename with zeros in calls to MeshWrVTK()
   if ( (p_FAST%n_VTKTime>0) .and. (p_FAST%n_TMax_m1+1>0) ) then
      Twidth = CEILING( log10( real(p_FAST%n_TMax_m1+1, ReKi) / p_FAST%n_VTKTime ) ) + 1
   else
      Twidth = 1
   endif
   
   
   NumBl = 0
   if (allocated(ED%Output(1)%BladeRootMotion)) then
      NumBl = SIZE(ED%Output(1)%BladeRootMotion)    
   end if

   VTK_path = get_vtkroot_path( p_FAST%OutFileRoot )

! Ground (written at initialization)
   
! Wave elevation
   if ( allocated( p_FAST%VTK_Surface%WaveElev ) ) call WrVTK_WaveElev( t_global, p_FAST, y_FAST, HD)
   
   
! Nacelle
   call MeshWrVTK_PointSurface (p_FAST%TurbinePos, ED%Output(1)%NacelleMotion, trim(VTK_path)//'.NacelleSurface', &
                                y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, Twidth , verts = p_FAST%VTK_Surface%NacelleBox, Sib=ED%Input(1)%NacelleLoads )
   
   
! Hub
   call MeshWrVTK_PointSurface (p_FAST%TurbinePos, ED%Output(1)%HubPtMotion, trim(VTK_path)//'.HubSurface', &
                                y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, Twidth , &
                                NumSegments=p_FAST%VTK_Surface%NumSectors, radius=p_FAST%VTK_Surface%HubRad, Sib=ED%Input(1)%HubPtLoad )
   
! Blades
   IF ( p_FAST%CompAero == Module_AD ) THEN  ! These meshes may have airfoil data associated with nodes...
      DO K=1,NumBl
         call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, AD%Input(1)%BladeMotion(K), trim(VTK_path)//'.Blade'//trim(num2lstr(k))//'Surface', &
                                    y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, Twidth , verts=p_FAST%VTK_Surface%BladeShape(K)%AirfoilCoords &
                                    ,Sib=AD%y%BladeLoad(k) )
      END DO                  
   ELSE IF ( p_FAST%CompElast == Module_BD ) THEN
      DO K=1,NumBl                 
         call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, BD%y(k)%BldMotion, trim(VTK_path)//'.Blade'//trim(num2lstr(k))//'Surface', &
                                    y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, Twidth , verts=p_FAST%VTK_Surface%BladeShape(K)%AirfoilCoords )
      END DO  
   ELSE
      DO K=1,NumBl        
         call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, ED%Output(1)%BladeLn2Mesh(K), trim(VTK_path)//'.Blade'//trim(num2lstr(k))//'Surface', &
                                    y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, Twidth , verts=p_FAST%VTK_Surface%BladeShape(K)%AirfoilCoords )
      END DO  
   END IF   
         
! Tower motions
   call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, ED%Output(1)%TowerLn2Mesh, trim(VTK_path)//'.TowerSurface', &
                              y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, Twidth, p_FAST%VTK_Surface%NumSectors, p_FAST%VTK_Surface%TowerRad )
   
! Platform
! call MeshWrVTK_PointSurface (p_FAST%TurbinePos, ED%Output(1)%PlatformPtMesh, trim(p_FAST%OutFileRoot)//'.PlatformSurface', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, Radius = p_FAST%VTK_Surface%GroundRad )
   
   
! Substructure   
!   call MeshWrVTK(p_FAST%TurbinePos, ED%Output(1)%PlatformPtMesh, trim(p_FAST%OutFileRoot)//'.ED_PlatformPtMesh_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )     
!   IF ( p_FAST%CompSub == Module_SD ) THEN
!     call MeshWrVTK(p_FAST%TurbinePos, SD%Input(1)%TPMesh, trim(p_FAST%OutFileRoot)//'.SD_TPMesh_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )     
!      call MeshWrVTK(p_FAST%TurbinePos, SD%y%y2Mesh, trim(p_FAST%OutFileRoot)//'.SD_y2Mesh_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )        
!   END IF     
      
   IF ( HD%Input(1)%Morison%DistribMesh%Committed ) THEN 
      !if ( p_FAST%CompSub == Module_NONE ) then ! floating
      !   OutputFields = .false.
      !else
      !   OutputFields = p_FAST%VTK_fields
      !end if
         
      call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, HD%Input(1)%Morison%DistribMesh, trim(VTK_path)//'.MorisonSurface', &
                                 y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, Twidth, p_FAST%VTK_Surface%NumSectors, &
                                 p_FAST%VTK_Surface%MorisonRad, Sib=HD%y%Morison%DistribMesh )
   END IF
   
   
! Mooring Lines?            
!   IF ( p_FAST%CompMooring == Module_MAP ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, MAPp%Input(1)%PtFairDisplacement, trim(p_FAST%OutFileRoot)//'.MAP_PtFair_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )        
!   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, MD%Input(1)%PtFairleadDisplacement, trim(p_FAST%OutFileRoot)//'.MD_PtFair_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )        
!   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, FEAM%Input(1)%PtFairleadDisplacement, trim(p_FAST%OutFileRoot)//'FEAM_PtFair_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2   )        
!   END IF
         
   
   if (p_FAST%VTK_fields) then
      call WrVTK_BasicMeshes(p_FAST, y_FAST, MeshMapData, ED, BD, AD14, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
   end if
   
   
END SUBROUTINE WrVTK_Surfaces 
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine writes the wave elevation data for a given time step
SUBROUTINE WrVTK_WaveElev(t_global, p_FAST, y_FAST, HD)

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code

   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data

   ! local variables
   INTEGER(IntKi)                        :: Un                    ! fortran unit number
   INTEGER(IntKi)                        :: n, iy, ix             ! loop counters
   REAL(SiKi)                            :: t
   CHARACTER(1024)                       :: FileName
   INTEGER(IntKi)                        :: NumberOfPoints 
   INTEGER(IntKi), parameter             :: NumberOfLines = 0
   INTEGER(IntKi)                        :: NumberOfPolys
   INTEGER(IntKi)                        :: Twidth
   CHARACTER(1024)                       :: VTK_path, Tstr = ''
   INTEGER(IntKi)                        :: ErrStat2 
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*),PARAMETER                :: RoutineName = 'WrVTK_WaveElev'

   
   NumberOfPoints = size(p_FAST%VTK_surface%WaveElevXY,2)
      ! I'm going to make triangles for now. we should probably just make this a structured file at some point
   NumberOfPolys  = ( p_FAST%VTK_surface%NWaveElevPts(1) - 1 ) * &
                    ( p_FAST%VTK_surface%NWaveElevPts(2) - 1 ) * 2
   
   !.................................................................
   ! write the data that potentially changes each time step:
   !.................................................................
   ! Calculate the number of digits for the maximum number of output steps to be written.
   ! This will be used to pad the write-out step in the VTK filename with zeros in calls to MeshWrVTK()
   if ( (p_FAST%n_VTKTime>0) .and. (p_FAST%n_TMax_m1+1>0) ) then
      Twidth = CEILING( log10( real(p_FAST%n_TMax_m1+1, ReKi) / p_FAST%n_VTKTime ) ) + 1
   else
      Twidth = 1
   endif

   VTK_path = get_vtkroot_path( p_FAST%OutFileRoot )

   ! construct the string for the zero-padded VTK write-out step
   write(Tstr(1 : Twidth), '(i' // trim(Num2LStr(Twidth)) //'.'// trim(Num2LStr(Twidth)) // ')') y_FAST%VTK_count
      
   ! PolyData (.vtp) - Serial vtkPolyData (unstructured) file
   FileName = TRIM(VTK_path)//'.WaveSurface.'//TRIM(Tstr)//'.vtp'
      
   call WrVTK_header( FileName, NumberOfPoints, NumberOfLines, NumberOfPolys, Un, ErrStat2, ErrMsg2 )    
      if (ErrStat2 >= AbortErrLev) return
         
! points (nodes, augmented with NumSegments):   
      WRITE(Un,'(A)')         '      <Points>'
      WRITE(Un,'(A)')         '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'

      ! I'm not going to interpolate in time; I'm just going to get the index of the closest wave time value
      t = REAL(t_global,SiKi)
      call GetWaveElevIndx( t, HD%p%WaveTime, y_FAST%VTK_LastWaveIndx )
      
      n = 1
      do ix=1,p_FAST%VTK_surface%NWaveElevPts(1)
         do iy=1,p_FAST%VTK_surface%NWaveElevPts(2)            
            WRITE(Un,VTK_AryFmt) p_FAST%VTK_surface%WaveElevXY(:,n), p_FAST%VTK_surface%WaveElev(y_FAST%VTK_LastWaveIndx,n) 
            n = n+1
         end do
      end do
                     
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Points>'
  
                  
      WRITE(Un,'(A)')         '      <Polys>'      
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'         
      
      do ix=1,p_FAST%VTK_surface%NWaveElevPts(1)-1
         do iy=1,p_FAST%VTK_surface%NWaveElevPts(2)-1
            n = p_FAST%VTK_surface%NWaveElevPts(1)*(ix-1)+iy - 1 ! points start at 0
            
            WRITE(Un,'(3(i7))') n,   n+1,                                    n+p_FAST%VTK_surface%NWaveElevPts(2)
            WRITE(Un,'(3(i7))') n+1, n+1+p_FAST%VTK_surface%NWaveElevPts(2), n+p_FAST%VTK_surface%NWaveElevPts(2)
            
         end do
      end do            
      WRITE(Un,'(A)')         '        </DataArray>'      
      
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="offsets" format="ascii">'                  
      do n=1,NumberOfPolys
         WRITE(Un,'(i7)') 3*n
      end do      
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Polys>'      
                  
      call WrVTK_footer( Un )       
      
END SUBROUTINE WrVTK_WaveElev
!----------------------------------------------------------------------------------------------------------------------------------
!> This function returns the index, Ind, of the XAry closest to XValIn, where XAry is assumed to be periodic. It starts
!! searching at the value of Ind from a previous step.
SUBROUTINE GetWaveElevIndx( XValIn, XAry, Ind )

      ! Argument declarations.

   INTEGER, INTENT(INOUT)       :: Ind                ! Initial and final index into the arrays.

   REAL(SiKi), INTENT(IN)       :: XAry    (:)        !< Array of X values to be interpolated.
   REAL(SiKi), INTENT(IN)       :: XValIn             !< X value to be found

   
   INTEGER                      :: AryLen             ! Length of the arrays.
   REAL(SiKi)                   :: XVal               !< X to be found (wrapped/periodic)
   
   
   AryLen = size(XAry)
   
      ! Wrap XValIn into the range XAry(1) to XAry(AryLen)
   XVal = MOD(XValIn, XAry(AryLen))

   
   
        ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      Ind = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      Ind = AryLen
      RETURN
   ELSE
      ! Set the Ind to the first index if we are at the beginning of XAry
      IF ( XVal <= XAry(2) )  THEN  
         Ind = 1
      END IF      
   END IF


     ! Let's interpolate!

   Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

   DO

      IF ( XVal < XAry(Ind) )  THEN

         Ind = Ind - 1

      ELSE IF ( XVal >= XAry(Ind+1) )  THEN

         Ind = Ind + 1

      ELSE
         
         ! XAry(Ind) <= XVal < XAry(Ind+1)
         ! this would make it the "closest" node, but I'm not going to worry about that for visualization purposes
         !if ( XVal > (XAry(Ind+1) + XAry(Ind))/2.0_SiKi ) Ind = Ind + 1

         RETURN

      END IF

   END DO

   RETURN
END SUBROUTINE GetWaveElevIndx

END MODULE