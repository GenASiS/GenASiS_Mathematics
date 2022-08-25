program Slope_DFV_DP__Form_Test

  !-- Slope_DivergenceFiniteVolume_DivergencePart__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Slopes

  implicit none

  integer ( KDI ) :: &
    nCompute
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    Sm
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( CurrentSetForm ), allocatable :: &
    CS
  type ( DivergencePart_CS_Form ), allocatable :: &
    DP
  type ( RiemannSolver_HLL_Form ), allocatable :: &
    RS
  type ( Slope_DFV_DP_Form ), allocatable :: &
    S

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Slope_DFV_DP__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( Sm )
  call Sm % Initialize ( A, GIS )

  allocate ( G )
  call G % Initialize ( A )
  call G % SetStream ( Sm )

  allocate ( CS )
  call CS % Initialize ( G )
  call CS % SetStream ( Sm )

  allocate ( DP )
  call DP % Initialize ( CS )

  allocate ( RS )
  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call RS % Initialize ( CS )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  allocate ( S )
  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call S % Initialize ( RS, DP, IgnorabilityOption = A % IGNORABILITY )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  call S % SetStream ( Sm )

  call  A % Show ( )
  call  G % Show ( )
  call CS % Show ( )
  call DP % Show ( )
  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call RS % Show ( )
  call  S % Show ( )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  call Sm % Show ( )

  call SetWave ( CS, G )

  nCompute  =  1000
  call PROGRAM_HEADER % GetParameter ( nCompute, 'nCompute' )

!  call CONSOLE % SetVerbosity ( 'INFO_4' )
  call TestSlope ( S, Sm )
!  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  deallocate ( S )
  deallocate ( RS )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  deallocate ( DP )
  deallocate ( CS )
  deallocate ( G )
  deallocate ( Sm )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( CS, G )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G

    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Offset, &
      Amplitude, &
      Speed
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber

    select type ( A  =>  CS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    nWavelengths  =  0
    nWavelengths ( 1 : C % nDimensions )  =  1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    Offset     =  2.0_KDR
    Amplitude  =  1.0_KDR
    Speed      =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( Amplitude, 'Amplitude' )
    call PROGRAM_HEADER % GetParameter ( Speed, 'Speed' )

    associate ( BoxSize  =>  C % MaxCoordinate  -  C % MinCoordinate )
    where ( BoxSize  >  0.0_KDR )
      Wavenumber  =  nWavelengths / BoxSize
    elsewhere
      Wavenumber  =  0.0_KDR
    end where
    end associate !-- BoxSize

    associate &
      (  GV  =>   G % Storage_GS % Value, &
        CSV  =>  CS % Storage_GS % Value )
    associate &
      (     X  =>   GV ( :,  G % CENTER_U_1 ), &
            Y  =>   GV ( :,  G % CENTER_U_2 ), &
            Z  =>   GV ( :,  G % CENTER_U_3 ), &
          Rho  =>  CSV ( :, CS % DENSITY_CS ), &
          V_1  =>  CSV ( :, CS % VELOCITY_CS_U_1 ), &
          V_2  =>  CSV ( :, CS % VELOCITY_CS_U_2 ), &
          V_3  =>  CSV ( :, CS % VELOCITY_CS_U_3 ), &
            K  =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    Rho  =  Offset  &
            +  Amplitude  &
               *  sin ( TwoPi * (    K ( 1 ) * X  &
                                  +  K ( 2 ) * Y  &
                                  +  K ( 3 ) * Z  ) )

    V_1  =  Speed  *  K ( 1 )  /  Abs_K
    V_2  =  Speed  *  K ( 2 )  /  Abs_K
    V_3  =  Speed  *  K ( 3 )  /  Abs_K
    
    end associate !-- Rho, etc.
    end associate !-- GV, etc.
    end associate !-- C
    end select !-- A

  end subroutine SetWave


  subroutine TestSlope ( S, Sm )

    class ( Slope_DFV_DP_Form ), intent ( inout ) :: &
      S
    class ( StreamForm ), intent ( inout ) :: &
      Sm

    integer ( KDI ) :: &
      iC, &  !-- iCompute
      iD     !-- iDimension
    type ( TimerForm ), pointer :: &
      T, &
      T_RPP

    call Show ( 'Slope computation' )
    call Show ( S % Name, 'Slope' )
    call Show ( nCompute, 'nCompute' )

    associate ( C  =>  S % Atlas % Chart ( 1 ) % Element )

    T      =>  S % Timer ( Level = 1 )
    T_RPP  =>  RS % Timer_P ( Level = T % Level + 1 )
    call T % Start ( )
    do iC  =  1,  nCompute
      do iD  =  1,  C % nDimensions
        call T_RPP % Start ( )
        call RS % Prepare ( iC = 1, iD = iD, T_Option = T_RPP )
        call T_RPP % Stop ( )
        call S % ComputePartialDerivative ( iC = 1, iD = iD, T_Option = T )
      end do !-- iD
      call S % ComputeConnectionFlat ( iC = 1, T_Option = T )
      call S % AddComponents ( T_Option = T )
    end do !-- iC
    call T % Stop ( )

    end associate !-- C

    T  =>  Sm % TimerWrite ( Level = 1 )
    call T % Start ( )
    call GIS % Open ( GIS % ACCESS_CREATE )
    call Sm % Write ( )
    call GIS % Close ( )
    call T % Stop ( )

  end subroutine TestSlope


end program Slope_DFV_DP__Form_Test
