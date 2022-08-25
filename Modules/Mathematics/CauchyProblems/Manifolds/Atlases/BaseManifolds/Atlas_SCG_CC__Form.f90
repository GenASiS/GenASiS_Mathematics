module Atlas_SCG_CC__Form

  !-- Atlas_SingleChartGrid_CentralCore_Form

  use Basics
  use Charts
  use Atlas_SCG_C__Form

  implicit none
  private

  type, public, extends ( Atlas_SCG_C_Form ) :: Atlas_SCG_CC_Form
    class ( Chart_GS_CC_Form ), pointer :: &
      Chart_GS_CC => null ( )
  contains
    procedure, private, pass :: &
      Initialize_SCG
    procedure, private, pass :: &
      Initialize_SCG_CC
    procedure, private, pass :: &
      Initialize_SCG_CC_A  !-- Average
    generic, public :: &
      Initialize => Initialize_SCG_CC, Initialize_SCG_CC_A
    final :: &
      Finalize
  end type Atlas_SCG_CC_Form


contains


  subroutine Initialize_SCG &
               ( A, CommunicatorOption, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, NameOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nGhostLayersOption, nBricksOption, &
                 nBricksCompatibleOption, IgnorabilityOption, &
                 nDimensionsOption, nEqualOption, iDimensionalityOption )

    class ( Atlas_SCG_CC_Form ), intent ( inout ), target :: &
      A
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption, &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption, &
      NameOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption, &
      RatioOption, &
      ScaleOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nGhostLayersOption, &
      nBricksOption, &
      nBricksCompatibleOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      nDimensionsOption, &
      nEqualOption, &
      iDimensionalityOption

      call Show ( 'The method Initialize_SCG is not appropriate for ' &
                  // 'this class.', CONSOLE % ERROR )
      call Show ( 'Please use a different Initialize interface.', &
                  CONSOLE % ERROR )
      call Show ( 'Atlas_SCG_CC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize_SCG', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )

  end subroutine Initialize_SCG


  subroutine Initialize_SCG_CC &
               ( A, RadiusMax, RadiusCore, CommunicatorOption, NameOption, &
                 CoordinateUnitOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsPolarOption, nEqualOption, nDimensionsOption )

    class ( Atlas_SCG_CC_Form ), intent ( inout ), target :: &
      A
    real ( KDR ), intent ( in ) :: &
      RadiusMax, &
      RadiusCore
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nEqualOption, &
      nDimensionsOption

    logical :: &
      PreviouslyAllocated

    if ( A % Type  ==  '' ) &
      A % Type  =  'an Atlas_SCG_CC'

    if ( .not. allocated ( A % Chart ) ) &
      allocate ( A % Chart ( 1 ) )
    
    if ( allocated ( A % Chart ( 1 ) % Element ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( Chart_GS_CC_Form :: A % Chart ( 1 ) % Element )
    end if

    call A % Atlas_SCG_Form % Initialize ( NameOption = NameOption )

    if ( .not. PreviouslyAllocated ) then

      select type ( C  =>  A % Chart ( 1 ) % Element )
      class is ( Chart_GS_CC_Form )

      call C % Initialize &
             ( RadiusMax, RadiusCore, &
               CommunicatorOption = CommunicatorOption, &
               NameOption = NameOption, &
               CoordinateUnitOption = CoordinateUnitOption, &
               RadialRatioOption = RadialRatioOption, &
               nGhostLayersOption = nGhostLayersOption, &
               nCellsPolarOption = nCellsPolarOption, &
               nEqualOption = nEqualOption, &
               nDimensionsOption = nDimensionsOption )

      end select !--  C

    end if !-- PreviouslyAllocated  

    select type ( C  =>  A % Chart ( 1 ) % Element )
    class is ( Chart_GS_CC_Form )
      A % Chart_GS_C  =>  C
      A % Chart_GS_CC  =>  C
    end select !-- C
      
  end subroutine Initialize_SCG_CC


  subroutine Initialize_SCG_CC_A ( A, A_S, nDimensions )

    class ( Atlas_SCG_CC_Form ), intent ( inout ) :: &
      A
    class ( Atlas_SCG_CC_Form ), intent ( in ) :: &
      A_S  !-- Source
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    character ( 3 ) :: &
      Suffix

    associate ( C_S  =>  A_S % Chart_GS_CC )

    select case ( nDimensions )
    case ( 1 )
      Suffix = '_SA'
    case ( 2 )
      Suffix = '_AA'
    case default
      call Show ( 'Expecting nDimensions = 1 or 2', CONSOLE % ERROR )
      call Show ( '(spherical or azimuthal average)', CONSOLE % ERROR )
      call Show ( nDimensions, 'nDimensions', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    call A % Initialize &
           ( RadiusMax = C_S % RadiusMax, &
             RadiusCore = C_S % RadiusScale, &
             CommunicatorOption = C_S % Communicator, &
             NameOption = trim ( A_S % Name ) // Suffix, &
             CoordinateUnitOption = C_S % CoordinateUnit, &
             RadialRatioOption = C_S % RadialRatio, &
             nGhostLayersOption = C_S % nGhostLayers, &
             nCellsPolarOption = C_S % nCellsPolar, &
             nEqualOption = C_S % nEqual, &
             nDimensionsOption = nDimensions )

    end associate !-- C_S

  end subroutine Initialize_SCG_CC_A


  impure elemental subroutine Finalize ( A )

    type ( Atlas_SCG_CC_Form ), intent ( inout ) :: &
      A

    nullify ( A % Chart_GS_CC )

  end subroutine Finalize


end module Atlas_SCG_CC__Form
