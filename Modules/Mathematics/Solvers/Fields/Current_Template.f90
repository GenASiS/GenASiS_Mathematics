!-- Current is a template for a 4-current, of the type appearing in
!   conservation laws.

module Current_Template

  use Basics
  use Manifolds

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_TEMPLATE = 0, &
      N_CONSERVED_TEMPLATE = 0, &
      N_FIELDS_TEMPLATE    = 6, &
      N_VECTORS_TEMPLATE   = 2

  type, public, extends ( VariableGroupForm ), abstract :: CurrentTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      N_PRIMITIVE  = 0, &
      N_CONSERVED  = 0, &
      N_FIELDS     = 0, &
      N_VECTORS    = 0, &
      N_PRIMITIVE_TEMPLATE = N_PRIMITIVE_TEMPLATE, &
      N_CONSERVED_TEMPLATE = N_CONSERVED_TEMPLATE, &
      N_FIELDS_TEMPLATE    = N_FIELDS_TEMPLATE, &
      N_VECTORS_TEMPLATE   = N_VECTORS_TEMPLATE
    integer ( KDI ), dimension ( 3 ) :: &
      FAST_EIGENSPEED_PLUS  = 0, &
      FAST_EIGENSPEED_MINUS = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPrimitive, &
      iaConserved
    character ( LDL ) :: &
      Type = ''
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, private, pass :: &
      ComputeFromPrimitiveSelf
    procedure, private, pass ( C ) :: &
      ComputeFromPrimitiveOther
    generic, public :: &
      ComputeFromPrimitive &
        => ComputeFromPrimitiveSelf, ComputeFromPrimitiveOther
    procedure, private, pass :: &
      ComputeFromConservedSelf
    procedure, private, pass ( C ) :: &
      ComputeFromConservedOther
    generic, public :: &
      ComputeFromConserved &
        => ComputeFromConservedSelf, ComputeFromConservedOther
    procedure ( CRF ), public, pass ( C ), deferred :: &
      ComputeRawFluxes
    procedure, public, pass ( C ) :: &
      ComputeRiemannSolverInput
    procedure, public, pass ( C ) :: &
      ComputeDiffusionFactor
    procedure, public, pass :: &
      FinalizeTemplate
    procedure ( CFPC ), public, pass ( C ), deferred :: &
      ComputeFromPrimitiveCommon
    procedure ( CFCC ), public, pass ( C ), deferred :: &
      ComputeFromConservedCommon
    procedure, public, nopass :: &
      SetDiffusionFactorUnity
  end type CurrentTemplate

  abstract interface

    subroutine CRF ( RawFlux, C, G, Value_C, Value_G, iDimension, &
                     nValuesOption, oValueOption )
      use Basics
      use Manifolds
      import CurrentTemplate
      real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
        RawFlux
      class ( CurrentTemplate ), intent ( in ) :: &
        C
      class ( GeometryFlatForm ), intent ( in ) :: &
        G
      real ( KDR ), dimension ( :, : ), intent ( in ) :: &
        Value_C, &
        Value_G
      integer ( KDI ), intent ( in ) :: &
        iDimension
      integer ( KDI ), intent ( in ), optional :: &
        nValuesOption, &
        oValueOption
    end subroutine CRF
    
    subroutine CFPC ( Value_C, C, G, Value_G, nValuesOption, oValueOption )
      use Basics
      use Manifolds
      import CurrentTemplate
      real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
        Value_C
      class ( CurrentTemplate ), intent ( in ) :: &
        C
      class ( GeometryFlatForm ), intent ( in ) :: &
        G
      real ( KDR ), dimension ( :, : ), intent ( in ) :: &
        Value_G
      integer ( KDI ), intent ( in ), optional :: &
        nValuesOption, &
        oValueOption
    end subroutine CFPC

    subroutine CFCC ( C, G, Value, nValuesOption, oValueOption )
      use Basics
      use Manifolds
      import CurrentTemplate
      class ( CurrentTemplate ), intent ( in ) :: &
        C
      class ( GeometryFlatForm ), intent ( in ) :: &
        G
      real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
        Value
      integer ( KDI ), intent ( in ), optional :: &
        nValuesOption, &
        oValueOption
    end subroutine CFCC

  end interface

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeRiemannSolverInputKernel

contains


  subroutine InitializeTemplate &
               ( C, VelocityUnit, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), &
      optional :: &
        VectorIndicesOption

    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDF ) :: &
      Name 
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector
    logical ( KDL ) :: &
      Clear

    call InitializeBasics &
           ( C, Variable, Vector, Name, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call SetUnits ( VariableUnit, C, VelocityUnit )

    Clear = .true.
    if ( present ( ClearOption ) ) Clear = ClearOption

    call C % VariableGroupForm % Initialize &
           ( [ nValues, C % N_FIELDS ], &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = Name, ClearOption = Clear, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

  end subroutine InitializeTemplate


  subroutine ComputeFromPrimitiveSelf ( C, G, nValuesOption, oValueOption )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeFromPrimitiveCommon &
           ( C % Value, G, G % Value, nValuesOption, oValueOption )
    
  end subroutine ComputeFromPrimitiveSelf


  subroutine ComputeFromPrimitiveOther &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value_C
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeFromPrimitiveCommon &
           ( Value_C, G, Value_G, nValuesOption, oValueOption )
    
  end subroutine ComputeFromPrimitiveOther


  subroutine ComputeFromConservedSelf ( C, G, nValuesOption, oValueOption )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeFromConservedCommon &
           ( G, C % Value, nValuesOption, oValueOption )
    
  end subroutine ComputeFromConservedSelf


  subroutine ComputeFromConservedOther &
               ( Value_C, C, G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value_C
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeFromConservedCommon &
           ( G, Value_C, nValuesOption, oValueOption )
    
  end subroutine ComputeFromConservedOther


  subroutine ComputeRiemannSolverInput &
               ( Increment, DiffusionFactor_I, AlphaPlus_I, AlphaMinus_I, &
                 C, Grid, C_IL, C_IR, iDimension )
    
    class ( * ), intent ( inout ) :: &
      Increment
    type ( VariableGroupForm ), intent ( inout ) :: &
      DiffusionFactor_I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      AlphaPlus_I, &
      AlphaMinus_I
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( * ), intent ( in ) :: &
      Grid
    type ( VariableGroupForm ), intent ( in ) :: &
      C_IL, &
      C_IR
    integer ( KDI ), intent ( in ) :: &
      iDimension
    
    call ComputeRiemannSolverInputKernel &
           ( AlphaPlus_I, AlphaMinus_I, &
             C_IL % Value ( :, C % FAST_EIGENSPEED_PLUS ( iDimension ) ), &
             C_IR % Value ( :, C % FAST_EIGENSPEED_PLUS ( iDimension ) ), &
             C_IL % Value ( :, C % FAST_EIGENSPEED_MINUS ( iDimension ) ), &
             C_IR % Value ( :, C % FAST_EIGENSPEED_MINUS ( iDimension ) ) )

    call C % ComputeDiffusionFactor ( DiffusionFactor_I, Grid, iDimension )

  end subroutine ComputeRiemannSolverInput


  subroutine ComputeDiffusionFactor ( DF_I, Grid, C, iDimension )

    type ( VariableGroupForm ), intent ( inout ) :: &
      DF_I
    class ( * ), intent ( in ), target :: &
      Grid
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iDimension

    call C % SetDiffusionFactorUnity ( DF_I % Value )

  end subroutine ComputeDiffusionFactor


  impure elemental subroutine FinalizeTemplate ( C )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C

    if ( allocated ( C % iaConserved ) ) &
      deallocate ( C % iaConserved )
    if ( allocated ( C % iaPrimitive ) ) &
      deallocate ( C % iaPrimitive )
 
    call Show ( 'Finalizing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )
   
  end subroutine FinalizeTemplate


  subroutine SetDiffusionFactorUnity ( DFV_I )

    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      DFV_I

    integer ( KDI ) :: &
      iV, jV, &
      nValues

    nValues = size ( DFV_I, dim = 1 )

    do jV = 1, size ( DFV_I, dim = 2 )
      !$OMP parallel do private ( iV ) 
      do iV = 1, nValues
        DFV_I ( iV, jV )  =  1.0_KDR
      end do !-- iV
      !$OMP end parallel do
    end do !-- jV

  end subroutine SetDiffusionFactorUnity


  subroutine InitializeBasics &
               ( C, Variable, Vector, Name, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, NameOption, VariableUnitOption, &
                 VectorIndicesOption )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable, &
      Vector
    character ( LDF ), intent ( out ) :: &
      Name
    type ( MeasuredValueForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    !-- FIXME: intent(out) here caused ICE with Intel Compiler 15
    !          Temporarily set to intent(inout)
    !type ( Integer_1D_Form ), dimension ( : ), allocatable, &
    !  intent ( out ) :: &
    type ( Integer_1D_Form ), dimension ( : ), allocatable, &
      intent ( inout ) :: &
        VectorIndices
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption

    integer ( KDI ) :: &
      iV, &  !-- iVector
      iF  !-- iField
    character ( LDF ), dimension ( : ), allocatable :: &
      PrimitiveName, &
      ConservedName

    if ( C % Type == '' ) &
      C % Type = 'a Current'

    Name = 'Current'
    if ( present ( NameOption ) ) &
      Name = NameOption

    C % IGNORABILITY = CONSOLE % INFO_4 
    call Show ( 'Initializing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( Name, 'Name', C % IGNORABILITY )

    !-- variable indices

    if ( C % N_FIELDS == 0 ) &
      C % N_FIELDS = C % N_FIELDS_TEMPLATE

    C % FAST_EIGENSPEED_PLUS  = [ 1, 2, 3 ]
    C % FAST_EIGENSPEED_MINUS = [ 4, 5, 6 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( C % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : C % N_FIELDS_TEMPLATE ) &
      = [ 'FastEigenspeedPlus_1           ', &
          'FastEigenspeedPlus_2           ', &
          'FastEigenspeedPlus_3           ', &
          'FastEigenspeedMinus_1          ', &
          'FastEigenspeedMinus_2          ', &
          'FastEigenspeedMinus_3          ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( C % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    if ( C % N_VECTORS == 0 ) C % N_VECTORS = C % N_VECTORS_TEMPLATE

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( C % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( 1 : C % N_VECTORS_TEMPLATE ) &
      = [ 'FastEigenspeedPlus             ', &
          'FastEigenspeedMinus            ' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = C % N_VECTORS_TEMPLATE + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( C % N_VECTORS ) )
    end if

    call VectorIndices ( 1 ) % Initialize ( C % FAST_EIGENSPEED_PLUS )
    call VectorIndices ( 2 ) % Initialize ( C % FAST_EIGENSPEED_MINUS )

    !-- show primitive, conserved

    allocate ( PrimitiveName ( C % N_PRIMITIVE ) )
    do iF = 1, C % N_PRIMITIVE
      PrimitiveName ( iF ) = Variable ( C % iaPrimitive ( iF ) )
    end do !-- iF
    call Show ( PrimitiveName, 'Primitive', C % IGNORABILITY )
 
    allocate ( ConservedName ( C % N_CONSERVED ) )
    do iF = 1, C % N_CONSERVED
      ConservedName ( iF ) = Variable ( C % iaConserved ( iF ) )
    end do !-- iF
    call Show ( ConservedName, 'Conserved', C % IGNORABILITY )
 
  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, C, VelocityUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit

    integer ( KDI ) :: &
      iD

    do iD = 1, 3
      VariableUnit ( C % FAST_EIGENSPEED_PLUS ( iD ) )  = VelocityUnit ( iD )
      VariableUnit ( C % FAST_EIGENSPEED_MINUS ( iD ) ) = VelocityUnit ( iD )
    end do

  end subroutine SetUnits


  subroutine ComputeRiemannSolverInputKernel &
               ( AP_I, AM_I, LP_IL, LP_IR, LM_IL, LM_IR )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      AP_I, &
      AM_I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      LP_IL, LP_IR, &
      LM_IL, LM_IR

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( AP_I )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      AP_I ( iV ) = max ( 0.0_KDR, + LP_IL ( iV ), + LP_IR ( iV ) )
      AM_I ( iV ) = max ( 0.0_KDR, - LM_IL ( iV ), - LM_IR ( iV ) )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeRiemannSolverInputKernel


end module Current_Template