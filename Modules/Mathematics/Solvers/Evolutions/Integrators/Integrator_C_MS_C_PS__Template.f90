!-- Integrator_C_MS_C_PS is a template for time evolution of a conserved 
!   currents on position space.

module Integrator_C_MS_C_PS__Template

  !-- Integrator_Current_MomentumSpace_C_PositionSpace__Template

  use Basics
  use Manifolds
  use Fields
  use Steps
  use Integrator_C_1D_PS__Template

  implicit none
  private

  type, public, extends ( Integrator_C_1D_PS_Template ), abstract :: &
    Integrator_C_MS_C_PS_Template
      logical ( KDL ) :: &
        UseLimiterParameter_S, &
        UseLimiterParameter_F
      class ( Current_BSLL_ASC_CSLD_Template ), allocatable :: &
        Current_BSLL_ASC_CSLD
  contains
    procedure, public, pass :: &  !-- 1
      InitializeTemplate_C_MS_C_PS
    procedure, public, pass :: &  !-- 1
      FinalizeTemplate_C_MS_C_PS
    procedure, private, pass :: &  !-- 2
      ComputeCycle
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, private, pass :: &
      ComputeCycle_BSLL_ASC_CSLD
    procedure, private, pass :: &
      ComputeTimeStepLocal
  end type Integrator_C_MS_C_PS_Template

      integer ( KDI ), private, parameter :: &
        iMS_S = 1, &
        iPS   = 2, &
        N_PS  = 2

contains


  subroutine InitializeTemplate_C_MS_C_PS &
               ( I, Name, UseLimiterParameter_S_Option, &
                 UseLimiterParameter_F_Option, TimeUnitOption, &
                 FinishTimeOption, nWriteOption )

    class ( Integrator_C_MS_C_PS_Template ), intent ( inout ) :: &
      I
    character ( * ), intent ( in )  :: &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterParameter_S_Option, &
      UseLimiterParameter_F_Option
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    if ( I % Type == '' ) &
      I % Type = 'an Integrator_C_MS_C_PS'

    if ( .not. allocated ( I % Current_BSLL_ASC_CSLD ) ) then
      call Show ( 'Current_BSLL_ASC_CSLD not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_MS_C_PS__Template', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_MS_C_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    call I % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, nWriteOption = nWriteOption )

    I % UseLimiterParameter_S = .false.
    I % UseLimiterParameter_F = .false.
    if ( present ( UseLimiterParameter_S_Option ) ) &
      I % UseLimiterParameter_S = UseLimiterParameter_S_Option
    if ( present ( UseLimiterParameter_F_Option ) ) &
      I % UseLimiterParameter_F = UseLimiterParameter_F_Option
    call PROGRAM_HEADER % GetParameter &
           ( I % UseLimiterParameter_S, 'UseLimiterParameter_S' )
    call PROGRAM_HEADER % GetParameter &
           ( I % UseLimiterParameter_F, 'UseLimiterParameter_F' )

  end subroutine InitializeTemplate_C_MS_C_PS


  subroutine FinalizeTemplate_C_MS_C_PS ( I )

    class ( Integrator_C_MS_C_PS_Template ), intent ( inout ) :: &
      I

   if ( allocated ( I % Current_BSLL_ASC_CSLD ) ) &
     deallocate ( I % Current_BSLL_ASC_CSLD )

    call I % FinalizeTemplate_C_PS ( )

  end subroutine FinalizeTemplate_C_MS_C_PS


  subroutine ComputeCycle ( I )

    class ( Integrator_C_MS_C_PS_Template ), intent ( inout ) :: &
      I

    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerCycle )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( MS => I % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
      call I % ComputeCycle_BSLL_ASC_CSLD ( MS )
    end select

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeCycle


  subroutine ComputeTally ( I, ComputeChangeOption, IgnorabilityOption )

    class ( Integrator_C_MS_C_PS_Template ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    type ( TimerForm ), pointer :: &
      Timer

    if ( .not. allocated ( I % Current_BSLL_ASC_CSLD ) ) &
      return

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerTally )
    if ( associated ( Timer ) ) call Timer % Start ( )

    associate &
      ( CB => I % Current_BSLL_ASC_CSLD, &
        CA => I % Current_ASC )

    call CB % ComputeTally &
           ( ComputeChangeOption = ComputeChangeOption, &
             IgnorabilityOption  = IgnorabilityOption )
    call CA % ComputeTally &
           ( ComputeChangeOption = ComputeChangeOption, &
             IgnorabilityOption  = IgnorabilityOption )

    end associate !-- CB, etc.

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeTally


  subroutine ComputeCycle_BSLL_ASC_CSLD ( I, MS )

    class ( Integrator_C_MS_C_PS_Template ), intent ( inout ) :: &
      I
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
      MS

    integer ( KDI ) :: &
      iS  !-- iSection
    real ( KDR ) :: &
      TimeNew

    call I % PrepareCycle ( )
    call I % ComputeNewTime ( TimeNew )

    select type ( S => I % Step )
    class is ( Step_RK_C_BSLL_ASC_CSLD_C_ASC_Template )

    associate ( TimeStep => TimeNew - I % Time )    

    call S % Compute ( I % Time, TimeStep )

    I % iCycle = I % iCycle + 1
    I % Time = I % Time + TimeStep

    if ( I % WriteTimeExact ) then
      if ( I % Time == I % WriteTime ) &
        I % IsCheckpointTime = .true.
    else 
      if ( I % Time  >  I % WriteTime &
           .or. abs ( I % Time - I % WriteTime )  <  0.5_KDR * TimeStep ) &
        I % IsCheckpointTime = .true.
    end if

    end associate !-- TimeStep
    end select !-- S

  end subroutine ComputeCycle_BSLL_ASC_CSLD


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( Integrator_C_MS_C_PS_Template ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    integer ( KDI ) :: &
      iC  !-- iCurrent
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( CurrentTemplate ), pointer :: &
      C
    class ( Current_ASC_Template ), pointer :: &
      CA

    if ( .not. allocated ( I % Current_BSLL_ASC_CSLD ) ) &
      return

    associate ( CB  => I % Current_BSLL_ASC_CSLD )
    associate ( CSL => CB % Bundle_SLL_ASC_CSLD % Base_CSLD )

    G => CSL % Geometry ( )

    do iC = 1, N_PS

      select case ( iC )
      case ( iMS_S )
        select type ( CAE => CB % Section % Atlas ( 1 ) % Element )
        class is ( Current_ASC_Template )
          CA => CAE
        end select !-- CAE
      case ( iPS )
        CA => I % Current_ASC
      end select !-- iC
      C => CA % Current ( )

      call I % ComputeTimeStepKernel_CSL &
             ( CSL % IsProperCell, &
               C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
               G % Value ( :, G % WIDTH ( 1 ) ), &
               G % Value ( :, G % WIDTH ( 2 ) ), & 
               G % Value ( :, G % WIDTH ( 3 ) ), &
               CSL % nDimensions, TimeStepCandidate ( iC ) )

    end do !-- iC

    end associate !-- CSL
    end associate !-- CB

    TimeStepCandidate = I % CourantFactor * TimeStepCandidate

    nullify ( G, C, CA )

  end subroutine ComputeTimeStepLocal


end module Integrator_C_MS_C_PS__Template
