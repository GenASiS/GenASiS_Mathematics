!-- Storage_BSLL_ASC_CSLD is a template for storage on a 
!   Bundle_SLL_ASC_CSLD.

module Storage_BSLL_ASC_CSLD__Form

  use Basics
  use Atlases
  use Bundle_SLL_ASC_CSLD__Form
  use Field_BSLL_ASC_CSLD__Template

  !-- Storage_BundleSingleLevelLocal_AtlasSingleChart
  !   _ChartSingleLevelDistributed_Form

  implicit none
  private

  type, public, extends ( Field_BSLL_ASC_CSLD_Template ) :: &
    Storage_BSLL_ASC_CSLD_Form
      integer ( KDI ) :: &
        nFields
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Storage_BSLL_ASC_CSLD_Form

contains


  subroutine Initialize ( SB, B, nFields, NameOutputOption )

    class ( Storage_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      SB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    integer ( KDI ), intent ( in ) :: &
      nFields
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    SB % nFields = nFields

    call SB % InitializeTemplate ( B, NameOutputOption = NameOutputOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( SB )

    type ( Storage_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      SB

    call SB % FinalizeTemplate_BSLL_ASC_CSLD ( )

  end subroutine Finalize


  subroutine SetField ( FB, NameOutputOption )

    class ( Storage_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      FB
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iS     !-- iSection

    associate ( B => FB % Bundle_SLL_ASC_CSLD )

    !-- Fibers

    allocate ( FB % Fiber )
    associate ( FBF => FB % Fiber )
    call FBF % Initialize ( FB % nFibers )

    do iF = 1, FB % nFibers
      allocate ( Storage_ASC_Form :: FBF % Atlas ( iF ) % Element )
      select type ( SA => FBF % Atlas ( iF ) % Element )
      class is ( Storage_ASC_Form )
        select type ( AF => B % Fiber % Atlas ( iF ) % Element )
        class is ( Atlas_SC_Form )
          call SA % Initialize ( AF, FB % nFields )
        end select !-- AF
      end select !-- SA
    end do !-- iF

    end associate !-- FBF

    !-- Sections

    allocate ( FB % Section )
    associate ( FBS => FB % Section )
    call FBS % Initialize ( FB % nSections )

    do iS = 1, FB % nSections
      allocate ( Storage_ASC_Form :: FBS % Atlas ( iS ) % Element )
      select type ( SA => FBS % Atlas ( iS ) % Element )
      class is ( Storage_ASC_Form )
        call SA % Initialize ( B % Base_ASC, FB % nFields ) 
      end select !-- SA
    end do !-- iS

    end associate !-- FBF

    end associate !-- B

  end subroutine SetField


end module Storage_BSLL_ASC_CSLD__Form
