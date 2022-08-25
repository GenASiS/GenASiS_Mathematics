program Chart_GS_CE__Form_Test

  !-- Chart_GridStructured_CentralCore_Form_Test

  use Basics
  use StructuredGrids

  implicit none

  type ( Chart_GS_CE_Form ), allocatable :: &
    C

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Chart_GS_CE__Form_Test', DimensionalityOption = '2D' )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  allocate ( C )
  call C % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusExcision = 0.45_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  call C % Show ( )

  deallocate ( C )
  deallocate ( PROGRAM_HEADER )

end program Chart_GS_CE__Form_Test
