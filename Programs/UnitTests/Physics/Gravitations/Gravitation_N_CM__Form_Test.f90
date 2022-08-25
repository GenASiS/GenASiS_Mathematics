program Gravitation_N_CM__Form_Test

  !-- Gravitation_Newton_CentralMass__Form_Test

  use Basics
  use Mathematics
  use Gravitations

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CE_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    S
  type ( FieldSetForm ), allocatable :: &
    F  !-- Dummy field
  type ( Gravitation_N_CM_Form ), allocatable :: &
    G

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Gravitation_N_CM__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusExcision = 0.45_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'PositionSpace' )

  allocate ( S )
  call S % Initialize ( A, GIS )

  allocate ( F )
  call F % Initialize ( A )

  allocate ( G )
  call G % Initialize ( A, GravitationalConstant = 1.0_KDR, Mass = 1.0_KDR )
  call G % SetStream ( S )

  call A % Show ( )
  call G % Show ( )
  call S % Show ( )

  call G % Solve ( F, iBaryonMass = 0, iBaryonDensity = 0 )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  deallocate ( G )
  deallocate ( F )
  deallocate ( S )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Gravitation_N_CM__Form_Test
