#include "Preprocessor"

submodule ( Slope_DFV_PD__Form ) Slope_DFV_PD__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iF, &
      iV, jV, kV, &
      nF
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVP, &
      lV, uV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    lV  =  1
    where ( shape ( S ( :, :, :, 1 ) )  >  1 )
      lV  =  oV + 1
    end where
    
    uV  =  1
    where ( shape ( S ( :, :, :, 1 ) )  >  1 )
      uV  =  shape ( S ( :, :, :, 1 ) )  -  oV
    end where
      
    iaS  =  0
    iaS ( iD )  =  1
    
    nF  =  size ( S, dim = 4 )

    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iaVP )
      do iF  =  1,  nF
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iaVP  =  [ iV, jV, kV ]  +  iaS

              S ( iV, jV, kV, iF )  &
                =  S ( iV, jV, kV, iF )  &
                   -  (    A_I ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  &
                           *  F_I ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )  &
                        -  A_I ( iV, jV, kV )  &
                           *  F_I ( iV, jV, kV, iF ) ) &
                       /  V ( iV, jV, kV )

            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iF
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else !-- use host
              
      !$OMP parallel do collapse ( 4 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iaVP )
      do iF  =  1,  nF
        do kV  =  lV ( 3 ),  uV ( 3 ) 
          do jV  =  lV ( 2 ),  uV ( 2 )
            do iV  =  lV ( 1 ),  uV ( 1 )

              iaVP  =  [ iV, jV, kV ]  +  iaS

              S ( iV, jV, kV, iF )  &
                =  S ( iV, jV, kV, iF )  &
                   -  (    A_I ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ) )  &
                           *  F_I ( iaVP ( 1 ), iaVP ( 2 ), iaVP ( 3 ), iF )  &
                        -  A_I ( iV, jV, kV )  &
                           *  F_I ( iV, jV, kV, iF ) ) &
                       /  V ( iV, jV, kV )

            end do !-- iV
          end do !-- jV
        end do !-- kV
      end do !-- iF
      !$OMP end parallel do
      
    end if !-- UseDevice
        
  end procedure ComputeKernel


  module procedure RecordBoundaryFlux_SCG_Kernel

    integer ( KDI ) :: &
      iV, jV, kV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then 
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
            BF ( iV, jV, kV ) &
              =  F ( oB ( 1 ) + iV, oB ( 2 ) + jV, oB ( 3 ) + kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else 
    
      !$OMP parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
            BF ( iV, jV, kV ) &
              =  F ( oB ( 1 ) + iV, oB ( 2 ) + jV, oB ( 3 ) + kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
    
    end if

  end procedure RecordBoundaryFlux_SCG_Kernel


end submodule Slope_DFV_PD__Kernel
