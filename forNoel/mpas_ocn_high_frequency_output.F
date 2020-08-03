!Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_high_frequency_output
!
!> \brief MPAS ocean analysis mode member: high_frequency_output
!> \author Todd Ringler and Phillip Wolfram
!> \date   2015/06/12
!> \details
!>  MPAS ocean analysis mode member: high_frequency_output
!>
!-----------------------------------------------------------------------

module ocn_high_frequency_output

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_init_high_frequency_output, &
             ocn_compute_high_frequency_output, &
             ocn_restart_high_frequency_output, &
             ocn_finalize_high_frequency_output

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_high_frequency_output
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler and Phillip Wolfram
!> \date    2015/06/12
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_high_frequency_output(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_init_high_frequency_output!}}}

!***********************************************************************
!
!  routine ocn_compute_high_frequency_output
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Todd Ringler and Phillip Wolfram
!> \date    2015/06/12
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_high_frequency_output(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: highFrequencyOutputAMPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: scratchPool

      integer :: iLevel, iLevelTarget, iCell, iEdge, i, cell1, cell2, k, eoe
      integer :: iLevel0100, iLevel0250, iLevel0700, iLevel2000
      integer, dimension(2) :: iLevel0005, iLevel0010, iLevel0015, iLevel0020, iLevel0025,   &
                               iLevel0030, iLevel0040, iLevel0050, iLevel0070, iLevel0095
      real (kind=RKIND) :: sumLayerThickness
      integer, pointer :: nVertLevels, nCells, nEdges
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelEdgeTop, maxLevelEdgeBot, maxLevelCell, nEdgesOnEdge
      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnEdge, edgesOnEdge

      real (kind=RKIND) :: invAreaCell1, layerThicknessEdge1, coeff, weightedNormalVel, cellArea
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, kineticEnergyAt250m, kineticEnergyAtSurface, relativeVorticityAt250m
      real (kind=RKIND), dimension(:), pointer :: divergenceAt250m, relativeVorticityVertexAt250m
      real (kind=RKIND), dimension(:), pointer :: divergenceAtBottom,relativeVorticityAtBottom,kineticEnergyAtBottom
      real (kind=RKIND), dimension(:), pointer :: vertVelAt250m
      real (kind=RKIND), dimension(:), pointer :: BruntVaisalaFreqTopAtSFC, BruntVaisalaFreqTopAt250m, BruntVaisalaFreqTopAtBottom

      real (kind=RKIND), dimension(:), pointer :: normalVelAtSFC, normalVelAt250m, normalVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: tangentialVelAtSFC, tangentialVelAt250m, tangentialVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: zonalVelAtSFC, zonalVelAt250m, zonalVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: meridionalVelAtSFC, meridionalVelAt250m, meridionalVelAtBottom

      real (kind=RKIND), dimension(:), pointer :: zonalAreaWeightedCellVelAtSFC, zonalAreaWeightedCellVelAt250m, zonalAreaWeightedCellVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: meridionalAreaWeightedCellVelAtSFC, meridionalAreaWeightedCellVelAt250m, meridionalAreaWeightedCellVelAtBottom

      real (kind=RKIND), dimension(:), pointer :: normalBarotropicVel
      real (kind=RKIND), dimension(:), pointer :: tangentialBarotropicVel
      real (kind=RKIND), dimension(:), pointer :: zonalBarotropicVel
      real (kind=RKIND), dimension(:), pointer :: meridionalBarotropicVel

      real (kind=RKIND), dimension(:), pointer :: normalBaroclinicVelAtSFC, normalBaroclinicVelAt250m, normalBaroclinicVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: tangentialBaroclinicVelAtSFC, tangentialBaroclinicVelAt250m, tangentialBaroclinicVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: zonalBaroclinicVelAtSFC, zonalBaroclinicVelAt250m, zonalBaroclinicVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: meridionalBaroclinicVelAtSFC, meridionalBaroclinicVelAt250m, meridionalBaroclinicVelAtBottom

      real (kind=RKIND), dimension(:), pointer :: normalGMBolusVelAtSFC, normalGMBolusVelAt250m, normalGMBolusVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: tangentialGMBolusVelAtSFC, tangentialGMBolusVelAt250m, tangentialGMBolusVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: zonalGMBolusVelAtSFC, zonalGMBolusVelAt250m, zonalGMBolusVelAtBottom
      real (kind=RKIND), dimension(:), pointer :: meridionalGMBolusVelAtSFC, meridionalGMBolusVelAt250m, meridionalGMBolusVelAtBottom

      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge
      real (kind=RKIND), dimension(:), pointer :: vertTransportVelocityAt250m, vertGMvelocityAt250m
      real (kind=RKIND), dimension(:), pointer :: vertVelSFC, vertTransportVelocitySFC, vertGMvelocitySFC
      real (kind=RKIND), dimension(:), pointer :: barotropicSpeed, columnIntegratedSpeed, dvEdge, dcEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell, relativeVorticityCell
      real (kind=RKIND), dimension(:,:), pointer :: activeTracersAt250m, activeTracersAvgTopto0100, activeTracersAvg0100to0250
      real (kind=RKIND), dimension(:,:), pointer :: activeTracersAtSurface, activeTracersAtBottom
      real (kind=RKIND), dimension(:,:), pointer :: activeTracersAvg0250to0700, activeTracersAvg0700to2000, activeTracersAvg2000toBottom
      real (kind=RKIND), dimension(:,:), pointer :: relativeVorticity, divergence, layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, tangentialVelocity, BruntVaisalaFreqTop
      real (kind=RKIND), dimension(:,:), pointer :: vertVelocityTop, vertGMBolusVelocityTop, vertTransportVelocityTop
      real (kind=RKIND), dimension(:,:), pointer :: normalGMBolusVelocity
      real (kind=RKIND), dimension(:,:), pointer :: normalBaroclinicVelocity ! baroclinic/barotropic pressume split explicit
      real (kind=RKIND), dimension(:), pointer :: normalBarotropicVelocity ! baroclinic/barotropic pressume split explicit

!maltrud 10 samples in top 100m
      real (kind=RKIND), dimension(:,:), pointer :: activeTracersAt5m, activeTracersAt10m, &
           activeTracersAt15m, activeTracersAt20m, activeTracersAt25m, activeTracersAt30m, &
           activeTracersAt40m, activeTracersAt50m, activeTracersAt70m, activeTracersAt95m
      real (kind=RKIND), dimension(:), pointer :: zonalVelCellAt5m, zonalVelCellAt10m, &
           zonalVelCellAt15m, zonalVelCellAt20m, zonalVelCellAt25m, zonalVelCellAt30m, &
           zonalVelCellAt40m, zonalVelCellAt50m, zonalVelCellAt70m, zonalVelCellAt95m
      real (kind=RKIND), dimension(:), pointer :: meridionalVelCellAt5m, meridionalVelCellAt10m, &
           meridionalVelCellAt15m, meridionalVelCellAt20m, meridionalVelCellAt25m, meridionalVelCellAt30m, &
           meridionalVelCellAt40m, meridionalVelCellAt50m, meridionalVelCellAt70m, meridionalVelCellAt95m
      real (kind=RKIND), dimension(2) ::     zedWeightsAt5m,  zedWeightsAt10m, &
           zedWeightsAt15m, zedWeightsAt20m, zedWeightsAt25m, zedWeightsAt30m, &
           zedWeightsAt40m, zedWeightsAt50m, zedWeightsAt70m, zedWeightsAt95m
      real (kind=RKIND), dimension(:), allocatable :: dZed, zedMid, dZedBot
      real (kind=RKIND) :: refDepth

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      character (len=StrKIND), pointer :: config_time_integrator

      ! scratch space
      type (field1DReal), pointer :: normalThicknessFluxSumField, layerThicknessSumEdgeField
      real (kind=RKIND), dimension(:), pointer :: normalThicknessFluxSum, layerThicknessSumEdge

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         ! get dimensions
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(domain % blocklist % structs, 'highFrequencyOutputAM', highFrequencyOutputAMPool)

         ! get static data from mesh pool
         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

         ! get arrays that will be 'sliced' and put into high frequency output
         call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticityCell', relativeVorticityCell, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticity', relativeVorticity, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'divergence', divergence, timeLevel)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'vertGMBolusVelocityTop', vertGMBolusVelocityTop)
         call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)
         call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)
         call mpas_pool_get_array(diagnosticsPool, 'tangentialVelocity', tangentialVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
         call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)

         ! get arrays that can be written to output at high freqency
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'kineticEnergyAt250m', kineticEnergyAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'kineticEnergyAtSurface', kineticEnergyAtSurface)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'relativeVorticityAt250m', relativeVorticityAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'divergenceAt250m', divergenceAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'relativeVorticityAtBottom', relativeVorticityAtBottom)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'divergenceAtBottom', divergenceAtBottom)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'kineticEnergyAtBottom', kineticEnergyAtBottom)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'relativeVorticityVertexAt250m', relativeVorticityVertexAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAtSurface', activeTracersAtSurface)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt250m', activeTracersAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAtBottom', activeTracersAtBottom)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAvgTopto0100', activeTracersAvgTopto0100)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAvg0100to0250', activeTracersAvg0100to0250)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAvg0250to0700', activeTracersAvg0250to0700)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAvg0700to2000', activeTracersAvg0700to2000)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAvg2000toBottom', activeTracersAvg2000toBottom)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertGMvelocitySFC', vertGMvelocitySFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertTransportVelocitySFC', vertTransportVelocitySFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertVelAt250m', vertVelAt250m)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalVelAtSFC', normalVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalVelAt250m', normalVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalVelAtBottom', normalVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialVelAtSFC', tangentialVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialVelAt250m', tangentialVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialVelAtBottom', tangentialVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelAtSFC', zonalVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelAt250m', zonalVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelAtBottom', zonalVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelAtSFC', meridionalVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelAt250m', meridionalVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelAtBottom', meridionalVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalAreaWeightedCellVelAtSFC', zonalAreaWeightedCellVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalAreaWeightedCellVelAt250m', zonalAreaWeightedCellVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalAreaWeightedCellVelAtBottom', zonalAreaWeightedCellVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalAreaWeightedCellVelAtSFC', meridionalAreaWeightedCellVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalAreaWeightedCellVelAt250m', meridionalAreaWeightedCellVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalAreaWeightedCellVelAtBottom', meridionalAreaWeightedCellVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalBarotropicVel', normalBarotropicVel)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialBarotropicVel', tangentialBarotropicVel)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalBarotropicVel', zonalBarotropicVel)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalBarotropicVel', meridionalBarotropicVel)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalBaroclinicVelAtSFC', normalBaroclinicVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalBaroclinicVelAt250m', normalBaroclinicVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalBaroclinicVelAtBottom', normalBaroclinicVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialBaroclinicVelAtSFC', tangentialBaroclinicVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialBaroclinicVelAt250m', tangentialBaroclinicVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialBaroclinicVelAtBottom', tangentialBaroclinicVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalBaroclinicVelAtSFC', zonalBaroclinicVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalBaroclinicVelAt250m', zonalBaroclinicVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalBaroclinicVelAtBottom', zonalBaroclinicVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalBaroclinicVelAtSFC', meridionalBaroclinicVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalBaroclinicVelAt250m', meridionalBaroclinicVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalBaroclinicVelAtBottom', meridionalBaroclinicVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalGMBolusVelAtSFC', normalGMBolusVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalGMBolusVelAt250m', normalGMBolusVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'normalGMBolusVelAtBottom', normalGMBolusVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialGMBolusVelAtSFC', tangentialGMBolusVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialGMBolusVelAt250m', tangentialGMBolusVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'tangentialGMBolusVelAtBottom', tangentialGMBolusVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalGMBolusVelAtSFC', zonalGMBolusVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalGMBolusVelAt250m', zonalGMBolusVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalGMBolusVelAtBottom', zonalGMBolusVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalGMBolusVelAtSFC', meridionalGMBolusVelAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalGMBolusVelAt250m', meridionalGMBolusVelAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalGMBolusVelAtBottom', meridionalGMBolusVelAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'BruntVaisalaFreqTopAtSFC', BruntVaisalaFreqTopAtSFC)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'BruntVaisalaFreqTopAt250m', BruntVaisalaFreqTopAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'BruntVaisalaFreqTopAtBottom', BruntVaisalaFreqTopAtBottom)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertGMvelocityAt250m', vertGMvelocityAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertTransportVelocityAt250m', vertTransportVelocityAt250m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'vertVelSFC', vertVelSFC)

!maltrud 10 samples in top 100m
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt5m',  activeTracersAt5m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt10m', activeTracersAt10m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt15m', activeTracersAt15m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt20m', activeTracersAt20m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt25m', activeTracersAt25m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt30m', activeTracersAt30m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt40m', activeTracersAt40m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt50m', activeTracersAt50m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt70m', activeTracersAt70m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'activeTracersAt95m', activeTracersAt95m)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelCellAt5m',  zonalVelCellAt5m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelCellAt10m', zonalVelCellAt10m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelCellAt15m', zonalVelCellAt15m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelCellAt20m', zonalVelCellAt20m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelCellAt25m', zonalVelCellAt25m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelCellAt30m', zonalVelCellAt30m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelCellAt40m', zonalVelCellAt40m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelCellAt50m', zonalVelCellAt50m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelCellAt70m', zonalVelCellAt70m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'zonalVelCellAt95m', zonalVelCellAt95m)

         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelCellAt5m',  meridionalVelCellAt5m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelCellAt10m', meridionalVelCellAt10m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelCellAt15m', meridionalVelCellAt15m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelCellAt20m', meridionalVelCellAt20m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelCellAt25m', meridionalVelCellAt25m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelCellAt30m', meridionalVelCellAt30m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelCellAt40m', meridionalVelCellAt40m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelCellAt50m', meridionalVelCellAt50m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelCellAt70m', meridionalVelCellAt70m)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'meridionalVelCellAt95m', meridionalVelCellAt95m)

         ! split explicit specific arrays
         call mpas_pool_get_config(ocnConfigs, 'config_time_integrator', config_time_integrator)
         if ( config_time_integrator == trim('split_explicit')) then
           call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocity, 1)
           call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocity, 1)
         endif

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
         call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
         call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
         call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'barotropicSpeed', barotropicSpeed)
         call mpas_pool_get_array(highFrequencyOutputAMPool, 'columnIntegratedSpeed', columnIntegratedSpeed)

         ! find vertical level that is just above the 100 m reference level
         iLevel0100 = 1
         do iLevel=2,nVertLevels
           if(refBottomDepth(iLevel) > 100.0_RKIND) then
              iLevel0100 = iLevel-1
              exit
           endif
         enddo

         ! find vertical level that is just above the 250 m reference level
         iLevel0250 = 1
         do iLevel=iLevel0100,nVertLevels
           if(refBottomDepth(iLevel) > 250.0_RKIND) then
              iLevel0250 = iLevel-1
              exit
           endif
         enddo

         ! find vertical level that is just above the 700 m reference level
         iLevel0700 = 1
         do iLevel=iLevel0250,nVertLevels
           if(refBottomDepth(iLevel) > 700.0_RKIND) then
              iLevel0700 = iLevel-1
              exit
           endif
         enddo

         ! find vertical level that is just above the 2000 m reference level
         iLevel2000 = 1
         do iLevel=iLevel0700,nVertLevels
           if(refBottomDepth(iLevel) > 2000.0_RKIND) then
              iLevel2000 = iLevel-1
              exit
           endif
         enddo

!maltrud 10 samples in top 100m
         allocate (dZed(nVertLevels), zedMid(nVertLevels), dZedBot(nVertLevels))

         dZed(1) = refBottomDepth(1)
         zedMid(1) = 0.5_RKIND*refBottomDepth(1)
         do iLevel=2,nVertLevels
           dZed(iLevel) = refBottomDepth(iLevel) - refBottomDepth(iLevel-1)
           zedMid(iLevel) = refBottomDepth(iLevel-1) + 0.5_RKIND*dZed(iLevel)
           dZedBot(iLevel-1) = 0.5_RKIND*(dZed(iLevel) + dZed(iLevel-1))
         enddo
         dZedBot(nVertLevels) = 0.5_RKIND*dZed(nVertLevels)
!maltrud debug
!write(*,*)' dZed, etc: ',dZed,zedMid,dZedBot

         ! find vertical level that is just above the 5 m reference level
         refDepth = 5.0_RKIND
         do iLevel=1, iLevel0100
           if(zedMid(iLevel) > refDepth) then
              iLevel0005(2) = iLevel
              exit
           endif
         enddo
         iLevel0005(1) = iLevel0005(2) - 1
         if (iLevel0005(1) == 0) then
           iLevel0005(1) = 1
           zedWeightsAt5m(1) = 0.0_RKIND
           zedWeightsAt5m(2) = 1.0_RKIND
         else
           zedWeightsAt5m(1) =  (zedMid(iLevel0005(2)) - refDepth)/dZedBot(iLevel0005(1))
           zedWeightsAt5m(2) = -(zedMid(iLevel0005(1)) - refDepth)/dZedBot(iLevel0005(1))
         endif
!maltrud debug
!write(*,*)' zedWeightsAt5m: ',iLevel0005,zedWeightsAt5m

         ! find vertical level that is just above the 10 m reference level
         refDepth = 10.0_RKIND
         do iLevel=1, iLevel0100
           if(zedMid(iLevel) > refDepth) then
              iLevel0010(2) = iLevel
              exit
           endif
         enddo
         iLevel0010(1) = iLevel0010(2) - 1
         if (iLevel0010(1) == 0) then
           iLevel0010(1) = 1
           zedWeightsAt10m(1) = 0.0_RKIND
           zedWeightsAt10m(2) = 1.0_RKIND
         else
           zedWeightsAt10m(1) =  (zedMid(iLevel0010(2)) - refDepth)/dZedBot(iLevel0010(1))
           zedWeightsAt10m(2) = -(zedMid(iLevel0010(1)) - refDepth)/dZedBot(iLevel0010(1))
         endif
!maltrud debug
!write(*,*)' zedWeightsAt10m: ',iLevel0010,zedWeightsAt10m

!maltrud from here on assume we are outside of the top level so no need to test for iLevel(1) = 0

         ! find vertical level that is just above the 15 m reference level
         refDepth = 15.0_RKIND
         do iLevel=1, iLevel0100
           if(zedMid(iLevel) > refDepth) then
              iLevel0015(2) = iLevel
              exit
           endif
         enddo
         iLevel0015(1) = iLevel0015(2) - 1
         zedWeightsAt15m(1) =  (zedMid(iLevel0015(2)) - refDepth)/dZedBot(iLevel0015(1))
         zedWeightsAt15m(2) = -(zedMid(iLevel0015(1)) - refDepth)/dZedBot(iLevel0015(1))
!maltrud debug
!write(*,*)' zedWeightsAt15m: ',iLevel0015,zedWeightsAt15m

         ! find vertical level that is just above the 20 m reference level
         refDepth = 20.0_RKIND
         do iLevel=1, iLevel0100
           if(zedMid(iLevel) > refDepth) then
              iLevel0020(2) = iLevel
              exit
           endif
         enddo
         iLevel0020(1) = iLevel0020(2) - 1
         zedWeightsAt20m(1) =  (zedMid(iLevel0020(2)) - refDepth)/dZedBot(iLevel0020(1))
         zedWeightsAt20m(2) = -(zedMid(iLevel0020(1)) - refDepth)/dZedBot(iLevel0020(1))
!maltrud debug
!write(*,*)' zedWeightsAt20m: ',iLevel0020,zedWeightsAt20m

         ! find vertical level that is just above the 25 m reference level
         refDepth = 25.0_RKIND
         do iLevel=1, iLevel0100
           if(zedMid(iLevel) > refDepth) then
              iLevel0025(2) = iLevel
              exit
           endif
         enddo
         iLevel0025(1) = iLevel0025(2) - 1
         zedWeightsAt25m(1) =  (zedMid(iLevel0025(2)) - refDepth)/dZedBot(iLevel0025(1))
         zedWeightsAt25m(2) = -(zedMid(iLevel0025(1)) - refDepth)/dZedBot(iLevel0025(1))
!maltrud debug
!write(*,*)' zedWeightsAt25m: ',iLevel0025,zedWeightsAt25m

         ! find vertical level that is just above the 30 m reference level
         refDepth = 30.0_RKIND
         do iLevel=1, iLevel0100
           if(zedMid(iLevel) > refDepth) then
              iLevel0030(2) = iLevel
              exit
           endif
         enddo
         iLevel0030(1) = iLevel0030(2) - 1
         zedWeightsAt30m(1) =  (zedMid(iLevel0030(2)) - refDepth)/dZedBot(iLevel0030(1))
         zedWeightsAt30m(2) = -(zedMid(iLevel0030(1)) - refDepth)/dZedBot(iLevel0030(1))
!maltrud debug
!write(*,*)' zedWeightsAt30m: ',iLevel0030,zedWeightsAt30m

         ! find vertical level that is just above the 40 m reference level
         refDepth = 40.0_RKIND
         do iLevel=1, iLevel0100
           if(zedMid(iLevel) > refDepth) then
              iLevel0040(2) = iLevel
              exit
           endif
         enddo
         iLevel0040(1) = iLevel0040(2) - 1
         zedWeightsAt40m(1) =  (zedMid(iLevel0040(2)) - refDepth)/dZedBot(iLevel0040(1))
         zedWeightsAt40m(2) = -(zedMid(iLevel0040(1)) - refDepth)/dZedBot(iLevel0040(1))
!maltrud debug
!write(*,*)' zedWeightsAt40m: ',iLevel0040,zedWeightsAt40m

         ! find vertical level that is just above the 50 m reference level
         refDepth = 50.0_RKIND
         do iLevel=1, iLevel0100
           if(zedMid(iLevel) > refDepth) then
              iLevel0050(2) = iLevel
              exit
           endif
         enddo
         iLevel0050(1) = iLevel0050(2) - 1
         zedWeightsAt50m(1) =  (zedMid(iLevel0050(2)) - refDepth)/dZedBot(iLevel0050(1))
         zedWeightsAt50m(2) = -(zedMid(iLevel0050(1)) - refDepth)/dZedBot(iLevel0050(1))
!maltrud debug
!write(*,*)' zedWeightsAt50m: ',iLevel0050,zedWeightsAt50m

         ! find vertical level that is just above the 70 m reference level
         refDepth = 70.0_RKIND
         do iLevel=1, iLevel0100
           if(zedMid(iLevel) > refDepth) then
              iLevel0070(2) = iLevel
              exit
           endif
         enddo
         iLevel0070(1) = iLevel0070(2) - 1
         zedWeightsAt70m(1) =  (zedMid(iLevel0070(2)) - refDepth)/dZedBot(iLevel0070(1))
         zedWeightsAt70m(2) = -(zedMid(iLevel0070(1)) - refDepth)/dZedBot(iLevel0070(1))
!maltrud debug
!write(*,*)' zedWeightsAt70m: ',iLevel0070,zedWeightsAt70m

         ! find vertical level that is just above the 95 m reference level
         refDepth = 95.0_RKIND
         do iLevel=1, iLevel0250
           if(zedMid(iLevel) > refDepth) then
              iLevel0095(2) = iLevel
              exit
           endif
         enddo
         iLevel0095(1) = iLevel0095(2) - 1
         zedWeightsAt95m(1) =  (zedMid(iLevel0095(2)) - refDepth)/dZedBot(iLevel0095(1))
         zedWeightsAt95m(2) = -(zedMid(iLevel0095(1)) - refDepth)/dZedBot(iLevel0095(1))
!maltrud debug
!write(*,*)' zedWeightsAt95m: ',iLevel0095,zedWeightsAt95m

         ! copy data into high frequency output fields
         kineticEnergyAt250m(:) = kineticEnergyCell(iLevel0250,:)
         relativeVorticityAt250m(:) = relativeVorticityCell(iLevel0250,:)
         divergenceAt250m(:) = divergence(iLevel0250,:)
         relativeVorticityVertexAt250m(:) = relativeVorticity(iLevel0250,:)
         kineticEnergyAtSurface(:) = kineticEnergyCell(1,:)
         activeTracersAtSurface(1,:) = activeTracers(1,1,:)
         activeTracersAtSurface(2,:) = activeTracers(2,1,:)
         activeTracersAt250m(1,:) = activeTracers(1,iLevel0250,:)
         activeTracersAt250m(2,:) = activeTracers(2,iLevel0250,:)
         vertGMvelocitySFC(:) = vertGMBolusVelocityTop(1,:)
         vertTransportVelocitySFC(:) = vertTransportVelocityTop(1,:)
         vertVelSFC(:) = vertVelocityTop(1,:)
         vertGMvelocityAt250m(:) = vertGMBolusVelocityTop(iLevel0250,:)
         vertTransportVelocityAt250m(:) = vertTransportVelocityTop(iLevel0250,:)
         vertVelAt250m(:) = vertVelocityTop(iLevel0250,:)

!maltrud 10 samples in top 100m
         activeTracersAt5m(:,:)   = activeTracers(:,iLevel0005(1),:)*zedWeightsAt5m(1) +  &
                                    activeTracers(:,iLevel0005(2),:)*zedWeightsAt5m(2)
         where (iLevel0005(1) >  maxLevelCell) activeTracersAt5m(1,:) = 0.0_RKIND
         where (iLevel0005(1) == maxLevelCell) activeTracersAt5m(1,:) = activeTracers(1,iLevel0005(1),:)
         where (iLevel0005(1) >  maxLevelCell) activeTracersAt5m(2,:) = 0.0_RKIND
         where (iLevel0005(1) == maxLevelCell) activeTracersAt5m(2,:) = activeTracers(2,iLevel0005(1),:)

         activeTracersAt10m(:,:)  = activeTracers(:,iLevel0010(1),:)*zedWeightsAt10m(1) +  &
                                    activeTracers(:,iLevel0010(2),:)*zedWeightsAt10m(2)
         where (iLevel0010(1) >  maxLevelCell) activeTracersAt10m(1,:) = 0.0_RKIND
         where (iLevel0010(1) == maxLevelCell) activeTracersAt10m(1,:) = activeTracers(1,iLevel0010(1),:)
         where (iLevel0010(1) >  maxLevelCell) activeTracersAt10m(2,:) = 0.0_RKIND
         where (iLevel0010(1) == maxLevelCell) activeTracersAt10m(2,:) = activeTracers(2,iLevel0010(1),:)

         activeTracersAt15m(:,:)  = activeTracers(:,iLevel0015(1),:)*zedWeightsAt15m(1) +  &
                                    activeTracers(:,iLevel0015(2),:)*zedWeightsAt15m(2)
         where (iLevel0015(1) >  maxLevelCell) activeTracersAt15m(1,:) = 0.0_RKIND
         where (iLevel0015(1) == maxLevelCell) activeTracersAt15m(1,:) = activeTracers(1,iLevel0015(1),:)
         where (iLevel0015(1) >  maxLevelCell) activeTracersAt15m(2,:) = 0.0_RKIND
         where (iLevel0015(1) == maxLevelCell) activeTracersAt15m(2,:) = activeTracers(2,iLevel0015(1),:)

         activeTracersAt20m(:,:)  = activeTracers(:,iLevel0020(1),:)*zedWeightsAt20m(1) +  &
                                    activeTracers(:,iLevel0020(2),:)*zedWeightsAt20m(2)
         where (iLevel0020(1) >  maxLevelCell) activeTracersAt20m(1,:) = 0.0_RKIND
         where (iLevel0020(1) == maxLevelCell) activeTracersAt20m(1,:) = activeTracers(1,iLevel0020(1),:)
         where (iLevel0020(1) >  maxLevelCell) activeTracersAt20m(2,:) = 0.0_RKIND
         where (iLevel0020(1) == maxLevelCell) activeTracersAt20m(2,:) = activeTracers(2,iLevel0020(1),:)

         activeTracersAt25m(:,:)  = activeTracers(:,iLevel0025(1),:)*zedWeightsAt25m(1) +  &
                                    activeTracers(:,iLevel0025(2),:)*zedWeightsAt25m(2)
         where (iLevel0025(1) >  maxLevelCell) activeTracersAt25m(1,:) = 0.0_RKIND
         where (iLevel0025(1) == maxLevelCell) activeTracersAt25m(1,:) = activeTracers(1,iLevel0025(1),:)
         where (iLevel0025(1) >  maxLevelCell) activeTracersAt25m(2,:) = 0.0_RKIND
         where (iLevel0025(1) == maxLevelCell) activeTracersAt25m(2,:) = activeTracers(2,iLevel0025(1),:)

         activeTracersAt30m(:,:)  = activeTracers(:,iLevel0030(1),:)*zedWeightsAt30m(1) +  &
                                    activeTracers(:,iLevel0030(2),:)*zedWeightsAt30m(2)
         where (iLevel0030(1) >  maxLevelCell) activeTracersAt30m(1,:) = 0.0_RKIND
         where (iLevel0030(1) == maxLevelCell) activeTracersAt30m(1,:) = activeTracers(1,iLevel0030(1),:)
         where (iLevel0030(1) >  maxLevelCell) activeTracersAt30m(2,:) = 0.0_RKIND
         where (iLevel0030(1) == maxLevelCell) activeTracersAt30m(2,:) = activeTracers(2,iLevel0030(1),:)

         activeTracersAt40m(:,:)  = activeTracers(:,iLevel0040(1),:)*zedWeightsAt40m(1) +  &
                                    activeTracers(:,iLevel0040(2),:)*zedWeightsAt40m(2)
         where (iLevel0040(1) >  maxLevelCell) activeTracersAt40m(1,:) = 0.0_RKIND
         where (iLevel0040(1) == maxLevelCell) activeTracersAt40m(1,:) = activeTracers(1,iLevel0040(1),:)
         where (iLevel0040(1) >  maxLevelCell) activeTracersAt40m(2,:) = 0.0_RKIND
         where (iLevel0040(1) == maxLevelCell) activeTracersAt40m(2,:) = activeTracers(2,iLevel0040(1),:)

         activeTracersAt50m(:,:)  = activeTracers(:,iLevel0050(1),:)*zedWeightsAt50m(1) +  &
                                    activeTracers(:,iLevel0050(2),:)*zedWeightsAt50m(2)
         where (iLevel0050(1) >  maxLevelCell) activeTracersAt50m(1,:) = 0.0_RKIND
         where (iLevel0050(1) == maxLevelCell) activeTracersAt50m(1,:) = activeTracers(1,iLevel0050(1),:)
         where (iLevel0050(1) >  maxLevelCell) activeTracersAt50m(2,:) = 0.0_RKIND
         where (iLevel0050(1) == maxLevelCell) activeTracersAt50m(2,:) = activeTracers(2,iLevel0050(1),:)

         activeTracersAt70m(:,:)  = activeTracers(:,iLevel0070(1),:)*zedWeightsAt70m(1) +  &
                                    activeTracers(:,iLevel0070(2),:)*zedWeightsAt70m(2)
         where (iLevel0070(1) >  maxLevelCell) activeTracersAt70m(1,:) = 0.0_RKIND
         where (iLevel0070(1) == maxLevelCell) activeTracersAt70m(1,:) = activeTracers(1,iLevel0070(1),:)
         where (iLevel0070(1) >  maxLevelCell) activeTracersAt70m(2,:) = 0.0_RKIND
         where (iLevel0070(1) == maxLevelCell) activeTracersAt70m(2,:) = activeTracers(2,iLevel0070(1),:)

         activeTracersAt95m(:,:)  = activeTracers(:,iLevel0095(1),:)*zedWeightsAt95m(1) +  &
                                    activeTracers(:,iLevel0095(2),:)*zedWeightsAt95m(2)
         where (iLevel0095(1) >  maxLevelCell) activeTracersAt95m(1,:) = 0.0_RKIND
         where (iLevel0095(1) == maxLevelCell) activeTracersAt95m(1,:) = activeTracers(1,iLevel0095(1),:)
         where (iLevel0095(1) >  maxLevelCell) activeTracersAt95m(2,:) = 0.0_RKIND
         where (iLevel0095(1) == maxLevelCell) activeTracersAt95m(2,:) = activeTracers(2,iLevel0095(1),:)

         do iEdge = 1, nEdges
           normalGMBolusVelAtSFC(iEdge) = normalGMBolusVelocity(1,iEdge)
           normalGMBolusVelAt250m(iEdge) = normalGMBolusVelocity(iLevel0250,iEdge)
           normalGMBolusVelAtBottom(iEdge) = normalGMBolusVelocity(maxLevelEdgeBot(iEdge),iEdge)

           tangentialVelAtSFC(iEdge) = tangentialVelocity(1,iEdge)
           tangentialVelAt250m(iEdge) = tangentialVelocity(iLevel0250,iEdge)
           tangentialVelAtBottom(iEdge) = tangentialVelocity(maxLevelEdgeBot(iEdge),iEdge)

           normalVelAtSFC(iEdge) = normalVelocity(1,iEdge)
           normalVelAt250m(iEdge) = normalVelocity(iLevel0250,iEdge)
           normalVelAtBottom(iEdge) = normalVelocity(maxLevelEdgeBot(iEdge),iEdge)

           zonalVelAtSFC(iEdge) = normalVelAtSFC(iEdge)*cos(angleEdge(iEdge)) &
                                  - tangentialVelAtSFC(iEdge)*sin(angleEdge(iEdge))
           meridionalVelAtSFC(iEdge) = normalVelAtSFC(iEdge)*sin(angleEdge(iEdge)) &
                                       + tangentialVelAtSFC(iEdge)*cos(angleEdge(iEdge))

           zonalVelAt250m(iEdge) = normalVelAt250m(iEdge)*cos(angleEdge(iEdge)) &
                                  - tangentialVelAt250m(iEdge)*sin(angleEdge(iEdge))
           meridionalVelAt250m(iEdge) = normalVelAt250m(iEdge)*sin(angleEdge(iEdge)) &
                                       + tangentialVelAt250m(iEdge)*cos(angleEdge(iEdge))

           zonalVelAtBottom(iEdge) = normalVelAtBottom(iEdge)*cos(angleEdge(iEdge)) &
                                  - tangentialVelAtBottom(iEdge)*sin(angleEdge(iEdge))
           meridionalVelAtBottom(iEdge) = normalVelAtBottom(iEdge)*sin(angleEdge(iEdge)) &
                                       + tangentialVelAtBottom(iEdge)*cos(angleEdge(iEdge))
         end do

         do iEdge = 1, nEdges
           tangentialGMBolusVelAtSFC(iEdge) = 0.0_RKIND
           tangentialGMBolusVelAt250m(iEdge) = 0.0_RKIND
           tangentialGMBolusVelAtBottom(iEdge) = 0.0_RKIND
           do i = 1, nEdgesOnEdge(iEdge)
             eoe = edgesOnEdge(i,iEdge)
             tangentialGMBolusVelAtSFC(iEdge) = tangentialGMBolusVelAtSFC(iEdge) + weightsOnEdge(i,iEdge) * normalGMBolusVelAtSFC(eoe)
             tangentialGMBolusVelAt250m(iEdge) = tangentialGMBolusVelAt250m(iEdge) + weightsOnEdge(i,iEdge) * normalGMBolusVelAt250m(eoe)
             tangentialGMBolusVelAtBottom(iEdge) = tangentialGMBolusVelAtBottom(iEdge) + weightsOnEdge(i,iEdge) * normalGMBolusVelAtBottom(eoe)
           end do

           zonalGMBolusVelAtSFC(iEdge) = normalGMBolusVelAtSFC(iEdge)*cos(angleEdge(iEdge)) &
             - tangentialGMBolusVelAtSFC(iEdge)*sin(angleEdge(iEdge))
           meridionalGMBolusVelAtSFC(iEdge) = normalGMBolusVelAtSFC(iEdge)*sin(angleEdge(iEdge)) &
             + tangentialGMBolusVelAtSFC(iEdge)*cos(angleEdge(iEdge))

           zonalGMBolusVelAt250m(iEdge) = normalGMBolusVelAt250m(iEdge)*cos(angleEdge(iEdge)) &
             - tangentialGMBolusVelAt250m(iEdge)*sin(angleEdge(iEdge))
           meridionalGMBolusVelAt250m(iEdge) = normalGMBolusVelAt250m(iEdge)*sin(angleEdge(iEdge)) &
             + tangentialGMBolusVelAt250m(iEdge)*cos(angleEdge(iEdge))

           zonalGMBolusVelAtBottom(iEdge) = normalGMBolusVelAtBottom(iEdge)*cos(angleEdge(iEdge)) &
             - tangentialGMBolusVelAtBottom(iEdge)*sin(angleEdge(iEdge))
           meridionalGMBolusVelAtBottom(iEdge) = normalGMBolusVelAtBottom(iEdge)*sin(angleEdge(iEdge)) &
             + tangentialGMBolusVelAtBottom(iEdge)*cos(angleEdge(iEdge))
         end do

         do iCell=1,nCells
           relativeVorticityAtBottom(iCell) = relativeVorticityCell(maxLevelCell(iCell),iCell)
           divergenceAtBottom(iCell) = divergence(maxLevelCell(iCell),iCell)
           kineticEnergyAtBottom(iCell) = kineticEnergyCell(maxLevelCell(iCell),iCell)
           activeTracersAtBottom(1,iCell) = activeTracers(1,maxLevelCell(iCell),iCell)
           activeTracersAtBottom(2,iCell) = activeTracers(2,maxLevelCell(iCell),iCell)

           BruntVaisalaFreqTopAtSFC(iCell) = BruntVaisalaFreqTop(1,iCell)
           BruntVaisalaFreqTopAt250m(iCell) = BruntVaisalaFreqTop(iLevel0250,iCell)
           BruntVaisalaFreqTopAtBottom(iCell) = BruntVaisalaFreqTop(maxLevelCell(iCell),iCell)
         end do

         do iCell = 1, nCells
           zonalAreaWeightedCellVelAtSFC(iCell) = 0.0_RKIND
           zonalAreaWeightedCellVelAt250m(iCell) = 0.0_RKIND
           zonalAreaWeightedCellVelAtBottom(iCell) = 0.0_RKIND

           meridionalAreaWeightedCellVelAtSFC(iCell) = 0.0_RKIND
           meridionalAreaWeightedCellVelAt250m(iCell) = 0.0_RKIND
           meridionalAreaWeightedCellVelAtBottom(iCell) = 0.0_RKIND

           zonalVelCellAt5m(iCell)       = 0.0_RKIND
           meridionalVelCellAt5m(iCell)  = 0.0_RKIND
           zonalVelCellAt10m(iCell)      = 0.0_RKIND
           meridionalVelCellAt10m(iCell) = 0.0_RKIND
           zonalVelCellAt15m(iCell)      = 0.0_RKIND
           meridionalVelCellAt10m(iCell) = 0.0_RKIND
           zonalVelCellAt20m(iCell)      = 0.0_RKIND
           meridionalVelCellAt20m(iCell) = 0.0_RKIND
           zonalVelCellAt25m(iCell)      = 0.0_RKIND
           meridionalVelCellAt25m(iCell) = 0.0_RKIND
           zonalVelCellAt30m(iCell)      = 0.0_RKIND
           meridionalVelCellAt30m(iCell) = 0.0_RKIND
           zonalVelCellAt40m(iCell)      = 0.0_RKIND
           meridionalVelCellAt40m(iCell) = 0.0_RKIND
           zonalVelCellAt50m(iCell)      = 0.0_RKIND
           meridionalVelCellAt50m(iCell) = 0.0_RKIND
           zonalVelCellAt70m(iCell)      = 0.0_RKIND
           meridionalVelCellAt70m(iCell) = 0.0_RKIND
           zonalVelCellAt95m(iCell)      = 0.0_RKIND
           meridionalVelCellAt95m(iCell) = 0.0_RKIND

           cellArea = 0.0_RKIND
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               ! note, we don't necessarily have d_{i,j} computed (Fringer, 2006), so we allow averaging between cells via 0.5_RKIND*dcEdge(iEdge) to compute
               ! segment of cell area d_{i,j}-- the key is that the edge-normal velocities are a convex combination
               cellArea = cellArea + 0.5_RKIND*dcEdge(iEdge)*dvEdge(iEdge)

               weightedNormalVel = normalVelAtSFC(iEdge)*0.5_RKIND*dcEdge(iEdge)*dvEdge(iEdge)
               zonalAreaWeightedCellVelAtSFC(iCell) = zonalAreaWeightedCellVelAtSFC(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalAreaWeightedCellVelAtSFC(iCell) = meridionalAreaWeightedCellVelAtSFC(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel

               weightedNormalVel = normalVelAt250m(iEdge)*0.5_RKIND*dcEdge(iEdge)*dvEdge(iEdge)
               zonalAreaWeightedCellVelAt250m(iCell) = zonalAreaWeightedCellVelAt250m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalAreaWeightedCellVelAt250m(iCell) = meridionalAreaWeightedCellVelAt250m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel

               weightedNormalVel = normalVelAtBottom(iEdge)*0.5_RKIND*dcEdge(iEdge)*dvEdge(iEdge)
               zonalAreaWeightedCellVelAtBottom(iCell) = zonalAreaWeightedCellVelAtBottom(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalAreaWeightedCellVelAtBottom(iCell) = meridionalAreaWeightedCellVelAtBottom(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel

!maltrud 10 samples in top 100m
! do each level separately so we can reuse weightedNormalVel

! 5m
               weightedNormalVel = normalVelocity(iLevel0005(1),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt5m(1)
               if (iLevel0005(1) >  maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               if (iLevel0005(1) == maxLevelCell(iCell)) weightedNormalVel = weightedNormalVel/zedWeightsAt5m(1)
               zonalVelCellAt5m(iCell)       =      zonalVelCellAt5m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt5m(iCell)  = meridionalVelCellAt5m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
               weightedNormalVel = normalVelocity(iLevel0005(2),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt5m(2)
               if (iLevel0005(1) >= maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               zonalVelCellAt5m(iCell)       =      zonalVelCellAt5m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt5m(iCell)  = meridionalVelCellAt5m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
! 10m
               weightedNormalVel = normalVelocity(iLevel0010(1),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt10m(1)
               if (iLevel0010(1) >  maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               if (iLevel0010(1) == maxLevelCell(iCell)) weightedNormalVel = weightedNormalVel/zedWeightsAt10m(1)
               zonalVelCellAt10m(iCell)      =      zonalVelCellAt10m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt10m(iCell) = meridionalVelCellAt10m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
               weightedNormalVel = normalVelocity(iLevel0010(2),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt10m(2)
               if (iLevel0010(1) >= maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               zonalVelCellAt10m(iCell)      =      zonalVelCellAt10m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt10m(iCell) = meridionalVelCellAt10m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
! 15m
               weightedNormalVel = normalVelocity(iLevel0015(1),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt15m(1)
               if (iLevel0015(1) >  maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               if (iLevel0015(1) == maxLevelCell(iCell)) weightedNormalVel = weightedNormalVel/zedWeightsAt15m(1)
               zonalVelCellAt15m(iCell)      =      zonalVelCellAt15m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt15m(iCell) = meridionalVelCellAt15m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
               weightedNormalVel = normalVelocity(iLevel0015(2),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt15m(2)
               if (iLevel0015(1) >= maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               zonalVelCellAt15m(iCell)      =      zonalVelCellAt15m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt15m(iCell) = meridionalVelCellAt15m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
! 20m
               weightedNormalVel = normalVelocity(iLevel0020(1),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt20m(1)
               if (iLevel0020(1) >  maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               if (iLevel0020(1) == maxLevelCell(iCell)) weightedNormalVel = weightedNormalVel/zedWeightsAt20m(1)
               zonalVelCellAt20m(iCell)      =      zonalVelCellAt20m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt20m(iCell) = meridionalVelCellAt20m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
               weightedNormalVel = normalVelocity(iLevel0020(2),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt20m(2)
               if (iLevel0020(1) >= maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               zonalVelCellAt20m(iCell)      =      zonalVelCellAt20m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt20m(iCell) = meridionalVelCellAt20m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
! 25m
               weightedNormalVel = normalVelocity(iLevel0025(1),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt25m(1)
               if (iLevel0025(1) >  maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               if (iLevel0025(1) == maxLevelCell(iCell)) weightedNormalVel = weightedNormalVel/zedWeightsAt25m(1)
               zonalVelCellAt25m(iCell)      =      zonalVelCellAt25m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt25m(iCell) = meridionalVelCellAt25m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
               weightedNormalVel = normalVelocity(iLevel0025(2),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt25m(2)
               if (iLevel0025(1) >= maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               zonalVelCellAt25m(iCell)      =      zonalVelCellAt25m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt25m(iCell) = meridionalVelCellAt25m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
! 30m
               weightedNormalVel = normalVelocity(iLevel0030(1),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt30m(1)
               if (iLevel0030(1) >  maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               if (iLevel0030(1) == maxLevelCell(iCell)) weightedNormalVel = weightedNormalVel/zedWeightsAt30m(1)
               zonalVelCellAt30m(iCell)      =      zonalVelCellAt30m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt30m(iCell) = meridionalVelCellAt30m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
               weightedNormalVel = normalVelocity(iLevel0030(2),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt30m(2)
               if (iLevel0030(1) >= maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               zonalVelCellAt30m(iCell)      =      zonalVelCellAt30m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt30m(iCell) = meridionalVelCellAt30m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
! 40m
               weightedNormalVel = normalVelocity(iLevel0040(1),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt40m(1)
               if (iLevel0040(1) >  maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               if (iLevel0040(1) == maxLevelCell(iCell)) weightedNormalVel = weightedNormalVel/zedWeightsAt40m(1)
               zonalVelCellAt40m(iCell)      =      zonalVelCellAt40m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt40m(iCell) = meridionalVelCellAt40m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
               weightedNormalVel = normalVelocity(iLevel0040(2),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt40m(2)
               if (iLevel0040(1) >= maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               zonalVelCellAt40m(iCell)      =      zonalVelCellAt40m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt40m(iCell) = meridionalVelCellAt40m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
! 50m
               weightedNormalVel = normalVelocity(iLevel0050(1),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt50m(1)
               if (iLevel0050(1) >  maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               if (iLevel0050(1) == maxLevelCell(iCell)) weightedNormalVel = weightedNormalVel/zedWeightsAt50m(1)
               zonalVelCellAt50m(iCell)      =      zonalVelCellAt50m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt50m(iCell) = meridionalVelCellAt50m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
               weightedNormalVel = normalVelocity(iLevel0050(2),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt50m(2)
               if (iLevel0050(1) >= maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               zonalVelCellAt50m(iCell)      =      zonalVelCellAt50m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt50m(iCell) = meridionalVelCellAt50m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
! 70m
               weightedNormalVel = normalVelocity(iLevel0070(1),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt70m(1)
               if (iLevel0070(1) >  maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               if (iLevel0070(1) == maxLevelCell(iCell)) weightedNormalVel = weightedNormalVel/zedWeightsAt70m(1)
               zonalVelCellAt70m(iCell)      =      zonalVelCellAt70m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt70m(iCell) = meridionalVelCellAt70m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
               weightedNormalVel = normalVelocity(iLevel0070(2),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt70m(2)
               if (iLevel0070(1) >= maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               zonalVelCellAt70m(iCell)      =      zonalVelCellAt70m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt70m(iCell) = meridionalVelCellAt70m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
! 95m
               weightedNormalVel = normalVelocity(iLevel0095(1),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt95m(1)
               if (iLevel0095(1) >  maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               if (iLevel0095(1) == maxLevelCell(iCell)) weightedNormalVel = weightedNormalVel/zedWeightsAt95m(1)
               zonalVelCellAt95m(iCell)      =      zonalVelCellAt95m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt95m(iCell) = meridionalVelCellAt95m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel
               weightedNormalVel = normalVelocity(iLevel0095(2),iEdge)*0.5_RKIND*dcEdge(iEdge)  &
                                   *dvEdge(iEdge)*zedWeightsAt95m(2)
               if (iLevel0095(1) >= maxLevelCell(iCell)) weightedNormalVel = 0.0_RKIND
               zonalVelCellAt95m(iCell)      =      zonalVelCellAt95m(iCell) + cos(angleEdge(iEdge)) * weightedNormalVel
               meridionalVelCellAt95m(iCell) = meridionalVelCellAt95m(iCell) + sin(angleEdge(iEdge)) * weightedNormalVel

           end do

           zonalAreaWeightedCellVelAtSFC(iCell) = zonalAreaWeightedCellVelAtSFC(iCell)/cellArea
           meridionalAreaWeightedCellVelAtSFC(iCell) = meridionalAreaWeightedCellVelAtSFC(iCell)/cellArea

           zonalAreaWeightedCellVelAt250m(iCell) = zonalAreaWeightedCellVelAt250m(iCell)/cellArea
           meridionalAreaWeightedCellVelAt250m(iCell) = meridionalAreaWeightedCellVelAt250m(iCell)/cellArea

           zonalAreaWeightedCellVelAtBottom(iCell) = zonalAreaWeightedCellVelAtBottom(iCell)/cellArea
           meridionalAreaWeightedCellVelAtBottom(iCell) = meridionalAreaWeightedCellVelAtBottom(iCell)/cellArea

!maltrud 10 samples in top 100m
           zonalVelCellAt5m(iCell)       = zonalVelCellAt5m(iCell)/cellArea
           meridionalVelCellAt5m(iCell)  = meridionalVelCellAt5m(iCell)/cellArea

           zonalVelCellAt10m(iCell)      = zonalVelCellAt10m(iCell)/cellArea
           meridionalVelCellAt10m(iCell) = meridionalVelCellAt10m(iCell)/cellArea

           zonalVelCellAt15m(iCell)      = zonalVelCellAt15m(iCell)/cellArea
           meridionalVelCellAt15m(iCell) = meridionalVelCellAt15m(iCell)/cellArea

           zonalVelCellAt20m(iCell)      = zonalVelCellAt20m(iCell)/cellArea
           meridionalVelCellAt20m(iCell) = meridionalVelCellAt20m(iCell)/cellArea

           zonalVelCellAt25m(iCell)      = zonalVelCellAt25m(iCell)/cellArea
           meridionalVelCellAt25m(iCell) = meridionalVelCellAt25m(iCell)/cellArea

           zonalVelCellAt30m(iCell)      = zonalVelCellAt30m(iCell)/cellArea
           meridionalVelCellAt30m(iCell) = meridionalVelCellAt30m(iCell)/cellArea

           zonalVelCellAt40m(iCell)      = zonalVelCellAt40m(iCell)/cellArea
           meridionalVelCellAt40m(iCell) = meridionalVelCellAt40m(iCell)/cellArea

           zonalVelCellAt50m(iCell)      = zonalVelCellAt50m(iCell)/cellArea
           meridionalVelCellAt50m(iCell) = meridionalVelCellAt50m(iCell)/cellArea

           zonalVelCellAt70m(iCell)      = zonalVelCellAt70m(iCell)/cellArea
           meridionalVelCellAt70m(iCell) = meridionalVelCellAt70m(iCell)/cellArea

           zonalVelCellAt95m(iCell)      = zonalVelCellAt95m(iCell)/cellArea
           meridionalVelCellAt95m(iCell) = meridionalVelCellAt95m(iCell)/cellArea

         end do

         if ( config_time_integrator == trim('split_explicit')) then

           do iEdge = 1, nEdges
             normalBarotropicVel(iEdge) = normalBarotropicVelocity(iEdge)
             normalBaroclinicVelAtSFC(iEdge) = normalBaroclinicVelocity(1,iEdge)
             normalBaroclinicVelAt250m(iEdge) = normalBaroclinicVelocity(iLevel0250,iEdge)
             normalBaroclinicVelAtBottom(iEdge) = normalBaroclinicVelocity(maxLevelEdgeBot(iEdge),iEdge)
           end do

           do iEdge = 1, nEdges
             tangentialBarotropicVel(iEdge) = 0.0_RKIND
             tangentialBaroclinicVelAtSFC(iEdge) = 0.0_RKIND
             tangentialBaroclinicVelAt250m(iEdge) = 0.0_RKIND
             tangentialBaroclinicVelAtBottom(iEdge) = 0.0_RKIND
             do i = 1, nEdgesOnEdge(iEdge)
               eoe = edgesOnEdge(i,iEdge)
               tangentialBarotropicVel(iEdge) = tangentialBarotropicVel(iEdge) + weightsOnEdge(i,iEdge) * normalBarotropicVel(eoe)
               tangentialBaroclinicVelAtSFC(iEdge) = tangentialBaroclinicVelAtSFC(iEdge) + weightsOnEdge(i,iEdge) * normalBaroclinicVelAtSFC(eoe)
               tangentialBaroclinicVelAt250m(iEdge) = tangentialBaroclinicVelAt250m(iEdge) + weightsOnEdge(i,iEdge) * normalBaroclinicVelAt250m(eoe)
               tangentialBaroclinicVelAtBottom(iEdge) = tangentialBaroclinicVelAtBottom(iEdge) + weightsOnEdge(i,iEdge) * normalBaroclinicVelAtBottom(eoe)
             end do

             zonalBarotropicVel(iEdge) = normalBarotropicVel(iEdge)*cos(angleEdge(iEdge)) &
               - tangentialBarotropicVel(iEdge)*sin(angleEdge(iEdge))
             meridionalBarotropicVel(iEdge) = normalBarotropicVel(iEdge)*sin(angleEdge(iEdge)) &
               + tangentialBarotropicVel(iEdge)*cos(angleEdge(iEdge))

             zonalBaroclinicVelAtSFC(iEdge) = normalBaroclinicVelAtSFC(iEdge)*cos(angleEdge(iEdge)) &
               - tangentialBaroclinicVelAtSFC(iEdge)*sin(angleEdge(iEdge))
             meridionalBaroclinicVelAtSFC(iEdge) = normalBaroclinicVelAtSFC(iEdge)*sin(angleEdge(iEdge)) &
               + tangentialBaroclinicVelAtSFC(iEdge)*cos(angleEdge(iEdge))

             zonalBaroclinicVelAt250m(iEdge) = normalBaroclinicVelAt250m(iEdge)*cos(angleEdge(iEdge)) &
               - tangentialBaroclinicVelAt250m(iEdge)*sin(angleEdge(iEdge))
             meridionalBaroclinicVelAt250m(iEdge) = normalBaroclinicVelAt250m(iEdge)*sin(angleEdge(iEdge)) &
               + tangentialBaroclinicVelAt250m(iEdge)*cos(angleEdge(iEdge))

             zonalBaroclinicVelAtBottom(iEdge) = normalBaroclinicVelAtBottom(iEdge)*cos(angleEdge(iEdge)) &
               - tangentialBaroclinicVelAtBottom(iEdge)*sin(angleEdge(iEdge))
             meridionalBaroclinicVelAtBottom(iEdge) = normalBaroclinicVelAtBottom(iEdge)*sin(angleEdge(iEdge)) &
               + tangentialBaroclinicVelAtBottom(iEdge)*cos(angleEdge(iEdge))
           end do

         endif

         !
         ! compute layer averaged tracers
         !
         activeTracersAvgTopto0100(:,:) = -1.0e34_RKIND
         activeTracersAvg0100to0250(:,:) = -1.0e34_RKIND
         activeTracersAvg0250to0700(:,:) = -1.0e34_RKIND
         activeTracersAvg0700to2000(:,:) = -1.0e34_RKIND
         activeTracersAvg2000toBottom(:,:) = -1.0e34_RKIND
         do iCell = 1, nCells
           sumLayerThickness = layerThickness(1,iCell)
           activeTracersAvgTopto0100(:,iCell) = activeTracers(:,1,iCell)*layerThickness(1,iCell)
           do k=2, min(maxLevelCell(iCell),iLevel0100)
              sumLayerThickness = sumLayerThickness + layerThickness(k,iCell)
              activeTracersAvgTopto0100(:,iCell) = activeTracersAvgTopto0100(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
           enddo
           activeTracersAvgTopto0100(:,iCell) = activeTracersAvgTopto0100(:,iCell) / max(sumLayerThickness,1.0_RKIND)

           if (iLevel0100+1.le.maxLevelCell(iCell)) then
              sumLayerThickness = layerThickness(iLevel0100+1,iCell)
              activeTracersAvg0100to0250(:,iCell) = activeTracers(:,iLevel0100+1,iCell)*layerThickness(iLevel0100,iCell)
              do k=iLevel0100+2, min(maxLevelCell(iCell),iLevel0250)
                 sumLayerThickness = sumLayerThickness + layerThickness(k,iCell)
                 activeTracersAvg0100to0250(:,iCell) = activeTracersAvg0100to0250(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
              enddo
              activeTracersAvg0100to0250(:,iCell) = activeTracersAvg0100to0250(:,iCell) / max(sumLayerThickness,1.0_RKIND)
           endif

           if (iLevel0250+1.le.maxLevelCell(iCell)) then
              sumLayerThickness = layerThickness(iLevel0250+1,iCell)
              activeTracersAvg0250to0700(:,iCell) = activeTracers(:,iLevel0250+1,iCell)*layerThickness(iLevel0250+1,iCell)
              do k=iLevel0250+2, min(maxLevelCell(iCell),iLevel0700)
                 sumLayerThickness = sumLayerThickness + layerThickness(k,iCell)
                 activeTracersAvg0250to0700(:,iCell) = activeTracersAvg0250to0700(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
              enddo
              activeTracersAvg0250to0700(:,iCell) = activeTracersAvg0250to0700(:,iCell) / max(sumLayerThickness,1.0_RKIND)
           endif

           if (iLevel0700+1.le.maxLevelCell(iCell)) then
              sumLayerThickness = layerThickness(iLevel0700+1,iCell)
              activeTracersAvg0700to2000(:,iCell) = activeTracers(:,iLevel0700+1,iCell)*layerThickness(iLevel0700+1,iCell)
              do k=iLevel0700+2, min(maxLevelCell(iCell),iLevel2000)
                 sumLayerThickness = sumLayerThickness + layerThickness(k,iCell)
                 activeTracersAvg0700to2000(:,iCell) = activeTracersAvg0700to2000(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
              enddo
              activeTracersAvg0700to2000(:,iCell) = activeTracersAvg0700to2000(:,iCell) / max(sumLayerThickness,1.0_RKIND)
           endif

           if (iLevel2000+1.le.maxLevelCell(iCell)) then
              sumLayerThickness = layerThickness(iLevel2000+1,iCell)
              activeTracersAvg2000toBottom(:,iCell) = activeTracers(:,iLevel2000+1,iCell)*layerThickness(iLevel2000+1,iCell)
              do k=iLevel2000+2, maxLevelCell(iCell)
                 sumLayerThickness = sumLayerThickness + layerThickness(k,iCell)
                 activeTracersAvg2000toBottom(:,iCell) = activeTracersAvg2000toBottom(:,iCell) + activeTracers(:,k,iCell)*layerThickness(k,iCell)
              enddo
              activeTracersAvg2000toBottom(:,iCell) = activeTracersAvg2000toBottom(:,iCell) / max(sumLayerThickness,1.0_RKIND)
           endif
         enddo

         !
         ! compute barotropic kinetic energy
         !

         ! scratch arrays
         call mpas_pool_get_field(scratchPool, 'normalThicknessFluxSum', normalThicknessFluxSumField)
         call mpas_allocate_scratch_field(normalThicknessFluxSumField, .true.)
         normalThicknessFluxSum => normalThicknessFluxSumField % array

         call mpas_pool_get_field(scratchPool, 'layerThicknessSumEdge', layerThicknessSumEdgeField)
         call mpas_allocate_scratch_field(layerThicknessSumEdgeField, .true.)
         layerThicknessSumEdge => layerThicknessSumEdgeField % array

         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            ! normal Barotropic Velocity = sum(h*u)/sum(h) on each edge
            layerThicknessEdge1 = 0.5_RKIND*( layerThickness(1,cell1) + layerThickness(1,cell2) )
            normalThicknessFluxSum(iEdge) = layerThicknessEdge1 * normalVelocity(1,iEdge)
            layerThicknessSumEdge(iEdge) = layerThicknessEdge1

            do k=2, maxLevelEdgeTop(iEdge)
               layerThicknessEdge1 = 0.5_RKIND*( layerThickness(k,cell1) + layerThickness(k,cell2) )

               normalThicknessFluxSum(iEdge) = normalThicknessFluxSum(iEdge) + layerThicknessEdge1 * normalVelocity(k,iEdge)
               layerThicknessSumEdge(iEdge) = layerThicknessSumEdge(iEdge) + layerThicknessEdge1
            enddo
         enddo

         do iCell = 1, nCells
            invAreaCell1 = 1.0_RKIND / areaCell(iCell)
            barotropicSpeed(iCell) = 0.0_RKIND
            !columnIntegratedSpeed(iCell) = 0.0_RKIND
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)

               coeff = 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge) * invAreaCell1
               ! this is kinetic energy, in units of m^2/sec^2
               barotropicSpeed(iCell)    = barotropicSpeed(iCell) &
                    + coeff * (normalThicknessFluxSum(iEdge) / layerThicknessSumEdge(iEdge))**2

            end do
            barotropicSpeed(iCell)    = sqrt(2.0_RKIND*barotropicSpeed(iCell))

            ! columnIntegratedSpeed = sum(h*sqrt(2*ke)), where ke is kineticEnergyCell
            !   and the sum is over the full column at cell centers.
            columnIntegratedSpeed(iCell) = layerThickness(1,iCell)*sqrt( 2.0_RKIND * kineticEnergyCell(1,iCell) )
            do k=2, maxLevelCell(iCell)
               columnIntegratedSpeed(iCell) = columnIntegratedSpeed(iCell) &
                    + layerThickness(k,iCell)*sqrt( 2.0_RKIND * kineticEnergyCell(k,iCell) )
            enddo

         end do

         block => block % next
      end do

   end subroutine ocn_compute_high_frequency_output!}}}

!***********************************************************************
!
!  routine ocn_restart_high_frequency_output
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler and Phillip Wolfram
!> \date    2015/06/12
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_high_frequency_output(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_high_frequency_output!}}}

!***********************************************************************
!
!  routine ocn_finalize_high_frequency_output
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Todd Ringler and Phillip Wolfram
!> \date    2015/06/12
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_high_frequency_output(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_finalize_high_frequency_output!}}}

end module ocn_high_frequency_output

! vim: foldmethod=marker
