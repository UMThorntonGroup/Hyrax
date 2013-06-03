# This file is to test/duplicate the results from Simmons 2001 concurrent nucleation and growth paper
# to duplicate the site saturation results.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 15 #480
  ny = 15 #480
  nz = 0
  xmin = 0
  xmax = 153.6 
  ymin = 0
  ymax = 153.6 
  zmin = 0
  zmax = 0
  elem_type = QUAD4
  uniform_refine = 5
[]

[Variables]
 [./concentration]
    order = THIRD
    family = HERMITE
   [./InitialCondition]
     #type = PolySpecifiedSmoothCircleIC
     #outvalue = 0.05
     #invalue = 0.6
     #int_width = 1.0
     #radius = 2.0
     #x_positions = ''
     #y_positions = ''
     #z_positions = ''
     #n_seeds = 20
     type = RandomIC
     min = 0.0561
     max = 0.0563
   [../]
  [../]

  [./n1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
     type = PolySpecifiedSmoothCircleIC
     outvalue = 0
     invalue = 1.6
     radius = 2.0
     int_width = 1.0
     x_positions = '39  65  73 27 51 124 151 130  32 129 139 119  82 106 76  35 113 58 79 112'
     y_positions = '64 121 112 67 27  65  44 133 122 141 142  41 105 139  4 104  94 55 71 103'
     z_positions = '  0   0   0  0   0  0   0   0  0   0   0  0  0  0  0  0  0   0   0   0'
     n_seeds = 20
    [../]
  [../]
[]

[Kernels]
  [./dcdt]
    type = TimeDerivative
    variable = concentration
  [../]

  [./dn1dt]
    type = TimeDerivative
    variable = n1
  [../]

  [./CHSolid]
    type = CHBulkCoupled
    variable = concentration
    mob_name = M
    coupled_OP_var = n1
  [../]

  [./CHInterface]
    type = CHInterface
    variable = concentration
    kappa_name = kappa_c
    mob_name = M
    grad_mob_name = grad_M
  [../]

  [./ACSolid]
    type = ACBulkCoupled
    variable = n1
    coupled_CH_var = concentration
    mob_name = L
  [../]

  [./ACInterface]
    type = ACInterface
    variable = n1
    mob_name = L
    kappa_name = kappa_n
  [../]
[]

[BCs]
  [./Neumann_C]
    type = NeumannBC
    variable = concentration 
    boundary = '0 1 2 3'
  [../]

  [./Neumann_n1]
    type = NeumannBC
    variable = n1
    boundary = '0 1 2 3'
  []
[]

[Materials]
  [./constant]
    type = PFMobilityLandau
    block = 0
    mob_CH = 0.4
    mob_AC = 0.4
    kappa_CH = 0.8
    kappa_AC = 0.8
    A1 = 18.5
    A2 = -8.5
    A3 = 11.5
    A4 = 4.5
    C1 = 0.006
    C2 = 0.59
  [../]
[]

[Postprocessors]
  [./ElementIntegral_n1]
    output = file
    type = ElementIntegralVariablePostprocessor
    variable = n1
  [../]

  [./ElementIntegral_c]
    output = file
    type = ElementIntegralVariablePostprocessor
    variable = concentration
  [../]

  [./NodalMaxValue_n1]
    output = file
    type = NodalMaxValue
    variable = n1
  [../]

  [./NodalMaxValue_c]
    output = file
    type = NodalMaxValue
    variable = concentration
  [../]

  [./VolumeFraction]
    type = NodalVolumeFraction
    output = file
    bubble_volume_file = SS_size_distribution_10.csv
    Avrami_file = Avrami_10.csv
    threshold = 0.75
    variable = n1
    mesh_volume = Volume
    equil_fraction = 0.67369140625
  [../]

  [./Volume]
    type = VolumePostprocessor
    execute_on = initial
  [../]
[]

[Executioner]
  type = Transient
   scheme = 'crank-nicolson'
   start_time = 0.0
   dtmin = 0.0001
   dtmax = 0.3
   end_time = 50

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 1e-2
  [../]

  petsc_options = '-snes_mf_operator'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 101'

  [./Adaptivity]
   coarsen_fraction = 0.05
   refine_fraction = 0.75
   max_h_level = 5
#   initial_adaptivity = 5
  [../]
[]

[Output]
  file_base = testSiteSaturationGuoSS_10
  postprocessor_csv = true
  output_initial = true
  interval = 50
  exodus = true
  perf_log = true
[]
