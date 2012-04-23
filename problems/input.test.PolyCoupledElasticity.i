[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 25
  ny = 25
  nz = 0
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 100
  zmin = 0
  zmax = 0
  elem_type = QUAD9
  uniform_refine = 2
[]

[Variables]
  [./concentration]
    order = THIRD
    family = HERMITE
    [./InitialCondition]
      type = PolySpecifiedSmoothCircleIC
      invalue = 0.6
      outvalue = 0.1
      radius = 2.0
      int_width = 1.5
      x_positions = '10 50 80'
      y_positions = '10 50 80'
      z_positions = '0 0 0'
      n_seeds = 3
    [../]
  [../]

  [./n1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = PolySpecifiedSmoothCircleIC
      invalue = 1.6
      outvalue = 0.0
      radius = 2.0
      int_width = 1.5
      x_positions = '10'
      y_positions = '10'
      z_positions = '0'
      n_seeds = 1
    [../]
  [../]

  [./n2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = PolySpecifiedSmoothCircleIC
      invalue = 1.6
      outvalue = 0.0
      radius = 2.0
      int_width = 1.5
      x_positions = '50'
      y_positions = '50'
      z_positions = '0'
      n_seeds = 1
    [../]
  [../]

  [./n3]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = PolySpecifiedSmoothCircleIC
      invalue = 1.6
      outvalue = 0.0
      radius = 2.0
      int_width = 1.5
      x_positions = '80'
      y_positions = '80'
      z_positions = '0'
      n_seeds = 1
    [../]
  [../]

  [./dispx]
    order = FIRST
    family = LAGRANGE
  [../]

  [./dispy]
    order = FIRST
    family = LAGRANGE
  [../]
[../]

[SolidMechanics]
  [./solid]
    disp_x = dispx
    disp_y = dispy
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

  [./dn2dt]
    type = TimeDerivative
    variable = n2
  [../]

  [./dn3dt]
    type = TimeDerivative
    variable = n3
  [../]

  [./CHSolid]
    type = CHBulkPolyCoupled
    variable = concentration
    mob_name = M
    n_OP_variables = 3
    OP_variable_names = 'n1 n2 n3'
  [../]

  [./CHInterface]
    type = CHInterface
    variable = concentration
    kappa_name = kappa_c
    mob_name = M
    grad_mob_name = grad_M
  [../]

  [./ACSolid1]
    type = ACBulkPolyCoupled
    variable = n1
    mob_name = L
    coupled_CH_var = concentration
    n_OP_vars = 3
    OP_var_names = 'n1 n2 n3'
    OP_number = 1
  [../]

  [./ACInterface1]
    type = ACInterface
    variable = n1
    mob_name = L
    kappa_name = kappa_n
  [../]

  [./ACTransformElasticDF1]
    type = ACTransformElasticDF
    variable = n1
    OP_var_names = 'n1 n2 n3'
    n_OP_vars = 3
    OP_number = 1
  [../]

  [./ACSolid2]
    type = ACBulkPolyCoupled
    variable = n2
    mob_name = L
    coupled_CH_var = concentration
    n_OP_vars = 3
    OP_var_names = 'n1 n2 n3'
    OP_number = 2
  [../]

  [./ACInterface2]
    type = ACInterface
    variable = n2
    mob_name = L
    kappa_name = kappa_n
  [../]

  [./ACTransformElasticDF2]
    type = ACTransformElasticDF
    variable = n2
    OP_var_names = 'n1 n2 n3'
    n_OP_vars = 3
    OP_number = 2
  [../]

  [./ACSolid3]
    type = ACBulkPolyCoupled
    variable = n3
    mob_name = L
    coupled_CH_var = concentration
    n_OP_vars = 3
    OP_var_names = 'n1 n2 n3'
    OP_number = 3
  [../]

  [./ACInterface3]
    type = ACInterface
    variable = n3
    mob_name = L
    kappa_name = kappa_n
  [../]

  [./ACTransformElasticDF3]
    type = ACTransformElasticDF
    variable = n3
    OP_var_names = 'n1 n2 n3'
    n_OP_vars = 3
    OP_number = 3
  [../]

[]

[BCs]
active = 'Periodic'
  [./Periodic]
    [./left_right]
      primary = 0
      secondary = 2
      translation = '0 100 0'
    [../]

    [./top_bottom]
      primary = 1
      secondary = 3
      translation = '-100 0 0'
    [../]
  [../]
[]

[Materials]
  [./constant]
    type = PFMobilityLandau
    block = 0
    mob_CH = 0.4
    mob_AC = 0.4
    kappa_CH = 1.5
    kappa_AC = 1.5
    A1 = 18.5
    A2 = -8.5
    A3 = 11.5
    A4 = 4.5
    A5 = 1.0
    A6 = 1.0
    A7 = 1.0
    C1 = 0.006
    C2 = 0.59
  [../]

  [./test_material]
    type = LinearSingleCrystalPrecipitateMaterial
    block = 0
    disp_x = dispx
    disp_y = dispy
    C_matrix = '155.4 68.03 64.60 155.4 64.6 172.51 36.31 36.31 44.09'
#    C_precipitate = '100.4 40.03 42.60 100.4 42.6 150.51 26.31 26.31 32.09'
    C_precipitate = '155.4 68.03 64.60 155.4 64.6 172.51 36.31 36.31 44.09'
    e_precipitate = '0.00551 0.0564 0.0570 0.0 0.0 0.0'
    n_variants = 3
    variable_names = 'n1 n2 n3'
    all_21 = false
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'crank-nicolson'
  petsc_options = '-snes_mf_operator -ksp_monitor'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 101'

  l_max_its = 15
  nl_max_its = 50
  nl_abs_tol = 5.0e-5

  start_time = 0.0
  num_steps = 1
  dt = 0.003
[]

[Output]
  file_base = testPolyCoupledElasticityFalse
  output_initial = true
  interval = 1
  exodus = true
  perf_log = true
[]
