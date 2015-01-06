[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
  nz = 0
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
  active = 'u'

  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  active = 'diff'

  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  active = 'left right'

  [./left]
    type = DirichletBC
    variable = u
    boundary = 'left'
    value = 0
  [../]

  [./right]
    type = DirichletBC
    variable = u
    boundary = 'right'
    value = 1
  [../]
[]

[Postprocessors]
  [./change_variable]
    type = ChangeVariableData
    variable = u
  [../]
[]

[Executioner]
  type = Steady
  petsc_options = '-snes_mf_operator'
[]

[Output]
  output_initial = true
  interval = 1
  exodus = true
[]
