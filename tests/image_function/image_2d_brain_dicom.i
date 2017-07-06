[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 64
  ny = 64
  nz=16
[]

[Variables]
  [./u]
  [../]
[]

[Functions]
  [./image_func]
    type = ImageFunctionItk
    file_base = new_stack/shoulder
    file_suffix = dcm
    file_range = '001' # 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020'  # file_range is a vector input, a single entry means "read only 1 file"
   # upper_value=1.0
   # lower_value=0.9
    component=0
  [../]
[]

[ICs]
  [./u_ic]
    type = FunctionIC
    function = image_func
    variable = u
  [../]
[]

[Problem]
  type = FEProblem
  solve = false
[../]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  num_steps = 1
  dt = 0.1
[]

[Outputs]
  exodus = true
[]
