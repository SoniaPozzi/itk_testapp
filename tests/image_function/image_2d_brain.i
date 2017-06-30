[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 391
  ny = 391
[]

[Variables]
  [./u]
  [../]
[]

[Functions]
  [./image_func]
    type = ImageFunction
    file_base = new_stack/test
    file_suffix = png
    file_range = '0' # file_range is a vector input, a single entry means "read only 1 file"
    upper_value=1.0
    lower_value=0.9
    
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