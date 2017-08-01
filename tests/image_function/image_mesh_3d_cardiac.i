[GlobalParams]
  dicomDirectory=/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/CRT001/SE000000/
  cells_per_pixel_vector='0.25 0.25 8.0'  
  filtering_params='0.0 0.125'

  file_base =   MR0000   
  seed_index='102 102 0'
  lower_upper_threshold_values='26 55'

[]



[Mesh]
  type = ImageMeshItk
  dim = 3
  scale_to_one=false
[]

[MeshModifiers]
  [./meshModifier]
    type = ImageThresholdMesh
    lower_threshold=13
    upper_threshold=120
  [../]
[]



[Variables]
  [./u]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./image_func]
    type = ImageFunctionItk
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