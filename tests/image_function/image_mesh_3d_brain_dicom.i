[GlobalParams]
  dicomDirectory=/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/
  cells_per_pixel_vector='0.5 0.50 1.0'  
  filtering_params='5.0 0.125'

#  file_base =brain
#  seed_index='60 116 8'
#  lower_upper_threshold_values='90 112'


  file_base =brain
  seed_index='60 116 8'
  lower_upper_threshold_values='90 112'


[]



[Mesh]
  type = ImageMeshItk
  dim = 3
  scale_to_one=false
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
    dicomDirectory=/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/
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
