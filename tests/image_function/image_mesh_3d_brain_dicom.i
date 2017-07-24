[GlobalParams]
  dicomDirectory=/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/CRT001/
  cells_per_pixel_vector='1.0 1.0 2.0'  
  filtering_params='5.0 0.125'

# file_base =shoulder
# seed_index='261 257 8  '
#  lower_upper_threshold_values='80 255'


 file_base =brain
  seed_index='128 54 19   177 90 19   87 76 19   85 65 18   174 67 18   67 107 18    91 48 17  61 98 17   179 59 17  114 30 16   69 67 16   178 51 16   178 48 16 198 73 15 129 26 15   75 50 15     182 48 15   127 22 14   86 30 13    97 18 12  124 15 12  133 13 11    111 39 7   167 56 7   '
  lower_upper_threshold_values='123 255'

  file_base =bran
  seed_index='61 86 4'
  lower_upper_threshold_values='23 46'

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
