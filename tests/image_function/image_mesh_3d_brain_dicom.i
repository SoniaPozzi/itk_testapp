[GlobalParams]
  dicomDirectory=/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/    # CRT001/SE000000
  cells_per_pixel_vector='0.25 0.25 8.0'  
  filtering_params='0.0 0.125'

# file_base =shoulder
# seed_index='261 257 8  '
#  lower_upper_threshold_values='80 255'


 file_base =brain
  seed_index='128 54 19   177 90 19   87 76 19   85 65 18   174 67 18   67 107 18    91 48 17  61 98 17   179 59 17  114 30 16   69 67 16   178 51 16   178 48 16 198 73 15 129 26 15   75 50 15     182 48 15   127 22 14  137 19 13 127 17 13 86 30 13 202 69 13 93 31 14 181 42 14 202 76 14 190 52 14  54 88 15 88 40 15 167 33 15 148 28 15  181 33 12 112 14 12 167 21 12   97 18 12  124 15 12  133 13 11    111 39 7   167 56 7  50 105 15 77 48 15 189 56 15 58 92 16 54 113 16 80 52 16 142 31 16 92 41 16 201 100 16  115 30 16 123 29 16 165 39 16 115 30 17 123 29 17 165 39 17 129 36 17 194 90 17 106 40 17 113 38 17 132 44 18 190 103 18 133 45 18 '
  lower_upper_threshold_values='122.6 255'

#  file_base =   MR0000    # MR0000
#  seed_index='102 102 0'
#  lower_upper_threshold_values='26 55'

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
