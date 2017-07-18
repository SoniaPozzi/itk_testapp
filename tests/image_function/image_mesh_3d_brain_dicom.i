[Mesh]
  type = ImageMeshItk
  dim = 3
  file_base =brain
  dicomDirectory=/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/
  cells_per_pixel_vector='1.0 1.0 1.0'   
  scale_to_one=false

[]

[Variables]
  [./u]
  [../]
[]

[Functions]
  [./image_func]
    # ImageFunction gets its file range parameters from ImageMesh,
    # when it is present.  This prevents duplicating information in
    # input files.
    type = ImageFunctionItk
    dicomDirectory=/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/
    file_base =brain
    scale_to_one=false

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
