#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  project_name:
    type: string
    doc: "Name of your project (projectID)"

  infile:
    type: File
    doc: "Input data file"

  slice_size:
    type: int
    default: 1000

outputs:

  mi.final:
    type: File
    outputSource: mergeAndnorm/MI_h5

steps: 
  prep:
    run: 11_prep_MIE.cwl
    in:
      out_name: project_name
      input_file: infile
      n_slice: slice_size
    out: [mat_pair_array]

  calc_pairwise_MI:
    run: 12_calc_MIpair.cwl
    in:
      h5tmp_file: prep/mat_pair_array
    scatter: h5tmp_file
    out: [MI_files]

  mergeAndnorm:
    run: 13_merge_and_norm.cwl # merge requires big memory
    in:
      mi_pairs: calc_pairwise_MI/MI_files
      out_name: project_name
    out: [MI_h5]

