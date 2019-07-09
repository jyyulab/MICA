#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

baseCommand: 01_prep.py

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1

  out_name:
    type: string
    inputBinding:
      position: 2

  n_slice:
    type: int
    default: 1000
    inputBinding:
      position: 3

outputs:
  mat_pair_array:
    type: File[]
    outputBinding:
      glob: "*.h5.tmp"

