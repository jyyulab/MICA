#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

baseCommand: prep.py

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
      prefix: -i

  out_name:
    type: string
    inputBinding:
      position: 2
      prefix: -o

  n_slice:
    type: int
    default: 1000
    inputBinding:
      position: 3
      prefix: -s

outputs:
  mat_pair_array:
    type: File[]
    outputBinding:
      glob: "*.h5.tmp"
