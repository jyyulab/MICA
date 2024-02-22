#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
baseCommand: merge_and_norm.py

inputs:

  out_name:
    type: string
    inputBinding:
      position: 1
      prefix: -o

  method:
    type: string
    inputBinding:
      position: 2
      prefix: -m

  mi_pairs:
    type: File[]
    inputBinding:
      position: 3
      prefix: -i

outputs:
  MI_h5:
    type: File
    outputBinding:
      glob: "*_dist.h5"
