#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
baseCommand: 03_merge_and_norm.py

inputs:

  out_name:
    type: string
    inputBinding:
      position: 1

  method:
    type: string
    inputBinding:
      position: 2

  mi_pairs:
    type: File[]
    inputBinding:
      position: 3

outputs:
  MI_h5:
    type: File
    outputBinding:
      glob: "*_dist.h5"
