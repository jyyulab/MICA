#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

baseCommand: 02_calc_scatter.py

inputs:
  h5tmp_file:
    type: File
    inputBinding:
      position: 1

  method:
    type: string
    default: "mi"
    inputBinding:
      position: 2

outputs:
  MI_files:
    type: File
    outputBinding:
      glob: $(inputs.h5tmp_file.nameroot)
