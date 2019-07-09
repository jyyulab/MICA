#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
 - class: InlineJavascriptRequirement

baseCommand: 04_transform.py

inputs:
  out_name:
    type: string
    inputBinding:
      position: 1

  in_file:
    type: File
    inputBinding:
      position: 2

  trans:
    type: string
    default: "MDS"
    inputBinding:
      position: 3

  max_d:
    type: int
    default: 19
    inputBinding:
      position: 4

  dist_method:
    type: string
    inputBinding:
      position: 5

outputs:
  reduced_mi:
    type: File
    outputBinding:
      glob: "*_reduced.h5"

  preview_plots:
    type: File
    outputBinding:
      glob: "*.pdf"
