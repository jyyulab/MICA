#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
 - class: InlineJavascriptRequirement

baseCommand: transform.py

inputs:
  out_name:
    type: string
    inputBinding:
      position: 1
      prefix: -o

  in_file:
    type: File
    inputBinding:
      position: 2
      prefix: -i

  trans:
    type: string
    default: "MDS"
    inputBinding:
      position: 3
      prefix: -t

  max_d:
    type: int
    default: 19
    inputBinding:
      position: 4
      prefix: -d

  dist_method:
    type: string
    inputBinding:
      position: 5
      prefix: -m

outputs:
  reduced_mi:
    type: File
    outputBinding:
      glob: "*_reduced.h5"

  preview_plots:
    type: File[]
    outputBinding:
      glob: "*.pdf"
