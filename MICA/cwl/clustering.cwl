#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}

baseCommand: clustering.py

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
      prefix: -i

  trans:
    type: string
    default: "MDS"
    inputBinding:
      position: 2
      prefix: -dr

  n_cluster:
    type: int
    inputBinding:
      position: 3
      prefix: -k

  n_bootstrap:
    type: int
    inputBinding:
      position: 4
      prefix: -n

  out_name:
    type: string
    inputBinding:
      position: 5
      prefix: -o

  plot:
    type: string
    default: "tsne"
    inputBinding:
      position: 6
      prefix: -m

  umap_min_dist:
    type: float
    default: 0.1
    inputBinding:
      position: 7
      prefix: -d

  tsne_perplexity:
    type: int
    default: 30
    inputBinding:
      position: 8
      prefix: -p

  dim_plot:
    type: int
    default: 19
    inputBinding:
      position: 9
      prefix: -dim

  n_thread:
    type: int
    default: 5
    inputBinding:
      position: 10
      prefix: -t

  dim_km:
    type: int[]
    default: [19]
    inputBinding:
      position: 11
      prefix: -km


outputs:
  out_fig:
    type: File[]
    outputBinding:
      glob: "*.pdf"

  out_txt:
    type: File
    outputBinding:
      glob: "clustering_*.txt"