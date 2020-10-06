#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}

baseCommand: graph_clustering.py

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

  n_neighbor:
    type: int
    inputBinding:
      position: 3
      prefix: -k

  out_name:
    type: string
    inputBinding:
      position: 4
      prefix: -o

  plot:
    type: string
    default: "tsne"
    inputBinding:
      position: 5
      prefix: -m

  umap_min_dist:
    type: float
    default: 0.1
    inputBinding:
      position: 6
      prefix: -d

  tsne_perplexity:
    type: int
    default: 30
    inputBinding:
      position: 7
      prefix: -p

  dim_plot:
    type: int
    default: 19
    inputBinding:
      position: 8
      prefix: -dim

  dims:
    type: int
    default: 19
    inputBinding:
      position: 9
      prefix: -ds


outputs:
  out_fig:
    type: File
    outputBinding:
      glob: "*.png"

  out_txt:
    type: File
    outputBinding:
      glob: "*_ClusterMem.txt"

