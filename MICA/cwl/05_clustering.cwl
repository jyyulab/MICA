#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
 - class: InlineJavascriptRequirement

baseCommand: 05_Clustering.py

inputs:

  input_file:
    type: File
    inputBinding:
      position: 1

  trans:
    type: string
    default: "MDS"
    inputBinding:
      position: 2

  n_cluster:
    type: int
    inputBinding:
      position: 3

  n_bootstrap:
    type: int
    inputBinding:
      position: 4

  out_name:
    type: string
    inputBinding:
      position: 5

  plot:
    type: string
    default: "tsne"
    inputBinding:
      position: 6

  umap_min_dist:
    type: float
    default: 0.1
    inputBinding:
      position: 7

  tsne_perplexity:
    type: int
    default: 30
    inputBinding:
      position: 8

  dim_plot:
    type: int
    default: 19
    inputBinding:
      position: 9

  n_thread:
    type: int
    default: 10
    inputBinding:
      position: 10

  dim_km:
    type: int[]
    default: [19]
    inputBinding:
      position: 11


outputs:
  out_fig:
    type: File
    outputBinding:
      glob: "*.png"

  out_txt:
    type: File
    outputBinding:
      glob: "*_ClusterMem.txt"

