#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

requirements:
 - class: ScatterFeatureRequirement
 - class: SubworkflowFeatureRequirement
 - class: InlineJavascriptRequirement

inputs:

  project_id:
    type: string
    doc: "Name of your project (projectID)"

  in_mi:
    type: File
    doc: "Input pairwise MI file"

  k:
    type: int[]
    default: [4, 5, 6]
    doc: "Number of clusters"

  visualize:
    type: string
    default: "tsne"
    doc: "Two visualization options: uamp and tsne"

  perplexity:
    type: int
    default: 30
    doc: "[TSNE] Visualization parameter determining how dense the clusters are"

  min_dist:
    type: float
    default: 0.1

  transformation:
    type: string
    default: "MDS"
    doc: "Transformation method used for dimension reduction [MDS | PCA | LPL | LPCA] (default: MDS)"

  bootstrap:
    type: int
    default: 10
    doc: "Maximum number of iterations per dimension (default: 10)"

  dims:
    type: int[]
    default: [19]
    doc: "Dimensions used in clustering (default: 19)"

  plot_dim:
    type: int
    default: 19
    doc: "Number of dimensions used in visualization (default: 19)"

  thread_number:
    type: int
    default: 6
    doc: "Number of thread used in kmeans (default: 6)"

outputs:

  initial_plot:
    type: File
    outputSource: dimension_reduce/preview_plots

  reduced_mi:
    type: File
    outputSource: dimension_reduce/reduced_mi

  cluster_mem:
    type: File[]
    outputSource: clustering/out_txt

  final_plot:
    type: File[]
    outputSource: clustering/out_fig


steps:

  dimension_reduce:
    run: 21_dim_reduce.cwl
    in:
      out_name: project_id
      in_file: in_mi
      trans: transformation
      max_d: plot_dim
    out: [reduced_mi, preview_plots]

  clustering:
    run: 22_clustering.cwl
    in:
      input_file: dimension_reduce/reduced_mi
      out_name: project_id
      n_cluster: k
      trans: transformation
      n_bootstrap: bootstrap
      dim_km: dims
      dim_plot: plot_dim
      plot: visualize
      umap_min_dist: min_dist
      tsne_perplexity: perplexity
      n_thread: thread_number
    scatter: n_cluster
    out: [out_fig, out_txt]
