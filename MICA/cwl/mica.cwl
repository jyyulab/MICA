#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

requirements:
 - class: ScatterFeatureRequirement
 - class: SubworkflowFeatureRequirement
 - class: StepInputExpressionRequirement

inputs:
  project_name:
    type: string
    doc: "Name of your project"

  infile:
    type: File
    doc: "Input data file"

  k:
    type: int[]
    default: [4, 5, 6]
    doc: "Number of clusters, array or int"

  visualization:
    type: string
    default: "tsne"
    doc: "Two visualization options: umap and tsne"

  dim_reduction:
    type: string
    default: "MDS"
    doc: "Transformation method used for dimension reduction [MDS | PCA | LPL | LPCA] (default: MDS)"

  dist_metrics:
    type: string
    default: "mi"
    doc: "Method for distance matrix calculation [mi | euclidean | spearman | pearson] (default:mi)"

  iterations_km:
    type: int
    default: 10
    doc: "Maximum number of iterations per dimension (default: 10)"

  dims_km:
    type: int[]
    default: [19]
    doc: "Dimensions used in clustering (default: 19)"

  dims_plot:
    type: int
    default: 19
    doc: "Number of dimensions used in visualization (default: 19)"

  perplexity:
    type: int
    default: 30
    doc: "[TSNE] Visualization parameter determining how dense the clusters are (default: 30)"

  min_dist:
    type: float
    default: 0.1
    doc: "[UMAP] Visualization parameter determining how dense the cluster are (default: 0.1)"

  slice_size:
    type: int
    default: 1000
    doc: "Number of cells in each MI sub-matrix (default: 1000)"

  thread_number:
    type: int
    default: 10
    doc: "Number of pooling used for multiple kmeans iterations, usually equals to iterations_km (default: 10)"

outputs:
  mi.h5:
    type: File
    outputSource: mergeAndnorm/MI_h5

  mi.reduced:
    type: File
    outputSource: dimension_reduce/reduced_mi

  initial_plots:
    type: File[]
    outputSource: dimension_reduce/preview_plots

  clustering_res:
    type: File[]
    outputSource: clustering/out_txt

  clustering_plots:
    type:
      type: array
      items:
        type: array
        items: File
    outputSource: clustering/out_fig

steps:
  prep:
    run: prep.cwl
    in:
      out_name: project_name
      input_file: infile
      n_slice: slice_size
    out: [mat_pair_array]

  calc_pairwise_MI:
    run: calc_MI_pair.cwl
    in:
      h5tmp_file: prep/mat_pair_array
      method: dist_metrics

    scatter: h5tmp_file
    out: [MI_files]

  mergeAndnorm:
    run: merge_and_norm.cwl     # merge requires large memory
    in:
      mi_pairs: calc_pairwise_MI/MI_files
      method: dist_metrics
      out_name: project_name

    out: [MI_h5]

  dimension_reduce:
    run: dim_reduce.cwl
    in:
      out_name: project_name
      in_file: mergeAndnorm/MI_h5
      trans: dim_reduction
      max_d: dims_plot
      dist_method: dist_metrics
    out: [reduced_mi, preview_plots]

  clustering:
    run: clustering.cwl
    in:
      input_file: dimension_reduce/reduced_mi
      out_name: project_name
      n_cluster: k
      trans: dim_reduction
      n_bootstrap: iterations_km
      dim_km: dims_km
      dim_plot: dims_plot
      plot: visualization
      umap_min_dist: min_dist
      tsne_perplexity: perplexity
      n_thread: thread_number
    scatter: n_cluster
    out: [out_fig, out_txt]
