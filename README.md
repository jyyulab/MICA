# MICA-Py
Mutual informtion-based Consensus Clustering for single-cell RNA-Seq data<br/>
To use this package download it and run python MICA.py from the downloaded directory.<br />
Running MICA without any arguments will show the following help message:<br /><br />

* -m OR --mode&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Running mode of the application, default: clustering.<br/>
  * Possible options: [clustering, validate, b-quality, c-quality, opt-k, overlay, opt-k-2]<br/>
* -i OR --infile&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Input File]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Path to the input file, must be symmetric simmilarity matrix (required)<br/>
  * Additional options can be added using '?' after the input file and putting each additional option into '{name=value}' format<br/>
  * Available additional options:<br/>
    * dmin&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Number >= 2], minimum eigen vector dimension, default: 16.<br/>
      * clustering, b-quality and c-quality mode: Used as the minimum number of components to do consensus clustering.<br/>
      * validate mode: Ignored.<br/>
    * dmax&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Number <= input file dimension], maximum eigen vector dimension, default: 19.<br/>
      * clustering, b-quality and c-quality mode: Used as the maximum number of components to do consensus clustering.<br/>
      * validate mode: Ignored.</br>
    * kn&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Range or list of numbers], number of clusters, default: [3:15]. (ignored in validate running mode)<br/>
    * B&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Number], number of bootstraps for k-means clustering, default: 10.<br/>
      * clustering and c-quality mode: Used as the number of k-means clustering.<br/>
      * b-quality mode: Used as the number of iterations.<br/>
      * validate and opt-k mode: Ignored.<br/>
    * delimiter&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;["," or "\t"], determines the delimiter in the txt input file.<br/>
  * Example:<br/>
    * <code>--infile /path_to_file/infile.txt?{delimiter=","}{kn="[3,4]"}{dmin=18}{dmax=19}{B=5}</code><br/>
* -d OR --decomposition&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Decomposition Method]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Name of the decomposition method for dimension reduction: [PCA, MDS, LPL, LPCA, LPCA2], default: MDS. (optional)<br/>
* -J OR --Job-Name&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Job Name]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Name that is used for all jobs related to the process. (optional)<br/>
* -T OR --true-labels&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Path to true label file], must be in the tab-separated format. (required if mode is validate, b-quality, or c-quality)<br/>
* -L OR --MICA-labels&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Path to MICA consensus clustering label file. (required if mode is validate or overlay)<br/>
* -e OR --expression\tExpression matrix path. (required if mode overlay)<br/>
  * Additional options can added using '?' after the expression file and putting each additional option into '{name=value}' format<br/>
  * Available additional options:<br/>
    * ncol&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Number> 0], determines the first column that includes a value, default: 2.<br/>
    * delimiter&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;["," or "\t"], determines the delimiter in the txt expression file.<br/>
  * Example:<br/>
    * <code>--expression /path_to_file/infile.txt?{ncol=3}{delimiter=","}</code><br/>
* -bm OR --biomarkers&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;List of biomarkers which are used to overlay gene expressions on cells. (required if mode overlay)<br/>
  * Example:<br/>
    * <code>--bimarkers "[gene1, gene2, gene3]"</code><br/>
		
## Download
<code>git clone https://github.com/jyyulab/MICA-Py</code>

## Run
<code>python MICA.py</code>

### Optimal K Analysis
<code>python MICA.py --infile your_input_file?{delimiter=","}{kn="[3:6]"} -m opt-k-2 </code>

### Clustering
<code>python MICA.py --infile your_input_file?{delimiter=","}{kn="[3,4]"}{dmin=18}{dmax=19}{B=5}</code>

### Validate
<code>python MICA.py --infile your_input_file?{delimiter=","} -m validate -T path_to_ground_truth -L path_to_clustering_labels </code>

### Component-by-Component Quality
<code>python MICA.py --infile your_input_file?{delimiter=","}{kn=[3]"}{dmin=1}{dmax=50}{B=10} -m c-quality -T path_to_ground_truth </code>

### Bootstrap Quality
<code>python MICA.py --infile your_input_file?{delimiter=","}{kn="[3]"}{B=50} -m b-quality -T path_to_ground_truth </code>

### Gene Expression Analysis
<code>python MICA.py --infile your_input_file?{delimiter=","} -m overlay -e path_to_expression_matrix?{delimiter="\t"}{ncol=1} -bm "[gene1, gene2, gene3, ...]" -L path_to_clustering_labels </code>

