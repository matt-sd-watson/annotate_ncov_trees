# annotate_ncov_trees

Scripts and workflow elements to automate the annotation of newick trees for COVID-19

## Installation

The package requires conda to be installed in the current environment. 

```
git clone https://github.com/matt-sd-watson/annotate_ncov_trees.git
cd annotate_ncov_trees
conda env create -f workflow/environment.yml
conda activate ncov_phylo

```

## Running the workflow

To execute the pipeline: 

```
conda activate ncov_phylo
snakemake --cores all -s annotate_ncov_trees/workflow/Snakefile \
          --configfile annotate_ncov_trees/workflow/config.yaml

```

The arguments specified in the config.yaml should be as follows: 

```
# directory to hold all lineage folders and files
output_dir: replace this path with the desired output directory

# lineages to evaluate
lineages: A python-=like list of lineages to investigate. Example: [B.1.1.7, B.1.351, P.1]

# number of background Gisaid sequences to include in the tree
subset_number: Number of background Gisaid sequences to include for each lineage. Example: 300

# path to Gisaid files
gisaid_dir: the directory where master Gisaid seqwuences are held.

# path to PHO master files
master_dir: the directory where master Gisaid seqwuences are held.

# path to reference genome
reference_genome: The path to the ncov reference in gb format.

#location for the QC list required to filter the final annotated tree
qc_list: path_to_qc_list.csv

```
