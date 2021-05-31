#/bin/bash/
## Author	Matthew Watson	January 25, 2021


# this script will execute as a smoothed pipeline to perform the following commands: 
# - Subset a gisaid metadata file and PHO lineage report to get sample names from a specific lineage
# - Sub-sample from master multi-fasta files to get the subsets of Gisaid and PHO samples for lineage analysis


while getopts ":g:l:i:r:s:n:p:o:e" opt; do
  case $opt in
    g) gisaid_metadata="$OPTARG"
    ;;
    l) lineage_report="$OPTARG"
    ;;
    i) lineage_name="$OPTARG"
    ;;
    r) aug_reference="$OPTARG"
    ;;
    s) gisaid_sequences="$OPTARG"
    ;;
    n) subset_num="$OPTARG"
    ;;
    p) pho_sequences="$OPTARG"
    ;;
    o) output_folder="$OPTARG"
    ;;
    e) empty_reference="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

mkdir -p $output_folder

mkdir -p $output_folder/civet_data

# make subsets of the sample names based on lineage using an R script
# this assumes that the rscript is in the samew directory as this script

# specify current directory when searching for R script
script_full_path=$(dirname "$0")

Rscript $script_full_path/lineage_report_make_subsets.R --input_gisaid $gisaid_metadata --input_lineage $lineage_report --output_filenames $output_folder/pho_filenames.txt \
--output_gisaid_names $output_folder/gisaid_names.txt \
--lineage_id $lineage_name --subset_number $subset_num --metadata_output $output_folder/civet_data/gisaid_subset_metadata.csv

echo "Data subsets created using R"
	
# extract the samples fastas using the subset sample names
# this assumes that the user has downloaded and given permissions to faSomerecords and made a copy in /usr/local/bin/
/usr/local/bin/./faSomeRecords $gisaid_sequences $output_folder/gisaid_names.txt $output_folder/gisaid_subset.fa

/usr/local/bin/./faSomeRecords $pho_sequences $output_folder/pho_filenames.txt $output_folder/pho_subset.fa


# if the fasta has duplicate headers, remove the duplicate entries
awk '/^>/{f=!d[$1];d[$1]=1}f' $output_folder/pho_subset.fa > $output_folder/pho_subset_filtered.fa
awk '/^>/{f=!d[$1];d[$1]=1}f' $output_folder/gisaid_subset.fa > $output_folder/gisaid_subset_filtered.fa

rm $output_folder/pho_subset.fa
rm $output_folder/gisaid_subset.fa

gisaid_length=$(grep ">" $output_folder/gisaid_subset_filtered.fa | wc -l)
pho_length=$(grep ">" $output_folder/pho_subset_filtered.fa | wc -l)

echo "Fasta sequence subsets created- Gisaid size: $gisaid_length, PHO size: $pho_length"

# if number of gisaid sequences for the lineage is less than the requested subset number, print a warning
if [ "$subset_num" -gt "$gisaid_length" ]; then
               echo "WARNING: Subset number is greater than total Gisaid sequences with requested lineage. Using all sequences in subset"
            fi
# if number of retained sequences is greater than requested, give warning
if [ "$gisaid_length" -gt "$subset_num" ]; then
               echo "WARNING: Number of Canadian sequences greater than requested subset number. Using all Canadian sequences for context"
            fi

# if there are no PHO samples with the lineage, print a warning
if [ $pho_length -eq 0 ]; then
               echo "WARNING: No PHO samples found with requested lineage: $lineage_name"
               echo "Check that the lineage report names match the FASTA header names"
            fi

# concatenate the gisaid and PHO samples together for alignment and tree
cat $output_folder/gisaid_subset_filtered.fa $output_folder/pho_subset_filtered.fa > $output_folder/merged.fa

rm $output_folder/pho_filenames.txt
rm $output_folder/gisaid_names.txt

















