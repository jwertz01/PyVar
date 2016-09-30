# PyVar

Compares variant annotation files (e.g., from different annotators) and
outputs an HTML file containing annotator comparison plots,
and a text file containing annotations of each individual variant by each
annotator.

##Requirements

matplotlib (http://matplotlib.org/)

matplotlib_venn (https://pypi.python.org/pypi/matplotlib-venn)

numpy (http://www.numpy.org/)

##Usage

  --anv_multianno_txt_filenames: Tab-separated text files output by table_annovar.pl.
  
  --anv_var_func_filenames: Tab-separated .variant_function files output by annotate_variation.pl.
  
  --anv_exonic_var_func_filenames Tab-separated .exonic_variant_function files output by annotate_variation.pl.
  
  --vep_txt_filenames Tab-separated VEP-format files output by Ensembl Variant Effect Predictor.
  
  --filegroup_1: Filenames included in one filegroup, for which variants will be combined during analysis. Files in  groups should also be specified under ***_filenames args. Exactly two files or file groups required to run script.
    
  --filegroup_2: Filenames included in one filegroup, for which variants will be combined during analysis. Files in groups should also be specified under ***_filenames args. Exactly two files or file groups required to run script.
    
  --image_dir: Directory in which to store individual images shown in HTML output file. Default: None.
  
  --out_txt_filepath: Path at which to store output text file containing annotations of individual variants by each annotator. Default: annotators_out.txt
    
  --out_html_filepath: Path at which to store output HTML file containing plots comparing annotators overall. Default:annotators_out_plots.html


##Examples
python pyvar_compare.py --anv_var_func_filenames annovar/gatk3.3HCensAnnovar.variant_function --anv_exonic_var_func_filenames annovar/1463_15_HC_ens/gatk3.3HCensAnnovar.exonic_variant_function --vep_txt_filenames vep/1463_15_gatk3.3HC_ens.vep --filegroup_1 gatk3.3HCensAnnovar.exonic_variant_function gatk3.3HCensAnnovar.variant_function

python pyvar_compare.py --anv_multianno_txt_filenames annovar/gatk3.3HCrefseqAnnovar.hg19_multianno.txt --vep_txt_filenames vep/1463_15_gatk3.3HC_refseq.vep --out_txt_filepath out.txt --out_html_filepath out.html
