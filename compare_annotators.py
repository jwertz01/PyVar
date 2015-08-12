"""Compares variant annotation files (e.g., from different annotators) and
outputs an HTML file containing annotator comparison plots,
and a text file containing annotations of each individual variant by each
annotator.
"""

import sys
import os
import operator
import re
import math
import argparse
import itertools
import pylab as plt
import numpy as np
import matplotlib
import matplotlib_venn
import shutil
import base64


class Variant(object):
    def __init__(
        self, string=None, chrom=None, start_pos=None, end_pos=None,
        consequence=None, normalized_consequence=None, transcript=None
    ):
        self.string = string
        self.chrom = chrom
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.consequence = consequence
        self.normalized_consequence = normalized_consequence
        self.transcript = transcript

    def __str__(self):
        return str(self.__dict__)


class AnnotationFile(object):
    """General class: Output file from annotator."""
    def __init__(
        self, annotator=None, delimiter=None, file_obj=None,
        header_list = None, curr_variant=None, next_variant=None,
        has_header=None
    ):
        self.annotator = annotator
        self.delimiter = delimiter
        self.file_obj = file_obj
        self.header_list = header_list
        self.eof = False
        self.curr_variant = curr_variant
        self.next_variant = next_variant
        self.other_consequence = 'other'
        self.has_header = has_header

    def __str__(self):
        return str(self.__dict__)

    def process_next_variant(self):
        """Extract variant chrom, pos, consequence, etc. from file."""
        pass


class AnvMultiannoTxt(AnnotationFile):
    """Annovar output multianno text file."""
    def __init__(self, file_obj=None,):
        super(AnvMultiannoTxt, self).__init__(
            annotator='annovar', delimiter='\t', file_obj=file_obj,
            has_header=True
        )

    def process_next_variant(self):
        v = self.next_variant
        h = self.header_list
        var_list = v.string.strip().split(self.delimiter)
        v.chrom = var_list[h.index('Chr')]
        v.start_pos = int(var_list[h.index('Start')])
        v.end_pos = int(var_list[h.index('End')])
        for i, categ in enumerate(h):
            if re.findall('^Func.*gene', categ):
                func_index = i
            if re.findall('ExonicFunc.*gene', categ):
                exonic_func_index = i
            if re.findall('AAChange.*gene', categ):
                aa_change_index = i
        func = var_list[func_index].strip('"')
        exonic_func = var_list[exonic_func_index].strip('"')
        dict_key = exonic_func if (func == 'exonic') else func
        v.consequence = [dict_key]
        v.normalized_consequence = get_norm_consequence(
            self.annotator, v.consequence, self.other_consequence
        )
        if var_list[aa_change_index] not in ['.', 'UNKNOWN']:
            v.transcript = [
                z.split(':')[1] for z in var_list[aa_change_index].split(',')
            ]


class AnvVariantFunction(AnnotationFile):
    """Annovar output variant_function text file."""
    def __init__(self, file_obj=None):
        super(AnvVariantFunction, self).__init__(
            annotator='annovar', delimiter='\t', file_obj=file_obj,
            has_header=False
        )

    def process_next_variant(self):
        v = self.next_variant
        var_list = v.string.strip().split(self.delimiter)
        v.chrom = var_list[2]
        v.start_pos = int(var_list[3])
        v.end_pos = int(var_list[4])
        v.consequence = [var_list[0]]
        v.normalized_consequence = get_norm_consequence(
            self.annotator, v.consequence, self.other_consequence
        )
        v.transcript = (
            # remove everything in parentheses
            re.sub(r'\([^\)]*\)', '', var_list[1]).split(',')
        )


class AnvExonicVariantFunction(AnnotationFile):
    """Annovar output exonic_variant_function text file."""
    def __init__(self, file_obj=None):
        super(AnvExonicVariantFunction, self).__init__(
            annotator='annovar', delimiter='\t', file_obj=file_obj,
            has_header=False
        )

    def process_next_variant(self):
        v = self.next_variant
        var_list = v.string.strip().split(self.delimiter)
        v.chrom = var_list[3]
        v.start_pos = int(var_list[4])
        v.end_pos = int(var_list[5])
        v.consequence = [var_list[1]]
        v.normalized_consequence = get_norm_consequence(
            self.annotator, v.consequence, self.other_consequence
        )
        if var_list[2] not in ['.', 'UNKNOWN']:
            v.transcript = [
                z.split(':')[1] for z in var_list[2].strip(',').split(',')
            ]

class VEPTxt(AnnotationFile):
    """VEP output text file."""
    def __init__(self, file_obj=None):
        super(VEPTxt, self).__init__(
            annotator='vep', delimiter='\t', file_obj=file_obj,
            has_header=True
        )

    def process_next_variant(self):
        v = self.next_variant
        h = self.header_list
        var_list = v.string.strip().split(self.delimiter)
        loc = re.split('[:-]', var_list[h.index('Location')])
        v.chrom = 'chr%s' % loc[0]
        v.start_pos = int(loc[1])
        v.end_pos = int(loc[2]) if (len(loc) > 2) else v.start_pos
        v.consequence = var_list[h.index('Consequence')].split(',')
        v.normalized_consequence = get_norm_consequence(
            self.annotator, v.consequence, self.other_consequence
        )
        if (
            'Feature_type' in h and 'Feature' in h and
            var_list[h.index('Feature_type')] == 'Transcript'
        ):
            v.transcript = [var_list[h.index('Feature')]]


consequence_names = {
    'annovar': {
        'frameshift deletion' : 'frameshift',
        'frameshift insertion' : 'frameshift',
        'frameshift_block_subsitution' : 'frameshift',
        'stopgain' : 'stopgain',
        'stoploss' : 'stoploss',
        'splicing' : 'splicing',
        'exonic;splicing' : 'splicing',
        'nonframeshift deletion' : 'inframe',
        'nonframeshift insertion' : 'inframe',
        'nonsynonymous SNV' : 'nonsynonymous',
        'nonframeshift_block_substitution' : 'nonsynonymous',
        'synonymous SNV' : 'synonymous',
        'UTR3' : 'UTR3',
        'UTR5' : 'UTR5',
        'UTR5;UTR3' : 'UTR5',
        'upstream' : 'upstream',
        'upstream;downstream' : 'upstream',
        'downstream' : 'downstream',
        'intronic' : 'intron',
        'intergenic' : 'intergenic',
        'ncRNA_exonic' : 'nc_exon',
        'ncRNA_intronic' : 'nc_intron',
        'ncRNA_splicing' : 'nc_splicing',
        'ncRNA_exonic;splicing' : 'nc_splicing',
        'ncRNA_UTR3' : 'nc_UTR3',
        'ncRNA_UTR5' : 'nc_UTR5'
    },
    'vep': {
        'frameshift_variant' : 'frameshift',
        'stop_gained' : 'stopgain',
        'stop_lost' : 'stoploss',
        'splice_donor_variant' : 'splicing',
        'splice_acceptor_variant' : 'splicing',
        'splice_region_variant' : 'splicing',
        'inframe_insertion' : 'inframe',
        'inframe_deletion' : 'inframe',
        'initiator_codon_variant' : 'nonsynonymous',
        'missense_variant' : 'nonsynonymous',
        'incomplete_terminal_codon_variant' : 'nonsynonymous',
        'synonymous_variant' : 'synonymous',
        'stop_retained_variant' : 'synonymous',
        '3_prime_UTR_variant' : 'UTR3',
        '5_prime_UTR_variant' : 'UTR5',
        'upstream_gene_variant' : 'upstream',
        'downstream_gene_variant' : 'downstream',
        'intron_variant' : 'intron',
        'intergenic_variant' : 'intergenic',
        'non_coding_transcript_exon_variant': 'nc_exon',
        'non_coding_transcript_variant': 'nc_intron',
        # not included: 'nc_splicing', 'nc_UTR3', 'nc_UTR5'
    }
}

severity_ranking = [
    'frameshift', 'stopgain', 'stoploss', 'splicing', 'inframe',
    'nonsynonymous', 'synonymous', 'UTR5', 'UTR3', 'nc_exon', 'nc_splicing',
    'intron', 'nc_intron',  'upstream', 'downstream', 'nc_UTR5', 'nc_UTR3',
    'intergenic', 'None'
]


def main():
    args, fg1_paths, fg2_paths = parse_arguments(argparse.ArgumentParser())
    files = open_files(args)
    file_groups = []
    if len(fg1_paths) > 1:
        file_groups.append(tuple(fg1_paths))
    if len(fg2_paths) > 1:
        file_groups.append(tuple(fg2_paths))

    (
        out_file, empty_files, file_combo_counts, grouped_files,
        all_norm_consq
    ) = init_files(files, args.out_txt_filepath, file_groups)

    process_files(files, out_file, file_combo_counts)
    close_files(files, out_file)

    image_dir = args.image_dir
    if not os.path.exists(image_dir):
        os.mkdir(image_dir)
    else:
        for f_name in os.listdir(image_dir):
            os.remove(os.path.join(image_dir, f_name))

    # output
    filegroup_count = len(grouped_files)
    set_labels = tuple(sorted(grouped_files))
    set_labels_wrapped = tuple(
        [wrap_str(';\n'.join(x.split(';')), 30) for x in set_labels]
    )

    subsets_transcript = find_subsets(
        filegroup_count, file_combo_counts['transcripts'], set_labels
    )
    subsets_consq = find_subsets(
        filegroup_count, file_combo_counts['norm_consq'], set_labels
    )

    piechart_colors = [
        'skyblue', 'lightcoral', 'beige', 'violet', 'cornflowerblue',
        'goldenrod', 'lightgreen', 'sandybrown', 'lightpink', 'teal',
        'thistle', 'darkseagreen', 'lightgray', 'bisque', 'aquamarine',
        'coral', 'greenyellow', 'hotpink', 'rosybrown', 'mediumpurple'
    ]

    #print file_combo_counts
    create_pie_chart(
        'Common Transcripts',
        'Percentage of variants that have a common transcript\n(a transcript '
        'that all annotators used for annotation).',
        os.path.join(image_dir, '01_transcr_pie.png'),
        ['Common transcript', 'No common transcript'],
        [subsets_transcript[-1], sum(subsets_transcript[:-1])],piechart_colors
    )
    create_venn(
        'Common Transcripts',
        'Middle: Variants annotated with >=1 transcript that both annotators '
        'used.\nLeft and right: Variants annotated only with transcripts '
        'other annotator did not use. ',
        os.path.join(image_dir, '02_transcr_venn.png'),
        set_labels_wrapped, subsets_transcript, filegroup_count
    )
    create_venn(
        'Normalized Consequences',
        'Middle: Variants annotated with >=1 consequence that both '
        'annotators identified.\nLeft and right: Variants annotated only '
        'with consequences other annotator did not identify.\n'
        'Includes only consequences using transcripts common to both '
        'annotators.',
        os.path.join(image_dir, '03_norm_consq_venn.png'),
        set_labels_wrapped, subsets_consq, filegroup_count
    )

    # make color map
    consq_color_map = {}
    all_norm_consq_list = list(all_norm_consq)
    for i in xrange(len(all_norm_consq)):
        consq_color_map[all_norm_consq_list[i]] = piechart_colors[i]
    consq_color_map['other'] = piechart_colors[i + 1]

    # make normalized consequence pie charts
    for i, x in enumerate(grouped_files):
        subsets_consq_percateg = []
        for y in sorted (all_norm_consq):
            subsets_consq_percateg.append(
                file_combo_counts['norm_consq_per_type'][(x,)][y]
            )
        create_pie_chart(
            'Normalized Consequences Per File Group',
            wrap_str(';\n'.join(x.split(';')), 80),
            os.path.join(image_dir, '0%d_consq_pie' % (i + 4)),
            sorted(all_norm_consq), subsets_consq_percateg,
            piechart_colors, consq_color_map
        )
        curr_image = i + 4

    # make heatmaps
    heatmap_data = []
    for x in severity_ranking:
        row = []
        for y in severity_ranking:
            if (x, y) in file_combo_counts['consq_mismatch_pairs']:
                row.append(
                    file_combo_counts['consq_mismatch_pairs'][(x, y)]
                )
            else:
                row.append(0)
        heatmap_data.append(row)

    row_norm_heatmap_data = normalize_2d_list(
        np.array(heatmap_data).transpose()
    )
    col_norm_heatmap_data = normalize_2d_list(heatmap_data)

    wrapped_severity = [
        wrap_str(';\n'.join(z.split(';')), 50) for z in
        sorted(grouped_files)
    ]
    heatmap_description_template = (
        'Counts of most severe normalized consequence identified '
        'by each annotator.\n(Colors correspond to logarithm of raw counts.) '
        'Heatmap is normalized per '
    )
    create_heatmap(
        row_norm_heatmap_data, 'Row-Normalized Consequence Heatmap',
        heatmap_description_template + 'row.', os.path.join(
            image_dir, '0%d_heatmap.png' % (curr_image + 1)
        ), severity_ranking, *wrapped_severity
    )
    create_heatmap(
        col_norm_heatmap_data, 'Column-Normalized Consequence Heatmap',
        heatmap_description_template + 'column.', os.path.join(
            image_dir, '0%d_heatmap.png' % (curr_image + 2)
        ), severity_ranking, *wrapped_severity
    )
    create_html(
        args.out_html_filepath, image_dir,
        re.sub(';', ';\n', '\n\nvs.\n\n'.join(grouped_files))
    )
    if image_dir == '.image_files_temp':
        shutil.rmtree(image_dir)


def get_norm_consequence(annotator, consequence, other_consequence):
    return [
        (
            consequence_names[annotator][z] if
            (z in consequence_names[annotator]) else other_consequence
        ) for z in consequence
    ]


def parse_arguments(parser):
    parser.add_argument(
        '--anv_multianno_txt_filenames', nargs='*',
        help='Tab-separated text files output by table_annovar.pl.'
    )
    parser.add_argument(
        '--anv_var_func_filenames', nargs='*',
        help=(
            'Tab-separated .variant_function files output by '
            'annotate_variation.pl.'
        )
    )
    parser.add_argument(
        '--anv_exonic_var_func_filenames', nargs='*',
        help=(
            'Tab-separated .exonic_variant_function files output by '
            'annotate_variation.pl.'
        )
    )
    parser.add_argument(
        '--vep_txt_filenames', nargs='*',
        help=(
            'Tab-separated VEP-format files output by '
            'Ensembl Variant Effect Predictor.'
        )
    )
    parser.add_argument(
        '--filegroup_1', nargs='*', help=(
            'Filenames included in one filegroup, for which variants will '
            'be combined during analysis. Files in groups should '
            'also be specified under ***_filenames args. '
            'Exactly two files or file groups required to run script.'
        )
    )
    parser.add_argument(
        '--filegroup_2', nargs='*', help=(
            'Filenames included in one filegroup, for which variants will '
            'be combined during analysis. Files in groups should '
            'also be specified under ***_filenames args. '
            'Exactly two files or file groups required to run script.'
        )
    )
    parser.add_argument(
        '--image_dir', help=(
            'Directory in which to store individual images shown in HTML '
            'output file. Default: None.'
        ), default='.image_files_temp'
    )
    parser.add_argument(
        '--out_txt_filepath', help=(
            'Path at which to store output text file containing annotations '
            'of individual variants by each annotator. '
            'Default: annotators_out.txt'
        ), default='annotators_out.txt'
    )
    parser.add_argument(
        '--out_html_filepath', help=(
            'Path at which to store output HTML file containing plots '
            'comparing annotators overall. '
            'Default: annotators_out_plots.html'
        ), default='annotators_out_plots.html'
    )

    # validate args
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    if not args.out_txt_filepath.endswith('.txt'):
        parser.print_help()
        print 'Error: Argument out_txt_filepath requires extension .txt.'
        sys.exit(1)
    if not args.out_html_filepath.endswith('.html'):
        parser.print_help()
        print 'Error: Argument out_html_filepath requires extension .html.'
        sys.exit(1)

    args_dict_old = args.__dict__
    args_dict = {}
    for x in args_dict_old:
        if args_dict_old[x]:
            args_dict[x] = args_dict_old[x]
    all_filepaths = [
        x for y in [
            args_dict[z] for z in args_dict if ('filenames' in z)
        ] for x in y
    ]
    if not all_filepaths:
        parser.print_help()
        print 'Error: All ***_filenames args are empty.'
        sys.exit(1)

    files_in_groups = []
    count_filegroups_used = 0
    if args.filegroup_1:
        check_file_grouping(
            args.filegroup_1, all_filepaths, 'filegroup_1', parser
        )
        files_in_groups += args.filegroup_1
        count_filegroups_used += 1
    if args.filegroup_2:
        check_file_grouping(
            args.filegroup_2, all_filepaths, 'filegroup_2', parser
        )
        files_in_groups += args.filegroup_2
        count_filegroups_used += 1

    files_not_in_groups = [
        z for z in all_filepaths if (
            (os.path.basename(z) not in files_in_groups) and
            (z not in files_in_groups)
        )
    ]
    filegroup_count = count_filegroups_used + len(files_not_in_groups)

    fg1_paths = []  # (path, filetype)
    fg2_paths = []
    other_paths = []
    if args.filegroup_1:
        for x in args.filegroup_1:
            for y in all_filepaths:
                filetype = [z for z in args_dict if y in args_dict[z]][0]
                if x == y or os.path.basename(y) == x:
                    fg1_paths.append((y, filetype))
                elif os.path.basename(x) == y:
                    fg1_paths.append((x, filetype))
    if args.filegroup_2:
        for x in args.filegroup_2:
            for y in all_filepaths:
                if x == y or os.path.basename(y) == x:
                    fg2_paths.append(y)
                elif os.path.basename(x) == y:
                    fg2_paths.append(x)
    for x in files_not_in_groups:
        filetype = [z for z in args_dict if x in args_dict[z]][0]
        if not fg1_paths:
            fg1_paths.append((x, filetype))
        elif not fg2_paths:
            fg2_paths.append((x, filetype))
        else:
            other_paths.append((x, filetype))

    detected_fgroups_str = '\nDetected filegroups\n-------------------\n'
    if fg1_paths:
        detected_fgroups_str += '-Filegroup 1: %s\n' % '; '.join([
            '%s (%s)' % (fname, ftype[:ftype.index('_filenames')])
            for (fname, ftype) in fg1_paths
        ])
    if fg2_paths:
        detected_fgroups_str += '-Filegroup 2: %s\n' % '; '.join([
            '%s (%s)' % (fname, ftype[:ftype.index('_filenames')])
            for (fname, ftype) in fg2_paths
        ])
    if other_paths:
        detected_fgroups_str += '-Other files: %s\n' % '; '.join([
            '%s (%s)' % (fname, ftype[:ftype.index('_filenames')])
            for (fname, ftype) in other_paths
        ])

    if filegroup_count != 2:
        parser.print_help()
        print (
            '\nError: Exactly two files or file groups required. '
            'All files must be specified under ***_filenames args. '
            'Filegroups may optionally be specified, in which case variants '
            'from multiple files will be combined during analysis.'
        )
        print detected_fgroups_str
        sys.exit(1)

    print detected_fgroups_str
    return (
        args,
        [os.path.basename(a) for (a, b) in fg1_paths],
        [os.path.basename(c) for (c, d) in fg2_paths]
    )


def check_file_grouping(file_grouping, all_filepaths, arg_label, parser):
    for x in file_grouping:
        if (
            (x not in all_filepaths) and
            (x not in [os.path.basename(z) for z in all_filepaths]) and
            (os.path.basename(x) not in all_filepaths)
        ):
            parser.print_help()
            print (
                'Error: File present in %s that is not present '
                'in ***_filenames args.' % arg_label
            )
            sys.exit(1)


def open_files(args):
    files = []
    open_filetype(args.anv_multianno_txt_filenames, AnvMultiannoTxt, files)
    open_filetype(args.anv_var_func_filenames, AnvVariantFunction, files)
    open_filetype(
        args.anv_exonic_var_func_filenames, AnvExonicVariantFunction, files
    )
    open_filetype(args.vep_txt_filenames, VEPTxt, files)
    return files


def open_filetype(filenames, filetype, files):
    """For each file in filenames, open file, create file obj of given
    filetype, and append to files list.
    """
    if filenames:
        for f_name in filenames:
            files.append(filetype(file_obj=open(f_name)))


def init_files(files, out_filename, file_groups):
    empty_files = []
    for i, f in enumerate(files):
        if f.has_header:
            for row in f.file_obj:
                if not row.strip().startswith('##'):
                    break
            f.header_list = row.strip().split(f.delimiter)
        try:
            f.next_variant = Variant(string=f.file_obj.next())
        except StopIteration:
            f.next_variant = Variant(string='')
        if not f.next_variant.string.strip():
            empty_files.append(f)
            del files[i]
            continue
        f.process_next_variant()

    categs = [
        'transcripts', 'norm_consq', 'norm_consq_per_type',
        'consq_mismatch_pairs'
    ]
    # set up agreement counts for each combination of file groups
    file_combo_counts = {}
    files_in_groups = set([])
    # list of file groups (semicolon-delimited strings of filenames)
    grouped_files = []
    for x in file_groups:
        for y in x:
            files_in_groups.add(y)
        grouped_files.append(';'.join(x))
    for x in files:
        f_name = os.path.basename(x.file_obj.name)
        if f_name not in files_in_groups:
            grouped_files.append(f_name)

    filename_combinations = []
    for i in range(1, len(files) + 1):
        filename_combinations.append([
            list(y) for y in itertools.combinations(grouped_files, i)
        ])
    filename_combinations = [x for y in filename_combinations for x in y]
    for x in categs:
        file_combo_counts[x] = {}
        for y in filename_combinations:
            if x == 'norm_consq_per_type':
                y = tuple(sorted(y))
                file_combo_counts[x][y] = {}
                all_norm_consq = set([
                    a for b in
                    [consequence_names[w].values() for w in consequence_names]
                    for a in b
                ])
                for z in all_norm_consq:
                    file_combo_counts[x][y][z] = 0
            elif x == 'consq_mismatch_pairs':
                if len(y) == 1:
                    file_combo_counts[x] = {}
            else:
                file_combo_counts[x][tuple(sorted(y))] = 0

    out_f = open(out_filename, 'w')
    out_f.write(
        'Annotator\tFilename\tChrom\tStart_pos\tEnd_pos\tConsequence\t'
        'Normalized_consequence\tTranscript\tCommon_transcript_exists\t'
        'Consequences_match'
    )
    return (
        out_f, empty_files, file_combo_counts, grouped_files, all_norm_consq
    )


def process_files(files, out_file, file_combo_counts):
    min_chrom = -1
    min_pos = -1
    variants_curr_pos = []

    while not(all([z.eof for z in files])):
        old_min_chrom = min_chrom
        min_chrom = min_chromosome([
            z.next_variant.chrom for z in files if not z.eof
        ])
        if old_min_chrom != min_chrom:
            # if old_min_chrom > 0:
            #     sys.exit(1)
            print 'Processing %s...' % min_chrom
        old_min_pos = min_pos
        min_pos = min([
            z.next_variant.start_pos for z in files if (
                (
                    parse_chromosome(z.next_variant.chrom) ==
                    parse_chromosome(min_chrom)
                ) and not z.eof
            )
        ])
        if old_min_pos != min_pos:
            process_curr_pos(variants_curr_pos, out_file, file_combo_counts)
            variants_curr_pos = []

        for f in files:
            f.curr_variant = f.next_variant
            if (
                (
                    parse_chromosome(f.curr_variant.chrom) ==
                    parse_chromosome(min_chrom)
                )
                and f.curr_variant.start_pos == min_pos and not f.eof
            ):
                # output
                v = f.curr_variant
                variants_curr_pos.append((v, f))
                try:
                    f.next_variant = Variant(string=f.file_obj.next())
                except StopIteration:
                    f.eof = True
                    continue
                f.process_next_variant()
    # process last variant
    if variants_curr_pos:
        process_curr_pos(
            variants_curr_pos, out_file, file_combo_counts
        )


def close_files(files, out_file):
    for f in files:
        f.file_obj.close()
    out_file.close()


def min_chromosome(chroms):
    """Args: Chromosome names ['chrN', 'chrP', ...]. Return name of
    chromosome which would appear first in ordered file.
    """
    alpha_chrom_order = ['X', 'Y', 'U', 'M']
    chroms_list = list(chroms)
    for i in xrange(len(chroms_list)):
        chroms_list[i] = parse_chromosome(chroms_list[i])
    if any([isinstance(z, int) for z in chroms_list]):
        return 'chr%s' % min(
            [y for y in chroms_list if isinstance(y, int)]
        )
    else:
        first_letter = alpha_chrom_order[
            min([alpha_chrom_order.index(z[0]) for z in chroms_list])
        ]
        return 'chr%s' % min(
            [y for y in chroms_list if y[0] == first_letter]
        )


def parse_chromosome(chrom):
    """Parse chromosome name."""
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    try:
        chrom = int(chrom)
    except ValueError:
        pass
    return chrom


def process_curr_pos(variants_curr_pos, out_file, file_combo_counts):
    """Process variant at current position in each annotation file.
    Arg: list of (variant obj, file obj) at a single position,
    over different files.
    """
    common_transcr_exists = update_file_combo_counts(
        file_combo_counts, variants_curr_pos, 'transcripts', 'transcript'
    )
    matching_norm_consq = update_file_combo_counts(
        file_combo_counts, variants_curr_pos, 'norm_consq',
        'normalized_consequence'
    )
    for v, f in variants_curr_pos:
        out_file.write(
            '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n' % (
                f.annotator, os.path.basename(f.file_obj.name),
                v.chrom, v.start_pos, v.end_pos, ', '.join(v.consequence),
                ', '.join(v.normalized_consequence),
                ', '.join(v.transcript) if v.transcript else 'unknown',
                common_transcr_exists, matching_norm_consq
            )
        )
    out_file.write('\n')


def update_file_combo_counts(
    file_combo_counts, variants_curr_pos, file_combo_counts_categ,
    variant_attribute
):
    """Create dict {file : set of transcripts/annotations/etc.}
    For each possible file combination, find whether all keys (files)
    have a value in common.
    If so, increment file_combo_counts[file_combo_counts_categ][files].
    Return whether common value for file_combo_counts_categ exists.
    """
    file_categ = {}  # {file: set of file_combo_counts_categ values}
    # {file: list of (transcript, normalized consequence) tuples}
    file_trans_consq = {}
    for grp in file_combo_counts['norm_consq']:
        if len(grp) == 1:
            file_trans_consq[grp[0]] = []
    for v, f in variants_curr_pos:
        filename = os.path.basename(f.file_obj.name)
        if filename not in file_categ:
            file_categ[filename] = set([])
        if v.__dict__[variant_attribute]:
            for y in v.__dict__[variant_attribute]:
                if y != 'other':
                    file_categ[filename].add(y)
        if v.transcript and v.normalized_consequence:
            for x in file_trans_consq:
                if filename in x.split(';'):
                    file_trans_consq[x].append(
                        (v.transcript, v.normalized_consequence)
                    )
    if variants_curr_pos:
        chrom = variants_curr_pos[0][0].chrom
        pos = variants_curr_pos[0][0].start_pos

    curr_values = set([]) # set of dict values for all keys
    for x in file_categ.values():
        curr_values = curr_values.union(x)

    filegroup_combos = file_combo_counts[file_combo_counts_categ].keys()
    filegroup_combos.sort(key=len, reverse=True)

    # update heat map counts (consq_mismatch_pairs)
    if file_combo_counts_categ == 'norm_consq':
        most_sev_fgroup_consq = []
        for f_group in sorted(
            [z[0] for z in filegroup_combos if len(z) == 1]
        ):
            filegroup_consq = set([])
            for w in f_group.split(';'):
                if w in file_categ:
                    filegroup_consq = filegroup_consq.union(file_categ[w])
            if filegroup_consq:
                most_sev_fgroup_consq.append(severity_ranking[min(
                    [severity_ranking.index(z) for z in filegroup_consq]
                )])
            else:
                most_sev_fgroup_consq.append('None')
        dict_key = tuple(most_sev_fgroup_consq)
        if dict_key not in file_combo_counts['consq_mismatch_pairs']:
            file_combo_counts['consq_mismatch_pairs'][dict_key] = 0
        file_combo_counts['consq_mismatch_pairs'][dict_key] += 1

    incremented = False
    for i, filegroup_combo in enumerate(filegroup_combos):
        # only increment count for filegroup combinations with the most files
        if incremented and (
            len(filegroup_combo) < len(filegroup_combos[i - 1])
        ):
            all_agree = len(filegroup_combos[i - 1]) == max(
                [len(j) for j in filegroup_combos]
            )
            return all_agree

        for x in curr_values:
            curr_value_present_in_all_filegroups = True
            for filegroup in filegroup_combo:
                if not any([
                    (z in file_categ and x in file_categ[z]) for z in
                    filegroup.split(';')
                ]):
                    curr_value_present_in_all_filegroups = False
                    break

            if curr_value_present_in_all_filegroups:
                if file_combo_counts_categ == 'norm_consq':
                    # transcripts corresponding to consequence x
                    file_transcr_sets = []
                    for y in file_trans_consq:
                        if y in filegroup_combo:
                            curr_transcr = [
                                a for b in [
                                    [u for u in w[0] if x in w[1]]
                                    for w in file_trans_consq[y]
                                ] for a in b
                            ]
                            file_transcr_sets.append(set(curr_transcr))
                    # transcript exists that is present in all files
                    # in filegroup_combo
                    if set.intersection(*file_transcr_sets):
                        file_combo_counts['norm_consq_per_type'][
                            filegroup_combo
                        ][x] += 1
                        incremented = True
                else:
                    incremented = True
                if incremented:
                    file_combo_counts[file_combo_counts_categ][
                        filegroup_combo
                    ] += 1
                    break
    return False


def wrap_str(str_, wrap_len):
    if wrap_len < 1:
        return str_
    str_wrapped = ''
    i = 0  # number of chars since \n
    for j, ch in enumerate(str_):
        i += 1
        str_wrapped += ch
        if ch == '\n':
            i = 0
        elif (
                (i == wrap_len) and (j < len(str_) - 1) and
                (str_[j + 1] != '\n')
        ):
            str_wrapped += '\n'
            i = 0
    return str_wrapped


def find_subsets(filegroup_count, file_combo_counts, set_labels):
    if filegroup_count == 2:
        subsets = (
            file_combo_counts[(set_labels[0],)],
            file_combo_counts[(set_labels[1],)],
            file_combo_counts[set_labels]
        )
    elif filegroup_count == 3:
        subsets = (
            file_combo_counts[(set_labels[0],)],
            file_combo_counts[(set_labels[1],)],
            file_combo_counts[(set_labels[0], set_labels[1])],
            file_combo_counts[(set_labels[2],)],
            file_combo_counts[(set_labels[0], set_labels[2])],
            file_combo_counts[(set_labels[1], set_labels[2])],
            file_combo_counts[set_labels]
        )
    else:
        print 'Charts only supported for 2 or 3 file groups.'
        return None

    return subsets


def create_venn(
    plot_title, description, outfile_path, set_labels, subsets,
    filegroup_count
):
    """Write Venn diagram to outfile_path"""
    plt.clf()
    if filegroup_count == 2:
        v = matplotlib_venn.venn2_unweighted(
            set_labels=set_labels, subsets=subsets
        )
    elif filegroup_count == 3:
        v = matplotlib_venn.venn3_unweighted(
            set_labels=set_labels, subsets=subsets
        )
    else:
        print 'Charts only supported for 2 or 3 files.'
        return

    for x in 'ABC'[:filegroup_count]:
        v.get_label_by_id(x).set_size(10)
    plt.suptitle(
        '%s\n\n' % plot_title, fontdict={'weight':'bold', 'size':'large'}
    )
    plt.title(description, fontdict={'size':'small'}, y=.96)
    plt.savefig(outfile_path, bbox_inches='tight', pad_inches=.5)


def create_pie_chart(
        plot_title, description, outfile_path, labels, subsets,
        default_colors, color_map=None
    ):
    plt.clf()
    plt.axis('equal')
    plt.suptitle(
        '%s\n\n' % plot_title, fontdict={'weight':'bold', 'size':'large'}
    )
    plt.title(description, fontdict={'size':'small'}, y=.97)
    most_freq = sorted(zip(labels, subsets), key=operator.itemgetter(1))[::-1]
    categ_limit = 6
    min_percent = 3
    if len(most_freq) > categ_limit:
        percentages = [
            100.0 * z[1] / sum([w[1] for w in most_freq]) for z in most_freq
        ]
        other_percent = sum([
            x for i, x in enumerate(percentages)
            if (x < min_percent or i >= categ_limit)
        ])
        percentages = [
            x for x in percentages[:categ_limit] if x >= min_percent
        ]
        subsets = [z[1] for z in most_freq[:len(percentages)]]
        labels = [z[0] for z in most_freq[:len(percentages)]]
        if other_percent > 0:
            subsets.append(sum([z[1] for z in most_freq[len(percentages):]]))
            percentages.append(other_percent)
            labels.append('other')
        for i in xrange(len(labels)):
            labels[i] = '%s: %.1f%%' % (labels[i], percentages[i])

    if color_map:
        try:
            colors = [color_map[x] for x in [y.split(':')[0] for y in labels]]
        except KeyError:
            colors = default_colors
    else:
        colors = default_colors
    patches, texts, autotexts = plt.pie(
        subsets, autopct='%1.1f%%', startangle=90, colors=colors,
        labels=[x.split(':')[0] for x in labels], labeldistance=1.05
    )
    for t in texts + autotexts:
        t.set_size('x-small')
    #plt.legend(
    #    patches, labels, loc='best', fontsize='x-small', ncol=2,
    #    framealpha=.2
    #)
    plt.savefig(outfile_path, bbox_inches='tight', pad_inches=.5)


def normalize_2d_list(list_):
    """Log of cell divided by sum(log of each cell in row/col)"""
    arr = np.array(list_)
    norm_list = []
    for row in arr:
        new_row = []
        for x in row:
            if x != 0:
                x = math.log(x)
            sum_row = sum([math.log(z) for z in row if z != 0])
            if sum_row:
                new_row.append(
                    float(x) / sum_row
                )
            else:
                new_row.append(0.0)
        norm_list.append(new_row)
    return norm_list


def create_heatmap(
    list_, plot_title, description, out_path, labels, xlabel, ylabel
):
    plt.clf()
    plt.suptitle(
        '%s\n\n' % plot_title, fontdict={'weight':'bold', 'size':'large'}
    )
    plt.title(description, fontdict={'size':'small'}, y=.98)
    plt.pcolor(np.array(list_), cmap=plt.cm.Reds, edgecolors='k')
    ax = plt.axes()
    # fix spines
    ax.spines['top'].set_bounds(0, len(list_))
    ax.spines['bottom'].set_bounds(0, len(list_))
    ax.spines['left'].set_bounds(0, len(list_))
    ax.spines['right'].set_bounds(0, len(list_))
    ax.spines['top'].set_position(('data', len(list_)))
    ax.spines['right'].set_position(('data', len(list_)))
    locs1, labels1 = plt.xticks(np.arange(len(list_[0])) + 0.5, labels)
    plt.yticks(np.arange(len(list_[0])) + 0.5, labels)
    plt.setp(labels1, rotation=90)
    plt.tick_params(
        axis='both', top='off', bottom='off', left='off', right='off',
        which='both'
    )
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10
    plt.xlabel(xlabel, style='oblique')
    plt.ylabel(ylabel, style='oblique')
    plt.savefig(out_path, bbox_inches='tight', pad_inches=.5)


def create_html(html_outfile, image_dir, page_subtitle):
    with open(html_outfile, 'w') as f:
        f.write(
            """
                <!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" 
                "http://www.w3.org/TR/html4/strict.dtd">
                <html><head>
                <meta http-equiv="Content-Type" content="text/html; 
                charset=UTF-8">
                <style type="text/css">
                body {
                    background-color: beige;
                    font-family: Verdana, Arial, sans-serif;
                }
                div{
                    text-align:center;
                }
                h1, h2, h3 {
                    text-align: center;
                }
                h1{
                    font-size:45px;
                }
                h2{
                    font-size:30px;
                }
                input{
                    height: 40px;
                    width: 70px;
                    font-size:150%;
                }
                #title_area{
                    margin: 50px;
                }
                .databox{
                    padding: 10px;
                    margin: 20px;
                }
                .databox1{
                    background-color: mistyrose;
                }
                .databox2{
                    background-color: lavender;
                }
                </style>
                <title>Annotator Comparison</title>
                </head><body>
            """
        )
        f.write(
            '<div id="title_area"><h1>Annotator Comparison</h1>'
            '<h3>%s</h3></div>' %
            re.sub('\n', '<br>', page_subtitle)
        )
        image_strs = {
            'heat' : [], 'consq' : [], 'other' : [], 'transcr' : []
        }
        image_fnames = sorted(os.listdir(image_dir))
        for image in image_fnames:
            with open(os.path.join(image_dir, image)) as img_file:
                image_str = (
                    '<img src="data:image.png;base64,%s"'
                    'alt="Annotator comparison plot">' %
                    base64.b64encode(img_file.read())
                )
                other_img_str = True
                for x in image_strs:
                    if x in image:
                        image_strs[x].append(image_str)
                        other_img_str = False
                if other_img_str:
                    image_strs['other'].append(image_str)

        # transcript similarity
        f.write(
            '<div class="databox databox2"><h2>Transcript Similarity</h2>'
        )
        create_js_slideshow(image_strs['transcr'], 'transcr', f)
        f.write('</div>')

        # consequence similarity
        f.write(
            '<div class="databox databox1"><h2>Consequence Similarity</h2>'
        )
        create_js_slideshow(image_strs['consq'], 'consq', f)
        f.write('</div>')

        # consequence heatmaps
        f.write('<div class="databox databox2"><h2>Consequence Heatmaps</h2>')
        create_js_slideshow(image_strs['heat'], 'heat', f)
        f.write('</div>')

        if image_strs['other']:
            f.write('<div class="databox databox2">')
            for image_str in image_strs['other']:
                f.write('<div>%s</div>' % image_str)
            f.write('</div>')

        f.write('</body></html>')


def create_js_slideshow(image_strs, slideshow_id, html_outfile):
    x = slideshow_id
    if not image_strs:
        return
    html_outfile.write(
        '<div id="%s">%s</div>' % (x, image_strs[0])
    )
    html_outfile.write(
        '<script type="text/javascript">'
        'var x%s = 0;' % x
    )
    for i, y in enumerate(image_strs):
        html_outfile.write(
            'var %s%d = "%s";' % (x, i, re.sub('"', '\'', y))
        )
    html_outfile.write(
        'function slideshowNext%s(){'
            'x%s = x%s + 1;'
            'if (x%s == %d) {x%s = 0};'
            'document.getElementById("%s").innerHTML = '
            'eval(\'%s\'+x%s)'
        '}' % (x, x, x, x, len(image_strs), x, x, x, x)
    )
    html_outfile.write(
        'function slideshowPrev%s(){'
            'x%s = x%s - 1;'
            'if (x%s < 0) {x%s = %d};'
            'document.getElementById("%s").innerHTML = '
            'eval(\'%s\'+x%s)'
        '}'
         % (x, x, x, x, x, len(image_strs) - 1, x, x, x)
    )
    html_outfile.write(
        '</script>'
        '<div>'
        '<input type="button" onclick="slideshowPrev%s()" value="<<">'
        '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
        '<input type="button" onclick="slideshowNext%s()" value=">>">'
        '</div>' % (x, x)
    )


if __name__ == '__main__':
    main()
