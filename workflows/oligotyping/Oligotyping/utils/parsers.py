#!/usr/bin/python
# -*- coding: utf-8

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import argparse
import Oligotyping as o


def decomposer():
    parser = argparse.ArgumentParser(description='Minimum Entropy Decomposition (version: %s)' % o.__version__)
    parser.add_argument('alignment', metavar = 'FILEPATH',
                        help = 'Alignment file that contains all samples and sequences in FASTA format')
    parser.add_argument('-m', '--min-entropy', type=float, default=0.0965, metavar="FLOAT",
                        help = 'Minimum entropy for a component to have in order to be picked as /path/to/input_dir
                                discriminant. Defeault: %(default)f')
    parser.add_argument('-X', '--skip-m-normalization', action = 'store_true', default = False,
                        help = 'Skip the m normalization heuristics. Normalization heuristics normalize /path/to/input_dir
                                parameter per node before decomposition depending on the node size, and the /path/to/input_dir
                                abundant sequence in the dataset. Please see the documentation for /path/to/input_dir
                                comprehensive description and best practices. Default: %(default)s')
    parser.add_argument('-d', '--number-of-discriminants', type=int, default=4, metavar="INTEGER",
                        help = 'Number of discriminant locations to be used for entropy /path/to/input_dir
                                discriminant. Defeault: %(default)d')
    parser.add_argument('-A', '--min-actual-abundance', type=int, default=0, metavar="INTEGER",
                        help = 'Minimum number of reads in a node for decomposition to continue. /path/to/input_dir
                                will continue for any node that has more reads than this number as far as /path/to/input_dir
                                present an entropy that is larger than --min-entropy. This number should /path/to/input_dir
                                chosen carefully depending on the size of the sample. Althought this /path/to/input_dir
                                is available to you for historical reasons, For noise filtering, you /path/to/input_dir
                                encouraged to use --min-substantive-abundance parameter instead.')
    parser.add_argument('-M', '--min-substantive-abundance', type=int, default=0, metavar = "INTEGER",
                        help = 'Unlike "actual" abundance, "substantive" abundance is interested in the /path/to/input_dir
                                of the most abundant read in a node. If the abundance of the most /path/to/input_dir
                                unique sequence in a node is smaller than the number given with this /path/to/input_dir
                                the node will be eliminated and not included in downstream analyses. This is /path/to/input_dir
                                most appropriate, and most cases, the only noise filtering parameter. If the /path/to/input_dir
                                does not set a value for minimum substantive abudnance, MED algorithm will /path/to/input_dir
                                one by default by dividing the number of reads in the input dataset by 5,000.')
    parser.add_argument('-V', '--maximum-variation-allowed', type=int, default=None, metavar = 'INTEGER',
                        help = 'This parameter is being used to remove "outliers" from nodes. The similarity of /path/to/input_dir
                                read in a node is less than --maximum-variation-allowed than the /path/to/input_dir
                                sequence of the node, it is identified as an outlier. If not set, this value is \
                                being computed depending on the average read length.')
    parser.add_argument('-t', '--sample-name-separator', type=str, default='_', metavar = "CHAR",
                        help = 'Character that separates sample name from unique info in the defline. For /path/to/input_dir
                                if the defline says >sample-1_GD7BRW402IVMZE, the separator should be set to "_"\
                                (which is the default character).')
    parser.add_argument('-o', '--output-directory', help = 'Output directory', default = None)
    parser.add_argument('-p', '--project', default = None, type=str, metavar = "STR",
                        help = 'When a project name is set, given name will be used in figures whenever possible.')
    parser.add_argument('-g', '--generate-frequency-curves', action = 'store_true', default = False,
                        help = 'When set, figure with frequency curve for unique reads and entropy /path/to/input_dir
                                will be generated for each node. Depending on the number of nodes, this /path/to/input_dir
                                be a time consuming step.')
    parser.add_argument('-S', '--skip-removing-outliers', action = 'store_true', default = False,
                        help = 'When set, outliers will not be removed from nodes.')
    parser.add_argument('-H', '--merge-homopolymer-splits', action = 'store_true', default = False,
                        help = 'When set, nodes that differ from each other by only one nucleotide that /path/to/input_dir
                                to be observed as an insertion at the upstream or downstream of a /path/to/input_dir
                                region will be merged.')
    parser.add_argument('-R', '--relocate-outliers', action = 'store_true', default = False,
                        help = 'Outliers are identified in two places: (1) during the raw topology /path/to/input_dir
                                and (2) during the refinement step where distant reads are removed from /path/to/input_dir
                                This parameter, when set, makes the pipeline go through each read identified /path/to/input_dir
                                an outlier and try to find the best nodes for them. Please read the /path/to/input_dir
                                for details. This step might take a long time. Default: %(default)s')
    parser.add_argument('-F', '--store-topology-dict', action = 'store_true', default = False,
                        help = 'When set, topology dict with read ids will be generated. This may take a very /path/to/input_dir
                                disk space and computation time for large data sets')
    parser.add_argument('-K', '--keep-tmp', action = 'store_true', default = False,
                        help = 'When set, directory with temporary BLAST results will not be deleted at the end of /path/to/input_dir
                                run. It may be necessary to debug the results')
    parser.add_argument('-T', '--no-threading', action = 'store_true', default = False,
                        help = 'When set, decomposer does not spawn multiple threads. Default behavior /path/to/input_dir
                                multi-threaded.')
    parser.add_argument('-N', '--number-of-threads', type=int, default = None, metavar = "INTEGER",
                        help = 'Number of threads to use. It is a good idea to keep this number smaller than the /path/to/input_dir
                                of CPU cores available. If not set, this number will be set to 90%% of available /path/to/input_dir
                                or (available cores - 1) if 10%% of the cores is a number smaller than 1')    
    parser.add_argument('-E', '--sample-mapping', metavar = 'FILEPATH', default = None,
                        help = 'TAB delimited categorical mapping of samples to be used for /path/to/input_dir
                                visualizations. Refer to the tutorial for the file format')
    parser.add_argument('--skip-gen-html', action = 'store_true', default = False,
                        help = 'By default, the process will end by generating a static HTML page to browse the /path/to/input_dir
                                results. The use of this flag will skip that')
    parser.add_argument('--skip-gen-figures', action = 'store_true', default = False,
                        help = 'When set, decomposer will not attempt to generate figures post analysis')
    parser.add_argument('--skip-check-input-file', action = 'store_true', default = False,
                        help = 'When set, input FASTA will not be checked for potential errors')
    parser.add_argument('--skip-gexf-files', action = 'store_true', default = False,
                        help = 'When set, GEXF files for network and topology will not be generated')
    parser.add_argument('--quick', action = 'store_true', default = False,
                        help = 'When set, the pipeline will do only the essential steps, skipping /path/to/input_dir
                                auxiliary, even if other parameters require otherwise. Please do not use it other /path/to/input_dir
                                benchmarking or testing purposes')
    parser.add_argument('--version', action = 'store_true', default = False,
                        help = 'Print version and exit.')    

    return parser


def oligotyping():
    parser = argparse.ArgumentParser(description='Oligotyping (version: %s)' % o.__version__)
    parser.add_argument('alignment', metavar = 'INPUT ALIGNMENT',
                        help = 'Alignment file that contains all samples and sequences in FASTA format')
    parser.add_argument('entropy', metavar = 'ENTROPY',
                        help = 'File that contains the columns and the entropy values computer previously')
    parser.add_argument('-o', '--output-directory', help = 'Output directory', default = None)
    parser.add_argument('-c', '--number-of-auto-components', type=int, default=None,
                        help = 'Number of components to use from alignment to generate oligotypes. /path/to/input_dir
                                is "5", which is a completely arbitrary value. Number of components /path/to/input_dir
                                be determined after a careful examination of entropy figure.')
    parser.add_argument('--qual-scores-file', metavar = 'QUAL SCORES FILE',
                        help = 'FASTA formatted file that contains PHRED base call /path/to/input_dir
                                for each read in the alignment file')
    parser.add_argument('--qual-scores-dict', metavar = 'QUAL SCORES DICT',
                        help = 'Previously computed and serialized dictionary that /path/to/input_dir
                                PHRED base call values for each read in the alignment file. If /path/to/input_dir
                                provide --qual-scores-file, that file will be used to recompute /path/to/input_dir
                                dictionary and the file you refer with this parameter /path/to/input_dir
                                not be ignored')
    parser.add_argument('--qual-stats-dict', metavar = 'QUAL STATS DICT',
                        help = 'Previously computed and serialized dictionary that /path/to/input_dir
                                PHRED base call quality score statistics for the alignment file. /path/to/input_dir
                                you provide --qual-scores-dict, it will be used to recompute /path/to/input_dir
                                dictionary and the file you refer to with this parameter /path/to/input_dir
                                actually not be used')
    parser.add_argument('-q', '--min-base-quality', type=int, default=15,
                        help = 'Minimum quality score for each base in locations of interest of a read to /path/to/input_dir
                                considered in an oligotype. When base quality score files are provided, /path/to/input_dir
                                value makes sure that low quality bases that are more likely to be the /path/to/input_dir
                                of random sequencing errors do not create artificial oligotypes. Any read that /path/to/input_dir
                                less quality score than the given value, will simply be discarded. This /path/to/input_dir
                                only in effect when --qual-scores-file or --qual-scores-dict parameters are used. \
                                Defeault is %(default)d.')
    parser.add_argument('-C', '--selected-components', type=str, default=None,
                        help = 'Comma separated entropy components to be used during the oligotyping process.')
    parser.add_argument('-s', '--min-number-of-samples', type=int, default=1,
                        help = 'Minimum number of samples oligotype expected to appear. The deafult is "5", /path/to/input_dir
                                is another completely arbitrary value. This parameter should be defined /path/to/input_dir
                                on the number of samples included in the analysis. If there are 10 /path/to/input_dir
                                3 might be a good choice, if there are 5 samples, 1 would be a better /path/to/input_dir
                                depending on the study. Default is %(default)d.')
    parser.add_argument('-a', '--min-percent-abundance', type=float, default=0.0,
                        help = 'Minimum percent abundance of an oligotype in at least one sample. The /path/to/input_dir
                                is "0.0". Just like --min-number-of-samples parameter, this parameter too /path/to/input_dir
                                to eliminate oligotypes that are formed by sequencing errors occured at /path/to/input_dir
                                component of interest. The value should be decided based on the average /path/to/input_dir
                                of sequences every sample has.')
    parser.add_argument('-A', '--min-actual-abundance', type=int, default=0,
                        help = 'Minimum total abundance of an oligotype in all datastes. The /path/to/input_dir
                                is "0". If the total abundance of an oligotype is smaller than the number /path/to/input_dir
                                with this parameter, oligotype would be eliminated and not included in /path/to/input_dir
                                analyses. Default is %(default)d.')
    parser.add_argument('-M', '--min-substantive-abundance', type=int, default=0,
                        help = 'Unlike "actual" abundance, "substantive" abundance is interested in the /path/to/input_dir
                                of the most abundant read in an oligotype. If the abundance of the most /path/to/input_dir
                                unique sequence in an oligotype smaller than the number given with this /path/to/input_dir
                                the oligotype will be eliminated and not included in downstream analyses. /path/to/input_dir
                                is %(default)d.')
    parser.add_argument('-t', '--sample-name-separator', type=str, default='_',
                        help = 'Character that separates sample name from unique info in the defline. For /path/to/input_dir
                                if the defline says >sample-1_GD7BRW402IVMZE, the separator should be set to "_"\
                                (which is the default character).')
    parser.add_argument('-l', '--limit-representative-sequences', type=int, default=None,
                        help = 'At the end of the oligotyping sequences that are being represented by the /path/to/input_dir
                                oligotype are being uniqued and stored in separate files. The number of /path/to/input_dir
                                to keep from the frequency ordered list can be defined with this parameter /path/to/input_dir
                                -l 10 would make it possible that only first 10 sequence would be stored). /path/to/input_dir
                                is 0, which stores everything, but when the sample size is too big, this /path/to/input_dir
                                take up disk space.')
    parser.add_argument('--limit-oligotypes-to', type = str, default = None,
                        help = 'Comma separated list of oligotypes to be taken into account during the /path/to/input_dir
                                All other oligotypes will be discarded if a list of oligotypes is being /path/to/input_dir
                                with this parameter.')
    parser.add_argument('-e', '--exclude-oligotypes', type = str, default = None,
                        help = 'Comma separated list of oligotypes to be excluded from the the analysis.')
    parser.add_argument('--quick', action = 'store_true', default = False,
                        help = 'Some relatively insignificant parts of the analysis may take a lot of time, such /path/to/input_dir
                                generating figures for representative sequences. When this parameter is set, /path/to/input_dir
                                trivial steps would be skipped to give results as soon as possible.')
    parser.add_argument('--no-figures', action = 'store_true', default = False,
                        help = 'When set, no figures will be generated or displayed.')
    parser.add_argument('--blast-ref-db', default = None, type=str,
                        help = 'When set, BLAST search will be done locally against the ref db (local BLAST /path/to/input_dir
                                requires NCBI+ tools)')
    parser.add_argument('--colors-list-file', default = None, type=str,
                        help = 'Optional file that contains HTML color codes in each line to color oligotypes. /path/to/input_dir
                                of colors in the file has to be equal or greater than the number of /path/to/input_dir
                                oligotypes, for which colors are going to be used for.')
    parser.add_argument('--do-blast-search', action = 'store_true', default = False,
                        help = 'When set, representative sequences will be searched on NCBI.')
    parser.add_argument('--no-display', action = 'store_true', default = False,
                        help = 'When set, no figures will be shown.')
    parser.add_argument('--skip-gen-html', action = 'store_true', default = False,
                        help = 'By default, the process will end by generating a static HTML page to browse the /path/to/input_dir
                                results. The use of this flag will skip that')
    parser.add_argument('--generate-sets', action = 'store_true', default = False,
                        help = 'Agglomerate oligotypes into oligotype sets when their frequency patterns /path/to/input_dir
                        samples are similar. Oligotype sets simply put oligotypes into the same set if they \
                        co-occur in samples consistenly.')
    parser.add_argument('-K', '--keep-tmp', action = 'store_true', default = False,
                        help = 'When set, directory with temporary results will not be deleted at the end of /path/to/input_dir
                                run. It may be necessary to debug the results')
    parser.add_argument('-S', '--cosine-similarity-threshold', default = 0.1, type=float, metavar='COS_SIM_TR'/path/to/input_dir
                        help = 'This value is used to agglomerate oligotypes into higher order groups. The /path/to/input_dir
                                the threshold is, the more oligotypes will be pulled together. Cosine /path/to/input_dir
                                would return 0 for perfectly similar two vectors. Default is %(default)f.')
    parser.add_argument('-E', '--sample-mapping', metavar = 'FILEPATH', default = None,
                        help = 'TAB delimited categorical mapping of samples to be used for /path/to/input_dir
                                visualizations. Refer to the tutorial for the file format')
    parser.add_argument('--project', default = None, type=str,
                        help = 'When a project name is set, given name will be used in figures whenever possible.')
    parser.add_argument('--skip-check-input-file', action = 'store_true', default = False,
                        help = 'When set, input FASTA will not be checked for potential errors')
    parser.add_argument('--skip-basic-analyses', action = 'store_true', default = False,
                        help = 'When set, basic analyses, such as basic NMDS plots and clustering, will /path/to/input_dir
                                skipped')
    parser.add_argument('--skip-gexf-network-file', action = 'store_true', default = False,
                        help = 'When set, GEXF network file will not be generated')
    parser.add_argument('-T', '--no-threading', action = 'store_true', default = False,
                        help = 'When set, oligotyping will not spawn multiple threads. Default behavior /path/to/input_dir
                                multi-threaded whenever possible.')
    parser.add_argument('-N', '--number-of-threads', type=int, default = None, metavar = "INTEGER",
                        help = 'Number of threads to use. It is a good idea to keep this number smaller than the /path/to/input_dir
                                of CPU cores available. If not set, this number will be set to 90%% of available /path/to/input_dir
                                or (available cores - 1) if 10%% of the cores is a number smaller than 1')    
    parser.add_argument('--version', action = 'store_true', default = False,
                        help = 'Print version and exit.')    


    return parser


def entropy():
    parser = argparse.ArgumentParser(description='Entropy Analysis (version: %s)' % o.__version__)
    parser.add_argument('alignment', metavar = 'ALIGNMENT', help = 'Alignment /path/to/input_dir
                         that contains all samples and sequences in FASTA format')
    parser.add_argument('--qual-scores-file', metavar = 'QUAL SCORES FILE',
                        help = 'FASTA formatted file that contains PHRED base call /path/to/input_dir
                         for each read in the alignment file')
    parser.add_argument('--qual-scores-dict', metavar = 'QUAL SCORES DICT',
                        help = 'Previously computed and serialized dictionary that /path/to/input_dir
                        PHRED base call values for each read in the alignment file. If /path/to/input_dir
                        provide --qual-scores-file, that file will be used to recompute /path/to/input_dir
                        dictionary and the file you refer with this parameter /path/to/input_dir
                        not be ignored')
    parser.add_argument('--qual-stats-dict', metavar = 'QUAL STATS DICT',
                        help = 'Previously computed and serialized dictionary that /path/to/input_dir
                        PHRED base call quality score statistics for the alignment file. /path/to/input_dir
                        you provide --qual-scores-dict, it will be used to recompute /path/to/input_dir
                        dictionary and the file you refer to with this parameter /path/to/input_dir
                        actually not be used')
    parser.add_argument('--uniqued', action = 'store_true', default = False,
                        help = 'When set, entropy computation will assume that the /path/to/input_dir
                        in FASTA file are unique. Frequency information of unique /path/to/input_dir
                        must be stored in the deflines. Every defline in the FASTA /path/to/input_dir
                        must present the frequency information in this /path/to/input_dir
                        "freq:NUMBER", e.g. ">Read_ID|X|Y|freq:42", or ">Read_ID|freq:42|X|Y"')
    parser.add_argument('--weighted', action = 'store_true', default = False,
                        help = 'When set, entropy computation per column will /path/to/input_dir
                        mean quality score for each column.')
    parser.add_argument('--amino-acid-sequences', action = 'store_true', default = False,
                        help = 'If sequences are composed of amino acids, instead /path/to/input_dir
                                nucleotides.')
    parser.add_argument('--quick', action = 'store_true', default = False,
                        help = 'When set, entropy values will be shown as fast /path/to/input_dir
                                possible (some visualization steps will be skipped).')
    parser.add_argument('--no-display', action = 'store_true', default = False,
                                help = 'When set, no figures will be shown.')
    parser.add_argument('--version', action = 'store_true', default = False,
                        help = 'Print version and exit.')    

    return parser
