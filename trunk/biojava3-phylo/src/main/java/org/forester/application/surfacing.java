// $Id: surfacing.java,v 1.161 2009/12/09 20:13:52 cmzmasek Exp $
//
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.application;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.go.GoId;
import org.forester.go.GoNameSpace;
import org.forester.go.GoTerm;
import org.forester.go.GoUtils;
import org.forester.go.OBOparser;
import org.forester.go.PfamToGoMapping;
import org.forester.go.PfamToGoParser;
import org.forester.io.parsers.HmmscanPerDomainTableParser;
import org.forester.io.parsers.HmmscanPerDomainTableParser.INDIVIDUAL_SCORE_CUTOFF;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogenyinference.DistanceMatrix;
import org.forester.phylogenyinference.NeighborJoining;
import org.forester.phylogenyinference.CharacterStateMatrix.Format;
import org.forester.surfacing.BasicDomainSimilarityCalculator;
import org.forester.surfacing.BasicGenomeWideCombinableDomains;
import org.forester.surfacing.BasicSpecies;
import org.forester.surfacing.BinaryDomainCombination;
import org.forester.surfacing.CombinationsBasedPairwiseDomainSimilarityCalculator;
import org.forester.surfacing.DomainCountsBasedPairwiseSimilarityCalculator;
import org.forester.surfacing.DomainCountsDifferenceUtil;
import org.forester.surfacing.DomainId;
import org.forester.surfacing.DomainLengthsTable;
import org.forester.surfacing.DomainParsimonyCalculator;
import org.forester.surfacing.DomainSimilarity;
import org.forester.surfacing.DomainSimilarityCalculator;
import org.forester.surfacing.GenomeWideCombinableDomains;
import org.forester.surfacing.MappingResults;
import org.forester.surfacing.PairwiseDomainSimilarityCalculator;
import org.forester.surfacing.PairwiseGenomeComparator;
import org.forester.surfacing.PrintableDomainSimilarity;
import org.forester.surfacing.Protein;
import org.forester.surfacing.ProteinCountsBasedPairwiseDomainSimilarityCalculator;
import org.forester.surfacing.Species;
import org.forester.surfacing.SurfacingUtil;
import org.forester.surfacing.DomainSimilarity.DomainSimilarityScoring;
import org.forester.surfacing.DomainSimilarity.DomainSimilaritySortField;
import org.forester.surfacing.DomainSimilarityCalculator.Detailedness;
import org.forester.surfacing.GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder;
import org.forester.surfacing.PrintableDomainSimilarity.PRINT_OPTION;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public class surfacing {

    public final static String                                DOMAIN_COMBINITONS_OUTPUT_OPTION_FOR_GRAPH_ANALYSIS                    = "graph_analysis_out";
    public final static String                                DOMAIN_COMBINITONS_OUTPUTFILE_SUFFIX_FOR_GRAPH_ANALYSIS                = "_dc.dot";
    public final static String                                PARSIMONY_OUTPUT_FITCH_PRESENT_BC_OUTPUTFILE_SUFFIX_FOR_GRAPH_ANALYSIS = "_fitch_present_dc.dot";
    public final static String                                DOMAIN_COMBINITON_COUNTS_OUTPUTFILE_SUFFIX                             = ".dcc";
    // gain/loss:
    public final static String                                PARSIMONY_OUTPUT_GL_SUFFIX_DOLLO_DOMAINS                               = "_dollo_gl_d";
    public final static String                                PARSIMONY_OUTPUT_GL_SUFFIX_DOLLO_BINARY_COMBINATIONS                   = "_dollo_gl_dc";
    public final static String                                PARSIMONY_OUTPUT_GL_SUFFIX_FITCH_DOMAINS                               = "_fitch_gl_d";
    public final static String                                PARSIMONY_OUTPUT_GL_SUFFIX_FITCH_BINARY_COMBINATIONS                   = "_fitch_gl_dc";
    // gain/loss counts:
    public final static String                                PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_DOLLO_DOMAINS                        = "_dollo_glc_d";
    public final static String                                PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_DOLLO_BINARY_COMBINATIONS            = "_dollo_glc_dc";
    public final static String                                PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_FITCH_DOMAINS                        = "_fitch_glc_d";
    public final static String                                PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_FITCH_BINARY_COMBINATIONS            = "_fitch_glc_dc";
    // tables:
    public final static String                                PARSIMONY_OUTPUT_FITCH_GAINS_BC                                        = "_fitch_gains_dc";
    public final static String                                PARSIMONY_OUTPUT_FITCH_GAINS_HTML_BC                                   = "_fitch_gains_dc.html";
    public final static String                                PARSIMONY_OUTPUT_FITCH_LOSSES_BC                                       = "_fitch_losses_dc";
    public final static String                                PARSIMONY_OUTPUT_FITCH_LOSSES_HTML_BC                                  = "_fitch_losses_dc.html";
    public final static String                                PARSIMONY_OUTPUT_FITCH_PRESENT_BC                                      = "_fitch_present_dc";
    public final static String                                PARSIMONY_OUTPUT_FITCH_PRESENT_HTML_BC                                 = "_fitch_present_dc.html";
    public final static String                                PARSIMONY_OUTPUT_DOLLO_GAINS_D                                         = "_dollo_gains_d";
    public final static String                                PARSIMONY_OUTPUT_DOLLO_GAINS_GOID_D                                    = "_dollo_gains_goid_d";
    public final static String                                PARSIMONY_OUTPUT_DOLLO_GAINS_HTML_D                                    = "_dollo_gains_d.html";
    public final static String                                PARSIMONY_OUTPUT_DOLLO_LOSSES_D                                        = "_dollo_losses_d";
    public final static String                                PARSIMONY_OUTPUT_DOLLO_LOSSES_HTML_D                                   = "_dollo_losses_d.html";
    public final static String                                PARSIMONY_OUTPUT_DOLLO_PRESENT_D                                       = "_dollo_present_d";
    public final static String                                PARSIMONY_OUTPUT_DOLLO_PRESENT_GOID_D                                  = "_dollo_present_goid_d";
    public final static String                                PARSIMONY_OUTPUT_DOLLO_PRESENT_HTML_D                                  = "_dollo_present_d.html";
    public final static String                                DOMAINS_PRESENT_NEXUS                                                  = "_dom.nex";
    public final static String                                BDC_PRESENT_NEXUS                                                      = "_dc.nex";
    // ---
    public final static String                                PRG_NAME                                                               = "surfacing";
    public static final String                                DOMAINS_PARSIMONY_TREE_OUTPUT_SUFFIX_DOLLO                             = "_d_dollo"
                                                                                                                                             + ForesterConstants.PHYLO_XML_SUFFIX;
    public static final String                                DOMAINS_PARSIMONY_TREE_OUTPUT_SUFFIX_FITCH                             = "_d_fitch"
                                                                                                                                             + ForesterConstants.PHYLO_XML_SUFFIX;
    public static final String                                BINARY_DOMAIN_COMBINATIONS_PARSIMONY_TREE_OUTPUT_SUFFIX_DOLLO          = "_dc_dollo"
                                                                                                                                             + ForesterConstants.PHYLO_XML_SUFFIX;
    public static final String                                BINARY_DOMAIN_COMBINATIONS_PARSIMONY_TREE_OUTPUT_SUFFIX_FITCH          = "_dc_fitch"
                                                                                                                                             + ForesterConstants.PHYLO_XML_SUFFIX;
    public static final String                                NEXUS_EXTERNAL_DOMAINS                                                 = "_dom.nex";
    public static final String                                NEXUS_EXTERNAL_DOMAIN_COMBINATIONS                                     = "_dc.nex";
    public static final String                                NEXUS_SECONDARY_FEATURES                                               = "_secondary_features.nex";
    public static final String                                PARSIMONY_OUTPUT_GL_SUFFIX_DOLLO_SECONDARY_FEATURES                    = "_dollo_gl_secondary_features";
    public static final String                                PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_DOLLO_SECONDARY_FEATURES             = "_dollo_glc_secondary_features";
    public static final String                                PARSIMONY_OUTPUT_DOLLO_GAINS_SECONDARY_FEATURES                        = "_dollo_gains_secondary_features";
    public static final String                                PARSIMONY_OUTPUT_DOLLO_LOSSES_SECONDARY_FEATURES                       = "_dollo_losses_secondary_features";
    public static final String                                PARSIMONY_OUTPUT_DOLLO_PRESENT_SECONDARY_FEATURES                      = "_dollo_present_secondary_features";
    public static final String                                SECONDARY_FEATURES_PARSIMONY_TREE_OUTPUT_SUFFIX_DOLLO                  = "_secondary_features_dollo"
                                                                                                                                             + ForesterConstants.PHYLO_XML_SUFFIX;
    public static final String                                PARSIMONY_OUTPUT_DOLLO_ALL_GOID_D_BIOLOGICAL_PROCESS                   = "_dollo_biol_proc_goid_d";
    public static final String                                PARSIMONY_OUTPUT_DOLLO_ALL_GOID_D_CELLULAR_COMPONENT                   = "_dollo_cell_comp_goid_d";
    public static final String                                PARSIMONY_OUTPUT_DOLLO_ALL_GOID_D_MOLECULAR_FUNCTION                   = "_dollo_mol_funct_goid_d";
    public static final String                                PARSIMONY_OUTPUT_DOLLO_ALL_GOID_D_ALL_NAMESPACES                       = "_dollo_goid_d";
    public static final String                                PARSIMONY_OUTPUT_FITCH_ALL_GOID_BC_BIOLOGICAL_PROCESS                  = "_fitch_biol_proc_goid_dc";
    public static final String                                PARSIMONY_OUTPUT_FITCH_ALL_GOID_BC_CELLULAR_COMPONENT                  = "_fitch_cell_comp_goid_dc";
    public static final String                                PARSIMONY_OUTPUT_FITCH_ALL_GOID_BC_MOLECULAR_FUNCTION                  = "_fitch_mol_funct_goid_dc";
    public static final String                                PARSIMONY_OUTPUT_FITCH_ALL_GOID_BC_ALL_NAMESPACES                      = "_fitch_goid_dc";
    final static private String                               HELP_OPTION_1                                                          = "help";
    final static private String                               HELP_OPTION_2                                                          = "h";
    final static private String                               OUTPUT_DIR_OPTION                                                      = "out_dir";
    final static private String                               SCORING_OPTION                                                         = "scoring";
    private static final DomainSimilarityScoring              SCORING_DEFAULT                                                        = DomainSimilarity.DomainSimilarityScoring.COMBINATIONS;
    final static private String                               SCORING_DOMAIN_COUNT_BASED                                             = "domains";
    final static private String                               SCORING_PROTEIN_COUNT_BASED                                            = "proteins";
    final static private String                               SCORING_COMBINATION_BASED                                              = "combinations";
    final static private String                               DETAILEDNESS_OPTION                                                    = "detail";
    private final static Detailedness                         DETAILEDNESS_DEFAULT                                                   = DomainSimilarityCalculator.Detailedness.PUNCTILIOUS;
    final static private String                               SPECIES_MATRIX_OPTION                                                  = "smatrix";
    final static private String                               DETAILEDNESS_BASIC                                                     = "basic";
    final static private String                               DETAILEDNESS_LIST_IDS                                                  = "list_ids";
    final static private String                               DETAILEDNESS_PUNCTILIOUS                                               = "punctilious";
    final static private String                               DOMAIN_SIMILARITY_SORT_OPTION                                          = "sort";
    private static final DomainSimilaritySortField            DOMAIN_SORT_FILD_DEFAULT                                               = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
    final static private String                               DOMAIN_SIMILARITY_SORT_MIN                                             = "min";
    final static private String                               DOMAIN_SIMILARITY_SORT_MAX                                             = "max";
    final static private String                               DOMAIN_SIMILARITY_SORT_SD                                              = "sd";
    final static private String                               DOMAIN_SIMILARITY_SORT_MEAN                                            = "mean";
    final static private String                               DOMAIN_SIMILARITY_SORT_DIFF                                            = "diff";
    final static private String                               DOMAIN_SIMILARITY_SORT_COUNTS_DIFF                                     = "count_diff";
    final static private String                               DOMAIN_SIMILARITY_SORT_ABS_COUNTS_DIFF                                 = "abs_count_diff";
    final static private String                               DOMAIN_SIMILARITY_SORT_SPECIES_COUNT                                   = "species";
    final static private String                               DOMAIN_SIMILARITY_SORT_ALPHA                                           = "alpha";
    final static private String                               DOMAIN_SIMILARITY_SORT_BY_SPECIES_COUNT_FIRST_OPTION                   = "species_first";
    final static private String                               DOMAIN_COUNT_SORT_OPTION                                               = "dc_sort";
    private static final GenomeWideCombinableDomainsSortOrder DOMAINS_SORT_ORDER_DEFAULT                                             = GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder.ALPHABETICAL_KEY_ID;
    final static private String                               DOMAIN_COUNT_SORT_ALPHA                                                = "alpha";
    final static private String                               DOMAIN_COUNT_SORT_KEY_DOMAIN_COUNT                                     = "dom";
    final static private String                               DOMAIN_COUNT_SORT_KEY_DOMAIN_PROTEINS_COUNT                            = "prot";
    final static private String                               DOMAIN_COUNT_SORT_COMBINATIONS_COUNT                                   = "comb";
    final static private String                               CUTOFF_SCORE_FILE_OPTION                                               = "cos";
    final static private String                               NOT_IGNORE_DUFS_OPTION                                                 = "dufs";
    final static private String                               MAX_E_VALUE_OPTION                                                     = "e";
    final static private String                               MAX_ALLOWED_OVERLAP_OPTION                                             = "mo";
    final static private String                               NO_ENGULFING_OVERLAP_OPTION                                            = "no_eo";
    final static private String                               IGNORE_COMBINATION_WITH_SAME_OPTION                                    = "ignore_self_comb";
    final static private String                               PAIRWISE_DOMAIN_COMPARISONS_PREFIX                                     = "pwc_";
    final static private String                               PAIRWISE_DOMAIN_COMPARISONS_OPTION                                     = "pwc";
    final static private String                               OUTPUT_FILE_OPTION                                                     = "o";
    final static private String                               PFAM_TO_GO_FILE_USE_OPTION                                             = "p2g";
    final static private String                               GO_OBO_FILE_USE_OPTION                                                 = "obo";
    final static private String                               GO_NAMESPACE_LIMIT_OPTION                                              = "go_namespace";
    final static private String                               GO_NAMESPACE_LIMIT_OPTION_MOLECULAR_FUNCTION                           = "molecular_function";
    final static private String                               GO_NAMESPACE_LIMIT_OPTION_BIOLOGICAL_PROCESS                           = "biological_process";
    final static private String                               GO_NAMESPACE_LIMIT_OPTION_CELLULAR_COMPONENT                           = "cellular_component";
    final static private String                               SECONDARY_FEATURES_PARSIMONY_MAP_FILE                                  = "secondary";
    final static private String                               DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_TAB_DELIMITED                    = "simple_tab";
    final static private String                               DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_HTML                             = "simple_html";
    final static private String                               DOMAIN_SIMILARITY_PRINT_OPTION_DETAILED_HTML                           = "detailed_html";
    final static private String                               DOMAIN_SIMILARITY_PRINT_OPTION                                         = "ds_output";
    private static final PRINT_OPTION                         DOMAIN_SIMILARITY_PRINT_OPTION_DEFAULT                                 = PrintableDomainSimilarity.PRINT_OPTION.HTML;
    final static private String                               IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_OPTION              = "ignore_singlet_domains";
    final static private String                               IGNORE_VIRAL_IDS                                                       = "ignore_viral_ids";
    final static private boolean                              IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_DEFAULT             = false;
    final static private String                               IGNORE_DOMAINS_SPECIFIC_TO_ONE_SPECIES_OPTION                          = "ignore_species_specific_domains";
    final static private boolean                              IGNORE_DOMAINS_SPECIFIC_TO_ONE_SPECIES_OPTION_DEFAULT                  = false;
    final static private String                               MATRIX_MEAN_SCORE_BASED_GENOME_DISTANCE_SUFFIX                         = "_mean_score.pwd";
    final static private String                               MATRIX_SHARED_DOMAINS_BASED_GENOME_DISTANCE_SUFFIX                     = "_domains.pwd";
    final static private String                               MATRIX_SHARED_BIN_COMBINATIONS_BASED_GENOME_DISTANCE_SUFFIX            = "_bin_combinations.pwd";
    final static private String                               NJ_TREE_MEAN_SCORE_BASED_GENOME_DISTANCE_SUFFIX                        = "_mean_score_NJ"
                                                                                                                                             + ForesterConstants.PHYLO_XML_SUFFIX;
    final static private String                               NJ_TREE_SHARED_DOMAINS_BASED_GENOME_DISTANCE_SUFFIX                    = "_domains_NJ"
                                                                                                                                             + ForesterConstants.PHYLO_XML_SUFFIX;
    final static private String                               NJ_TREE_SHARED_BIN_COMBINATIONS_BASED_GENOME_DISTANCE_SUFFIX           = "_bin_combinations_NJ"
                                                                                                                                             + ForesterConstants.PHYLO_XML_SUFFIX;
    final static private String                               DISPLAY_M_HISTOGRAMS_OPTION                                            = "mhisto";
    //  final static private boolean DISPLAY_M_HISTOGRAMS_OPTION_DEFAULT                                    = false;
    final static private String                               JACKNIFE_OPTION                                                        = "jack";
    final static private String                               JACKNIFE_RANDOM_SEED_OPTION                                            = "seed";
    final static private String                               JACKNIFE_RATIO_OPTION                                                  = "jack_ratio";
    private static final int                                  JACKNIFE_NUMBER_OF_RESAMPLINGS_DEFAULT                                 = 100;
    final static private long                                 JACKNIFE_RANDOM_SEED_DEFAULT                                           = 19;
    final static private double                               JACKNIFE_RATIO_DEFAULT                                                 = 0.5;
    //final static private String  INFER_SPECIES_TREES_OPTION                                             = "species_tree_inference";
    final static private String                               INFERRED_SD_BASED_NJ_SPECIES_TREE_SUFFIX                               = "_sd_nj.nh";
    final static private String                               INFERRED_SBC_BASED_NJ_SPECIES_TREE_SUFFIX                              = "_sbc_nj.nh";
    final static private String                               FILTER_POSITIVE_OPTION                                                 = "pos_filter";
    final static private String                               FILTER_NEGATIVE_OPTION                                                 = "neg_filter";
    final static private String                               FILTER_NEGATIVE_DOMAINS_OPTION                                         = "neg_dom_filter";
    final static private String                               INPUT_FILES_FROM_FILE_OPTION                                           = "input";
    final static private String                               INPUT_SPECIES_TREE_OPTION                                              = "species_tree";
    final static private String                               SEQ_EXTRACT_OPTION                                                     = "prot_extract";
    final static private char                                 SEPARATOR_FOR_INPUT_VALUES                                             = '#';
    final static private String                               PRG_VERSION                                                            = "2.00";
    final static private String                               PRG_DATE                                                               = "2009.11.10";
    final static private String                               E_MAIL                                                                 = "czmasek@burnham.org";
    final static private String                               WWW                                                                    = "www.phylosoft.org/forester/applications/surfacing";
    final static private boolean                              IGNORE_DUFS_DEFAULT                                                    = true;
    final static private boolean                              IGNORE_COMBINATION_WITH_SAME_DEFAULLT                                  = false;
    final static private double                               MAX_E_VALUE_DEFAULT                                                    = -1;
    final static private int                                  MAX_ALLOWED_OVERLAP_DEFAULT                                            = -1;
    final static private String                               DEFAULT_SEARCH_PARAMETER                                               = "ls";
    final private static boolean                              VERBOSE_DEFAULT                                                        = true;
    private static final String                               RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION                                 = "random_seed";
    private static final String                               CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS                               = "consider_bdc_direction";
    private static final String                               CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS_AND_ADJACENCY                 = "consider_bdc_adj";
    private static final String                               SEQ_EXTRACT_SUFFIX                                                     = ".prot";
    private static final String                               PLUS_MINUS_ANALYSIS_OPTION                                             = "plus_minus";
    private static final String                               PLUS_MINUS_DOM_SUFFIX                                                  = "_plus_minus_dom.txt";
    private static final String                               PLUS_MINUS_DOM_SUFFIX_HTML                                             = "_plus_minus_dom.html";
    private static final String                               PLUS_MINUS_DC_SUFFIX_HTML                                              = "_plus_minus_dc.html";
    private static final int                                  PLUS_MINUS_ANALYSIS_MIN_DIFF_DEFAULT                                   = 0;
    private static final double                               PLUS_MINUS_ANALYSIS_FACTOR_DEFAULT                                     = 1.0;
    private static final String                               PLUS_MINUS_ALL_GO_IDS_DOM_SUFFIX                                       = "_plus_minus_go_ids_all.txt";
    private static final String                               PLUS_MINUS_PASSING_GO_IDS_DOM_SUFFIX                                   = "_plus_minus_go_ids_passing.txt";
    private static final String                               OUTPUT_LIST_OF_ALL_PROTEINS_OPTIONS                                    = "all_prot";
    private static final boolean                              VERBOSE                                                                = false;
    private static final String                               OUTPUT_DOMAIN_COMBINATIONS_GAINED_MORE_THAN_ONCE_ANALYSIS_SUFFIX       = "_fitch_dc_gains_counts";
    private static final String                               OUTPUT_DOMAIN_COMBINATIONS_LOST_MORE_THAN_ONCE_ANALYSIS_SUFFIX         = "_fitch_dc_losses_counts";
    private static final String                               DOMAIN_LENGTHS_ANALYSIS_SUFFIX                                         = "_domain_lengths_analysis";
    private static final boolean                              PERFORM_DOMAIN_LENGTH_ANALYSIS                                         = true;
    public static final String                                ALL_PFAMS_ENCOUNTERED_SUFFIX                                           = "_all_encountered_pfams";
    public static final String                                ALL_PFAMS_ENCOUNTERED_WITH_GO_ANNOTATION_SUFFIX                        = "_all_encountered_pfams_with_go_annotation";
    public static final String                                ENCOUNTERED_PFAMS_SUMMARY_SUFFIX                                       = "_encountered_pfams_summary";
    public static final String                                ALL_PFAMS_GAINED_AS_DOMAINS_SUFFIX                                     = "_all_pfams_gained_as_domains";
    public static final String                                ALL_PFAMS_LOST_AS_DOMAINS_SUFFIX                                       = "_all_pfams_lost_as_domains";
    public static final String                                ALL_PFAMS_GAINED_AS_DC_SUFFIX                                          = "_all_pfams_gained_as_dc";
    public static final String                                ALL_PFAMS_LOST_AS_DC_SUFFIX                                            = "_all_pfams_lost_as_dc";
    public static final String                                BASE_DIRECTORY_PER_NODE_DOMAIN_GAIN_LOSS_FILES                         = "PER_NODE_EVENTS";
    public static final String                                BASE_DIRECTORY_PER_SUBTREE_DOMAIN_GAIN_LOSS_FILES                      = "PER_SUBTREE_EVENTS";
    public static final String                                D_PROMISCUITY_FILE_SUFFIX                                              = "_domain_promiscuities";
    private static final String                               LOG_FILE_SUFFIX                                                        = "_log.txt";
    private static final INDIVIDUAL_SCORE_CUTOFF              INDIVIDUAL_SCORE_CUTOFF_DEFAULT                                        = INDIVIDUAL_SCORE_CUTOFF.FULL_SEQUENCE;

    // final String error = ForesterUtil.isReadableFile( new File(
    // input_file_properties[ i ][ 0 ] ) );
    // if ( !ForesterUtil.isEmpty( error ) ) {
    // ForesterUtil.fatalError( surfacing.PRG_NAME, error );
    // }
    private static void checkWriteabilityForPairwiseComparisons( final PrintableDomainSimilarity.PRINT_OPTION domain_similarity_print_option,
                                                                 final String[][] input_file_properties,
                                                                 final String automated_pairwise_comparison_suffix,
                                                                 final File outdir ) {
        for( int i = 0; i < input_file_properties.length; ++i ) {
            for( int j = 0; j < i; ++j ) {
                final String species_i = input_file_properties[ i ][ 1 ];
                final String species_j = input_file_properties[ j ][ 1 ];
                String pairwise_similarities_output_file_str = PAIRWISE_DOMAIN_COMPARISONS_PREFIX + species_i + "_"
                        + species_j + automated_pairwise_comparison_suffix;
                switch ( domain_similarity_print_option ) {
                    case HTML:
                        if ( !pairwise_similarities_output_file_str.endsWith( ".html" ) ) {
                            pairwise_similarities_output_file_str += ".html";
                        }
                        break;
                }
                final String error = ForesterUtil
                        .isWritableFile( new File( outdir == null ? pairwise_similarities_output_file_str : outdir
                                + ForesterUtil.FILE_SEPARATOR + pairwise_similarities_output_file_str ) );
                if ( !ForesterUtil.isEmpty( error ) ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, error );
                }
            }
        }
    }

    private static StringBuilder createParametersAsString( final boolean ignore_dufs,
                                                           final double e_value_max,
                                                           final int max_allowed_overlap,
                                                           final boolean no_engulfing_overlaps,
                                                           final File cutoff_scores_file,
                                                           final BinaryDomainCombination.DomainCombinationType dc_type ) {
        final StringBuilder parameters_sb = new StringBuilder();
        parameters_sb.append( "E-value: " + e_value_max );
        if ( cutoff_scores_file != null ) {
            parameters_sb.append( ", Cutoff-scores-file: " + cutoff_scores_file );
        }
        else {
            parameters_sb.append( ", Cutoff-scores-file: not-set" );
        }
        if ( max_allowed_overlap != surfacing.MAX_ALLOWED_OVERLAP_DEFAULT ) {
            parameters_sb.append( ", Max-overlap: " + max_allowed_overlap );
        }
        else {
            parameters_sb.append( ", Max-overlap: not-set" );
        }
        if ( no_engulfing_overlaps ) {
            parameters_sb.append( ", Engulfing-overlaps: not-allowed" );
        }
        else {
            parameters_sb.append( ", Engulfing-overlaps: allowed" );
        }
        if ( ignore_dufs ) {
            parameters_sb.append( ", Ignore-dufs: true" );
        }
        else {
            parameters_sb.append( ", Ignore-dufs: false" );
        }
        parameters_sb.append( ", DC type (if applicable): " + dc_type );
        return parameters_sb;
    }

    /**
     * Warning: This sideeffects 'all_bin_domain_combinations_encountered'!
     * 
     * 
     * @param output_file
     * @param all_bin_domain_combinations_changed
     * @param sum_of_all_domains_encountered
     * @param all_bin_domain_combinations_encountered
     * @param is_gains_analysis
     * @throws IOException
     */
    private static void executeFitchGainsAnalysis( final File output_file,
                                                   final List<BinaryDomainCombination> all_bin_domain_combinations_changed,
                                                   final int sum_of_all_domains_encountered,
                                                   final SortedSet<BinaryDomainCombination> all_bin_domain_combinations_encountered,
                                                   final boolean is_gains_analysis ) throws IOException {
        SurfacingUtil.checkForOutputFileWriteability( output_file );
        final Writer out = ForesterUtil.createBufferedWriter( output_file );
        final SortedMap<Object, Integer> bdc_to_counts = ForesterUtil
                .listToSortedCountsMap( all_bin_domain_combinations_changed );
        final SortedSet<DomainId> all_domains_in_combination_changed_more_than_once = new TreeSet<DomainId>();
        final SortedSet<DomainId> all_domains_in_combination_changed_only_once = new TreeSet<DomainId>();
        int above_one = 0;
        int one = 0;
        for( final Object bdc_object : bdc_to_counts.keySet() ) {
            final BinaryDomainCombination bdc = ( BinaryDomainCombination ) bdc_object;
            final int count = bdc_to_counts.get( bdc_object );
            if ( count < 1 ) {
                ForesterUtil.unexpectedFatalError( PRG_NAME, "count < 1 " );
            }
            out.write( bdc + "\t" + count + ForesterUtil.LINE_SEPARATOR );
            if ( count > 1 ) {
                all_domains_in_combination_changed_more_than_once.add( bdc.getId0() );
                all_domains_in_combination_changed_more_than_once.add( bdc.getId1() );
                above_one++;
            }
            else if ( count == 1 ) {
                all_domains_in_combination_changed_only_once.add( bdc.getId0() );
                all_domains_in_combination_changed_only_once.add( bdc.getId1() );
                one++;
            }
        }
        final int all = all_bin_domain_combinations_encountered.size();
        int never_lost = -1;
        if ( !is_gains_analysis ) {
            all_bin_domain_combinations_encountered.removeAll( all_bin_domain_combinations_changed );
            never_lost = all_bin_domain_combinations_encountered.size();
            for( final BinaryDomainCombination bdc : all_bin_domain_combinations_encountered ) {
                out.write( bdc + "\t" + "0" + ForesterUtil.LINE_SEPARATOR );
            }
        }
        if ( is_gains_analysis ) {
            out.write( "Sum of all distinct domain combinations appearing once               : " + one
                    + ForesterUtil.LINE_SEPARATOR );
            out.write( "Sum of all distinct domain combinations appearing more than once     : " + above_one
                    + ForesterUtil.LINE_SEPARATOR );
            out.write( "Sum of all distinct domains in combinations apppearing only once     : "
                    + all_domains_in_combination_changed_only_once.size() + ForesterUtil.LINE_SEPARATOR );
            out.write( "Sum of all distinct domains in combinations apppearing more than once: "
                    + all_domains_in_combination_changed_more_than_once.size() + ForesterUtil.LINE_SEPARATOR );
        }
        else {
            out.write( "Sum of all distinct domain combinations never lost                   : " + never_lost
                    + ForesterUtil.LINE_SEPARATOR );
            out.write( "Sum of all distinct domain combinations lost once                    : " + one
                    + ForesterUtil.LINE_SEPARATOR );
            out.write( "Sum of all distinct domain combinations lost more than once          : " + above_one
                    + ForesterUtil.LINE_SEPARATOR );
            out.write( "Sum of all distinct domains in combinations lost only once           : "
                    + all_domains_in_combination_changed_only_once.size() + ForesterUtil.LINE_SEPARATOR );
            out.write( "Sum of all distinct domains in combinations lost more than once: "
                    + all_domains_in_combination_changed_more_than_once.size() + ForesterUtil.LINE_SEPARATOR );
        }
        out.write( "All binary combinations                                              : " + all
                + ForesterUtil.LINE_SEPARATOR );
        out.write( "All domains                                                          : "
                + sum_of_all_domains_encountered );
        out.close();
        ForesterUtil.programMessage( surfacing.PRG_NAME,
                                     "Wrote fitch domain combination dynamics counts analysis to \"" + output_file
                                             + "\"" );
    }

    private static void executePlusMinusAnalysis( final File output_file,
                                                  final List<String> plus_minus_analysis_high_copy_base,
                                                  final List<String> plus_minus_analysis_high_copy_target,
                                                  final List<String> plus_minus_analysis_low_copy,
                                                  final List<GenomeWideCombinableDomains> gwcd_list,
                                                  final SortedMap<Species, List<Protein>> protein_lists_per_species,
                                                  final Map<DomainId, List<GoId>> domain_id_to_go_ids_map,
                                                  final Map<GoId, GoTerm> go_id_to_term_map,
                                                  final List<Object> plus_minus_analysis_numbers ) {
        final Set<String> all_spec = new HashSet<String>();
        for( final GenomeWideCombinableDomains gwcd : gwcd_list ) {
            all_spec.add( gwcd.getSpecies().getSpeciesId() );
        }
        final File html_out_dom = new File( output_file + PLUS_MINUS_DOM_SUFFIX_HTML );
        final File plain_out_dom = new File( output_file + PLUS_MINUS_DOM_SUFFIX );
        final File html_out_dc = new File( output_file + PLUS_MINUS_DC_SUFFIX_HTML );
        final File all_domains_go_ids_out_dom = new File( output_file + PLUS_MINUS_ALL_GO_IDS_DOM_SUFFIX );
        final File passing_domains_go_ids_out_dom = new File( output_file + PLUS_MINUS_PASSING_GO_IDS_DOM_SUFFIX );
        final File proteins_file_base = new File( output_file + "" );
        final int min_diff = ( ( Integer ) plus_minus_analysis_numbers.get( 0 ) ).intValue();
        final double factor = ( ( Double ) plus_minus_analysis_numbers.get( 1 ) ).doubleValue();
        try {
            DomainCountsDifferenceUtil.calculateCopyNumberDifferences( gwcd_list,
                                                                       protein_lists_per_species,
                                                                       plus_minus_analysis_high_copy_base,
                                                                       plus_minus_analysis_high_copy_target,
                                                                       plus_minus_analysis_low_copy,
                                                                       min_diff,
                                                                       factor,
                                                                       plain_out_dom,
                                                                       html_out_dom,
                                                                       html_out_dc,
                                                                       domain_id_to_go_ids_map,
                                                                       go_id_to_term_map,
                                                                       all_domains_go_ids_out_dom,
                                                                       passing_domains_go_ids_out_dom,
                                                                       proteins_file_base );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getLocalizedMessage() );
        }
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote plus minus domain analysis results to \""
                + html_out_dom + "\"" );
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote plus minus domain analysis results to \""
                + plain_out_dom + "\"" );
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote plus minus domain analysis results to \"" + html_out_dc
                + "\"" );
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote plus minus domain analysis based passing GO ids to \""
                + passing_domains_go_ids_out_dom + "\"" );
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote plus minus domain analysis based all GO ids to \""
                + all_domains_go_ids_out_dom + "\"" );
    }

    private static Phylogeny[] getIntrees( final File[] intree_files,
                                           final int number_of_genomes,
                                           final String[][] input_file_properties ) {
        final Phylogeny[] intrees = new Phylogeny[ intree_files.length ];
        int i = 0;
        for( final File intree_file : intree_files ) {
            Phylogeny intree = null;
            final String error = ForesterUtil.isReadableFile( intree_file );
            if ( !ForesterUtil.isEmpty( error ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot read input tree file [" + intree_file + "]: "
                        + error );
            }
            try {
                final Phylogeny[] p_array = ParserBasedPhylogenyFactory.getInstance()
                        .create( intree_file, ForesterUtil.createParserDependingOnFileType( intree_file, true ) );
                if ( p_array.length < 1 ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "file [" + intree_file
                            + "] does not contain any phylogeny in phyloXML format" );
                }
                else if ( p_array.length > 1 ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "file [" + intree_file
                            + "] contains more than one phylogeny in phyloXML format" );
                }
                intree = p_array[ 0 ];
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "failed to read input tree from file [" + intree_file
                        + "]: " + error );
            }
            if ( ( intree == null ) || intree.isEmpty() ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "input tree [" + intree_file + "] is empty" );
            }
            if ( !intree.isRooted() ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "input tree [" + intree_file + "] is not rooted" );
            }
            if ( intree.getNumberOfExternalNodes() < number_of_genomes ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "number of external nodes ["
                        + intree.getNumberOfExternalNodes() + "] of input tree [" + intree_file
                        + "] is smaller than the number of genomes the be analyzed [" + number_of_genomes + "]" );
            }
            final StringBuilder parent_names = new StringBuilder();
            final int nodes_lacking_name = SurfacingUtil.getNumberOfNodesLackingName( intree, parent_names );
            if ( nodes_lacking_name > 0 ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "input tree [" + intree_file + "] has "
                        + nodes_lacking_name + " node(s) lacking a name [parent names:" + parent_names + "]" );
            }
            preparePhylogenyForParsimonyAnalyses( intree, input_file_properties );
            if ( !intree.isCompletelyBinary() ) {
                ForesterUtil.printWarningMessage( surfacing.PRG_NAME, "input tree [" + intree_file
                        + "] is not completely binary" );
            }
            intrees[ i++ ] = intree;
        }
        return intrees;
    }

    private static List<Phylogeny> inferSpeciesTrees( final File outfile, final List<DistanceMatrix> distances_list ) {
        final NeighborJoining nj = NeighborJoining.createInstance();
        final List<Phylogeny> phylogenies = nj.execute( distances_list );
        final PhylogenyWriter w = new PhylogenyWriter();
        try {
            w.toNewHampshire( phylogenies, true, true, outfile, ";" );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to outfile [" + outfile + "]: " + e.getMessage() );
        }
        return phylogenies;
    }

    private static void log( final String msg, final Writer w ) {
        try {
            w.write( msg );
            w.write( ForesterUtil.LINE_SEPARATOR );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getLocalizedMessage() );
        }
    }

    public static void main( final String args[] ) {
        final long start_time = new Date().getTime();
        // final StringBuffer log = new StringBuffer();
        final StringBuilder html_desc = new StringBuilder();
        ForesterUtil.printProgramInformation( surfacing.PRG_NAME,
                                              surfacing.PRG_VERSION,
                                              surfacing.PRG_DATE,
                                              surfacing.E_MAIL,
                                              surfacing.WWW );
        final String nl = ForesterUtil.LINE_SEPARATOR;
        html_desc.append( "<table>" + nl );
        html_desc.append( "<tr><td>Produced by:</td><td>" + surfacing.PRG_NAME + "</td></tr>" + nl );
        html_desc.append( "<tr><td>Version:</td><td>" + surfacing.PRG_VERSION + "</td></tr>" + nl );
        html_desc.append( "<tr><td>Release Date:</td><td>" + surfacing.PRG_DATE + "</td></tr>" + nl );
        html_desc.append( "<tr><td>Contact:</td><td>" + surfacing.E_MAIL + "</td></tr>" + nl );
        html_desc.append( "<tr><td>WWW:</td><td>" + surfacing.WWW + "</td></tr>" + nl );
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( surfacing.HELP_OPTION_1 ) || cla.isOptionSet( surfacing.HELP_OPTION_2 ) ) {
            surfacing.printHelp();
            System.exit( 0 );
        }
        if ( ( args.length < 1 ) ) {
            surfacing.printHelp();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( surfacing.NOT_IGNORE_DUFS_OPTION );
        allowed_options.add( surfacing.MAX_E_VALUE_OPTION );
        allowed_options.add( surfacing.DETAILEDNESS_OPTION );
        allowed_options.add( surfacing.OUTPUT_FILE_OPTION );
        allowed_options.add( surfacing.DOMAIN_SIMILARITY_SORT_OPTION );
        allowed_options.add( surfacing.SPECIES_MATRIX_OPTION );
        allowed_options.add( surfacing.SCORING_OPTION );
        allowed_options.add( surfacing.MAX_ALLOWED_OVERLAP_OPTION );
        allowed_options.add( surfacing.NO_ENGULFING_OVERLAP_OPTION );
        allowed_options.add( surfacing.DOMAIN_COUNT_SORT_OPTION );
        allowed_options.add( surfacing.CUTOFF_SCORE_FILE_OPTION );
        allowed_options.add( surfacing.DOMAIN_SIMILARITY_SORT_BY_SPECIES_COUNT_FIRST_OPTION );
        allowed_options.add( surfacing.OUTPUT_DIR_OPTION );
        allowed_options.add( surfacing.IGNORE_COMBINATION_WITH_SAME_OPTION );
        allowed_options.add( surfacing.PFAM_TO_GO_FILE_USE_OPTION );
        allowed_options.add( surfacing.GO_OBO_FILE_USE_OPTION );
        allowed_options.add( surfacing.DOMAIN_SIMILARITY_PRINT_OPTION );
        allowed_options.add( surfacing.GO_NAMESPACE_LIMIT_OPTION );
        allowed_options.add( surfacing.PAIRWISE_DOMAIN_COMPARISONS_OPTION );
        allowed_options.add( surfacing.IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_OPTION );
        allowed_options.add( surfacing.DISPLAY_M_HISTOGRAMS_OPTION );
        allowed_options.add( surfacing.CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS );
        allowed_options.add( JACKNIFE_OPTION );
        allowed_options.add( JACKNIFE_RANDOM_SEED_OPTION );
        allowed_options.add( JACKNIFE_RATIO_OPTION );
        allowed_options.add( INPUT_SPECIES_TREE_OPTION );
        //allowed_options.add( INFER_SPECIES_TREES_OPTION );
        allowed_options.add( FILTER_POSITIVE_OPTION );
        allowed_options.add( FILTER_NEGATIVE_OPTION );
        allowed_options.add( INPUT_FILES_FROM_FILE_OPTION );
        allowed_options.add( RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION );
        allowed_options.add( FILTER_NEGATIVE_DOMAINS_OPTION );
        allowed_options.add( IGNORE_VIRAL_IDS );
        allowed_options.add( SEQ_EXTRACT_OPTION );
        allowed_options.add( SECONDARY_FEATURES_PARSIMONY_MAP_FILE );
        allowed_options.add( PLUS_MINUS_ANALYSIS_OPTION );
        allowed_options.add( DOMAIN_COMBINITONS_OUTPUT_OPTION_FOR_GRAPH_ANALYSIS );
        allowed_options.add( OUTPUT_LIST_OF_ALL_PROTEINS_OPTIONS );
        allowed_options.add( CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS_AND_ADJACENCY );
        boolean ignore_dufs = surfacing.IGNORE_DUFS_DEFAULT;
        boolean ignore_combination_with_same = surfacing.IGNORE_COMBINATION_WITH_SAME_DEFAULLT;
        double e_value_max = surfacing.MAX_E_VALUE_DEFAULT;
        int max_allowed_overlap = surfacing.MAX_ALLOWED_OVERLAP_DEFAULT;
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        boolean output_binary_domain_combinationsfor_graph_analysis = false;
        if ( cla.isOptionSet( DOMAIN_COMBINITONS_OUTPUT_OPTION_FOR_GRAPH_ANALYSIS ) ) {
            output_binary_domain_combinationsfor_graph_analysis = true;
        }
        if ( cla.isOptionSet( surfacing.MAX_E_VALUE_OPTION ) ) {
            try {
                e_value_max = cla.getOptionValueAsDouble( surfacing.MAX_E_VALUE_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no acceptable value for E-value maximum" );
            }
        }
        if ( cla.isOptionSet( surfacing.MAX_ALLOWED_OVERLAP_OPTION ) ) {
            try {
                max_allowed_overlap = cla.getOptionValueAsInt( surfacing.MAX_ALLOWED_OVERLAP_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no acceptable value for maximal allowed domain overlap" );
            }
        }
        boolean no_engulfing_overlaps = false;
        if ( cla.isOptionSet( surfacing.NO_ENGULFING_OVERLAP_OPTION ) ) {
            no_engulfing_overlaps = true;
        }
        boolean ignore_virus_like_ids = false;
        if ( cla.isOptionSet( surfacing.IGNORE_VIRAL_IDS ) ) {
            ignore_virus_like_ids = true;
        }
        if ( cla.isOptionSet( surfacing.NOT_IGNORE_DUFS_OPTION ) ) {
            ignore_dufs = false;
        }
        if ( cla.isOptionSet( surfacing.IGNORE_COMBINATION_WITH_SAME_OPTION ) ) {
            ignore_combination_with_same = true;
        }
        boolean ignore_domains_without_combs_in_all_spec = IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_DEFAULT;
        if ( cla.isOptionSet( surfacing.IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_OPTION ) ) {
            ignore_domains_without_combs_in_all_spec = true;
        }
        boolean ignore_species_specific_domains = IGNORE_DOMAINS_SPECIFIC_TO_ONE_SPECIES_OPTION_DEFAULT;
        if ( cla.isOptionSet( surfacing.IGNORE_DOMAINS_SPECIFIC_TO_ONE_SPECIES_OPTION ) ) {
            ignore_species_specific_domains = true;
        }
        File output_file = null;
        if ( cla.isOptionSet( surfacing.OUTPUT_FILE_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.OUTPUT_FILE_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME,
                                         "no value for domain combinations similarities output file: -"
                                                 + surfacing.OUTPUT_FILE_OPTION + "=<file>" );
            }
            output_file = new File( cla.getOptionValue( surfacing.OUTPUT_FILE_OPTION ) );
            SurfacingUtil.checkForOutputFileWriteability( output_file );
        }
        File cutoff_scores_file = null;
        Map<String, Double> individual_score_cutoffs = null;
        if ( cla.isOptionSet( surfacing.CUTOFF_SCORE_FILE_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.CUTOFF_SCORE_FILE_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for individual domain score cutoffs file: -"
                        + surfacing.CUTOFF_SCORE_FILE_OPTION + "=<file>" );
            }
            cutoff_scores_file = new File( cla.getOptionValue( surfacing.CUTOFF_SCORE_FILE_OPTION ) );
            final String error = ForesterUtil.isReadableFile( cutoff_scores_file );
            if ( !ForesterUtil.isEmpty( error ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot read individual domain score cutoffs file: "
                        + error );
            }
            try {
                final BasicTable<String> scores_table = BasicTableParser.parse( cutoff_scores_file, " " );
                individual_score_cutoffs = scores_table.getColumnsAsMapDouble( 0, 1 );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot read from individual score cutoffs file: " + e );
            }
        }
        BinaryDomainCombination.DomainCombinationType dc_type = BinaryDomainCombination.DomainCombinationType.BASIC;
        if ( cla.isOptionSet( surfacing.CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS ) ) {
            dc_type = BinaryDomainCombination.DomainCombinationType.DIRECTED;
        }
        if ( cla.isOptionSet( surfacing.CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS_AND_ADJACENCY ) ) {
            dc_type = BinaryDomainCombination.DomainCombinationType.DIRECTED_ADJACTANT;
        }
        File out_dir = null;
        if ( cla.isOptionSet( surfacing.OUTPUT_DIR_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.OUTPUT_DIR_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for output directory: -"
                        + surfacing.OUTPUT_DIR_OPTION + "=<dir>" );
            }
            out_dir = new File( cla.getOptionValue( surfacing.OUTPUT_DIR_OPTION ) );
            if ( out_dir.exists() && ( out_dir.listFiles().length > 0 ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "\"" + out_dir + "\" aready exists and is not empty" );
            }
            if ( !out_dir.exists() ) {
                final boolean success = out_dir.mkdir();
                if ( !success || !out_dir.exists() ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "failed to create \"" + out_dir + "\"" );
                }
            }
            if ( !out_dir.canWrite() ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot write to \"" + out_dir + "\"" );
            }
        }
        File positive_filter_file = null;
        File negative_filter_file = null;
        File negative_domains_filter_file = null;
        if ( cla.isOptionSet( surfacing.FILTER_NEGATIVE_OPTION ) && cla.isOptionSet( surfacing.FILTER_POSITIVE_OPTION ) ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, "attempt to use both negative and positive protein filter" );
        }
        if ( cla.isOptionSet( surfacing.FILTER_NEGATIVE_DOMAINS_OPTION )
                && ( cla.isOptionSet( surfacing.FILTER_NEGATIVE_OPTION ) || cla
                        .isOptionSet( surfacing.FILTER_POSITIVE_OPTION ) ) ) {
            ForesterUtil
                    .fatalError( surfacing.PRG_NAME,
                                 "attempt to use both negative or positive protein filter together wirh a negative domains filter" );
        }
        if ( cla.isOptionSet( surfacing.FILTER_NEGATIVE_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.FILTER_NEGATIVE_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for negative filter: -"
                        + surfacing.FILTER_NEGATIVE_OPTION + "=<file>" );
            }
            negative_filter_file = new File( cla.getOptionValue( surfacing.FILTER_NEGATIVE_OPTION ) );
            final String msg = ForesterUtil.isReadableFile( negative_filter_file );
            if ( !ForesterUtil.isEmpty( msg ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "can not read from \"" + negative_filter_file + "\": "
                        + msg );
            }
        }
        else if ( cla.isOptionSet( surfacing.FILTER_POSITIVE_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.FILTER_POSITIVE_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for positive filter: -"
                        + surfacing.FILTER_POSITIVE_OPTION + "=<file>" );
            }
            positive_filter_file = new File( cla.getOptionValue( surfacing.FILTER_POSITIVE_OPTION ) );
            final String msg = ForesterUtil.isReadableFile( positive_filter_file );
            if ( !ForesterUtil.isEmpty( msg ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "can not read from \"" + positive_filter_file + "\": "
                        + msg );
            }
        }
        else if ( cla.isOptionSet( surfacing.FILTER_NEGATIVE_DOMAINS_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.FILTER_NEGATIVE_DOMAINS_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for negative domains filter: -"
                        + surfacing.FILTER_NEGATIVE_DOMAINS_OPTION + "=<file>" );
            }
            negative_domains_filter_file = new File( cla.getOptionValue( surfacing.FILTER_NEGATIVE_DOMAINS_OPTION ) );
            final String msg = ForesterUtil.isReadableFile( negative_domains_filter_file );
            if ( !ForesterUtil.isEmpty( msg ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "can not read from \"" + negative_domains_filter_file
                        + "\": " + msg );
            }
        }
        final List<String> plus_minus_analysis_high_copy_base_species = new ArrayList<String>();
        final List<String> plus_minus_analysis_high_copy_target_species = new ArrayList<String>();
        final List<String> plus_minus_analysis_high_low_copy_species = new ArrayList<String>();
        final List<Object> plus_minus_analysis_numbers = new ArrayList<Object>();
        processPlusMinusAnalysisOption( cla,
                                        plus_minus_analysis_high_copy_base_species,
                                        plus_minus_analysis_high_copy_target_species,
                                        plus_minus_analysis_high_low_copy_species,
                                        plus_minus_analysis_numbers );
        File input_files_file = null;
        String[] input_file_names_from_file = null;
        if ( cla.isOptionSet( surfacing.INPUT_FILES_FROM_FILE_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.INPUT_FILES_FROM_FILE_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for input files file: -"
                        + surfacing.INPUT_FILES_FROM_FILE_OPTION + "=<file>" );
            }
            input_files_file = new File( cla.getOptionValue( surfacing.INPUT_FILES_FROM_FILE_OPTION ) );
            final String msg = ForesterUtil.isReadableFile( input_files_file );
            if ( !ForesterUtil.isEmpty( msg ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "can not read from \"" + input_files_file + "\": " + msg );
            }
            try {
                input_file_names_from_file = ForesterUtil.file2array( input_files_file );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "failed to read from \"" + input_files_file + "\": " + e );
            }
        }
        if ( ( cla.getNumberOfNames() < 1 )
                && ( ( input_file_names_from_file == null ) || ( input_file_names_from_file.length < 1 ) ) ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME,
                                     "No hmmpfam output file indicated is input: use comand line directly or "
                                             + surfacing.INPUT_FILES_FROM_FILE_OPTION + "=<file>" );
        }
        DomainSimilarity.DomainSimilarityScoring scoring = SCORING_DEFAULT;
        if ( cla.isOptionSet( surfacing.SCORING_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.SCORING_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME,
                                         "no value for scoring method for domain combinations similarity calculation: -"
                                                 + surfacing.SCORING_OPTION + "=<"
                                                 + surfacing.SCORING_DOMAIN_COUNT_BASED + "|"
                                                 + surfacing.SCORING_PROTEIN_COUNT_BASED + "|"
                                                 + surfacing.SCORING_COMBINATION_BASED + ">\"" );
            }
            final String scoring_str = cla.getOptionValue( surfacing.SCORING_OPTION );
            if ( scoring_str.equals( surfacing.SCORING_DOMAIN_COUNT_BASED ) ) {
                scoring = DomainSimilarity.DomainSimilarityScoring.DOMAINS;
            }
            else if ( scoring_str.equals( surfacing.SCORING_COMBINATION_BASED ) ) {
                scoring = DomainSimilarity.DomainSimilarityScoring.COMBINATIONS;
            }
            else if ( scoring_str.equals( surfacing.SCORING_PROTEIN_COUNT_BASED ) ) {
                scoring = DomainSimilarity.DomainSimilarityScoring.PROTEINS;
            }
            else {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "unknown value \"" + scoring_str
                        + "\" for scoring method for domain combinations similarity calculation: \"-"
                        + surfacing.SCORING_OPTION + "=<" + surfacing.SCORING_DOMAIN_COUNT_BASED + "|"
                        + surfacing.SCORING_PROTEIN_COUNT_BASED + "|" + surfacing.SCORING_COMBINATION_BASED + ">\"" );
            }
        }
        boolean sort_by_species_count_first = false;
        if ( cla.isOptionSet( surfacing.DOMAIN_SIMILARITY_SORT_BY_SPECIES_COUNT_FIRST_OPTION ) ) {
            sort_by_species_count_first = true;
        }
        boolean species_matrix = false;
        if ( cla.isOptionSet( surfacing.SPECIES_MATRIX_OPTION ) ) {
            species_matrix = true;
        }
        boolean output_protein_lists_for_all_domains = false;
        if ( cla.isOptionSet( surfacing.OUTPUT_LIST_OF_ALL_PROTEINS_OPTIONS ) ) {
            output_protein_lists_for_all_domains = true;
        }
        Detailedness detailedness = DETAILEDNESS_DEFAULT;
        if ( cla.isOptionSet( surfacing.DETAILEDNESS_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.DETAILEDNESS_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for -" + surfacing.DETAILEDNESS_OPTION + "=<"
                        + surfacing.DETAILEDNESS_BASIC + "|" + surfacing.DETAILEDNESS_LIST_IDS + "|"
                        + surfacing.DETAILEDNESS_PUNCTILIOUS + ">\"" );
            }
            final String detness = cla.getOptionValue( surfacing.DETAILEDNESS_OPTION ).toLowerCase();
            if ( detness.equals( surfacing.DETAILEDNESS_BASIC ) ) {
                detailedness = DomainSimilarityCalculator.Detailedness.BASIC;
            }
            else if ( detness.equals( surfacing.DETAILEDNESS_LIST_IDS ) ) {
                detailedness = DomainSimilarityCalculator.Detailedness.LIST_COMBINING_DOMAIN_FOR_EACH_SPECIES;
            }
            else if ( detness.equals( surfacing.DETAILEDNESS_PUNCTILIOUS ) ) {
                detailedness = DomainSimilarityCalculator.Detailedness.PUNCTILIOUS;
            }
            else {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "unknown value \"" + detness + "\" for detailedness: \"-"
                        + surfacing.DETAILEDNESS_OPTION + "=<" + surfacing.DETAILEDNESS_BASIC + "|"
                        + surfacing.DETAILEDNESS_LIST_IDS + "|" + surfacing.DETAILEDNESS_PUNCTILIOUS + ">\"" );
            }
        }
        String automated_pairwise_comparison_suffix = null;
        boolean perform_pwc = false;
        boolean write_pwc_files = false;
        if ( cla.isOptionSet( surfacing.PAIRWISE_DOMAIN_COMPARISONS_OPTION ) ) {
            perform_pwc = true;
            if ( !cla.isOptionValueSet( surfacing.PAIRWISE_DOMAIN_COMPARISONS_OPTION ) ) {
                write_pwc_files = false;
            }
            else {
                write_pwc_files = true;
                automated_pairwise_comparison_suffix = "_"
                        + cla.getOptionValue( surfacing.PAIRWISE_DOMAIN_COMPARISONS_OPTION );
            }
        }
        String query_domain_ids = null;
        if ( cla.isOptionSet( surfacing.SEQ_EXTRACT_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.SEQ_EXTRACT_OPTION ) ) {
                ForesterUtil
                        .fatalError( surfacing.PRG_NAME,
                                     "no domain ids given for sequences with given domains to be extracted : -"
                                             + surfacing.SEQ_EXTRACT_OPTION
                                             + "=<ordered domain sequences, domain ids separated by '~', sequences separated by '#'>" );
            }
            query_domain_ids = cla.getOptionValue( surfacing.SEQ_EXTRACT_OPTION );
        }
        DomainSimilarity.DomainSimilaritySortField domain_similarity_sort_field = DOMAIN_SORT_FILD_DEFAULT;
        DomainSimilarity.DomainSimilaritySortField domain_similarity_sort_field_for_automated_pwc = DOMAIN_SORT_FILD_DEFAULT;
        if ( cla.isOptionSet( surfacing.DOMAIN_SIMILARITY_SORT_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.DOMAIN_SIMILARITY_SORT_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for domain combinations similarities sorting: -"
                        + surfacing.DOMAIN_SIMILARITY_SORT_OPTION + "=<" + surfacing.DOMAIN_SIMILARITY_SORT_ALPHA + "|"
                        + surfacing.DOMAIN_SIMILARITY_SORT_MAX + "|" + surfacing.DOMAIN_SIMILARITY_SORT_MIN + "|"
                        + surfacing.DOMAIN_SIMILARITY_SORT_MEAN + "|" + surfacing.DOMAIN_SIMILARITY_SORT_DIFF + "|"
                        + surfacing.DOMAIN_SIMILARITY_SORT_ABS_COUNTS_DIFF + "|"
                        + surfacing.DOMAIN_SIMILARITY_SORT_COUNTS_DIFF + "|"
                        + surfacing.DOMAIN_SIMILARITY_SORT_SPECIES_COUNT + "|" + surfacing.DOMAIN_SIMILARITY_SORT_SD
                        + ">\"" );
            }
            final String sort_str = cla.getOptionValue( surfacing.DOMAIN_SIMILARITY_SORT_OPTION ).toLowerCase();
            if ( sort_str.equals( surfacing.DOMAIN_SIMILARITY_SORT_ALPHA ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
            }
            else if ( sort_str.equals( surfacing.DOMAIN_SIMILARITY_SORT_MAX ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.MAX;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
            }
            else if ( sort_str.equals( surfacing.DOMAIN_SIMILARITY_SORT_MIN ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.MIN;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
            }
            else if ( sort_str.equals( surfacing.DOMAIN_SIMILARITY_SORT_MEAN ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.MEAN;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.MEAN;
            }
            else if ( sort_str.equals( surfacing.DOMAIN_SIMILARITY_SORT_SPECIES_COUNT ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.SPECIES_COUNT;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
            }
            else if ( sort_str.equals( surfacing.DOMAIN_SIMILARITY_SORT_SD ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.SD;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.DOMAIN_ID;
            }
            else if ( sort_str.equals( surfacing.DOMAIN_SIMILARITY_SORT_DIFF ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.MAX_DIFFERENCE;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.MAX_DIFFERENCE;
            }
            else if ( sort_str.equals( surfacing.DOMAIN_SIMILARITY_SORT_ABS_COUNTS_DIFF ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.ABS_MAX_COUNTS_DIFFERENCE;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.ABS_MAX_COUNTS_DIFFERENCE;
            }
            else if ( sort_str.equals( surfacing.DOMAIN_SIMILARITY_SORT_COUNTS_DIFF ) ) {
                domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.MAX_COUNTS_DIFFERENCE;
                domain_similarity_sort_field_for_automated_pwc = DomainSimilarity.DomainSimilaritySortField.MAX_COUNTS_DIFFERENCE;
            }
            else {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "unknown value \"" + sort_str
                        + "\" for domain combinations similarities sorting: \"-"
                        + surfacing.DOMAIN_SIMILARITY_SORT_OPTION + "=<" + surfacing.DOMAIN_SIMILARITY_SORT_ALPHA + "|"
                        + surfacing.DOMAIN_SIMILARITY_SORT_MAX + "|" + surfacing.DOMAIN_SIMILARITY_SORT_MIN + "|"
                        + surfacing.DOMAIN_SIMILARITY_SORT_MEAN + "|" + surfacing.DOMAIN_SIMILARITY_SORT_DIFF + "|"
                        + surfacing.DOMAIN_SIMILARITY_SORT_ABS_COUNTS_DIFF + "|"
                        + surfacing.DOMAIN_SIMILARITY_SORT_COUNTS_DIFF + "|" + "|"
                        + surfacing.DOMAIN_SIMILARITY_SORT_SPECIES_COUNT + "|" + surfacing.DOMAIN_SIMILARITY_SORT_SD
                        + ">\"" );
            }
        }
        PrintableDomainSimilarity.PRINT_OPTION domain_similarity_print_option = DOMAIN_SIMILARITY_PRINT_OPTION_DEFAULT;
        if ( cla.isOptionSet( surfacing.DOMAIN_SIMILARITY_PRINT_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.DOMAIN_SIMILARITY_PRINT_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for print option: -"
                        + surfacing.DOMAIN_SIMILARITY_PRINT_OPTION_DETAILED_HTML + "|"
                        + surfacing.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_HTML + "|"
                        + surfacing.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_TAB_DELIMITED + ">\"" );
            }
            final String sort = cla.getOptionValue( surfacing.DOMAIN_SIMILARITY_PRINT_OPTION ).toLowerCase();
            if ( sort.equals( surfacing.DOMAIN_SIMILARITY_PRINT_OPTION_DETAILED_HTML ) ) {
                domain_similarity_print_option = PrintableDomainSimilarity.PRINT_OPTION.HTML;
            }
            else if ( sort.equals( surfacing.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_HTML ) ) {
                // domain_similarity_print_option =
                // DomainSimilarity.PRINT_OPTION.SIMPLE_HTML;
                ForesterUtil.fatalError( surfacing.PRG_NAME, "simple HTML output not implemented yet :(" );
            }
            else if ( sort.equals( surfacing.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_TAB_DELIMITED ) ) {
                domain_similarity_print_option = PrintableDomainSimilarity.PRINT_OPTION.SIMPLE_TAB_DELIMITED;
            }
            else {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "unknown value \"" + sort + "\" for print option: -"
                        + surfacing.DOMAIN_SIMILARITY_PRINT_OPTION_DETAILED_HTML + "|"
                        + surfacing.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_HTML + "|"
                        + surfacing.DOMAIN_SIMILARITY_PRINT_OPTION_SIMPLE_TAB_DELIMITED + ">\"" );
            }
        }
        GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder dc_sort_order = DOMAINS_SORT_ORDER_DEFAULT;
        if ( cla.isOptionSet( surfacing.DOMAIN_COUNT_SORT_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.DOMAIN_COUNT_SORT_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for sorting of domain counts: -"
                        + surfacing.DOMAIN_COUNT_SORT_OPTION + "=<" + surfacing.DOMAIN_COUNT_SORT_ALPHA + "|"
                        + surfacing.DOMAIN_COUNT_SORT_KEY_DOMAIN_COUNT + "|"
                        + surfacing.DOMAIN_COUNT_SORT_KEY_DOMAIN_PROTEINS_COUNT + "|"
                        + surfacing.DOMAIN_COUNT_SORT_COMBINATIONS_COUNT + ">\"" );
            }
            final String sort = cla.getOptionValue( surfacing.DOMAIN_COUNT_SORT_OPTION ).toLowerCase();
            if ( sort.equals( surfacing.DOMAIN_COUNT_SORT_ALPHA ) ) {
                dc_sort_order = GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder.ALPHABETICAL_KEY_ID;
            }
            else if ( sort.equals( surfacing.DOMAIN_COUNT_SORT_KEY_DOMAIN_COUNT ) ) {
                dc_sort_order = GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder.KEY_DOMAIN_COUNT;
            }
            else if ( sort.equals( surfacing.DOMAIN_COUNT_SORT_KEY_DOMAIN_PROTEINS_COUNT ) ) {
                dc_sort_order = GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder.KEY_DOMAIN_PROTEINS_COUNT;
            }
            else if ( sort.equals( surfacing.DOMAIN_COUNT_SORT_COMBINATIONS_COUNT ) ) {
                dc_sort_order = GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder.COMBINATIONS_COUNT;
            }
            else {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "unknown value \"" + sort
                        + "\" for sorting of domain counts: \"-" + surfacing.DOMAIN_COUNT_SORT_OPTION + "=<"
                        + surfacing.DOMAIN_COUNT_SORT_ALPHA + "|" + surfacing.DOMAIN_COUNT_SORT_KEY_DOMAIN_COUNT + "|"
                        + surfacing.DOMAIN_COUNT_SORT_KEY_DOMAIN_PROTEINS_COUNT + "|"
                        + surfacing.DOMAIN_COUNT_SORT_COMBINATIONS_COUNT + ">\"" );
            }
        }
        String[][] input_file_properties = null;
        if ( input_file_names_from_file != null ) {
            input_file_properties = surfacing.processInputFileNames( input_file_names_from_file );
        }
        else {
            input_file_properties = surfacing.processInputFileNames( cla.getNames() );
        }
        final int number_of_genomes = input_file_properties.length;
        if ( number_of_genomes < 2 ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot analyze less than two files" );
        }
        if ( ( number_of_genomes < 3 ) && perform_pwc ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot use : -"
                    + surfacing.PAIRWISE_DOMAIN_COMPARISONS_OPTION
                    + "=<suffix> to turn on pairwise analyses with less than three input files" );
        }
        checkWriteabilityForPairwiseComparisons( domain_similarity_print_option,
                                                 input_file_properties,
                                                 automated_pairwise_comparison_suffix,
                                                 out_dir );
        for( int i = 0; i < number_of_genomes; i++ ) {
            File dcc_outfile = new File( input_file_properties[ i ][ 0 ]
                    + surfacing.DOMAIN_COMBINITON_COUNTS_OUTPUTFILE_SUFFIX );
            if ( out_dir != null ) {
                dcc_outfile = new File( out_dir + ForesterUtil.FILE_SEPARATOR + dcc_outfile );
            }
            SurfacingUtil.checkForOutputFileWriteability( dcc_outfile );
        }
        File pfam_to_go_file = null;
        Map<DomainId, List<GoId>> domain_id_to_go_ids_map = null;
        int domain_id_to_go_ids_count = 0;
        if ( cla.isOptionSet( surfacing.PFAM_TO_GO_FILE_USE_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.PFAM_TO_GO_FILE_USE_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for Pfam to GO mapping file: -"
                        + surfacing.PFAM_TO_GO_FILE_USE_OPTION + "=<file>" );
            }
            pfam_to_go_file = new File( cla.getOptionValue( surfacing.PFAM_TO_GO_FILE_USE_OPTION ) );
            final String error = ForesterUtil.isReadableFile( pfam_to_go_file );
            if ( !ForesterUtil.isEmpty( error ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot read Pfam to GO mapping file: " + error );
            }
            try {
                final PfamToGoParser parser = new PfamToGoParser( pfam_to_go_file );
                final List<PfamToGoMapping> pfam_to_go_mappings = parser.parse();
                domain_id_to_go_ids_map = SurfacingUtil.createDomainIdToGoIdMap( pfam_to_go_mappings );
                if ( parser.getMappingCount() < domain_id_to_go_ids_map.size() ) {
                    ForesterUtil.unexpectedFatalError( surfacing.PRG_NAME,
                                                       "parser.getMappingCount() < domain_id_to_go_ids_map.size()" );
                }
                domain_id_to_go_ids_count = parser.getMappingCount();
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot read from Pfam to GO mapping file: " + e );
            }
        }
        File go_obo_file = null;
        List<GoTerm> go_terms = null;
        if ( cla.isOptionSet( surfacing.GO_OBO_FILE_USE_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.GO_OBO_FILE_USE_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for GO OBO file: -"
                        + surfacing.GO_OBO_FILE_USE_OPTION + "=<file>" );
            }
            if ( ( domain_id_to_go_ids_map == null ) || ( domain_id_to_go_ids_map.size() < 1 ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot use GO OBO file (-"
                        + surfacing.GO_OBO_FILE_USE_OPTION + "=<file>) without Pfam to GO mapping file ("
                        + surfacing.PFAM_TO_GO_FILE_USE_OPTION + "=<file>)" );
            }
            go_obo_file = new File( cla.getOptionValue( surfacing.GO_OBO_FILE_USE_OPTION ) );
            final String error = ForesterUtil.isReadableFile( go_obo_file );
            if ( !ForesterUtil.isEmpty( error ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot read GO OBO file: " + error );
            }
            try {
                final OBOparser parser = new OBOparser( go_obo_file, OBOparser.ReturnType.BASIC_GO_TERM );
                go_terms = parser.parse();
                if ( parser.getGoTermCount() != go_terms.size() ) {
                    ForesterUtil
                            .unexpectedFatalError( surfacing.PRG_NAME, "parser.getGoTermCount() != go_terms.size()" );
                }
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot read from GO OBO file: " + e );
            }
        }
        Map<GoId, GoTerm> go_id_to_term_map = null;
        if ( ( ( domain_id_to_go_ids_map != null ) && ( domain_id_to_go_ids_map.size() > 0 ) )
                && ( ( go_terms != null ) && ( go_terms.size() > 0 ) ) ) {
            go_id_to_term_map = GoUtils.createGoIdToGoTermMap( go_terms );
        }
        GoNameSpace go_namespace_limit = null;
        if ( cla.isOptionSet( surfacing.GO_NAMESPACE_LIMIT_OPTION ) ) {
            if ( ( go_id_to_term_map == null ) || go_id_to_term_map.isEmpty() ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot use GO namespace limit (-"
                        + surfacing.GO_NAMESPACE_LIMIT_OPTION + "=<namespace>) without Pfam to GO mapping file ("
                        + surfacing.PFAM_TO_GO_FILE_USE_OPTION + "=<file>) and GO OBO file (-"
                        + surfacing.GO_OBO_FILE_USE_OPTION + "=<file>)" );
            }
            if ( !cla.isOptionValueSet( surfacing.GO_NAMESPACE_LIMIT_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for GO namespace limit: \"-"
                        + surfacing.GO_NAMESPACE_LIMIT_OPTION + "=<"
                        + surfacing.GO_NAMESPACE_LIMIT_OPTION_MOLECULAR_FUNCTION + "|"
                        + surfacing.GO_NAMESPACE_LIMIT_OPTION_BIOLOGICAL_PROCESS + "|"
                        + surfacing.GO_NAMESPACE_LIMIT_OPTION_CELLULAR_COMPONENT + ">\"" );
            }
            final String go_namespace_limit_str = cla.getOptionValue( surfacing.GO_NAMESPACE_LIMIT_OPTION )
                    .toLowerCase();
            if ( go_namespace_limit_str.equals( surfacing.GO_NAMESPACE_LIMIT_OPTION_MOLECULAR_FUNCTION ) ) {
                go_namespace_limit = GoNameSpace.createMolecularFunction();
            }
            else if ( go_namespace_limit_str.equals( surfacing.GO_NAMESPACE_LIMIT_OPTION_BIOLOGICAL_PROCESS ) ) {
                go_namespace_limit = GoNameSpace.createBiologicalProcess();
            }
            else if ( go_namespace_limit_str.equals( surfacing.GO_NAMESPACE_LIMIT_OPTION_CELLULAR_COMPONENT ) ) {
                go_namespace_limit = GoNameSpace.createCellularComponent();
            }
            else {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "unknown value \"" + go_namespace_limit_str
                        + "\" for GO namespace limit: \"-" + surfacing.GO_NAMESPACE_LIMIT_OPTION + "=<"
                        + surfacing.GO_NAMESPACE_LIMIT_OPTION_MOLECULAR_FUNCTION + "|"
                        + surfacing.GO_NAMESPACE_LIMIT_OPTION_BIOLOGICAL_PROCESS + "|"
                        + surfacing.GO_NAMESPACE_LIMIT_OPTION_CELLULAR_COMPONENT + ">\"" );
            }
        }
        if ( ( domain_similarity_sort_field == DomainSimilarity.DomainSimilaritySortField.MAX_COUNTS_DIFFERENCE )
                && ( number_of_genomes > 2 ) ) {
            domain_similarity_sort_field = DomainSimilarity.DomainSimilaritySortField.ABS_MAX_COUNTS_DIFFERENCE;
        }
        boolean jacknifed_distances = false;
        int jacknife_resamplings = JACKNIFE_NUMBER_OF_RESAMPLINGS_DEFAULT;
        double jacknife_ratio = JACKNIFE_RATIO_DEFAULT;
        long random_seed = JACKNIFE_RANDOM_SEED_DEFAULT;
        if ( cla.isOptionSet( surfacing.JACKNIFE_OPTION ) ) {
            if ( ( number_of_genomes < 3 ) || !perform_pwc ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot use jacknife resampling analysis (-"
                        + surfacing.JACKNIFE_OPTION + "[=<number of resamplings>]) without pairwise analyses ("
                        + surfacing.PAIRWISE_DOMAIN_COMPARISONS_OPTION
                        + "=<suffix for pairwise comparison output files>)" );
            }
            jacknifed_distances = true;
            if ( cla.isOptionHasAValue( surfacing.JACKNIFE_OPTION ) ) {
                try {
                    jacknife_resamplings = cla.getOptionValueAsInt( surfacing.JACKNIFE_OPTION );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "illegal format for number of resamplings" );
                }
                if ( jacknife_resamplings < 2 ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "attempt to use less than 2 resamplings" );
                }
            }
            if ( cla.isOptionSet( surfacing.JACKNIFE_RATIO_OPTION )
                    && cla.isOptionHasAValue( surfacing.JACKNIFE_RATIO_OPTION ) ) {
                try {
                    jacknife_ratio = cla.getOptionValueAsDouble( surfacing.JACKNIFE_RATIO_OPTION );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "illegal format for jacknife ratio" );
                }
                if ( ( jacknife_ratio <= 0.0 ) || ( jacknife_ratio >= 1.0 ) ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "attempt to use illegal value for jacknife ratio: "
                            + jacknife_ratio );
                }
            }
            if ( cla.isOptionSet( surfacing.JACKNIFE_RANDOM_SEED_OPTION )
                    && cla.isOptionHasAValue( surfacing.JACKNIFE_RANDOM_SEED_OPTION ) ) {
                try {
                    random_seed = cla.getOptionValueAsLong( surfacing.JACKNIFE_RANDOM_SEED_OPTION );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "illegal format for random generator seed" );
                }
            }
        }
        //        boolean infer_species_trees = false;
        //        if ( cla.isOptionSet( surfacing.INFER_SPECIES_TREES_OPTION ) ) {
        //            if ( ( output_file == null ) || ( number_of_genomes < 3 )
        //                    || ForesterUtil.isEmpty( automated_pairwise_comparison_suffix ) ) {
        //                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot infer species trees (-"
        //                        + surfacing.INFER_SPECIES_TREES_OPTION + " without pairwise analyses ("
        //                        + surfacing.PAIRWISE_DOMAIN_COMPARISONS_OPTION
        //                        + "=<suffix for pairwise comparison output files>)" );
        //            }
        //            infer_species_trees = true;
        //        }
        File[] intree_files = null;
        Phylogeny[] intrees = null;
        if ( cla.isOptionSet( surfacing.INPUT_SPECIES_TREE_OPTION ) ) {
            // TODO FIXME if jacknife.... maybe not
            if ( number_of_genomes < 3 ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot infer gains and losses on input species trees (-"
                        + surfacing.INPUT_SPECIES_TREE_OPTION + " without pairwise analyses ("
                        + surfacing.PAIRWISE_DOMAIN_COMPARISONS_OPTION
                        + "=<suffix for pairwise comparison output files>)" );
            }
            if ( !cla.isOptionValueSet( surfacing.INPUT_SPECIES_TREE_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for input tree: -"
                        + surfacing.INPUT_SPECIES_TREE_OPTION + "=<tree file in phyloXML format>" );
            }
            final String intrees_str = cla.getOptionValue( surfacing.INPUT_SPECIES_TREE_OPTION );
            if ( intrees_str.indexOf( "#" ) > 0 ) {
                final String[] intrees_strs = intrees_str.split( "#" );
                intree_files = new File[ intrees_strs.length ];
                int i = 0;
                for( final String s : intrees_strs ) {
                    intree_files[ i++ ] = new File( s.trim() );
                }
            }
            else {
                intree_files = new File[ 1 ];
                intree_files[ 0 ] = new File( intrees_str );
            }
            intrees = getIntrees( intree_files, number_of_genomes, input_file_properties );
        }
        long random_number_seed_for_fitch_parsimony = 0l;
        boolean radomize_fitch_parsimony = false;
        if ( cla.isOptionSet( surfacing.RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for random number seed: -"
                        + surfacing.RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION + "=<seed>" );
            }
            try {
                random_number_seed_for_fitch_parsimony = cla
                        .getOptionValueAsLong( RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
            }
            radomize_fitch_parsimony = true;
        }
        SortedSet<DomainId> filter = null;
        if ( ( positive_filter_file != null ) || ( negative_filter_file != null )
                || ( negative_domains_filter_file != null ) ) {
            filter = new TreeSet<DomainId>();
            if ( positive_filter_file != null ) {
                processFilter( positive_filter_file, filter );
            }
            else if ( negative_filter_file != null ) {
                processFilter( negative_filter_file, filter );
            }
            else if ( negative_domains_filter_file != null ) {
                processFilter( negative_domains_filter_file, filter );
            }
        }
        Map<DomainId, Set<String>>[] domain_id_to_secondary_features_maps = null;
        File[] secondary_features_map_files = null;
        final File domain_lengths_analysis_outfile = new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file
                + DOMAIN_LENGTHS_ANALYSIS_SUFFIX );
        if ( PERFORM_DOMAIN_LENGTH_ANALYSIS ) {
            SurfacingUtil.checkForOutputFileWriteability( domain_lengths_analysis_outfile );
        }
        if ( cla.isOptionSet( surfacing.SECONDARY_FEATURES_PARSIMONY_MAP_FILE ) ) {
            if ( !cla.isOptionValueSet( surfacing.SECONDARY_FEATURES_PARSIMONY_MAP_FILE ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for secondary features map file: -"
                        + surfacing.SECONDARY_FEATURES_PARSIMONY_MAP_FILE + "=<file>" );
            }
            final String[] secondary_features_map_files_strs = cla
                    .getOptionValue( surfacing.SECONDARY_FEATURES_PARSIMONY_MAP_FILE ).split( "#" );
            secondary_features_map_files = new File[ secondary_features_map_files_strs.length ];
            domain_id_to_secondary_features_maps = new Map[ secondary_features_map_files_strs.length ];
            int i = 0;
            for( final String secondary_features_map_files_str : secondary_features_map_files_strs ) {
                secondary_features_map_files[ i ] = new File( secondary_features_map_files_str );
                final String error = ForesterUtil.isReadableFile( secondary_features_map_files[ i ] );
                if ( !ForesterUtil.isEmpty( error ) ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot read secondary features map file: " + error );
                }
                try {
                    domain_id_to_secondary_features_maps[ i ] = SurfacingUtil
                            .createDomainIdToSecondaryFeaturesMap( secondary_features_map_files[ i ] );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "cannot read secondary features map file: "
                            + e.getMessage() );
                }
                catch ( final Exception e ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "problem with contents of features map file ["
                            + secondary_features_map_files[ i ] + "]: " + e.getMessage() );
                }
                i++;
            }
        }
        if ( out_dir == null ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, "no output directory indicated (-"
                    + surfacing.OUTPUT_DIR_OPTION + "=<dir>)" );
        }
        if ( output_file == null ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, "no name for (main) output file indicated (-"
                    + surfacing.OUTPUT_FILE_OPTION + "=<file>)" );
        }
        if ( ( domain_id_to_go_ids_map == null ) || domain_id_to_go_ids_map.isEmpty() ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME,
                                     "no (acceptable) Pfam to GO id mapping file provided ('pfam2go file') (-"
                                             + surfacing.PFAM_TO_GO_FILE_USE_OPTION + "=<file>)" );
        }
        if ( ( go_id_to_term_map == null ) || go_id_to_term_map.isEmpty() ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME,
                                     "no (acceptable) go id to term mapping file provided ('GO OBO file') (-"
                                             + surfacing.GO_OBO_FILE_USE_OPTION + "=<file>)" );
        }
        boolean display_histograms = false;
        if ( cla.isOptionSet( surfacing.DISPLAY_M_HISTOGRAMS_OPTION ) ) {
            display_histograms = true;
        }
        System.out.println( "Output directory            : " + out_dir );
        if ( input_file_names_from_file != null ) {
            System.out.println( "Input files names from      : " + input_files_file + " ["
                    + input_file_names_from_file.length + " input files]" );
            html_desc.append( "<tr><td>Input files names from:</td><td>" + input_files_file + " ["
                    + input_file_names_from_file.length + " input files]</td></tr>" + nl );
        }
        if ( positive_filter_file != null ) {
            final int filter_size = filter.size();
            System.out.println( "Positive protein filter     : " + positive_filter_file + " [" + filter_size
                    + " domain ids]" );
            html_desc.append( "<tr><td>Positive protein filter:</td><td>" + positive_filter_file + " [" + filter_size
                    + " domain ids]</td></tr>" + nl );
        }
        if ( negative_filter_file != null ) {
            final int filter_size = filter.size();
            System.out.println( "Negative protein filter     : " + negative_filter_file + " [" + filter_size
                    + " domain ids]" );
            html_desc.append( "<tr><td>Negative protein filter:</td><td>" + negative_filter_file + " [" + filter_size
                    + " domain ids]</td></tr>" + nl );
        }
        if ( negative_domains_filter_file != null ) {
            final int filter_size = filter.size();
            System.out.println( "Negative domain filter      : " + negative_domains_filter_file + " [" + filter_size
                    + " domain ids]" );
            html_desc.append( "<tr><td>Negative domain filter:</td><td>" + negative_domains_filter_file + " ["
                    + filter_size + " domain ids]</td></tr>" + nl );
        }
        if ( plus_minus_analysis_high_copy_base_species.size() > 0 ) {
            String plus0 = "";
            for( final String s : plus_minus_analysis_high_copy_base_species ) {
                plus0 += "+" + s + " ";
            }
            String plus1 = "";
            for( final String s : plus_minus_analysis_high_copy_target_species ) {
                plus1 += "*" + s + " ";
            }
            String minus = "";
            for( final String s : plus_minus_analysis_high_low_copy_species ) {
                minus += "-" + s + " ";
            }
            System.out.println( "Plus-minus analysis         : " + plus1 + "&& " + plus0 + "&& " + minus );
            html_desc.append( "<tr><td>Plus-minus analysis:</td><td>" + plus1 + "&& " + plus0 + "&& " + minus
                    + "</td></tr>" + nl );
        }
        if ( cutoff_scores_file != null ) {
            System.out.println( "Cutoff scores file          : " + cutoff_scores_file );
            html_desc.append( "<tr><td>Cutoff scores file:</td><td>" + cutoff_scores_file + "</td></tr>" + nl );
        }
        if ( e_value_max >= 0.0 ) {
            System.out.println( "E-value maximum (inclusive) : " + e_value_max );
            html_desc.append( "<tr><td>E-value maximum (inclusive):</td><td>" + e_value_max + "</td></tr>" + nl );
        }
        System.out.println( "Ignore DUFs                 : " + ignore_dufs );
        if ( ignore_virus_like_ids ) {
            System.out.println( "Ignore virus like ids       : " + ignore_virus_like_ids );
            html_desc.append( "<tr><td>Ignore virus, phage, transposition related ids:</td><td>"
                    + ignore_virus_like_ids + "</td></tr>" + nl );
        }
        html_desc.append( "<tr><td>Ignore DUFs:</td><td>" + ignore_dufs + "</td></tr>" + nl );
        if ( max_allowed_overlap != surfacing.MAX_ALLOWED_OVERLAP_DEFAULT ) {
            System.out.println( "Max allowed domain overlap  : " + max_allowed_overlap );
            html_desc.append( "<tr><td>Max allowed domain overlap:</td><td>" + max_allowed_overlap + "</td></tr>" + nl );
        }
        if ( no_engulfing_overlaps ) {
            System.out.println( "Ignore engulfed domains     : " + no_engulfing_overlaps );
            html_desc.append( "<tr><td>Ignore (lower confidence) engulfed domains:</td><td>" + no_engulfing_overlaps
                    + "</td></tr>" + nl );
        }
        System.out.println( "Ignore singlet domains      : " + ignore_domains_without_combs_in_all_spec );
        html_desc
                .append( "<tr><td>Ignore singlet domains for domain combination similarity analyses (not for parsimony analyses):</td><td>"
                        + ignore_domains_without_combs_in_all_spec + "</td></tr>" + nl );
        System.out.println( "Ignore species specific doms: " + ignore_species_specific_domains );
        html_desc
                .append( "<tr><td>Ignore species specific domains for domain combination similarity analyses (not for parsimony analyses):</td><td>"
                        + ignore_species_specific_domains + "</td></tr>" + nl );
        System.out.println( "Ignore combination with self: " + ignore_combination_with_same );
        html_desc.append( "<tr><td>Ignore combination with self for domain combination similarity analyses:</td><td>"
                + ignore_combination_with_same + "</td></tr>" + nl );
        ;
        System.out.println( "Consider directedness       : "
                + ( dc_type != BinaryDomainCombination.DomainCombinationType.BASIC ) );
        html_desc.append( "<tr><td>Consider directedness of binary domain combinations:</td><td>"
                + ( dc_type != BinaryDomainCombination.DomainCombinationType.BASIC ) + "</td></tr>" + nl );
        if ( dc_type != BinaryDomainCombination.DomainCombinationType.BASIC ) {
            System.out.println( "Consider adjacency          : "
                    + ( dc_type == BinaryDomainCombination.DomainCombinationType.DIRECTED_ADJACTANT ) );
            html_desc.append( "<tr><td>Consider djacency of binary domain combinations:</td><td>"
                    + ( dc_type == BinaryDomainCombination.DomainCombinationType.DIRECTED_ADJACTANT ) + "</td></tr>"
                    + nl );
        }
        System.out.print( "Domain counts sort order    : " );
        switch ( dc_sort_order ) {
            case ALPHABETICAL_KEY_ID:
                System.out.println( "alphabetical" );
                break;
            case KEY_DOMAIN_COUNT:
                System.out.println( "domain count" );
                break;
            case KEY_DOMAIN_PROTEINS_COUNT:
                System.out.println( "domain proteins count" );
                break;
            case COMBINATIONS_COUNT:
                System.out.println( "domain combinations count" );
                break;
            default:
                ForesterUtil.unexpectedFatalError( surfacing.PRG_NAME, "unknown value for dc sort order" );
        }
        if ( domain_id_to_go_ids_map != null ) {
            System.out.println( "Pfam to GO mappings from    : " + pfam_to_go_file + " [" + domain_id_to_go_ids_count
                    + " mappings]" );
            html_desc.append( "<tr><td>Pfam to GO mappings from:</td><td>" + pfam_to_go_file + " ["
                    + domain_id_to_go_ids_count + " mappings]" + "</td></tr>" + nl );
        }
        if ( go_terms != null ) {
            System.out.println( "GO terms from               : " + go_obo_file + " [" + go_terms.size() + " terms]" );
            html_desc.append( "<tr><td>GO terms from:</td><td>" + go_obo_file + " [" + go_terms.size() + " terms]"
                    + "</td></tr>" + nl );
        }
        if ( go_namespace_limit != null ) {
            System.out.println( "Limit GO terms to           : " + go_namespace_limit.toString() );
            html_desc.append( "<tr><td>Limit GO terms to</td><td>" + go_namespace_limit + "</td></tr>" + nl );
        }
        if ( perform_pwc ) {
            System.out.println( "Suffix for PWC files        : " + automated_pairwise_comparison_suffix );
            html_desc.append( "<tr><td>Suffix for PWC files</td><td>" + automated_pairwise_comparison_suffix
                    + "</td></tr>" + nl );
        }
        if ( out_dir != null ) {
            System.out.println( "Output directory            : " + out_dir );
        }
        if ( query_domain_ids != null ) {
            System.out.println( "Query domains (ordered)     : " + query_domain_ids );
            html_desc.append( "<tr><td></td><td>" + query_domain_ids + "</td></tr>" + nl );
        }
        System.out.println( "Write similarities to       : " + output_file );
        System.out.print( "  Scoring method            : " );
        html_desc.append( "<tr><td>Scoring method:</td><td>" );
        switch ( scoring ) {
            case COMBINATIONS:
                System.out.println( "domain combinations based" );
                html_desc.append( "domain combinations based" + "</td></tr>" + nl );
                break;
            case DOMAINS:
                System.out.println( "domain counts based" );
                html_desc.append( "domain counts based" + "</td></tr>" + nl );
                break;
            case PROTEINS:
                System.out.println( "domain proteins counts based" );
                html_desc.append( "domain proteins counts based" + "</td></tr>" + nl );
                break;
            default:
                ForesterUtil.unexpectedFatalError( surfacing.PRG_NAME, "unknown value for sorting for scoring" );
        }
        System.out.print( "  Sort by                   : " );
        html_desc.append( "<tr><td>Sort by:</td><td>" );
        switch ( domain_similarity_sort_field ) {
            case MIN:
                System.out.print( "score minimum" );
                html_desc.append( "score minimum" );
                break;
            case MAX:
                System.out.print( "score maximum" );
                html_desc.append( "score maximum" );
                break;
            case MEAN:
                System.out.print( "score mean" );
                html_desc.append( "score mean" );
                break;
            case SD:
                System.out.print( "score standard deviation" );
                html_desc.append( "score standard deviation" );
                break;
            case SPECIES_COUNT:
                System.out.print( "species number" );
                html_desc.append( "species number" );
                break;
            case DOMAIN_ID:
                System.out.print( "alphabetical domain identifier" );
                html_desc.append( "alphabetical domain identifier" );
                break;
            case MAX_DIFFERENCE:
                System.out.print( "(maximal) difference" );
                html_desc.append( "(maximal) difference" );
                break;
            case ABS_MAX_COUNTS_DIFFERENCE:
                System.out.print( "absolute (maximal) counts difference" );
                html_desc.append( "absolute (maximal) counts difference" );
                break;
            case MAX_COUNTS_DIFFERENCE:
                System.out.print( "(maximal) counts difference" );
                html_desc.append( "(maximal) counts  difference" );
                break;
            default:
                ForesterUtil.unexpectedFatalError( surfacing.PRG_NAME, "unknown value for sorting for similarities" );
        }
        if ( sort_by_species_count_first ) {
            System.out.println( " (sort by species count first)" );
            html_desc.append( " (sort by species count first)" );
        }
        else {
            System.out.println();
        }
        html_desc.append( "</td></tr>" + nl );
        System.out.print( "  Detailedness              : " );
        switch ( detailedness ) {
            case BASIC:
                System.out.println( "basic" );
                break;
            case LIST_COMBINING_DOMAIN_FOR_EACH_SPECIES:
                System.out.println( "list combining domains for each species" );
                break;
            case PUNCTILIOUS:
                System.out.println( "punctilious" );
                break;
            default:
                ForesterUtil.unexpectedFatalError( surfacing.PRG_NAME, "unknown value for sorting for detailedness" );
        }
        System.out.print( "  Print option              : " );
        switch ( domain_similarity_print_option ) {
            case HTML:
                System.out.println( "HTML" );
                break;
            case SIMPLE_TAB_DELIMITED:
                System.out.println( "simple tab delimited" );
                break;
            default:
                ForesterUtil.unexpectedFatalError( surfacing.PRG_NAME, "unknown value for print option" );
        }
        System.out.print( "  Species matrix            : " + species_matrix );
        System.out.println();
        if ( perform_pwc ) {
            System.out.println( "Pairwise comparisons: " );
            html_desc.append( "<tr><td>Pairwise comparisons:</td><td></td></tr>" );
            System.out.print( "  Sort by                   : " );
            html_desc.append( "<tr><td>Sort by:</td><td>" );
            switch ( domain_similarity_sort_field_for_automated_pwc ) {
                case MEAN:
                    System.out.print( "score mean" );
                    html_desc.append( "score mean" );
                    break;
                case DOMAIN_ID:
                    System.out.print( "alphabetical domain identifier" );
                    html_desc.append( "alphabetical domain identifier" );
                    break;
                case MAX_DIFFERENCE:
                    System.out.print( "difference" );
                    html_desc.append( "difference" );
                    break;
                case ABS_MAX_COUNTS_DIFFERENCE:
                    System.out.print( "absolute counts difference" );
                    html_desc.append( "absolute counts difference" );
                    break;
                case MAX_COUNTS_DIFFERENCE:
                    System.out.print( "counts difference" );
                    html_desc.append( "counts difference" );
                    break;
                default:
                    ForesterUtil
                            .unexpectedFatalError( surfacing.PRG_NAME, "unknown value for sorting for similarities" );
            }
            System.out.println();
            html_desc.append( "</td></tr>" + nl );
            if ( jacknifed_distances ) {
                html_desc.append( "<tr><td>Jacknife:</td><td>" + jacknife_resamplings + " resamplings</td></tr>" + nl );
                html_desc.append( "<tr><td>Jacknife ratio:</td><td>" + ForesterUtil.round( jacknife_ratio, 2 )
                        + "</td></tr>" + nl );
                html_desc.append( "<tr><td>Jacknife random number seed:</td><td>" + random_seed + "</td></tr>" + nl );
                System.out.println( "  Jacknife                  : " + jacknife_resamplings + " resamplings" );
                System.out.println( "    Ratio                   : " + ForesterUtil.round( jacknife_ratio, 2 ) );
                System.out.println( "    Random number seed      : " + random_seed );
            }
            //                if ( infer_species_trees ) {
            //                    html_desc.append( "<tr><td>Infer species trees:</td><td>true</td></tr>" + nl );
            //                    System.out.println( "  Infer species trees       : true" );
            //                }
            if ( ( intrees != null ) && ( intrees.length > 0 ) ) {
                for( final File intree_file : intree_files ) {
                    html_desc.append( "<tr><td>Intree for gain/loss parsimony analysis:</td><td>" + intree_file
                            + "</td></tr>" + nl );
                    System.out.println( "  Intree for gain/loss pars.: " + intree_file );
                }
            }
            if ( radomize_fitch_parsimony ) {
                html_desc.append( "<tr><td>    Random number seed for Fitch parsimony analysis:</td><td>"
                        + random_number_seed_for_fitch_parsimony + "</td></tr>" + nl );
                System.out.println( "    Random number seed      : " + random_number_seed_for_fitch_parsimony );
            }
            if ( ( domain_id_to_secondary_features_maps != null ) && ( domain_id_to_secondary_features_maps.length > 0 ) ) {
                for( int i = 0; i < secondary_features_map_files.length; i++ ) {
                    html_desc.append( "<tr><td>Secondary features map file:</td><td>"
                            + secondary_features_map_files[ i ] + "</td></tr>" + nl );
                    System.out.println( "Secondary features map file : " + secondary_features_map_files[ i ]
                            + " [mappings for " + domain_id_to_secondary_features_maps[ i ].size() + " domain ids]" );
                    if ( VERBOSE ) {
                        System.out.println();
                        System.out.println( "Domain ids to secondary features map:" );
                        for( final DomainId domain_id : domain_id_to_secondary_features_maps[ i ].keySet() ) {
                            System.out.print( domain_id.getId() );
                            System.out.print( " => " );
                            for( final String sec : domain_id_to_secondary_features_maps[ i ].get( domain_id ) ) {
                                System.out.print( sec );
                                System.out.print( " " );
                            }
                            System.out.println();
                        }
                    }
                }
            }
        } // if ( perform_pwc ) {
        System.out.println();
        html_desc.append( "<tr><td>Command line:</td><td>" + cla.getCommandLineArgsAsString() + "</td></tr>" + nl );
        System.out.println( "Command line                : " + cla.getCommandLineArgsAsString() );
        BufferedWriter[] query_domains_writer_ary = null;
        List<DomainId>[] query_domain_ids_array = null;
        if ( query_domain_ids != null ) {
            final String[] query_domain_ids_str_array = query_domain_ids.split( "#" );
            query_domain_ids_array = new ArrayList[ query_domain_ids_str_array.length ];
            query_domains_writer_ary = new BufferedWriter[ query_domain_ids_str_array.length ];
            for( int i = 0; i < query_domain_ids_str_array.length; i++ ) {
                String query_domain_ids_str = query_domain_ids_str_array[ i ];
                final String[] query_domain_ids_str_ary = query_domain_ids_str.split( "~" );
                final List<DomainId> query = new ArrayList<DomainId>();
                for( final String element : query_domain_ids_str_ary ) {
                    query.add( new DomainId( element ) );
                }
                query_domain_ids_array[ i ] = query;
                query_domain_ids_str = query_domain_ids_str.replace( '~', '_' );
                String protein_names_writer_str = query_domain_ids_str + surfacing.SEQ_EXTRACT_SUFFIX;
                if ( out_dir != null ) {
                    protein_names_writer_str = out_dir + ForesterUtil.FILE_SEPARATOR + protein_names_writer_str;
                }
                try {
                    query_domains_writer_ary[ i ] = new BufferedWriter( new FileWriter( protein_names_writer_str ) );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, "Could not open [" + protein_names_writer_str + "]: "
                            + e.getLocalizedMessage() );
                }
            }
        }
        SortedMap<Species, List<Protein>> protein_lists_per_species = null; //This will only be created if neede.
        boolean need_protein_lists_per_species = false;
        if ( ( plus_minus_analysis_high_copy_base_species.size() > 0 ) || output_protein_lists_for_all_domains ) {
            need_protein_lists_per_species = true;
        }
        if ( need_protein_lists_per_species ) {
            protein_lists_per_species = new TreeMap<Species, List<Protein>>();
        }
        final List<GenomeWideCombinableDomains> gwcd_list = new ArrayList<GenomeWideCombinableDomains>( number_of_genomes );
        final SortedSet<DomainId> all_domains_encountered = new TreeSet<DomainId>();
        final SortedSet<BinaryDomainCombination> all_bin_domain_combinations_encountered = new TreeSet<BinaryDomainCombination>();
        List<BinaryDomainCombination> all_bin_domain_combinations_gained_fitch = null;
        List<BinaryDomainCombination> all_bin_domain_combinations_lost_fitch = null;
        if ( ( intrees != null ) && ( intrees.length == 1 ) ) {
            all_bin_domain_combinations_gained_fitch = new ArrayList<BinaryDomainCombination>();
            all_bin_domain_combinations_lost_fitch = new ArrayList<BinaryDomainCombination>();
        }
        final DomainLengthsTable domain_lengths_table = new DomainLengthsTable();
        final File per_genome_domain_promiscuity_statistics_file = new File( out_dir + ForesterUtil.FILE_SEPARATOR
                + output_file + D_PROMISCUITY_FILE_SUFFIX );
        BufferedWriter per_genome_domain_promiscuity_statistics_writer = null;
        try {
            per_genome_domain_promiscuity_statistics_writer = new BufferedWriter( new FileWriter( per_genome_domain_promiscuity_statistics_file ) );
            per_genome_domain_promiscuity_statistics_writer.write( "Species:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "Mean:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "SD:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "Median:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "Min:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "Max:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "N:\t" );
            per_genome_domain_promiscuity_statistics_writer.write( "Max Promiscuous Domains:"
                    + ForesterUtil.LINE_SEPARATOR );
        }
        catch ( final IOException e2 ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e2.getMessage() );
        }
        final File log_file = new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file + LOG_FILE_SUFFIX );
        BufferedWriter log_writer = null;
        try {
            log_writer = new BufferedWriter( new FileWriter( log_file ) );
        }
        catch ( final IOException e2 ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e2.getMessage() );
        }
        for( int i = 0; i < number_of_genomes; ++i ) {
            System.out.println();
            System.out.println( ( i + 1 ) + "/" + number_of_genomes );
            log( ( i + 1 ) + "/" + number_of_genomes, log_writer );
            System.out.println( "Processing                                     : " + input_file_properties[ i ][ 0 ] );
            log( "Genome                                         : " + input_file_properties[ i ][ 0 ], log_writer );
            HmmscanPerDomainTableParser parser = null;
            INDIVIDUAL_SCORE_CUTOFF ind_score_cutoff = INDIVIDUAL_SCORE_CUTOFF.NONE;
            if ( individual_score_cutoffs != null ) {
                ind_score_cutoff = INDIVIDUAL_SCORE_CUTOFF_DEFAULT;
            }
            if ( ( positive_filter_file != null ) || ( negative_filter_file != null )
                    || ( negative_domains_filter_file != null ) ) {
                HmmscanPerDomainTableParser.FilterType filter_type = HmmscanPerDomainTableParser.FilterType.NONE;
                if ( positive_filter_file != null ) {
                    filter_type = HmmscanPerDomainTableParser.FilterType.POSITIVE_PROTEIN;
                }
                else if ( negative_filter_file != null ) {
                    filter_type = HmmscanPerDomainTableParser.FilterType.NEGATIVE_PROTEIN;
                }
                else if ( negative_domains_filter_file != null ) {
                    filter_type = HmmscanPerDomainTableParser.FilterType.NEGATIVE_DOMAIN;
                }
                parser = new HmmscanPerDomainTableParser( new File( input_file_properties[ i ][ 0 ] ),
                                                          input_file_properties[ i ][ 1 ],
                                                          filter,
                                                          filter_type,
                                                          ind_score_cutoff );
            }
            else {
                parser = new HmmscanPerDomainTableParser( new File( input_file_properties[ i ][ 0 ] ),
                                                          input_file_properties[ i ][ 1 ],
                                                          ind_score_cutoff );
            }
            if ( e_value_max >= 0.0 ) {
                parser.setEValueMaximum( e_value_max );
            }
            parser.setIgnoreDufs( ignore_dufs );
            parser.setIgnoreVirusLikeIds( ignore_virus_like_ids );
            parser.setIgnoreEngulfedDomains( no_engulfing_overlaps );
            if ( max_allowed_overlap != surfacing.MAX_ALLOWED_OVERLAP_DEFAULT ) {
                parser.setMaxAllowedOverlap( max_allowed_overlap );
            }
            parser
                    .setReturnType( HmmscanPerDomainTableParser.ReturnType.UNORDERED_PROTEIN_DOMAIN_COLLECTION_PER_PROTEIN );
            if ( individual_score_cutoffs != null ) {
                parser.setIndividualScoreCutoffs( individual_score_cutoffs );
            }
            List<Protein> protein_list = null;
            try {
                protein_list = parser.parse();
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
            }
            catch ( final Exception e ) {
                ForesterUtil.unexpectedFatalError( surfacing.PRG_NAME, e.getMessage(), e );
            }
            if ( VERBOSE ) {
                System.out.println( "Domains ignored due to negative domain filter: " );
                ForesterUtil.printCountingMap( parser.getDomainsIgnoredDueToNegativeDomainFilterCountsMap() );
                System.out.println( "Domains ignored due to virus like id: " );
                ForesterUtil.printCountingMap( parser.getDomainsIgnoredDueToVirusLikeIdCountsMap() );
            }
            System.out.println( "Number of proteins encountered                 : " + parser.getProteinsEncountered() );
            log( "Number of proteins encountered                 : " + parser.getProteinsEncountered(), log_writer );
            System.out.println( "Number of proteins stored                      : " + protein_list.size() );
            log( "Number of proteins stored                      : " + protein_list.size(), log_writer );
            System.out.println( "Domains encountered                            : " + parser.getDomainsEncountered() );
            log( "Domains encountered                            : " + parser.getDomainsEncountered(), log_writer );
            System.out.println( "Domains stored                                 : " + parser.getDomainsStored() );
            log( "Domains stored                                 : " + parser.getDomainsStored(), log_writer );
            System.out.println( "Distinct domains stored                        : "
                    + parser.getDomainsStoredSet().size() );
            log( "Distinct domains stored                        : " + parser.getDomainsStoredSet().size(), log_writer );
            System.out.println( "Domains ignored due to individual score cutoffs: "
                    + parser.getDomainsIgnoredDueToIndividualScoreCutoff() );
            log( "Domains ignored due to individual score cutoffs: "
                    + parser.getDomainsIgnoredDueToIndividualScoreCutoff(), log_writer );
            System.out.println( "Domains ignored due to E-value                 : "
                    + parser.getDomainsIgnoredDueToEval() );
            log( "Domains ignored due to E-value                 : " + parser.getDomainsIgnoredDueToEval(), log_writer );
            System.out.println( "Domains ignored due to DUF designation         : "
                    + parser.getDomainsIgnoredDueToDuf() );
            log( "Domains ignored due to DUF designation         : " + parser.getDomainsIgnoredDueToDuf(), log_writer );
            if ( ignore_virus_like_ids ) {
                System.out.println( "Domains ignored due virus like ids             : "
                        + parser.getDomainsIgnoredDueToVirusLikeIds() );
                log( "Domains ignored due virus like ids             : " + parser.getDomainsIgnoredDueToVirusLikeIds(),
                     log_writer );
            }
            System.out.println( "Domains ignored due negative domain filter     : "
                    + parser.getDomainsIgnoredDueToNegativeDomainFilter() );
            log( "Domains ignored due negative domain filter     : "
                    + parser.getDomainsIgnoredDueToNegativeDomainFilter(), log_writer );
            System.out.println( "Domains ignored due to overlap                 : "
                    + parser.getDomainsIgnoredDueToOverlap() );
            log( "Domains ignored due to overlap                 : " + parser.getDomainsIgnoredDueToOverlap(),
                 log_writer );
            if ( negative_filter_file != null ) {
                System.out.println( "Proteins ignored due to negative filter        : "
                        + parser.getProteinsIgnoredDueToFilter() );
                log( "Proteins ignored due to negative filter        : " + parser.getProteinsIgnoredDueToFilter(),
                     log_writer );
            }
            if ( positive_filter_file != null ) {
                System.out.println( "Proteins ignored due to positive filter        : "
                        + parser.getProteinsIgnoredDueToFilter() );
                log( "Proteins ignored due to positive filter        : " + parser.getProteinsIgnoredDueToFilter(),
                     log_writer );
            }
            System.out.println( "Time for processing                            : " + parser.getTime() + "ms" );
            log( "", log_writer );
            html_desc.append( "<tr><td>" + input_file_properties[ i ][ 0 ] + " [species: "
                    + input_file_properties[ i ][ 1 ] + "]" + ":</td><td>domains analyzed: "
                    + parser.getDomainsStored() + "; domains ignored: [ind score cutoffs: "
                    + parser.getDomainsIgnoredDueToIndividualScoreCutoff() + "] [E-value cutoff: "
                    + parser.getDomainsIgnoredDueToEval() + "] [DUF: " + parser.getDomainsIgnoredDueToDuf()
                    + "] [virus like ids: " + parser.getDomainsIgnoredDueToVirusLikeIds()
                    + "] [negative domain filter: " + parser.getDomainsIgnoredDueToNegativeDomainFilter()
                    + "] [overlap: " + parser.getDomainsIgnoredDueToOverlap() + "]" );
            if ( negative_filter_file != null ) {
                html_desc.append( "; proteins ignored due to negative filter: "
                        + parser.getProteinsIgnoredDueToFilter() );
            }
            if ( positive_filter_file != null ) {
                html_desc.append( "; proteins ignored due to positive filter: "
                        + parser.getProteinsIgnoredDueToFilter() );
            }
            html_desc.append( "</td></tr>" + nl );
            // domain_partner_counts_array[ i ] =
            // Methods.getDomainPartnerCounts( protein_domain_collections_array[
            // i ],
            // false, input_file_properties[ i ][ 1 ] );
            gwcd_list.add( BasicGenomeWideCombinableDomains
                    .createInstance( protein_list,
                                     ignore_combination_with_same,
                                     new BasicSpecies( input_file_properties[ i ][ 1 ] ),
                                     domain_id_to_go_ids_map,
                                     dc_type ) );
            domain_lengths_table.addLengths( protein_list );
            if ( gwcd_list.get( i ).getSize() > 0 ) {
                SurfacingUtil.writeDomainCombinationsCountsFile( input_file_properties,
                                                                 out_dir,
                                                                 per_genome_domain_promiscuity_statistics_writer,
                                                                 gwcd_list.get( i ),
                                                                 i,
                                                                 dc_sort_order );
                if ( output_binary_domain_combinationsfor_graph_analysis ) {
                    SurfacingUtil.writeBinaryDomainCombinationsFileForGraphAnalysis( input_file_properties,
                                                                                     out_dir,
                                                                                     gwcd_list.get( i ),
                                                                                     i,
                                                                                     dc_sort_order );
                }
                SurfacingUtil.addAllDomainIdsToSet( gwcd_list.get( i ), all_domains_encountered );
                SurfacingUtil.addAllBinaryDomainCombinationToSet( gwcd_list.get( i ),
                                                                  all_bin_domain_combinations_encountered );
            }
            if ( query_domains_writer_ary != null ) {
                for( int j = 0; j < query_domain_ids_array.length; j++ ) {
                    try {
                        SurfacingUtil.extractProteinNames( protein_list,
                                                           query_domain_ids_array[ j ],
                                                           query_domains_writer_ary[ j ],
                                                           "\t" );
                        query_domains_writer_ary[ j ].flush();
                    }
                    catch ( final IOException e ) {
                        e.printStackTrace();
                    }
                }
            }
            if ( need_protein_lists_per_species ) {
                protein_lists_per_species.put( new BasicSpecies( input_file_properties[ i ][ 1 ] ), protein_list );
            }
            try {
                log_writer.flush();
            }
            catch ( final IOException e2 ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, e2.getLocalizedMessage() );
            }
            System.gc();
        } // for( int i = 0; i < number_of_hmmpfam_files_to_analyze; ++i ) {
        try {
            per_genome_domain_promiscuity_statistics_writer.flush();
            per_genome_domain_promiscuity_statistics_writer.close();
            log_writer.close();
        }
        catch ( final IOException e2 ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e2.getLocalizedMessage() );
        }
        ForesterUtil.programMessage( PRG_NAME, "Wrote domain promiscuities to: "
                + per_genome_domain_promiscuity_statistics_file );
        if ( query_domains_writer_ary != null ) {
            for( int j = 0; j < query_domain_ids_array.length; j++ ) {
                try {
                    query_domains_writer_ary[ j ].close();
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( surfacing.PRG_NAME, e.toString() );
                }
            }
        }
        if ( PERFORM_DOMAIN_LENGTH_ANALYSIS ) {
            try {
                SurfacingUtil.executeDomainLengthAnalysis( input_file_properties,
                                                           number_of_genomes,
                                                           domain_lengths_table,
                                                           domain_lengths_analysis_outfile );
            }
            catch ( final IOException e1 ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, e1.toString() );
            }
            System.out.println();
            ForesterUtil.programMessage( PRG_NAME, "Wrote domain length data to: " + domain_lengths_analysis_outfile );
            System.out.println();
        }
        final long analysis_start_time = new Date().getTime();
        PairwiseDomainSimilarityCalculator pw_calc = null;
        // double[] values_for_all_scores_histogram = null;
        final DomainSimilarityCalculator calc = new BasicDomainSimilarityCalculator( domain_similarity_sort_field,
                                                                                     sort_by_species_count_first,
                                                                                     number_of_genomes == 2 );
        switch ( scoring ) {
            case COMBINATIONS:
                pw_calc = new CombinationsBasedPairwiseDomainSimilarityCalculator();
                break;
            case DOMAINS:
                pw_calc = new DomainCountsBasedPairwiseSimilarityCalculator();
                break;
            case PROTEINS:
                pw_calc = new ProteinCountsBasedPairwiseDomainSimilarityCalculator();
                break;
            default:
                ForesterUtil.unexpectedFatalError( surfacing.PRG_NAME, "unknown value for sorting for scoring" );
        }
        DomainSimilarityCalculator.GoAnnotationOutput go_annotation_output = DomainSimilarityCalculator.GoAnnotationOutput.NONE;
        if ( domain_id_to_go_ids_map != null ) {
            go_annotation_output = DomainSimilarityCalculator.GoAnnotationOutput.ALL;
        }
        final SortedSet<DomainSimilarity> similarities = calc
                .calculateSimilarities( pw_calc,
                                        gwcd_list,
                                        ignore_domains_without_combs_in_all_spec,
                                        ignore_species_specific_domains );
        SurfacingUtil.decoratePrintableDomainSimilarities( similarities,
                                                           detailedness,
                                                           go_annotation_output,
                                                           go_id_to_term_map,
                                                           go_namespace_limit );
        DescriptiveStatistics pw_stats = null;
        try {
            String my_outfile = output_file.toString();
            if ( !my_outfile.endsWith( ".html" ) ) {
                my_outfile += ".html";
            }
            final Writer writer = new BufferedWriter( new FileWriter( out_dir == null ? my_outfile : out_dir
                    + ForesterUtil.FILE_SEPARATOR + my_outfile ) );
            List<Species> species_order = null;
            if ( species_matrix ) {
                species_order = new ArrayList<Species>();
                for( int i = 0; i < number_of_genomes; i++ ) {
                    species_order.add( new BasicSpecies( input_file_properties[ i ][ 1 ] ) );
                }
            }
            html_desc.append( "<tr><td>Sum of all distinct binary combinations:</td><td>"
                    + all_bin_domain_combinations_encountered.size() + "</td></tr>" + nl );
            html_desc.append( "<tr><td>Sum of all distinct domains:</td><td>" + all_domains_encountered.size()
                    + "</td></tr>" + nl );
            html_desc.append( "<tr><td>Analysis date/time:</td><td>"
                    + new java.text.SimpleDateFormat( "yyyy.MM.dd HH:mm:ss" ).format( new java.util.Date() )
                    + "</td></tr>" + nl );
            html_desc.append( "</table>" + nl );
            pw_stats = SurfacingUtil
                    .writeDomainSimilaritiesToFile( html_desc,
                                                    new StringBuilder( number_of_genomes + " genomes" ),
                                                    writer,
                                                    similarities,
                                                    number_of_genomes == 2,
                                                    species_order,
                                                    domain_similarity_print_option,
                                                    domain_similarity_sort_field,
                                                    scoring,
                                                    true );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote main output (includes domain similarities) to: \""
                    + ( out_dir == null ? my_outfile : out_dir + ForesterUtil.FILE_SEPARATOR + my_outfile ) + "\"" );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, "Failed to write similarites to: \"" + output_file + "\" ["
                    + e.getMessage() + "]" );
        }
        System.out.println();
        // values_for_all_scores_histogram = pw_stats.getDataAsDoubleArray();
        final Species[] species = new Species[ number_of_genomes ];
        for( int i = 0; i < number_of_genomes; ++i ) {
            species[ i ] = new BasicSpecies( input_file_properties[ i ][ 1 ] );
        }
        List<Phylogeny> inferred_trees = null;
        if ( ( number_of_genomes > 2 ) && perform_pwc ) {
            final PairwiseGenomeComparator pwgc = new PairwiseGenomeComparator();
            pwgc.performPairwiseComparisons( html_desc,
                                             sort_by_species_count_first,
                                             detailedness,
                                             ignore_domains_without_combs_in_all_spec,
                                             ignore_species_specific_domains,
                                             domain_similarity_sort_field_for_automated_pwc,
                                             domain_similarity_print_option,
                                             scoring,
                                             domain_id_to_go_ids_map,
                                             go_id_to_term_map,
                                             go_namespace_limit,
                                             species,
                                             number_of_genomes,
                                             gwcd_list,
                                             pw_calc,
                                             automated_pairwise_comparison_suffix,
                                             true,
                                             surfacing.PAIRWISE_DOMAIN_COMPARISONS_PREFIX,
                                             surfacing.PRG_NAME,
                                             display_histograms,
                                             out_dir,
                                             write_pwc_files );
            String matrix_output_file = new String( output_file.toString() );
            if ( matrix_output_file.indexOf( '.' ) > 1 ) {
                matrix_output_file = matrix_output_file.substring( 0, matrix_output_file.indexOf( '.' ) );
            }
            if ( out_dir != null ) {
                matrix_output_file = out_dir + ForesterUtil.FILE_SEPARATOR + matrix_output_file;
                output_file = new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file );
            }
            SurfacingUtil.writeMatrixToFile( new File( matrix_output_file
                    + surfacing.MATRIX_MEAN_SCORE_BASED_GENOME_DISTANCE_SUFFIX ), pwgc.getDomainDistanceScoresMeans() );
            SurfacingUtil.writeMatrixToFile( new File( matrix_output_file
                    + surfacing.MATRIX_SHARED_BIN_COMBINATIONS_BASED_GENOME_DISTANCE_SUFFIX ), pwgc
                    .getSharedBinaryCombinationsBasedDistances() );
            SurfacingUtil.writeMatrixToFile( new File( matrix_output_file
                    + surfacing.MATRIX_SHARED_DOMAINS_BASED_GENOME_DISTANCE_SUFFIX ), pwgc
                    .getSharedDomainsBasedDistances() );
            final Phylogeny nj_gd = SurfacingUtil.createNjTreeBasedOnMatrixToFile( new File( matrix_output_file
                    + surfacing.NJ_TREE_MEAN_SCORE_BASED_GENOME_DISTANCE_SUFFIX ), pwgc.getDomainDistanceScoresMeans()
                    .get( 0 ) );
            final Phylogeny nj_bc = SurfacingUtil.createNjTreeBasedOnMatrixToFile( new File( matrix_output_file
                    + surfacing.NJ_TREE_SHARED_BIN_COMBINATIONS_BASED_GENOME_DISTANCE_SUFFIX ), pwgc
                    .getSharedBinaryCombinationsBasedDistances().get( 0 ) );
            final Phylogeny nj_d = SurfacingUtil.createNjTreeBasedOnMatrixToFile( new File( matrix_output_file
                    + surfacing.NJ_TREE_SHARED_DOMAINS_BASED_GENOME_DISTANCE_SUFFIX ), pwgc
                    .getSharedDomainsBasedDistances().get( 0 ) );
            inferred_trees = new ArrayList<Phylogeny>();
            inferred_trees.add( nj_gd );
            inferred_trees.add( nj_bc );
            inferred_trees.add( nj_d );
            // final List<HistogramData> histogram_datas = pwgc.getHistogramDatas();
            //            if ( infer_species_trees ) {
            //                inferred_trees = new ArrayList<Phylogeny>();
            //                final List<Phylogeny> inferred_trees_bc =  inferSpeciesTrees( new File( output_file + INFERRED_SBC_BASED_NJ_SPECIES_TREE_SUFFIX ), pwgc
            //                        .getSharedBinaryCombinationsBasedDistances() );
            //                final List<Phylogeny> inferred_trees_d =  inferSpeciesTrees( new File( output_file + INFERRED_SD_BASED_NJ_SPECIES_TREE_SUFFIX ), pwgc
            //                        .getSharedDomainsBasedDistances() );
            //                inferred_trees.addAll( inferred_trees_bc );
            //                inferred_trees.addAll( inferred_trees_d );
            //            }
            if ( jacknifed_distances ) {
                pwgc.performPairwiseComparisonsJacknifed( species,
                                                          number_of_genomes,
                                                          gwcd_list,
                                                          true,
                                                          jacknife_resamplings,
                                                          jacknife_ratio,
                                                          random_seed );
                SurfacingUtil.writeMatrixToFile( new File( matrix_output_file + "_"
                        + ForesterUtil.round( jacknife_ratio, 2 ) + "_" + jacknife_resamplings
                        + surfacing.MATRIX_SHARED_BIN_COMBINATIONS_BASED_GENOME_DISTANCE_SUFFIX ), pwgc
                        .getSharedBinaryCombinationsBasedDistances() );
                SurfacingUtil.writeMatrixToFile( new File( matrix_output_file + "_"
                        + ForesterUtil.round( jacknife_ratio, 2 ) + "_" + jacknife_resamplings
                        + surfacing.MATRIX_SHARED_DOMAINS_BASED_GENOME_DISTANCE_SUFFIX ), pwgc
                        .getSharedDomainsBasedDistances() );
                //                if ( infer_species_trees ) {
                //                    inferSpeciesTrees( new File( output_file + "_" + jacknife_resamplings
                //                            + INFERRED_SBC_BASED_NJ_SPECIES_TREE_SUFFIX ), pwgc
                //                            .getSharedBinaryCombinationsBasedDistances() );
                //                    inferSpeciesTrees( new File( output_file + "_" + jacknife_resamplings
                //                            + INFERRED_SD_BASED_NJ_SPECIES_TREE_SUFFIX ), pwgc.getSharedDomainsBasedDistances() );
                //                }
            }
            if ( display_histograms ) {
                //                final List<HistogramData> histogram_datas_all = new ArrayList<HistogramData>();
                //                histogram_datas_all.add( new HistogramData( "all",
                //                                                            values_for_all_scores_histogram,
                //                                                            null,
                //                                                            20 ) );
                //                final HistogramsFrame hf_all = new HistogramsFrame( histogram_datas_all );
                //                final HistogramsFrame hf = new HistogramsFrame( histogram_datas );
                //                hf_all.setVisible( true );
                //                hf.setVisible( true );
            }
        } // if ( ( output_file != null ) && ( number_of_genomes > 2 ) && !isEmpty( automated_pairwise_comparison_suffix ) )
        if ( ( out_dir != null ) && ( !perform_pwc ) ) {
            output_file = new File( out_dir + ForesterUtil.FILE_SEPARATOR + output_file );
        }
        writePresentToNexus( output_file, positive_filter_file, filter, gwcd_list );
        if ( ( ( intrees != null ) && ( intrees.length > 0 ) ) && ( number_of_genomes > 2 ) ) {
            final StringBuilder parameters_sb = createParametersAsString( ignore_dufs,
                                                                          e_value_max,
                                                                          max_allowed_overlap,
                                                                          no_engulfing_overlaps,
                                                                          cutoff_scores_file,
                                                                          dc_type );
            String s = "_";
            if ( radomize_fitch_parsimony ) {
                s += random_number_seed_for_fitch_parsimony + "_";
            }
            int i = 0;
            for( final Phylogeny intree : intrees ) {
                final String outfile_name = ForesterUtil.removeSuffix( output_file.toString() ) + s
                        + ForesterUtil.removeSuffix( intree_files[ i ].toString() );
                final DomainParsimonyCalculator domain_parsimony = DomainParsimonyCalculator.createInstance( intree,
                                                                                                             gwcd_list );
                SurfacingUtil.executeParsimonyAnalysis( random_number_seed_for_fitch_parsimony,
                                                        radomize_fitch_parsimony,
                                                        outfile_name,
                                                        domain_parsimony,
                                                        intree,
                                                        domain_id_to_go_ids_map,
                                                        go_id_to_term_map,
                                                        go_namespace_limit,
                                                        parameters_sb.toString(),
                                                        domain_id_to_secondary_features_maps,
                                                        positive_filter_file == null ? null : filter,
                                                        output_binary_domain_combinationsfor_graph_analysis,
                                                        all_bin_domain_combinations_gained_fitch,
                                                        all_bin_domain_combinations_lost_fitch,
                                                        dc_type );
                // Listing of all domain combinations gained is only done if only one input tree is used. 
                if ( ( domain_id_to_secondary_features_maps != null )
                        && ( domain_id_to_secondary_features_maps.length > 0 ) ) {
                    int j = 0;
                    for( final Map<DomainId, Set<String>> domain_id_to_secondary_features_map : domain_id_to_secondary_features_maps ) {
                        final Map<Species, MappingResults> mapping_results_map = new TreeMap<Species, MappingResults>();
                        final DomainParsimonyCalculator secondary_features_parsimony = DomainParsimonyCalculator
                                .createInstance( intree, gwcd_list, domain_id_to_secondary_features_map );
                        SurfacingUtil
                                .executeParsimonyAnalysisForSecondaryFeatures( outfile_name
                                                                                       + "_"
                                                                                       + secondary_features_map_files[ j++ ],
                                                                               secondary_features_parsimony,
                                                                               intree,
                                                                               parameters_sb.toString(),
                                                                               mapping_results_map );
                        if ( i == 0 ) {
                            System.out.println();
                            System.out.println( "Mapping to secondary features:" );
                            for( final Species spec : mapping_results_map.keySet() ) {
                                final MappingResults mapping_results = mapping_results_map.get( spec );
                                final int total_domains = mapping_results.getSumOfFailures()
                                        + mapping_results.getSumOfSuccesses();
                                System.out.print( spec + ":" );
                                System.out.print( " mapped domains = " + mapping_results.getSumOfSuccesses() );
                                System.out.print( ", not mapped domains = " + mapping_results.getSumOfFailures() );
                                if ( total_domains > 0 ) {
                                    System.out.println( ", mapped ratio = "
                                            + ( 100 * mapping_results.getSumOfSuccesses() / total_domains ) + "%" );
                                }
                                else {
                                    System.out.println( ", mapped ratio = n/a (total domains = 0 )" );
                                }
                            }
                        }
                    }
                }
                i++;
            } // for( final Phylogeny intree : intrees ) {
        }
        if ( plus_minus_analysis_high_copy_base_species.size() > 0 ) {
            executePlusMinusAnalysis( output_file,
                                      plus_minus_analysis_high_copy_base_species,
                                      plus_minus_analysis_high_copy_target_species,
                                      plus_minus_analysis_high_low_copy_species,
                                      gwcd_list,
                                      protein_lists_per_species,
                                      domain_id_to_go_ids_map,
                                      go_id_to_term_map,
                                      plus_minus_analysis_numbers );
        }
        if ( output_protein_lists_for_all_domains ) {
            writeProteinListsForAllSpecies( out_dir, protein_lists_per_species, gwcd_list );
        }
        //        if ( ( intrees != null ) && ( intrees.length > 0 ) && ( inferred_trees != null ) && ( inferred_trees.size() > 0 ) ) {
        //            final StringBuilder parameters_sb = createParametersAsString( ignore_dufs,
        //                                                                          e_value_max,
        //                                                                          max_allowed_overlap,
        //                                                                          no_engulfing_overlaps,
        //                                                                          cutoff_scores_file );
        //            String s = "_";
        //            if ( radomize_fitch_parsimony ) {
        //                s += random_number_seed_for_fitch_parsimony + "_";
        //            }
        //            int i = 0;
        //            for( final Phylogeny inferred_tree : inferred_trees ) {
        //                if ( !inferred_tree.isRooted() ) { 
        //                    intrees[ 0 ].getRoot().getName();
        //                    inferred_tree.r
        //                }
        //                final String outfile_name = ForesterUtil.removeSuffix( inferred_tree.getName() ) + s;
        //                final DomainParsimonyCalculator domain_parsimony = DomainParsimonyCalculator
        //                        .createInstance( inferred_tree, gwcd_list );
        //                SurfacingUtil.executeParsimonyAnalysis( random_number_seed_for_fitch_parsimony,
        //                                                        radomize_fitch_parsimony,
        //                                                        outfile_name,
        //                                                        domain_parsimony,
        //                                                        inferred_tree,
        //                                                        domain_id_to_go_ids_map,
        //                                                        go_id_to_term_map,
        //                                                        go_namespace_limit,
        //                                                        parameters_sb.toString() );
        //                i++;
        //            }
        //        }
        if ( all_bin_domain_combinations_gained_fitch != null ) {
            try {
                executeFitchGainsAnalysis( new File( output_file
                                                   + surfacing.OUTPUT_DOMAIN_COMBINATIONS_GAINED_MORE_THAN_ONCE_ANALYSIS_SUFFIX ),
                                           all_bin_domain_combinations_gained_fitch,
                                           all_domains_encountered.size(),
                                           all_bin_domain_combinations_encountered,
                                           true );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
            }
        }
        if ( all_bin_domain_combinations_lost_fitch != null ) {
            try {
                executeFitchGainsAnalysis( new File( output_file
                                                   + surfacing.OUTPUT_DOMAIN_COMBINATIONS_LOST_MORE_THAN_ONCE_ANALYSIS_SUFFIX ),
                                           all_bin_domain_combinations_lost_fitch,
                                           all_domains_encountered.size(),
                                           all_bin_domain_combinations_encountered,
                                           false );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
            }
        }
        final Runtime rt = java.lang.Runtime.getRuntime();
        final long free_memory = rt.freeMemory() / 1000000;
        final long total_memory = rt.totalMemory() / 1000000;
        System.out.println();
        System.out.println( "Time for analysis : " + ( new Date().getTime() - analysis_start_time ) + "ms" );
        System.out.println( "Total running time: " + ( new Date().getTime() - start_time ) + "ms " );
        System.out.println( "Free memory       : " + free_memory + "MB, total memory: " + total_memory + "MB" );
        System.out.println();
        System.out.println( "If this application is useful to you, please cite:" );
        System.out.println( surfacing.WWW );
        System.out.println();
        ForesterUtil.programMessage( PRG_NAME, "OK" );
        System.out.println();
    }

    private static void preparePhylogenyForParsimonyAnalyses( final Phylogeny intree,
                                                              final String[][] input_file_properties ) {
        final String[] genomes = new String[ input_file_properties.length ];
        for( int i = 0; i < input_file_properties.length; ++i ) {
            if ( intree.getNodes( input_file_properties[ i ][ 1 ] ).size() > 1 ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "node named [" + input_file_properties[ i ][ 1 ]
                        + "] is not unique in input tree " + intree.getName() );
            }
            genomes[ i ] = input_file_properties[ i ][ 1 ];
        }
        PhylogenyMethods.deleteExternalNodesPositiveSelection( genomes, intree );
        for( int i = 0; i < input_file_properties.length; ++i ) {
            try {
                intree.getNode( input_file_properties[ i ][ 1 ] );
            }
            catch ( final IllegalArgumentException e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "node named [" + input_file_properties[ i ][ 1 ]
                        + "] not present/not unique in input tree" );
            }
        }
    }

    // public static StringBuffer stringCombinableDomainsMapToStringBuffer(
    // final SortedMap<String, CombinableDomains> map ) {
    // final StringBuffer sb = new StringBuffer();
    // for( final Iterator<String> iter = map.keySet().iterator();
    // iter.hasNext(); ) {
    // final Object key = iter.next();
    // sb.append( ForesterUtil.pad( new StringBuffer( key.toString() ), 18, ' ',
    // false ) );
    // final CombinableDomains domain_combination = map.get( key );
    // sb.append( ForesterUtil.pad( new StringBuffer( "" +
    // domain_combination.getNumberOfCombiningDomains() ), 8,
    // ' ', false ) );
    // sb.append( domain_combination.toStringBuffer() );
    // sb.append( ForesterUtil.getLineSeparator() );
    // }
    // return sb;
    // }
    private static void printHelp() {
        System.out.println();
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( "% java -Xms256m -Xmx512m -cp forester.jar org.forester.applications." + surfacing.PRG_NAME
                + " [options] <phylogen(y|ies) infile> [external node name 1] [name 2] ... [name n]" );
        System.out.println();
        System.out.println( " Note: This software might need a significant amount of memory (heap space);" );
        System.out
                .println( "       hence use \"-Xms128m -Xmx512m\" (or more) to prevent a \"java.lang.OutOfMemoryError\"." );
        System.out.println();
        System.out.println( " Options: " );
        System.out.println( surfacing.DETAILEDNESS_OPTION + ": level of detail for similarities output file (default:"
                + DETAILEDNESS_DEFAULT + ")" );
        System.out.println( surfacing.IGNORE_COMBINATION_WITH_SAME_OPTION
                + ": to ignore combinations with self (default: not to ignore)" );
        System.out
                .println( surfacing.IGNORE_DOMAINS_WITHOUT_COMBINATIONS_IN_ALL_SPECIES_OPTION
                        + ": to ignore domains without combinations in any species (for similarity calc purposes, not for parsimony analyses) (default: not to ignore)" );
        System.out
                .println( surfacing.IGNORE_DOMAINS_SPECIFIC_TO_ONE_SPECIES_OPTION
                        + ": to ignore domains specific to one species (for similarity calc purposes, not for parsimony analyses) (default: not to ignore)" );
        System.out.println( surfacing.NOT_IGNORE_DUFS_OPTION
                + ": to _not_ ignore DUFs (domains with unknown function) (default: ignore DUFs)" );
        System.out
                .println( surfacing.IGNORE_VIRAL_IDS
                        + ": to ignore domains with ids containing 'vir', 'retro', 'transpos', 'phage', or starting with 'rv' or 'gag_'" );
        System.out.println( surfacing.DOMAIN_SIMILARITY_SORT_OPTION + ": sorting for similarities (default: "
                + DOMAIN_SORT_FILD_DEFAULT + ")" );
        System.out.println( surfacing.OUTPUT_FILE_OPTION + ": name for (main) output file (mandatory)" );
        System.out.println( surfacing.MAX_E_VALUE_OPTION + ": max (inclusive) E-value" );
        System.out.println( surfacing.MAX_ALLOWED_OVERLAP_OPTION + ": maximal allowed domain overlap" );
        System.out.println( surfacing.NO_ENGULFING_OVERLAP_OPTION + ": to ignore engulfed lower confidence domains" );
        System.out.println( surfacing.SPECIES_MATRIX_OPTION + ": species matrix" );
        System.out.println( surfacing.SCORING_OPTION + ": scoring (default:" + SCORING_DEFAULT + ")" );
        System.out.println( surfacing.DOMAIN_COUNT_SORT_OPTION + ": sorting for domain counts (default:"
                + DOMAINS_SORT_ORDER_DEFAULT + ")" );
        System.out.println( surfacing.DOMAIN_SIMILARITY_PRINT_OPTION + ": domain similarity print option (default:"
                + DOMAIN_SIMILARITY_PRINT_OPTION_DEFAULT + ")" );
        System.out.println( surfacing.CUTOFF_SCORE_FILE_OPTION + ": cutoff score file" );
        System.out.println( surfacing.DOMAIN_SIMILARITY_SORT_BY_SPECIES_COUNT_FIRST_OPTION
                + ": sort by species count first" );
        System.out.println( surfacing.OUTPUT_DIR_OPTION + ": output directory" );
        System.out.println( surfacing.PFAM_TO_GO_FILE_USE_OPTION + ": Pfam to GO mapping file" );
        System.out.println( surfacing.GO_OBO_FILE_USE_OPTION + ": GO terms file (OBO format)" );
        System.out.println( surfacing.GO_NAMESPACE_LIMIT_OPTION + ": limit GO term to one GO namespace" );
        System.out.println( surfacing.PAIRWISE_DOMAIN_COMPARISONS_OPTION
                + "[=<suffix for pairwise comparison output files>]: to perform pairwise comparison based analyses" );
        System.out.println( surfacing.INPUT_SPECIES_TREE_OPTION
                + ": species tree, to perform (Dollo, Fitch) parismony analyses" );
        System.out
                .println( surfacing.DISPLAY_M_HISTOGRAMS_OPTION + ": to display multiple histograms (using fluorite)" );
        System.out
                .println( JACKNIFE_OPTION
                        + ": perform jacknife resampling for domain and binary domain combination based distance matrices [default resamplings: "
                        + JACKNIFE_NUMBER_OF_RESAMPLINGS_DEFAULT + "]" );
        System.out.println( JACKNIFE_RATIO_OPTION + ": ratio for jacknife resampling [default: "
                + JACKNIFE_RATIO_DEFAULT + "]" );
        System.out.println( JACKNIFE_RANDOM_SEED_OPTION
                + ": seed for random number generator for jacknife resampling [default: "
                + JACKNIFE_RANDOM_SEED_DEFAULT + "]" );
        //        System.out.println( surfacing.INFER_SPECIES_TREES_OPTION
        //                + ": to infer NJ species trees based on shared domains/binary domain combinations" );
        System.out
                .println( surfacing.INPUT_SPECIES_TREE_OPTION
                        + "=<treefiles in phyloXML format, separated by #>: to infer domain/binary domain combination gains/losses on given species trees" );
        System.out.println( surfacing.FILTER_POSITIVE_OPTION
                + "=<file>: to filter out proteins not containing at least one domain listed in <file>" );
        System.out.println( surfacing.FILTER_NEGATIVE_OPTION
                + "=<file>: to filter out proteins containing at least one domain listed in <file>" );
        System.out.println( surfacing.FILTER_NEGATIVE_DOMAINS_OPTION
                + "=<file>: to filter out (ignore) domains listed in <file>" );
        System.out.println( surfacing.INPUT_FILES_FROM_FILE_OPTION + "=<file>: to read input files from <file>" );
        System.out
                .println( surfacing.RANDOM_SEED_FOR_FITCH_PARSIMONY_OPTION
                        + "=<seed>: seed for random number generator for Fitch Parsimony analysis (type: long, default: no randomization - given a choice, prefer absence" );
        System.out.println( surfacing.CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS
                + ": to consider directedness in binary combinations: e.g. A-B != B-A" );
        System.out.println( surfacing.CONSIDER_DOMAIN_COMBINATION_DIRECTEDNESS_AND_ADJACENCY
                + ": to consider directedness and adjacency in binary combinations" );
        System.out
                .println( surfacing.SEQ_EXTRACT_OPTION
                        + "=<domain ids (Pfam names)>: to extract sequence names of sequences containing matching domains and/or domain-sequences (order N to C) (domain separator: '~', domain sequences speparator: '#', e.g. 'NACHT#BIR~CARD')" );
        System.out.println( surfacing.SECONDARY_FEATURES_PARSIMONY_MAP_FILE
                + "=<file>: to perfom parsimony analysis on secondary features" );
        System.out.println( surfacing.PLUS_MINUS_ANALYSIS_OPTION + "=<file>: to presence/absence genome analysis" );
        System.out.println( surfacing.DOMAIN_COMBINITONS_OUTPUT_OPTION_FOR_GRAPH_ANALYSIS
                + ": to output binary domain combinations for (downstream) graph analysis" );
        System.out.println( surfacing.OUTPUT_LIST_OF_ALL_PROTEINS_OPTIONS + ": to output all proteins per domain" );
        System.out.println();
        System.out.println();
        System.out.println( "Example: java -Xms128m -Xmx512m -cp path/to/forester.jar"
                + "org.forester.application.surfacing -detail=punctilious -o=TEST.html -pwc=TEST"
                + " -cos=Pfam_ls_22_TC2 -p2g=pfam2go -obo=gene_ontology_edit.obo "
                + "-dc_sort=dom -ignore_with_self -no_singles -e=0.001 -mo=1 -no_eo "
                + "-ds_output=detailed_html -scoring=domains -sort=alpha -" + JACKNIFE_OPTION
                + "=50 human mouse brafl strpu" );
        System.out.println();
    }

    private static void processFilter( final File filter_file, final SortedSet<DomainId> filter ) {
        SortedSet<String> filter_str = null;
        try {
            filter_str = ForesterUtil.file2set( filter_file );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        if ( filter_str != null ) {
            for( final String string : filter_str ) {
                filter.add( new DomainId( string ) );
            }
        }
        if ( VERBOSE ) {
            System.out.println( "Filter:" );
            for( final DomainId domainId : filter ) {
                System.out.println( domainId.getId() );
            }
        }
    }

    private static String[][] processInputFileNames( final String[] names ) {
        final String[][] input_file_properties = new String[ names.length ][];
        for( int i = 0; i < names.length; ++i ) {
            if ( names[ i ].indexOf( SEPARATOR_FOR_INPUT_VALUES ) < 0 ) {
                input_file_properties[ i ] = new String[ 2 ];
                input_file_properties[ i ][ 0 ] = names[ i ];
                input_file_properties[ i ][ 1 ] = names[ i ];
            }
            else {
                input_file_properties[ i ] = names[ i ].split( surfacing.SEPARATOR_FOR_INPUT_VALUES + "" );
                if ( input_file_properties[ i ].length != 3 ) {
                    ForesterUtil
                            .fatalError( surfacing.PRG_NAME,
                                         "properties for the input files (hmmpfam output) are expected "
                                                 + "to be in the following format \"<hmmpfam output file>#<species>\" (or just one word, which is both the filename and the species id), instead received \""
                                                 + names[ i ] + "\"" );
                }
            }
            final String error = ForesterUtil.isReadableFile( new File( input_file_properties[ i ][ 0 ] ) );
            if ( !ForesterUtil.isEmpty( error ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, error );
            }
        }
        return input_file_properties;
    }

    private static void processPlusMinusAnalysisOption( final CommandLineArguments cla,
                                                        final List<String> high_copy_base,
                                                        final List<String> high_copy_target,
                                                        final List<String> low_copy,
                                                        final List<Object> numbers ) {
        if ( cla.isOptionSet( surfacing.PLUS_MINUS_ANALYSIS_OPTION ) ) {
            if ( !cla.isOptionValueSet( surfacing.PLUS_MINUS_ANALYSIS_OPTION ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "no value for 'plus-minus' file: -"
                        + surfacing.PLUS_MINUS_ANALYSIS_OPTION + "=<file>" );
            }
            final File plus_minus_file = new File( cla.getOptionValue( surfacing.PLUS_MINUS_ANALYSIS_OPTION ) );
            final String msg = ForesterUtil.isReadableFile( plus_minus_file );
            if ( !ForesterUtil.isEmpty( msg ) ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, "can not read from \"" + plus_minus_file + "\": " + msg );
            }
            processPlusMinusFile( plus_minus_file, high_copy_base, high_copy_target, low_copy, numbers );
        }
    }

    // First numbers is minimal difference, second is factor.
    private static void processPlusMinusFile( final File plus_minus_file,
                                              final List<String> high_copy_base,
                                              final List<String> high_copy_target,
                                              final List<String> low_copy,
                                              final List<Object> numbers ) {
        Set<String> species_set = null;
        int min_diff = PLUS_MINUS_ANALYSIS_MIN_DIFF_DEFAULT;
        double factor = PLUS_MINUS_ANALYSIS_FACTOR_DEFAULT;
        try {
            species_set = ForesterUtil.file2set( plus_minus_file );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        if ( species_set != null ) {
            for( final String species : species_set ) {
                final String species_trimmed = species.substring( 1 );
                if ( species.startsWith( "+" ) ) {
                    if ( low_copy.contains( species_trimmed ) ) {
                        ForesterUtil.fatalError( surfacing.PRG_NAME,
                                                 "species/genome names can not appear with both '+' and '-' suffix, as appears the case for: \""
                                                         + species_trimmed + "\"" );
                    }
                    high_copy_base.add( species_trimmed );
                }
                else if ( species.startsWith( "*" ) ) {
                    if ( low_copy.contains( species_trimmed ) ) {
                        ForesterUtil.fatalError( surfacing.PRG_NAME,
                                                 "species/genome names can not appear with both '*' and '-' suffix, as appears the case for: \""
                                                         + species_trimmed + "\"" );
                    }
                    high_copy_target.add( species_trimmed );
                }
                else if ( species.startsWith( "-" ) ) {
                    if ( high_copy_base.contains( species_trimmed ) || high_copy_target.contains( species_trimmed ) ) {
                        ForesterUtil.fatalError( surfacing.PRG_NAME,
                                                 "species/genome names can not appear with both '+' or '*' and '-' suffix, as appears the case for: \""
                                                         + species_trimmed + "\"" );
                    }
                    low_copy.add( species_trimmed );
                }
                else if ( species.startsWith( "$D" ) ) {
                    try {
                        min_diff = Integer.parseInt( species.substring( 3 ) );
                    }
                    catch ( final NumberFormatException e ) {
                        ForesterUtil.fatalError( surfacing.PRG_NAME,
                                                 "could not parse integer value for minimal difference from: \""
                                                         + species.substring( 3 ) + "\"" );
                    }
                }
                else if ( species.startsWith( "$F" ) ) {
                    try {
                        factor = Double.parseDouble( species.substring( 3 ) );
                    }
                    catch ( final NumberFormatException e ) {
                        ForesterUtil.fatalError( surfacing.PRG_NAME, "could not parse double value for factor from: \""
                                + species.substring( 3 ) + "\"" );
                    }
                }
                else if ( species.startsWith( "#" ) ) {
                    // Comment, ignore.
                }
                else {
                    ForesterUtil
                            .fatalError( surfacing.PRG_NAME,
                                         "species/genome names in 'plus minus' file must begin with '*' (high copy target genome), '+' (high copy base genomes), '-' (low copy genomes), '$D=<integer>' minimal Difference (default is 1), '$F=<double>' factor (default is 1.0), double), or '#' (ignore) suffix, encountered: \""
                                                 + species + "\"" );
                }
                numbers.add( new Integer( min_diff + "" ) );
                numbers.add( new Double( factor + "" ) );
            }
        }
        else {
            ForesterUtil.fatalError( surfacing.PRG_NAME, "'plus minus' file [" + plus_minus_file + "] appears empty" );
        }
    }

    private static void writePresentToNexus( final File output_file,
                                             final File positive_filter_file,
                                             final SortedSet<DomainId> filter,
                                             final List<GenomeWideCombinableDomains> gwcd_list ) {
        try {
            SurfacingUtil
                    .writeMatrixToFile( DomainParsimonyCalculator
                            .createMatrixOfDomainPresenceOrAbsence( gwcd_list, positive_filter_file == null ? null
                                    : filter ), output_file + DOMAINS_PRESENT_NEXUS, Format.NEXUS_BINARY );
            SurfacingUtil.writeMatrixToFile( DomainParsimonyCalculator
                    .createMatrixOfBinaryDomainCombinationPresenceOrAbsence( gwcd_list ), output_file
                    + BDC_PRESENT_NEXUS, Format.NEXUS_BINARY );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getLocalizedMessage() );
        }
    }

    private static void writeProteinListsForAllSpecies( final File output_dir,
                                                        final SortedMap<Species, List<Protein>> protein_lists_per_species,
                                                        final List<GenomeWideCombinableDomains> gwcd_list ) {
        final SortedSet<DomainId> all_domains = new TreeSet<DomainId>();
        for( final GenomeWideCombinableDomains gwcd : gwcd_list ) {
            all_domains.addAll( gwcd.getAllDomainIds() );
        }
        for( final DomainId domain : all_domains ) {
            final File out = new File( output_dir + ForesterUtil.FILE_SEPARATOR + domain + SEQ_EXTRACT_SUFFIX );
            SurfacingUtil.checkForOutputFileWriteability( out );
            try {
                final Writer proteins_file_writer = new BufferedWriter( new FileWriter( out ) );
                SurfacingUtil.extractProteinNames( protein_lists_per_species, domain, proteins_file_writer, "\t" );
                proteins_file_writer.close();
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( surfacing.PRG_NAME, e.getLocalizedMessage() );
            }
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote proteins list to \"" + out + "\"" );
        }
    }
}
