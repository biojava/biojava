// $Id: DomainCountsDifferenceUtil.java,v 1.15 2008/08/16 00:41:13 cmzmasek Exp
// $
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

package org.forester.surfacing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.go.GoId;
import org.forester.go.GoTerm;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

/*
 * Poorly designed static class which essential has one method:
 * calculateCopyNumberDifferences.
 */
public final class DomainCountsDifferenceUtil {

    private final static NumberFormat          FORMATTER                                   = new DecimalFormat( "0.0E0" );
    private static final COPY_CALCULATION_MODE COPY_CALC_MODE_FOR_HIGH_COPY_TARGET_SPECIES = COPY_CALCULATION_MODE.MIN;
    private static final COPY_CALCULATION_MODE COPY_CALC_MODE_FOR_HIGH_COPY_BASE_SPECIES   = COPY_CALCULATION_MODE.MIN;
    private static final COPY_CALCULATION_MODE COPY_CALC_MODE_FOR_LOW_COPY_SPECIES         = COPY_CALCULATION_MODE.MAX;
    private static final String                PLUS_MINUS_PROTEINS_FILE_DOM_SUFFIX         = ".prot";

    //FIXME really needs to be tested! 
    private static void addCounts( final SortedMap<BinaryDomainCombination, List<Integer>> copy_counts,
                                   final BinaryDomainCombination dc,
                                   final GenomeWideCombinableDomains genome,
                                   final Set<BinaryDomainCombination> bdc ) {
        if ( !copy_counts.containsKey( dc ) ) {
            copy_counts.put( dc, new ArrayList<Integer>() );
        }
        if ( bdc.contains( dc )
                && ( ( ( BasicCombinableDomains ) genome.get( dc.getId0() ) ).getCombiningDomains().get( dc.getId1() ) != null ) ) {
            final int count = ( ( BasicCombinableDomains ) genome.get( dc.getId0() ) ).getCombiningDomains().get( dc
                    .getId1() );
            copy_counts.get( dc ).add( count );
        }
        else {
            copy_counts.get( dc ).add( 0 );
        }
    }

    private static void addCounts( final SortedMap<DomainId, List<Integer>> copy_counts,
                                   final DomainId domain,
                                   final GenomeWideCombinableDomains genome ) {
        if ( !copy_counts.containsKey( domain ) ) {
            copy_counts.put( domain, new ArrayList<Integer>() );
        }
        if ( genome.contains( domain ) ) {
            copy_counts.get( domain ).add( genome.get( domain ).getKeyDomainProteinsCount() );
        }
        else {
            copy_counts.get( domain ).add( 0 );
        }
    }

    private static StringBuilder addGoInformation( final DomainId d,
                                                   final Map<DomainId, List<GoId>> domain_id_to_go_ids_map,
                                                   final Map<GoId, GoTerm> go_id_to_term_map ) {
        final StringBuilder sb = new StringBuilder();
        if ( ( domain_id_to_go_ids_map == null ) || domain_id_to_go_ids_map.isEmpty()
                || !domain_id_to_go_ids_map.containsKey( d ) ) {
            return sb;
        }
        final List<GoId> go_ids = domain_id_to_go_ids_map.get( d );
        for( int i = 0; i < go_ids.size(); ++i ) {
            final GoId go_id = go_ids.get( i );
            if ( go_id_to_term_map.containsKey( go_id ) ) {
                appendGoTerm( sb, go_id_to_term_map.get( go_id ) );
                sb.append( "<br>" );
            }
            else {
                sb.append( "go id \"" + go_id + "\" not found [" + d.getId() + "]" );
            }
        }
        return sb;
    }

    private static void appendGoTerm( final StringBuilder sb, final GoTerm go_term ) {
        final GoId go_id = go_term.getGoId();
        sb.append( "<a href=\"" + SurfacingConstants.AMIGO_LINK + go_id + "\" target=\"amigo_window\">" + go_id
                + "</a>" );
        sb.append( ":" );
        sb.append( go_term.getName() );
        sb.append( " [" );
        sb.append( go_term.getGoNameSpace().toShortString() );
        sb.append( "]" );
    }

    public static void calculateCopyNumberDifferences( final List<GenomeWideCombinableDomains> genomes,
                                                       final SortedMap<Species, List<Protein>> protein_lists_per_species,
                                                       final List<String> high_copy_base_species,
                                                       final List<String> high_copy_target_species,
                                                       final List<String> low_copy_species,
                                                       final int min_diff,
                                                       final Double factor,
                                                       final File plain_output_dom,
                                                       final File html_output_dom,
                                                       final File html_output_dc,
                                                       final Map<DomainId, List<GoId>> domain_id_to_go_ids_map,
                                                       final Map<GoId, GoTerm> go_id_to_term_map,
                                                       final File all_domains_go_ids_out_dom,
                                                       final File passing_domains_go_ids_out_dom,
                                                       final File proteins_file_base ) throws IOException {
        if ( genomes.size() < 1 ) {
            throw new IllegalArgumentException( "attempt to use empty list of genomes for domain difference calculation" );
        }
        if ( ( high_copy_base_species.size() < 1 ) || ( low_copy_species.size() < 1 ) ) {
            throw new IllegalArgumentException( "attempt to use empty list of species for domain difference calculation" );
        }
        if ( high_copy_base_species.contains( high_copy_target_species )
                || low_copy_species.contains( high_copy_target_species ) ) {
            throw new IllegalArgumentException( "species [" + high_copy_target_species
                    + "] appears in other list as well" );
        }
        if ( min_diff < 0 ) {
            throw new IllegalArgumentException( "attempt to use negative addition [" + min_diff + "]" );
        }
        if ( factor <= 0.0 ) {
            throw new IllegalArgumentException( "attempt to use factor equal or smaller than 0.0 [" + factor + "]" );
        }
        SurfacingUtil.checkForOutputFileWriteability( plain_output_dom );
        SurfacingUtil.checkForOutputFileWriteability( html_output_dom );
        SurfacingUtil.checkForOutputFileWriteability( html_output_dc );
        SurfacingUtil.checkForOutputFileWriteability( all_domains_go_ids_out_dom );
        SurfacingUtil.checkForOutputFileWriteability( passing_domains_go_ids_out_dom );
        final Writer plain_writer = new BufferedWriter( new FileWriter( plain_output_dom ) );
        final Writer html_writer = new BufferedWriter( new FileWriter( html_output_dom ) );
        final Writer html_writer_dc = new BufferedWriter( new FileWriter( html_output_dc ) );
        final Writer all_gos_writer = new BufferedWriter( new FileWriter( all_domains_go_ids_out_dom ) );
        final Writer passing_gos_writer = new BufferedWriter( new FileWriter( passing_domains_go_ids_out_dom ) );
        final SortedMap<DomainId, Double> high_copy_base_values = new TreeMap<DomainId, Double>();
        final SortedMap<DomainId, Double> high_copy_target_values = new TreeMap<DomainId, Double>();
        final SortedMap<DomainId, Double> low_copy_values = new TreeMap<DomainId, Double>();
        final SortedMap<DomainId, List<Integer>> high_copy_base_copy_counts = new TreeMap<DomainId, List<Integer>>();
        final SortedMap<DomainId, List<Integer>> high_copy_target_copy_counts = new TreeMap<DomainId, List<Integer>>();
        final SortedMap<DomainId, List<Integer>> low_copy_copy_counts = new TreeMap<DomainId, List<Integer>>();
        final SortedSet<DomainId> all_domains = new TreeSet<DomainId>();
        final SortedMap<BinaryDomainCombination, Double> high_copy_base_values_dc = new TreeMap<BinaryDomainCombination, Double>();
        final SortedMap<BinaryDomainCombination, Double> high_copy_target_values_dc = new TreeMap<BinaryDomainCombination, Double>();
        final SortedMap<BinaryDomainCombination, Double> low_copy_values_dc = new TreeMap<BinaryDomainCombination, Double>();
        final SortedMap<BinaryDomainCombination, List<Integer>> high_copy_base_copy_counts_dc = new TreeMap<BinaryDomainCombination, List<Integer>>();
        final SortedMap<BinaryDomainCombination, List<Integer>> high_copy_target_copy_counts_dc = new TreeMap<BinaryDomainCombination, List<Integer>>();
        final SortedMap<BinaryDomainCombination, List<Integer>> low_copy_copy_counts_dc = new TreeMap<BinaryDomainCombination, List<Integer>>();
        final SortedSet<BinaryDomainCombination> all_dcs = new TreeSet<BinaryDomainCombination>();
        final Map<String, Set<BinaryDomainCombination>> bdcs_per_genome = new HashMap<String, Set<BinaryDomainCombination>>();
        final SortedSet<GoId> go_ids_of_passing_domains = new TreeSet<GoId>();
        final SortedSet<GoId> go_ids_all = new TreeSet<GoId>();
        for( final GenomeWideCombinableDomains genome : genomes ) {
            final SortedSet<DomainId> domains = genome.getAllDomainIds();
            final SortedSet<BinaryDomainCombination> dcs = genome.toBinaryDomainCombinations();
            final String species = genome.getSpecies().getSpeciesId();
            bdcs_per_genome.put( species, genome.toBinaryDomainCombinations() );
            for( final DomainId d : domains ) {
                all_domains.add( d );
                if ( domain_id_to_go_ids_map.containsKey( d ) ) {
                    go_ids_all.addAll( domain_id_to_go_ids_map.get( d ) );
                }
            }
            for( final BinaryDomainCombination dc : dcs ) {
                all_dcs.add( dc );
            }
        }
        for( final DomainId domain : all_domains ) {
            for( final GenomeWideCombinableDomains genome : genomes ) {
                final String species = genome.getSpecies().getSpeciesId();
                if ( high_copy_base_species.contains( species ) ) {
                    DomainCountsDifferenceUtil.addCounts( high_copy_base_copy_counts, domain, genome );
                }
                if ( high_copy_target_species.contains( species ) ) {
                    DomainCountsDifferenceUtil.addCounts( high_copy_target_copy_counts, domain, genome );
                }
                if ( low_copy_species.contains( species ) ) {
                    DomainCountsDifferenceUtil.addCounts( low_copy_copy_counts, domain, genome );
                }
            }
        }
        for( final BinaryDomainCombination dc : all_dcs ) {
            for( final GenomeWideCombinableDomains genome : genomes ) {
                final String species = genome.getSpecies().getSpeciesId();
                if ( high_copy_base_species.contains( species ) ) {
                    DomainCountsDifferenceUtil.addCounts( high_copy_base_copy_counts_dc, dc, genome, bdcs_per_genome
                            .get( species ) );
                }
                if ( high_copy_target_species.contains( species ) ) {
                    DomainCountsDifferenceUtil.addCounts( high_copy_target_copy_counts_dc, dc, genome, bdcs_per_genome
                            .get( species ) );
                }
                if ( low_copy_species.contains( species ) ) {
                    DomainCountsDifferenceUtil.addCounts( low_copy_copy_counts_dc, dc, genome, bdcs_per_genome
                            .get( species ) );
                }
            }
        }
        for( final DomainId domain : all_domains ) {
            calculateDomainCountsBasedValue( high_copy_target_values,
                                             high_copy_target_copy_counts,
                                             domain,
                                             COPY_CALC_MODE_FOR_HIGH_COPY_TARGET_SPECIES );
            calculateDomainCountsBasedValue( high_copy_base_values,
                                             high_copy_base_copy_counts,
                                             domain,
                                             COPY_CALC_MODE_FOR_HIGH_COPY_BASE_SPECIES );
            calculateDomainCountsBasedValue( low_copy_values,
                                             low_copy_copy_counts,
                                             domain,
                                             COPY_CALC_MODE_FOR_LOW_COPY_SPECIES );
        }
        for( final BinaryDomainCombination dc : all_dcs ) {
            calculateDomainCountsBasedValue( high_copy_target_values_dc,
                                             high_copy_target_copy_counts_dc,
                                             dc,
                                             COPY_CALC_MODE_FOR_HIGH_COPY_TARGET_SPECIES );
            calculateDomainCountsBasedValue( high_copy_base_values_dc,
                                             high_copy_base_copy_counts_dc,
                                             dc,
                                             COPY_CALC_MODE_FOR_HIGH_COPY_BASE_SPECIES );
            calculateDomainCountsBasedValue( low_copy_values_dc,
                                             low_copy_copy_counts_dc,
                                             dc,
                                             COPY_CALC_MODE_FOR_LOW_COPY_SPECIES );
        }
        writeDomainValuesToFiles( genomes,
                                  high_copy_base_species,
                                  high_copy_target_species,
                                  low_copy_species,
                                  min_diff,
                                  factor,
                                  domain_id_to_go_ids_map,
                                  go_id_to_term_map,
                                  plain_writer,
                                  html_writer,
                                  proteins_file_base,
                                  high_copy_base_values,
                                  high_copy_target_values,
                                  low_copy_values,
                                  all_domains,
                                  go_ids_of_passing_domains,
                                  protein_lists_per_species );
        writeDomainCombinationValuesToFiles( genomes,
                                             high_copy_base_species,
                                             high_copy_target_species,
                                             low_copy_species,
                                             min_diff,
                                             factor,
                                             html_writer_dc,
                                             high_copy_base_values_dc,
                                             high_copy_target_values_dc,
                                             low_copy_values_dc,
                                             all_dcs,
                                             bdcs_per_genome );
        writeGoIdsToFile( all_gos_writer, go_ids_all );
        writeGoIdsToFile( passing_gos_writer, go_ids_of_passing_domains );
    }

    private static void calculateDomainCountsBasedValue( final SortedMap<BinaryDomainCombination, Double> copy_values,
                                                         final SortedMap<BinaryDomainCombination, List<Integer>> copy_counts,
                                                         final BinaryDomainCombination bdc,
                                                         final COPY_CALCULATION_MODE copy_calc_mode ) {
        if ( copy_counts.containsKey( bdc ) ) {
            switch ( copy_calc_mode ) {
                case MAX:
                    DomainCountsDifferenceUtil.calculateMaxCount( copy_values, copy_counts, bdc );
                    break;
                case MIN:
                    DomainCountsDifferenceUtil.calculateMinCount( copy_values, copy_counts, bdc );
                    break;
                case MEAN:
                    DomainCountsDifferenceUtil.calculateMeanCount( copy_values, copy_counts, bdc );
                    break;
                case MEDIAN:
                    DomainCountsDifferenceUtil.calculateMedianCount( copy_values, copy_counts, bdc );
                    break;
                default:
                    throw new IllegalArgumentException();
            }
        }
        else {
            copy_values.put( bdc, Double.valueOf( 0.0 ) );
        }
    }

    private static void calculateDomainCountsBasedValue( final SortedMap<DomainId, Double> copy_values,
                                                         final SortedMap<DomainId, List<Integer>> copy_counts,
                                                         final DomainId domain,
                                                         final COPY_CALCULATION_MODE copy_calc_mode ) {
        if ( copy_counts.containsKey( domain ) ) {
            switch ( copy_calc_mode ) {
                case MAX:
                    DomainCountsDifferenceUtil.calculateMaxCount( copy_values, copy_counts, domain );
                    break;
                case MIN:
                    DomainCountsDifferenceUtil.calculateMinCount( copy_values, copy_counts, domain );
                    break;
                case MEAN:
                    DomainCountsDifferenceUtil.calculateMeanCount( copy_values, copy_counts, domain );
                    break;
                case MEDIAN:
                    DomainCountsDifferenceUtil.calculateMedianCount( copy_values, copy_counts, domain );
                    break;
                default:
                    throw new IllegalArgumentException();
            }
        }
        else {
            copy_values.put( domain, Double.valueOf( 0.0 ) );
        }
    }

    private static void calculateMaxCount( final SortedMap<BinaryDomainCombination, Double> results,
                                           final SortedMap<BinaryDomainCombination, List<Integer>> copy_counts,
                                           final BinaryDomainCombination bdc ) {
        final List<Integer> counts = copy_counts.get( bdc );
        int max = 0;
        for( final Integer count : counts ) {
            if ( count > max ) {
                max = count;
            }
        }
        results.put( bdc, ( double ) max );
    }

    private static void calculateMaxCount( final SortedMap<DomainId, Double> results,
                                           final SortedMap<DomainId, List<Integer>> copy_counts,
                                           final DomainId domain ) {
        final List<Integer> counts = copy_counts.get( domain );
        int max = 0;
        for( final Integer count : counts ) {
            if ( count > max ) {
                max = count;
            }
        }
        results.put( domain, ( double ) max );
    }

    private static void calculateMeanCount( final SortedMap<BinaryDomainCombination, Double> results,
                                            final SortedMap<BinaryDomainCombination, List<Integer>> copy_counts,
                                            final BinaryDomainCombination bdc ) {
        final List<Integer> counts = copy_counts.get( bdc );
        int sum = 0;
        for( final Integer count : counts ) {
            sum += count;
        }
        results.put( bdc, ( ( double ) sum ) / ( ( double ) counts.size() ) );
    }

    private static void calculateMeanCount( final SortedMap<DomainId, Double> results,
                                            final SortedMap<DomainId, List<Integer>> copy_counts,
                                            final DomainId domain ) {
        final List<Integer> counts = copy_counts.get( domain );
        int sum = 0;
        for( final Integer count : counts ) {
            sum += count;
        }
        results.put( domain, ( ( double ) sum ) / ( ( double ) counts.size() ) );
    }

    private static void calculateMedianCount( final SortedMap<BinaryDomainCombination, Double> results,
                                              final SortedMap<BinaryDomainCombination, List<Integer>> copy_counts,
                                              final BinaryDomainCombination bdc ) {
        final List<Integer> counts = copy_counts.get( bdc );
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final Integer count : counts ) {
            stats.addValue( count );
        }
        results.put( bdc, stats.median() );
    }

    private static void calculateMedianCount( final SortedMap<DomainId, Double> results,
                                              final SortedMap<DomainId, List<Integer>> copy_counts,
                                              final DomainId domain ) {
        final List<Integer> counts = copy_counts.get( domain );
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final Integer count : counts ) {
            stats.addValue( count );
        }
        results.put( domain, stats.median() );
    }

    private static void calculateMinCount( final SortedMap<BinaryDomainCombination, Double> results,
                                           final SortedMap<BinaryDomainCombination, List<Integer>> copy_counts,
                                           final BinaryDomainCombination bdc ) {
        final List<Integer> counts = copy_counts.get( bdc );
        int min = Integer.MAX_VALUE;
        for( final Integer count : counts ) {
            if ( count < min ) {
                min = count;
            }
        }
        results.put( bdc, ( double ) min );
    }

    private static void calculateMinCount( final SortedMap<DomainId, Double> results,
                                           final SortedMap<DomainId, List<Integer>> copy_counts,
                                           final DomainId domain ) {
        final List<Integer> counts = copy_counts.get( domain );
        int min = Integer.MAX_VALUE;
        for( final Integer count : counts ) {
            if ( count < min ) {
                min = count;
            }
        }
        results.put( domain, ( double ) min );
    }

    private static String combinableDomaindToString( final CombinableDomains cd ) {
        final StringBuilder sb = new StringBuilder();
        sb.append( cd.getKeyDomainProteinsCount() );
        sb.append( "\t[" );
        sb.append( FORMATTER.format( cd.getKeyDomainConfidenceDescriptiveStatistics().median() ) );
        sb.append( "]" );
        return sb.toString();
    }

    private static String combinableDomaindToStringHtml( final CombinableDomains cd ) {
        final StringBuilder sb = new StringBuilder();
        sb.append( "[" );
        sb.append( cd.getKeyDomainCount() );
        sb.append( ", <b>" );
        sb.append( cd.getKeyDomainProteinsCount() );
        sb.append( "</b>, " );
        sb.append( cd.getNumberOfCombinableDomains() );
        sb.append( "]</td><td>[" );
        sb.append( FORMATTER.format( cd.getKeyDomainConfidenceDescriptiveStatistics().median() ) );
        sb.append( "]</td><td>" );
        sb.append( cd.getCombiningDomainIdsAsStringBuilder() );
        return sb.toString();
    }

    private static void writeCopyNumberValues( final SortedMap<BinaryDomainCombination, Double> copy_means,
                                               final BinaryDomainCombination bdc,
                                               final GenomeWideCombinableDomains genome,
                                               final Map<String, Set<BinaryDomainCombination>> bdcs_per_genome,
                                               final String species,
                                               final Writer html_writer,
                                               final String color ) throws IOException {
        html_writer.write( "<td> " );
        if ( !ForesterUtil.isEmpty( color ) ) {
            html_writer.write( "<font color=\"" + color + "\">" );
        }
        html_writer.write( "<b>" + species + ":</b> " );
        if ( !ForesterUtil.isEmpty( color ) ) {
            html_writer.write( "</font>" );
        }
        html_writer.write( "</td><td>" );
        if ( bdcs_per_genome.get( species ).contains( bdc ) && ( copy_means.get( bdc ) > 0 ) ) {
            final int count = ( ( BasicCombinableDomains ) genome.get( bdc.getId0() ) ).getCombiningDomains().get( bdc
                    .getId1() );
            html_writer.write( count + "" );
        }
        else {
            html_writer.write( "0" );
        }
        html_writer.write( "</td>" );
    }

    private static void writeCopyNumberValues( final SortedMap<DomainId, Double> copy_means,
                                               final DomainId domain,
                                               final GenomeWideCombinableDomains genome,
                                               final String species,
                                               final Writer plain_writer,
                                               final Writer html_writer,
                                               final String color ) throws IOException {
        plain_writer.write( "  " + species + "\t" );
        html_writer.write( "<td> " );
        if ( !ForesterUtil.isEmpty( color ) ) {
            html_writer.write( "<font color=\"" + color + "\">" );
        }
        html_writer.write( "<b>" + species + ":</b> " );
        if ( !ForesterUtil.isEmpty( color ) ) {
            html_writer.write( "</font>" );
        }
        html_writer.write( "</td><td>" );
        if ( genome.contains( domain ) && ( copy_means.get( domain ) > 0 ) ) {
            plain_writer.write( DomainCountsDifferenceUtil.combinableDomaindToString( genome.get( domain ) ) );
            html_writer.write( DomainCountsDifferenceUtil.combinableDomaindToStringHtml( genome.get( domain ) ) );
        }
        else {
            plain_writer.write( "0" );
            html_writer.write( "0" );
        }
        html_writer.write( "</td>" );
        plain_writer.write( SurfacingConstants.NL );
    }

    private static void writeDomainCombinationValuesToFiles( final List<GenomeWideCombinableDomains> genomes,
                                                             final List<String> high_copy_base_species,
                                                             final List<String> high_copy_target_species,
                                                             final List<String> low_copy_species,
                                                             final int min_diff,
                                                             final Double factor,
                                                             final Writer html_writer,
                                                             final SortedMap<BinaryDomainCombination, Double> high_copy_base_values,
                                                             final SortedMap<BinaryDomainCombination, Double> high_copy_target_values,
                                                             final SortedMap<BinaryDomainCombination, Double> low_copy_values,
                                                             final SortedSet<BinaryDomainCombination> all_bdcs,
                                                             final Map<String, Set<BinaryDomainCombination>> bdcs_per_genome )
            throws IOException {
        int counter = 0;
        int total_absense_counter = 0;
        int not_total_absense_counter = 0;
        SurfacingUtil.addHtmlHead( html_writer, "Binary Domain Combination Copy Differences" );
        html_writer.write( "<body><table>" );
        for( final BinaryDomainCombination bdc : all_bdcs ) {
            if ( ( high_copy_base_values.get( bdc ) > 0 ) && ( high_copy_target_values.get( bdc ) > 0 )
                    && ( high_copy_base_values.get( bdc ) >= low_copy_values.get( bdc ) ) ) {
                if ( high_copy_target_values.get( bdc ) >= min_diff + ( factor * low_copy_values.get( bdc ) ) ) {
                    if ( low_copy_values.get( bdc ) <= 0.0 ) {
                        ++total_absense_counter;
                    }
                    else {
                        ++not_total_absense_counter;
                    }
                    ++counter;
                    html_writer.write( "<tr><td><a href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK + bdc.getId0()
                            + "\">" + bdc.getId0() + "</a> = <a href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK
                            + bdc.getId1() + "\">" + bdc.getId1() + "</a>" );
                    html_writer.write( "</td><td>" );
                    html_writer.write( "<table>" );
                    for( final GenomeWideCombinableDomains genome : genomes ) {
                        final String species = genome.getSpecies().getSpeciesId();
                        if ( high_copy_target_species.contains( species ) ) {
                            html_writer.write( "<tr>" );
                            writeCopyNumberValues( high_copy_target_values,
                                                   bdc,
                                                   genome,
                                                   bdcs_per_genome,
                                                   species,
                                                   html_writer,
                                                   "#0000FF" );
                            html_writer.write( "</tr>" );
                        }
                        else if ( low_copy_species.contains( species ) ) {
                            html_writer.write( "<tr>" );
                            writeCopyNumberValues( low_copy_values,
                                                   bdc,
                                                   genome,
                                                   bdcs_per_genome,
                                                   species,
                                                   html_writer,
                                                   "#A0A0A0" );
                            html_writer.write( "</tr>" );
                        }
                        else if ( high_copy_base_species.contains( species ) ) {
                            html_writer.write( "<tr>" );
                            writeCopyNumberValues( high_copy_base_values,
                                                   bdc,
                                                   genome,
                                                   bdcs_per_genome,
                                                   species,
                                                   html_writer,
                                                   "#404040" );
                            html_writer.write( "</tr>" );
                        }
                    }
                    html_writer.write( "</table>" );
                    html_writer.write( "</td></tr>" );
                    html_writer.write( SurfacingConstants.NL );
                }
            }
        }
        html_writer.write( "</table>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<hr>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "Rule 1: high-copy-base > 0 && high-copy-target > 0 && high-copy-base >= low-copy" );
        html_writer.write( "<br>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "Rule 2: high-copy-target >= minimal-difference + ( factor * low-copy )" );
        html_writer.write( "<br>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "Calculation mode for high copy target : " + COPY_CALC_MODE_FOR_HIGH_COPY_TARGET_SPECIES );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Calculation mode for high copy base : " + COPY_CALC_MODE_FOR_HIGH_COPY_BASE_SPECIES );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Calculation mode for low copy : " + COPY_CALC_MODE_FOR_LOW_COPY_SPECIES );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Minimal difference : " + min_diff );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Factor : " + factor );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Lower copy binary domain combinations : " + counter );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Total absence : " + total_absense_counter );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Not total absence : " + not_total_absense_counter );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Total binary domain combinations : " + all_bdcs.size() );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<hr>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "</body></html>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.close();
    }

    private static void writeDomainValuesToFiles( final List<GenomeWideCombinableDomains> genomes,
                                                  final List<String> high_copy_base_species,
                                                  final List<String> high_copy_target_species,
                                                  final List<String> low_copy_species,
                                                  final int min_diff,
                                                  final Double factor,
                                                  final Map<DomainId, List<GoId>> domain_id_to_go_ids_map,
                                                  final Map<GoId, GoTerm> go_id_to_term_map,
                                                  final Writer plain_writer,
                                                  final Writer html_writer,
                                                  final File proteins_file_base,
                                                  final SortedMap<DomainId, Double> high_copy_base_values,
                                                  final SortedMap<DomainId, Double> high_copy_target_values,
                                                  final SortedMap<DomainId, Double> low_copy_values,
                                                  final SortedSet<DomainId> all_domains,
                                                  final SortedSet<GoId> go_ids_of_passing_domains,
                                                  final SortedMap<Species, List<Protein>> protein_lists_per_species )
            throws IOException {
        int counter = 0;
        int total_absense_counter = 0;
        int not_total_absense_counter = 0;
        SurfacingUtil.addHtmlHead( html_writer, "Domain Copy Differences" );
        html_writer.write( "<body><table>" );
        for( final DomainId domain_id : all_domains ) {
            if ( ( high_copy_base_values.get( domain_id ) > 0 ) && ( high_copy_target_values.get( domain_id ) > 0 )
                    && ( high_copy_base_values.get( domain_id ) >= low_copy_values.get( domain_id ) ) ) {
                if ( high_copy_target_values.get( domain_id ) >= min_diff
                        + ( factor * low_copy_values.get( domain_id ) ) ) {
                    if ( low_copy_values.get( domain_id ) <= 0.0 ) {
                        ++total_absense_counter;
                    }
                    else {
                        ++not_total_absense_counter;
                    }
                    ++counter;
                    writeProteinsToFile( proteins_file_base, protein_lists_per_species, domain_id );
                    if ( domain_id_to_go_ids_map.containsKey( domain_id ) ) {
                        go_ids_of_passing_domains.addAll( domain_id_to_go_ids_map.get( domain_id ) );
                    }
                    plain_writer.write( domain_id.getId() );
                    plain_writer.write( SurfacingConstants.NL );
                    html_writer.write( "<tr><td><a href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK
                            + domain_id.getId() + "\">" + domain_id.getId() + "</a></td><td>" );
                    html_writer.write( addGoInformation( domain_id, domain_id_to_go_ids_map, go_id_to_term_map )
                            .toString() );
                    html_writer.write( "</td><td>" );
                    html_writer.write( "<table>" );
                    for( final GenomeWideCombinableDomains genome : genomes ) {
                        final String species = genome.getSpecies().getSpeciesId();
                        if ( high_copy_target_species.contains( species ) ) {
                            html_writer.write( "<tr>" );
                            writeCopyNumberValues( high_copy_target_values,
                                                   domain_id,
                                                   genome,
                                                   species,
                                                   plain_writer,
                                                   html_writer,
                                                   "#0000FF" );
                            html_writer.write( "</tr>" );
                        }
                        else if ( low_copy_species.contains( species ) ) {
                            html_writer.write( "<tr>" );
                            writeCopyNumberValues( low_copy_values,
                                                   domain_id,
                                                   genome,
                                                   species,
                                                   plain_writer,
                                                   html_writer,
                                                   "#A0A0A0" );
                            html_writer.write( "</tr>" );
                        }
                        else if ( high_copy_base_species.contains( species ) ) {
                            html_writer.write( "<tr>" );
                            writeCopyNumberValues( high_copy_base_values,
                                                   domain_id,
                                                   genome,
                                                   species,
                                                   plain_writer,
                                                   html_writer,
                                                   "#404040" );
                            html_writer.write( "</tr>" );
                        }
                    }
                    html_writer.write( "</table>" );
                    html_writer.write( "</td></tr>" );
                    html_writer.write( SurfacingConstants.NL );
                    plain_writer.write( SurfacingConstants.NL );
                }
            }
        }
        html_writer.write( "</table>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<hr>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "Rule 1: high-copy-base > 0 && high-copy-target > 0 && high-copy-base >= low-copy" );
        html_writer.write( "<br>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "Rule 2: high-copy-target >= minimal-difference + ( factor * low-copy )" );
        html_writer.write( "<br>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "Calculation mode for high copy target : " + COPY_CALC_MODE_FOR_HIGH_COPY_TARGET_SPECIES );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Calculation mode for high copy base : " + COPY_CALC_MODE_FOR_HIGH_COPY_BASE_SPECIES );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Calculation mode for low copy : " + COPY_CALC_MODE_FOR_LOW_COPY_SPECIES );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Minimal difference : " + min_diff );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Factor : " + factor );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Lower copy domains : " + counter );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Total absence : " + total_absense_counter );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Not total absence : " + not_total_absense_counter );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<br>" );
        html_writer.write( "Total domains : " + all_domains.size() );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "<hr>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.write( "</body></html>" );
        html_writer.write( SurfacingConstants.NL );
        html_writer.close();
        plain_writer.write( "# Rule 1: high-copy-base > 0 && high-copy-target > 0 && high-copy-base >= low-copy" );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.write( "# Rule 2: high-copy-target >= minimal-difference + ( factor * low-copy )" );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.write( "# Calculation mode for high copy target: " + COPY_CALC_MODE_FOR_HIGH_COPY_TARGET_SPECIES );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.write( "# Calculation mode for high copy base  : " + COPY_CALC_MODE_FOR_HIGH_COPY_BASE_SPECIES );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.write( "# Calculation mode for low copy        : " + COPY_CALC_MODE_FOR_LOW_COPY_SPECIES );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.write( "# Minimal difference: " + min_diff );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.write( "# Factor            : " + factor );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.write( "# Lower copy domains: " + counter );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.write( "# Total absence     : " + total_absense_counter );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.write( "# Not total absence : " + not_total_absense_counter );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.write( "# Total domains     : " + all_domains.size() );
        plain_writer.write( SurfacingConstants.NL );
        plain_writer.close();
    }

    private static void writeGoIdsToFile( final Writer writer, final SortedSet<GoId> gos ) throws IOException {
        for( final GoId go_id : gos ) {
            writer.write( go_id.toString() );
            writer.write( SurfacingConstants.NL );
        }
        writer.close();
    }

    private static void writeProteinsToFile( final File proteins_file_base,
                                             final SortedMap<Species, List<Protein>> protein_lists_per_species,
                                             final DomainId domain_id ) throws IOException {
        final File my_proteins_file = new File( proteins_file_base.getParentFile() + ForesterUtil.FILE_SEPARATOR
                + domain_id + PLUS_MINUS_PROTEINS_FILE_DOM_SUFFIX );
        SurfacingUtil.checkForOutputFileWriteability( my_proteins_file );
        final Writer proteins_file_writer = new BufferedWriter( new FileWriter( my_proteins_file ) );
        SurfacingUtil.extractProteinNames( protein_lists_per_species, domain_id, proteins_file_writer, "\t" );
        proteins_file_writer.close();
        System.out.println( "Wrote proteins list to \"" + my_proteins_file + "\"" );
    }

    public static enum COPY_CALCULATION_MODE {
        MEAN, MEDIAN, MAX, MIN
    }
}
