// $Id: BasicDomainSimilarityCalculator.java,v 1.15 2007/10/10 22:43:35 cmzmasek
// Exp $
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

import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;

public class BasicDomainSimilarityCalculator implements DomainSimilarityCalculator {

    final DomainSimilarity.DomainSimilaritySortField _sort;
    private final boolean                            _sort_by_species_count_first;
    private final boolean                            _treat_as_binary_comparison;

    public BasicDomainSimilarityCalculator( final DomainSimilarity.DomainSimilaritySortField sort,
                                            final boolean sort_by_species_count_first,
                                            final boolean treat_as_binary_comparison ) {
        _sort = sort;
        _sort_by_species_count_first = sort_by_species_count_first;
        _treat_as_binary_comparison = treat_as_binary_comparison;
    }

    public SortedSet<DomainSimilarity> calculateSimilarities( final PairwiseDomainSimilarityCalculator pairwise_calculator,
                                                              final List<GenomeWideCombinableDomains> cdc_list,
                                                              final boolean ignore_domains_without_combinations_in_any_genome,
                                                              final boolean ignore_domains_specific_to_one_genome ) {
        if ( cdc_list.size() < 2 ) {
            throw new IllegalArgumentException( "attempt to calculate multiple combinable domains similarity for less than two combinale domains collections" );
        }
        final SortedSet<DomainSimilarity> similarities = new TreeSet<DomainSimilarity>();
        final SortedSet<DomainId> keys = new TreeSet<DomainId>();
        for( final GenomeWideCombinableDomains cdc : cdc_list ) {
            keys.addAll( ( cdc ).getAllCombinableDomainsIds().keySet() );
        }
        for( final DomainId key : keys ) {
            final List<CombinableDomains> same_id_cd_list = new ArrayList<CombinableDomains>( cdc_list.size() );
            final List<Species> species_with_key_id_domain = new ArrayList<Species>();
            for( final GenomeWideCombinableDomains cdc : cdc_list ) {
                if ( cdc.contains( key ) ) {
                    same_id_cd_list.add( cdc.get( key ) );
                    species_with_key_id_domain.add( cdc.getSpecies() );
                }
            }
            if ( ignore_domains_without_combinations_in_any_genome ) { //TODO: test me..........................................<<<<<<<<<<<<<
                boolean without_combinations = true;
                for( final CombinableDomains cd : same_id_cd_list ) {
                    if ( cd.getNumberOfCombinableDomains() > 0 ) {
                        without_combinations = false;
                        break;
                    }
                }
                if ( without_combinations ) {
                    continue;
                }
            }
            // BIG CHANGE IN LOGIC: Tuesday July 08, 0;55
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // OLD: if ( same_id_cd_list.size() > 1 ) {
            if ( same_id_cd_list.size() > 0 ) {
                if ( !ignore_domains_specific_to_one_genome || ( same_id_cd_list.size() > 1 ) ) {
                    final DomainSimilarity s = calculateSimilarity( pairwise_calculator, same_id_cd_list );
                    if ( s != null ) {
                        similarities.add( s );
                    }
                    else {
                        throw new IllegalStateException( "similarity is null: this should not have happened" );
                    }
                }
            }
            // ~~~ NEW:
            else {
                throw new IllegalStateException( "this should not have happened" );
            }
            // ~~~ OLD:
            // else if ( same_id_cd_list.size() == 1 ) {
            // TODO need to go in file
            // System.out.println( "only in one species [" +
            // species_with_key_id_domain.get( 0 ) + "]: " + key_id );
            //}
            //else {
            //    throw new IllegalStateException( "this should not have happened" );
            // }
        }
        return similarities;
    }

    private DomainSimilarity calculateSimilarity( final PairwiseDomainSimilarityCalculator pairwise_calculator,
                                                  final List<CombinableDomains> domains_list ) {
        if ( domains_list.size() == 1 ) {
            // BIG CHANGE IN LOGIC: Tuesday July 08, 0;55
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // ~~~OLD:
            //throw new IllegalArgumentException( "attempt to calculate multiple combinable domains similarity for less than two combinable domains" );
            // ~~~new: 
            final SortedMap<Species, SpeciesSpecificDomainSimilariyData> species_data = new TreeMap<Species, SpeciesSpecificDomainSimilariyData>();
            species_data.put( domains_list.get( 0 ).getSpecies(),
                              createSpeciesSpecificDomainSimilariyData( domains_list.get( 0 ) ) );
            return new PrintableDomainSimilarity( domains_list.get( 0 ),
                                                  1.0,
                                                  1.0,
                                                  1.0,
                                                  1.0,
                                                  0.0,
                                                  0,
                                                  0,
                                                  0,
                                                  species_data,
                                                  getSort(),
                                                  isSortBySpeciesCountFirst(),
                                                  isTreatAsBinaryComparison() );
        }
        final DescriptiveStatistics stat = new BasicDescriptiveStatistics();
        final SortedMap<Species, SpeciesSpecificDomainSimilariyData> species_data = new TreeMap<Species, SpeciesSpecificDomainSimilariyData>();
        species_data.put( domains_list.get( 0 ).getSpecies(), createSpeciesSpecificDomainSimilariyData( domains_list
                .get( 0 ) ) );
        int max_difference_in_counts = 0;
        int max_difference = 0;
        final boolean is_domain_combination_based = pairwise_calculator instanceof CombinationsBasedPairwiseDomainSimilarityCalculator;
        for( int i = 1; i < domains_list.size(); ++i ) {
            species_data.put( domains_list.get( i ).getSpecies(),
                              createSpeciesSpecificDomainSimilariyData( domains_list.get( i ) ) );
            final CombinableDomains domains_i = domains_list.get( i );
            for( int j = 0; j < i; ++j ) {
                final PairwiseDomainSimilarity pairwise_similarity = pairwise_calculator
                        .calculateSimilarity( domains_i, domains_list.get( j ) );
                final int difference_in_counts = pairwise_similarity.getDifferenceInCounts();
                int difference = 0;
                if ( is_domain_combination_based ) {
                    difference = ( ( CombinationsBasedPairwiseDomainSimilarity ) pairwise_similarity )
                            .getNumberOfDifferentDomains();
                }
                else {
                    difference = difference_in_counts;
                }
                if ( Math.abs( difference_in_counts ) > Math.abs( max_difference_in_counts ) ) {
                    max_difference_in_counts = difference_in_counts;
                }
                if ( Math.abs( difference ) > Math.abs( max_difference ) ) {
                    max_difference = difference;
                }
                stat.addValue( pairwise_similarity.getSimilarityScore() );
            }
        }
        if ( stat.getN() < 1 ) {
            throw new AssertionError( "empty descriptive statistics: this should not have happened" );
        }
        if ( ( stat.getN() != 1 ) && isTreatAsBinaryComparison() ) {
            throw new IllegalArgumentException( "attmpt to treat similarity with N not equal to one as binary comparison" );
        }
        if ( ( /*stat.getN() != 1 ||*/!isTreatAsBinaryComparison() ) && ( max_difference_in_counts < 0 ) ) {
            max_difference_in_counts = Math.abs( max_difference_in_counts );
            if ( !is_domain_combination_based ) {
                max_difference = Math.abs( max_difference ); //=max_difference_in_counts for !is_domain_combination_based.
            }
        }
        DomainSimilarity similarity = null;
        if ( stat.getN() == 1 ) {
            similarity = new PrintableDomainSimilarity( domains_list.get( 0 ),
                                                        stat.getMin(),
                                                        stat.getMax(),
                                                        stat.arithmeticMean(),
                                                        stat.median(),
                                                        0.0,
                                                        stat.getN(),
                                                        max_difference_in_counts,
                                                        max_difference,
                                                        species_data,
                                                        getSort(),
                                                        isSortBySpeciesCountFirst(),
                                                        isTreatAsBinaryComparison() );
        }
        else {
            similarity = new PrintableDomainSimilarity( domains_list.get( 0 ),
                                                        stat.getMin(),
                                                        stat.getMax(),
                                                        stat.arithmeticMean(),
                                                        stat.median(),
                                                        stat.sampleStandardDeviation(),
                                                        stat.getN(),
                                                        max_difference_in_counts,
                                                        max_difference,
                                                        species_data,
                                                        getSort(),
                                                        isSortBySpeciesCountFirst(),
                                                        isTreatAsBinaryComparison() );
        }
        return similarity;
    }

    private DomainSimilarity.DomainSimilaritySortField getSort() {
        return _sort;
    }

    private boolean isSortBySpeciesCountFirst() {
        return _sort_by_species_count_first;
    }

    private boolean isTreatAsBinaryComparison() {
        return _treat_as_binary_comparison;
    }

    private static SpeciesSpecificDomainSimilariyData createSpeciesSpecificDomainSimilariyData( final CombinableDomains cd ) {
        final SpeciesSpecificDomainSimilariyData sd = new PrintableSpeciesSpecificDomainSimilariyData( cd
                .getKeyDomainProteinsCount(), cd.getKeyDomainCount(), cd.getNumberOfCombinableDomains(), cd
                .getKeyDomainConfidenceDescriptiveStatistics() );
        for( final DomainId domain : cd.getCombinableDomains() ) {
            sd.addProteinsExhibitingCombinationCount( domain, cd.getNumberOfProteinsExhibitingCombination( domain ) );
        }
        return sd;
    }
}
