// $Id: BasicExternalNodeBasedCoverageExtender.java,v 1.4 2008/09/08 22:07:16
// cmzmasek Exp $
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

package org.forester.pccx;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/*
 * @author Christian M. Zmasek
 */
public class BasicExternalNodeBasedCoverageExtender implements CoverageExtender {

    private String find( final CoverageCalculationOptions options,
                         final BranchCountingBasedScoringMethod scoring_method,
                         final List<SortedMap<PhylogenyNode, Double>> external_node_scores_list,
                         final List<SortedMap<PhylogenyNode, Double>> external_node_scores_list_temp,
                         final List<Phylogeny> phylogenies,
                         final Set<String> already_covered,
                         final PrintStream out,
                         final int i,
                         final double normalization_factor ) {
        final Phylogeny p = phylogenies.get( 0 );
        String best_name = null;
        double best_score = -Double.MAX_VALUE;
        for( final PhylogenyNodeIterator iter = p.iteratorExternalForward(); iter.hasNext(); ) {
            final String name = iter.next().getNodeName();
            if ( !already_covered.contains( name ) ) {
                final double score = BasicExternalNodeBasedCoverageExtender
                        .calculateCoverage( phylogenies,
                                            name,
                                            options,
                                            scoring_method,
                                            external_node_scores_list_temp,
                                            false );
                if ( score > best_score ) {
                    best_score = score;
                    best_name = name;
                }
            }
        }
        BasicExternalNodeBasedCoverageExtender.calculateCoverage( phylogenies,
                                                                  best_name,
                                                                  options,
                                                                  scoring_method,
                                                                  external_node_scores_list_temp,
                                                                  true );
        if ( out != null ) {
            out.println( i + "\t" + best_name + "\t" + ( best_score * normalization_factor ) );
        }
        return best_name;
    }

    /*
     * (non-Javadoc)
     * 
     * @see org.forester.tools.modeling.CoverageExtender#find(java.util.List,
     *      java.util.List, int,
     *      org.forester.tools.modeling.CoverageCalculationMethod,
     *      org.forester.tools.modeling.CoverageCalculationOptions,
     *      java.io.PrintStream)
     */
    public List<String> find( final List<Phylogeny> phylogenies,
                              final List<String> already_covered,
                              int number_names_to_find,
                              final CoverageCalculationOptions options,
                              final PrintStream out ) {
        final ExternalNodeBasedCoverageMethodOptions my_options = ( ExternalNodeBasedCoverageMethodOptions ) options;
        if ( ( my_options == null ) || ForesterUtil.isEmpty( my_options.getScoringMethod() ) ) {
            throw new IllegalArgumentException( "options for external node based coverage method appear to not have been set" );
        }
        BranchCountingBasedScoringMethod scoring_method;
        try {
            scoring_method = ( BranchCountingBasedScoringMethod ) ( Class.forName( my_options.getScoringMethod() ) )
                    .newInstance();
        }
        catch ( final Exception e ) {
            throw new IllegalArgumentException( "could not create scoring method class \""
                    + my_options.getScoringMethod() + "\"" );
        }
        final List<String> best_names = new ArrayList<String>();
        final Set<String> my_already_covered = new HashSet<String>();
        final List<SortedMap<PhylogenyNode, Double>> external_node_scores_list = new ArrayList<SortedMap<PhylogenyNode, Double>>();
        for( int i = 0; i < phylogenies.size(); ++i ) {
            external_node_scores_list.add( ModelingUtils.setUpExternalCoverageHashMap( phylogenies.get( i ) ) );
        }
        if ( already_covered != null ) {
            for( final String name : already_covered ) {
                my_already_covered.add( name );
                BasicExternalNodeBasedCoverageExtender.calculateCoverage( phylogenies,
                                                                          name,
                                                                          options,
                                                                          scoring_method,
                                                                          external_node_scores_list,
                                                                          true );
            }
        }
        if ( number_names_to_find < 1 ) {
            number_names_to_find = phylogenies.get( 0 ).getNumberOfExternalNodes() - my_already_covered.size();
        }
        final double normalization_factor = scoring_method.getNormalizationFactor( phylogenies.get( 0 ) );
        for( int i = 0; i < number_names_to_find; ++i ) {
            final String name = find( my_options,
                                      scoring_method,
                                      external_node_scores_list,
                                      external_node_scores_list,
                                      phylogenies,
                                      my_already_covered,
                                      out,
                                      i,
                                      normalization_factor );
            my_already_covered.add( name );
            best_names.add( name );
        }
        return best_names;
    }

    private static double calculateCoverage( final List<Phylogeny> phylogenies,
                                             final String name,
                                             final CoverageCalculationOptions options,
                                             final BranchCountingBasedScoringMethod scoring_method,
                                             final List<SortedMap<PhylogenyNode, Double>> external_node_scores_list,
                                             final boolean update_external_node_scores_list ) {
        int i = 0;
        double score_sum = 0.0;
        for( final Object element : phylogenies ) {
            SortedMap<PhylogenyNode, Double> external_node_scores;
            if ( update_external_node_scores_list ) {
                external_node_scores = external_node_scores_list.get( i++ );
            }
            else {
                external_node_scores = new TreeMap<PhylogenyNode, Double>( external_node_scores_list.get( i++ ) );
            }
            final Phylogeny phylogeny = ( Phylogeny ) element;
            scoring_method.calculateScoreForExternalNode( external_node_scores,
                                                          phylogeny,
                                                          phylogeny.getNode( name ),
                                                          options );
            for( final Object element2 : external_node_scores.values() ) {
                score_sum += ( ( Double ) element2 ).doubleValue();
            }
        }
        return score_sum / i;
    }
}
