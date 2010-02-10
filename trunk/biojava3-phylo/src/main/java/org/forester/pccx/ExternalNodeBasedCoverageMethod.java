// $Id: ExternalNodeBasedCoverageMethod.java,v 1.5 2008/03/19 01:57:01 cmzmasek
// Exp $
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

import java.awt.Color;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

/*
 * @author Christian M. Zmasek
 */
public class ExternalNodeBasedCoverageMethod implements CoverageCalculationMethod {

    private static final Color MEAN_COVERAGE_COLOR = new Color( 0, 0, 0 );
    private static final Color MAXIMAL_COV_COLOR   = new Color( 0, 255, 0 );
    private static final Color MINIMAL_COV_COLOR   = new Color( 255, 0, 0 );

    public Coverage calculateCoverage( final List<Phylogeny> phylogenies,
                                       final List<String> names,
                                       final CoverageCalculationOptions options,
                                       final boolean annotate_phylogenies ) {
        final DescriptiveStatistics normalized_score_stats = new BasicDescriptiveStatistics();
        final DescriptiveStatistics raw_score_stats = new BasicDescriptiveStatistics();
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
        final double normalization_factor = scoring_method.getNormalizationFactor( phylogenies.get( 0 ) );
        for( final Object element : phylogenies ) {
            final double raw_score = calculateCoverage( ( Phylogeny ) element,
                                                        names,
                                                        options,
                                                        scoring_method,
                                                        annotate_phylogenies,
                                                        normalization_factor );
            normalized_score_stats.addValue( raw_score * normalization_factor );
            raw_score_stats.addValue( raw_score );
        }
        return new ExternalNodeBasedCoverage( normalized_score_stats, raw_score_stats.arithmeticMean(), options );
    }

    private double calculateCoverage( final Phylogeny phylogeny,
                                      final List<String> names,
                                      final CoverageCalculationOptions options,
                                      final BranchCountingBasedScoringMethod scoring_method,
                                      final boolean annotate_phylogeny,
                                      final double normalization_factor ) {
        final SortedMap<PhylogenyNode, Double> external_node_scores = ModelingUtils
                .setUpExternalCoverageHashMap( phylogeny );
        for( final Object element : names ) {
            scoring_method.calculateScoreForExternalNode( external_node_scores, phylogeny, phylogeny
                    .getNode( ( String ) element ), options );
        }
        if ( annotate_phylogeny ) {
            colorizePhylogenyAccordingToCoverage( external_node_scores, phylogeny, normalization_factor );
        }
        double score = 0.0;
        for( final Object element : external_node_scores.values() ) {
            score += ( ( Double ) element ).doubleValue();
        }
        return score;
    }

    private void colorizePhylogenyAccordingToCoverage( final SortedMap<PhylogenyNode, Double> external_node_scores,
                                                       final Phylogeny phylogeny,
                                                       final double normalization_factor ) {
        final DescriptiveStatistics ds = new BasicDescriptiveStatistics();
        for( final Object element : external_node_scores.entrySet() ) {
            ds.addValue( ( Double ) ( ( Map.Entry ) element ).getValue() * normalization_factor );
        }
        final double min = ds.getMin();
        final double max = ds.getMax();
        final double median = ds.median();
        for( final Object element2 : external_node_scores.entrySet() ) {
            final Map.Entry element = ( Map.Entry ) element2;
            final PhylogenyNode node = ( PhylogenyNode ) element.getKey();
            final double normalized_value = ( Double ) element.getValue() * normalization_factor;
            PhylogenyMethods.setBranchColorValue( node, ForesterUtil
                    .calcColor( normalized_value,
                                min,
                                max,
                                median,
                                ExternalNodeBasedCoverageMethod.MINIMAL_COV_COLOR,
                                ExternalNodeBasedCoverageMethod.MAXIMAL_COV_COLOR,
                                ExternalNodeBasedCoverageMethod.MEAN_COVERAGE_COLOR ) );
        }
        PhylogenyMethods.postorderBranchColorAveragingExternalNodeBased( phylogeny );
    }
}
