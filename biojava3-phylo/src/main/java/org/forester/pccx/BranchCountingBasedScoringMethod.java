// $Id: BranchCountingBasedScoringMethod.java,v 1.3 2008/03/09 00:11:18 cmzmasek
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

import java.util.SortedMap;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/*
 * Scoring method according to an idea by Adam Godzik, PhD.
 * 
 * @author Christian M. Zmasek
 */
public class BranchCountingBasedScoringMethod implements ScoringMethodForExternalNode {

    double calculateScoreContributionPerExternalNode( final PhylogenyNode external_node,
                                                      final PhylogenyNode current_node ) {
        double score_contribution = 0.0;
        if ( current_node == external_node ) {
            score_contribution = 1.0;
        }
        else {
            score_contribution = 1.0 / ModelingUtils.calculateBranchSum( external_node, current_node );
        }
        return score_contribution;
    }

    public void calculateScoreForExternalNode( final SortedMap<PhylogenyNode, Double> external_node_scores,
                                               final Phylogeny phylogeny,
                                               final PhylogenyNode external_node,
                                               final CoverageCalculationOptions options ) {
        for( final Object element : external_node_scores.keySet() ) {
            final PhylogenyNode current_node = ( PhylogenyNode ) element;
            final double score_contribution = calculateScoreContributionPerExternalNode( external_node, current_node );
            final double prev_score_contribution = external_node_scores.get( current_node );
            if ( score_contribution > prev_score_contribution ) {
                external_node_scores.put( current_node, score_contribution );
            }
        }
    }

    public String getDesciption() {
        return "sum of 1/branch-segment-sum";
    }

    public double getNormalizationFactor( final Phylogeny phylogeny ) {
        return ( 1.0 / phylogeny.getNumberOfExternalNodes() );
    }
}
