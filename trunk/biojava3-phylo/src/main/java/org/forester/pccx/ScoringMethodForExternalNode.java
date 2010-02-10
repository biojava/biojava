// $Id: ScoringMethodForExternalNode.java,v 1.3 2008/03/09 00:11:19 cmzmasek Exp
// $
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
 * Interface providing implementations of scoring methods used by
 * ExternalNodeBasedCoverageMethod.
 * 
 * @author Christian M. Zmasek
 */
public interface ScoringMethodForExternalNode {

    /**
     * This calculates the coverage score for one external node.
     * 
     * 
     * @param external_node_scores
     *            SortedMap<PhylogenyNode, Double> in which the external node
     *            scores are stored (node->score)
     * @param phylogeny
     *            Phylogeny containing the external nodes to score
     * @param external_node
     *            PhylogenyNod for which to calculate the score
     * @param options
     *            CoverageCalculationOptions
     * @param annotate_phylogeny           
     *            
     */
    public void calculateScoreForExternalNode( final SortedMap<PhylogenyNode, Double> external_node_scores,
                                               final Phylogeny phylogeny,
                                               final PhylogenyNode external_node,
                                               final CoverageCalculationOptions options );

    /**
     * This returns a short description of this scoring method
     * 
     * @return short description of this scoring method
     */
    public String getDesciption();

    /**
     * This calculates a normalization factor, so that a normalized score of 1.0
     * means complete coverage.
     * 
     * 
     * @param phylogeny
     *            Phylogeny containing the external nodes to score
     * @return normalization factor
     */
    public double getNormalizationFactor( final Phylogeny phylogeny );
}
