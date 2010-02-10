// $Id: BranchLengthBasedScoringMethod.java,v 1.3 2008/03/09 00:11:19 cmzmasek
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

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/*
 * 
 * @author Christian M. Zmasek
 */
public class BranchLengthBasedScoringMethod extends BranchCountingBasedScoringMethod {

    public static final double MIN_ALLOWED_BL_VALUE = 0.001;

    @Override
    double calculateScoreContributionPerExternalNode( final PhylogenyNode external_node,
                                                      final PhylogenyNode current_node ) {
        double score_contribution = 0.0;
        if ( current_node == external_node ) {
            score_contribution = external_node.getDistanceToParent();
            // This, of course, is completely /ad hoc/.
        }
        else {
            score_contribution = ModelingUtils.calculateBranchLengthSum( external_node, current_node );
        }
        return 1.0 / ( score_contribution > BranchLengthBasedScoringMethod.MIN_ALLOWED_BL_VALUE ? score_contribution
                : BranchLengthBasedScoringMethod.MIN_ALLOWED_BL_VALUE );
    }

    @Override
    public String getDesciption() {
        return "sum of 1/branch-length-sum [for self: 1/branch-length] [min branch length: "
                + BranchLengthBasedScoringMethod.MIN_ALLOWED_BL_VALUE + "]";
    }

    @Override
    public double getNormalizationFactor( final Phylogeny phylogeny ) {
        double s = 0.0;
        double d = 0.0;
        for( final PhylogenyNodeIterator iter = phylogeny.iteratorExternalForward(); iter.hasNext(); ) {
            d = iter.next().getDistanceToParent();
            s += ( 1.0 / ( d > BranchLengthBasedScoringMethod.MIN_ALLOWED_BL_VALUE ? d
                    : BranchLengthBasedScoringMethod.MIN_ALLOWED_BL_VALUE ) );
        }
        return 1.0 / s;
    }
}
