// $Id: ModelingUtils.java,v 1.5 2009/10/26 23:29:39 cmzmasek Exp $
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
import java.util.TreeMap;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/*
 * @author Christian M. Zmasek
 */
public final class ModelingUtils {

    static double calculateBranchLengthSum( final PhylogenyNode n1, final PhylogenyNode n2 ) {
        final PhylogenyNode lca = PhylogenyMethods.getInstance().getLCA( n1, n2 );
        return ModelingUtils.calculateBranchLengthSumHelper( n1, lca )
                + ModelingUtils.calculateBranchLengthSumHelper( n2, lca );
    }

    private static double calculateBranchLengthSumHelper( final PhylogenyNode outer, final PhylogenyNode inner ) {
        PhylogenyNode my_outer = outer;
        double l = 0;
        while ( my_outer != inner ) {
            if ( my_outer.getDistanceToParent() > 0.0 ) {
                l += my_outer.getDistanceToParent();
            }
            my_outer = my_outer.getParent();
        }
        return l;
    }

    static int calculateBranchSum( final PhylogenyNode n1, final PhylogenyNode n2 ) {
        final PhylogenyNode lca = PhylogenyMethods.getInstance().getLCA( n1, n2 );
        return ModelingUtils.calculateBranchSumHelper( n1, lca ) + ModelingUtils.calculateBranchSumHelper( n2, lca );
    }

    private static int calculateBranchSumHelper( final PhylogenyNode outer, final PhylogenyNode inner ) {
        PhylogenyNode my_outer = outer;
        int s = 0;
        while ( my_outer != inner ) {
            s++;
            my_outer = my_outer.getParent();
        }
        return s;
    }

    static SortedMap<PhylogenyNode, Double> setUpExternalCoverageHashMap( final Phylogeny phylogeny ) {
        final SortedMap<PhylogenyNode, Double> external_node_coverage = new TreeMap<PhylogenyNode, Double>();
        for( final PhylogenyNodeIterator iter = phylogeny.iteratorExternalForward(); iter.hasNext(); ) {
            external_node_coverage.put( iter.next(), 0.0 );
        }
        return external_node_coverage;
    }
}
