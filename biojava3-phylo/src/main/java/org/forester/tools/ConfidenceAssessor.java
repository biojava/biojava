// $Id: ConfidenceAssessor.java,v 1.4 2009/12/17 02:28:12 cmzmasek Exp $
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

package org.forester.tools;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public final class ConfidenceAssessor {

    private ConfidenceAssessor() {
        // Hidden constructor.
    }

    public final static void evaluate( final String confidence_type,
                                       final Phylogeny[] evaluators,
                                       final Phylogeny target,
                                       final boolean strict,
                                       final double value ) {
        evaluate( confidence_type, evaluators, target, strict, value, 0, 0 );
    }

    public final static void evaluate( final String confidence_type,
                                       final Phylogeny[] evaluators,
                                       final Phylogeny target,
                                       final boolean strict,
                                       final double value,
                                       final int first,
                                       final int last ) {
        try {
            checkPreconditions( confidence_type, evaluators, target, value, first, last );
        }
        catch ( final IllegalArgumentException e ) {
            throw e;
        }
        boolean all = true;
        if ( ( first != 0 ) || ( last != 0 ) ) {
            all = false;
        }
        int counter = 0;
        final Map<PhylogenyNode, Set<PhylogenyNode>> node_to_ext_nodes_map = new HashMap<PhylogenyNode, Set<PhylogenyNode>>();
        for( final Phylogeny evaluator : evaluators ) {
            if ( all || ( ( counter >= first ) && ( counter <= last ) ) ) {
                if ( strict ) {
                    if ( evaluator.getNumberOfExternalNodes() != target.getNumberOfExternalNodes() ) {
                        throw new IllegalArgumentException( "evaluator #" + counter
                                + " does not have the same number of external nodes ["
                                + evaluator.getNumberOfExternalNodes() + "] than the corresponding target ["
                                + target.getNumberOfExternalNodes() + "]" );
                    }
                }
                final TreeSplitMatrix s = new TreeSplitMatrix( evaluator, strict, target );
                for( final PhylogenyNodeIterator it = target.iteratorPostorder(); it.hasNext(); ) {
                    final PhylogenyNode node = it.next();
                    if ( !node.isExternal() && !node.isRoot() ) {
                        if ( node.getParent().isRoot()
                                && ( target.getRoot().getNumberOfDescendants() == 2 )
                                && ( target.getRoot().getChildNode1().isExternal() || target.getRoot().getChildNode2()
                                        .isExternal() ) ) {
                            continue;
                        }
                        if ( !node_to_ext_nodes_map.containsKey( node ) ) {
                            addExternalNodesToMap( node_to_ext_nodes_map, node );
                        }
                        final Set<PhylogenyNode> ex_descs = node_to_ext_nodes_map.get( node );
                        if ( s.match( ex_descs ) ) {
                            final Confidence c = ConfidenceAssessor.getConfidence( node, confidence_type );
                            c.setValue( c.getValue() + value );
                        }
                    }
                }
            }
            ++counter;
        }
    }

    private final static void checkPreconditions( final String confidence_type,
                                                  final Phylogeny[] evaluators,
                                                  final Phylogeny target,
                                                  final double value,
                                                  final int first,
                                                  final int last ) {
        if ( ( first < 0 ) || ( last < 0 ) ) {
            throw new IllegalArgumentException( "attempt to set first or last evaluator topology to use to a number less than zero" );
        }
        if ( evaluators.length < 1 ) {
            throw new IllegalArgumentException( "need at least one evaluator topology" );
        }
        if ( ForesterUtil.isEmpty( confidence_type ) ) {
            throw new IllegalArgumentException( "attempt to use empty confidence type" );
        }
        if ( value <= 0 ) {
            throw new IllegalArgumentException( "attempt to use zero or negative \'count value\'" );
        }
        if ( ( first != 0 ) || ( last != 0 ) ) {
            if ( ( last >= evaluators.length ) || ( last <= first ) ) {
                throw new IllegalArgumentException( "illegal value for last evaluator topology to use" );
            }
        }
        for( final PhylogenyNodeIterator it = target.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            final List<Confidence> confidences = node.getBranchData().getConfidences();
            for( final Confidence confidence : confidences ) {
                if ( confidence.getType().equals( confidence_type ) ) {
                    throw new IllegalArgumentException( "confidence [" + confidence_type
                            + "] is already present in target" );
                }
            }
        }
    }

    private final static void addExternalNodesToMap( final Map<PhylogenyNode, Set<PhylogenyNode>> node_to_ext_nodes_map,
                                                     final PhylogenyNode node ) {
        final Set<PhylogenyNode> ex_descs = new HashSet<PhylogenyNode>();
        for( final PhylogenyNode n : node.getAllExternalDescendants() ) {
            if ( ex_descs.contains( n ) ) {
                throw new IllegalArgumentException( "node '" + n.toString() + "' of target is not unique" );
            }
            ex_descs.add( n );
        }
        node_to_ext_nodes_map.put( node, ex_descs );
    }

    private final static Confidence getConfidence( final PhylogenyNode n, final String confidence_type ) {
        final List<Confidence> confidences = n.getBranchData().getConfidences();
        Confidence match = null;
        for( final Confidence confidence : confidences ) {
            if ( confidence.getType().equals( confidence_type ) ) {
                if ( match != null ) {
                    throw new IllegalArgumentException( "confidence [" + confidence_type + "] is not unique" );
                }
                match = confidence;
            }
        }
        if ( match == null ) {
            match = new Confidence( 0, confidence_type );
            confidences.add( match );
        }
        return match;
    }
}
