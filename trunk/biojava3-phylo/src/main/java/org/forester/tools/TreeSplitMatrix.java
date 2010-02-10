// $Id: TreeSplitMatrix.java,v 1.6 2009/12/16 03:07:10 cmzmasek Exp $
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class TreeSplitMatrix {

    private final SortedMap<PhylogenyNode, List<Boolean>> _data;
    private final Map<Integer, Integer>                   _positive_counts;
    private final boolean                                 _strict;

    public TreeSplitMatrix( final Phylogeny evaluator, final boolean strict, final Phylogeny target ) {
        Set<PhylogenyNode> target_external_nodes = null;
        if ( !strict ) {
            if ( ( target == null ) || target.isEmpty() ) {
                throw new IllegalArgumentException( "target must not be null or empty if non-strict evalution is expected" );
            }
            target_external_nodes = new HashSet<PhylogenyNode>();
            for( final PhylogenyNodeIterator it = target.iteratorExternalForward(); it.hasNext(); ) {
                target_external_nodes.add( it.next() );
            }
        }
        _data = new TreeMap<PhylogenyNode, List<Boolean>>();
        _positive_counts = new HashMap<Integer, Integer>();
        _strict = strict;
        decompose( evaluator, target_external_nodes );
    }

    /**
     * If strict is true, target nodes (all external nodes of the phylogeny for
     * which support values are to be calculated) is not used for anything during construction.
     * 
     * 
     * @param target
     * @param evaluator
     * @param strict
     */
    public TreeSplitMatrix( final Phylogeny evaluator,
                            final boolean strict,
                            final Set<PhylogenyNode> target_external_nodes ) {
        if ( !strict && ( ( target_external_nodes == null ) || target_external_nodes.isEmpty() ) ) {
            throw new IllegalArgumentException( "target nodes list must not be null or empty if non-strict evalution is expected" );
        }
        _data = new TreeMap<PhylogenyNode, List<Boolean>>();
        _positive_counts = new HashMap<Integer, Integer>();
        _strict = strict;
        decompose( evaluator, target_external_nodes );
    }

    private boolean contains( final PhylogenyNode node ) {
        return _data.keySet().contains( node );
    }

    private void decompose( final Phylogeny phy, final Set<PhylogenyNode> target_external_nodes ) {
        setUpKeys( phy, target_external_nodes );
        setUpValues( phy, target_external_nodes );
        sanityCheck();
    }

    private int getNumberOfTrueValuesAt( final int index ) {
        if ( _positive_counts.containsKey( index ) ) {
            return _positive_counts.get( index );
        }
        return 0;
    }

    private boolean getValue( final PhylogenyNode node, final int index ) {
        if ( _data.containsKey( node ) ) {
            return _data.get( node ).get( index );
        }
        return false;
    }

    private char getValueAsChar( final PhylogenyNode node, final int index ) {
        if ( getValue( node, index ) ) {
            return '.';
        }
        else {
            return ' ';
        }
    }

    private Set<PhylogenyNode> keySet() {
        return _data.keySet();
    }

    public boolean match( final Set<PhylogenyNode> query_nodes ) {
        if ( _strict ) {
            if ( !keySet().containsAll( query_nodes ) ) {
                throw new IllegalArgumentException( "external nodes of target and evaluator do not match" );
            }
        }
        for( int i = 0; i < size(); ++i ) {
            if ( match( query_nodes, i ) ) {
                return true;
            }
        }
        return false;
    }

    private boolean match( final Set<PhylogenyNode> query_nodes, final int i ) {
        final int counts = getNumberOfTrueValuesAt( i );
        final int q_counts = query_nodes.size();
        boolean positive_matches = true;
        boolean negative_matches = true;
        if ( q_counts != counts ) {
            positive_matches = false;
        }
        if ( q_counts != keySet().size() - counts ) {
            negative_matches = false;
        }
        if ( !positive_matches && !negative_matches ) {
            return false;
        }
        for( final PhylogenyNode query_node : query_nodes ) {
            if ( !contains( query_node ) ) {
                if ( _strict ) {
                    //TODO remove me after testing
                    throw new IllegalStateException( "this should not have happened, for query " + query_node + ":\n"
                            + toString() );
                }
                else {
                    return false;
                }
            }
            if ( getValue( query_node, i ) ) {
                negative_matches = false;
            }
            else {
                positive_matches = false;
            }
            if ( !positive_matches && !negative_matches ) {
                return false;
            }
        }
        return true;
    }

    private void sanityCheck() {
        int size = -1;
        for( final PhylogenyNode key : keySet() ) {
            if ( size < 0 ) {
                size = size( key );
            }
            else if ( size != size( key ) ) {
                throw new IllegalStateException( "this should not have happened: failed to build split matrix" );
            }
        }
    }

    private void setUpKeys( final Phylogeny phy, final Set<PhylogenyNode> target_external_nodes ) {
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( _strict || target_external_nodes.contains( n ) ) {
                if ( _data.containsKey( n ) ) {
                    throw new IllegalArgumentException( "node '" + n.toString() + "' of evaluator is not unique" );
                }
                _data.put( n, new ArrayList<Boolean>() );
            }
        }
    }

    private void setUpValues( final Phylogeny phy, final Set<PhylogenyNode> target_external_nodes ) {
        int index = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            final List<PhylogenyNode> current_ext_descs = node.getAllExternalDescendants();
            for( final PhylogenyNode key : keySet() ) {
                if ( _strict || target_external_nodes.contains( key ) ) {
                    if ( current_ext_descs.contains( key ) ) {
                        _data.get( key ).add( index, true );
                        if ( !_positive_counts.containsKey( index ) ) {
                            _positive_counts.put( index, 1 );
                        }
                        else {
                            _positive_counts.put( index, _positive_counts.get( index ) + 1 );
                        }
                    }
                    else {
                        _data.get( key ).add( index, false );
                    }
                }
            }
            index++;
        }
    }

    private int size() {
        for( final PhylogenyNode key : keySet() ) {
            return size( key );
        }
        return 0;
    }

    private int size( final PhylogenyNode node ) {
        return _data.get( node ).size();
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        for( final PhylogenyNode key : keySet() ) {
            sb.append( key.getNodeName() );
            sb.append( ":" );
            for( int i = 0; i < size( key ); ++i ) {
                sb.append( " " );
                sb.append( getValueAsChar( key, i ) );
            }
            sb.append( "\n" );
        }
        return sb.toString();
    }
}
