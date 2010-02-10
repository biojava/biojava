// $Id: PhylogenyBranch.java,v 1.14 2009/01/13 19:49:31 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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

package org.forester.phylogeny;

import org.forester.phylogeny.data.PhylogenyData;

/*
 * @author Christian M. Zmasek
 */
public class PhylogenyBranch implements Edge {

    private final PhylogenyNode _node_1;
    private final PhylogenyNode _node_2;
    private PhylogenyData       _data;
    private final boolean       _is_directed;
    private boolean             _towards_1;

    public PhylogenyBranch( final PhylogenyNode first_node, final PhylogenyNode second_node ) {
        if ( ( first_node == null ) || ( second_node == null ) ) {
            throw new IllegalArgumentException( "Attempt to create a branch with a null node" );
        }
        _node_1 = first_node;
        _node_2 = second_node;
        _is_directed = false;
    }

    public PhylogenyBranch( final PhylogenyNode first_node,
                            final PhylogenyNode second_node,
                            final boolean direction_towards_first ) {
        if ( ( first_node == null ) || ( second_node == null ) ) {
            throw new IllegalArgumentException( "Attempt to create a branch with a null node" );
        }
        _node_1 = first_node;
        _node_2 = second_node;
        _is_directed = true;
        _towards_1 = direction_towards_first;
    }

    @Override
    public boolean equals( final Object obj ) {
        if ( this == obj ) {
            return true;
        }
        if ( obj == null ) {
            return false;
        }
        if ( getClass() != obj.getClass() ) {
            return false;
        }
        final PhylogenyBranch other = ( PhylogenyBranch ) obj;
        return hashCode() == other.hashCode();
    }

    public PhylogenyNode getConnectedNode( final PhylogenyNode node ) throws IllegalArgumentException {
        if ( node == _node_1 ) {
            return _node_2;
        }
        else if ( node == _node_2 ) {
            return _node_1;
        }
        else {
            throw new IllegalArgumentException( "Attempt to get " + "connected node on branch with node which is "
                    + "not connected by the branch" );
        }
    }

    public PhylogenyData getData() {
        return _data;
    }

    public PhylogenyNode getFirstNode() {
        return _node_1;
    }

    public PhylogenyNode getSecondNode() {
        return _node_2;
    }

    @Override
    public int hashCode() {
        final int PRIME = 31;
        int result = 1;
        final int node_1_hc = _node_1.hashCode();
        final int node_2_hc = _node_2.hashCode();
        int hc_1 = 0;
        int hc_2 = 0;
        if ( !_is_directed ) {
            if ( node_1_hc > node_2_hc ) {
                hc_1 = node_2_hc;
                hc_2 = node_1_hc;
            }
            else {
                hc_1 = node_1_hc;
                hc_2 = node_2_hc;
            }
        }
        else {
            if ( _towards_1 ) {
                hc_1 = node_2_hc;
                hc_2 = node_1_hc;
            }
            else {
                hc_1 = node_1_hc;
                hc_2 = node_2_hc;
            }
        }
        result = PRIME * result + ( ( _data == null ) ? 0 : _data.hashCode() );
        result = PRIME * result + ( _is_directed ? 1231 : 1237 );
        result = PRIME * result + hc_1;
        result = PRIME * result + hc_2;
        return result;
    }

    public boolean isDirected() {
        return _is_directed;
    }

    public boolean isDirectionTowards( final PhylogenyNode node ) throws IllegalStateException {
        if ( !isDirected() ) {
            throw new IllegalStateException( "Attempt to get direction of undirected branch" );
        }
        return ( ( node == _node_1 ) && _towards_1 );
    }

    public void setDirectionTowards( final PhylogenyNode node ) {
        _towards_1 = node == _node_1;
    }

    @Override
    public String toString() {
        if ( isDirected() ) {
            if ( isDirectionTowards( getFirstNode() ) ) {
                return ( getSecondNode().getNodeName() + " -> " + getFirstNode().getNodeName() );
            }
            else {
                return ( getFirstNode().getNodeName() + " -> " + getSecondNode().getNodeName() );
            }
        }
        else {
            return ( getFirstNode().getNodeName() + " -- " + getSecondNode().getNodeName() );
        }
    }
}
