// $Id: NeighborJoining.java,v 1.15 2009/10/26 23:29:39 cmzmasek Exp $
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

package org.forester.phylogenyinference;

import java.util.ArrayList;
import java.util.List;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyNodeI;
import org.forester.util.ForesterUtil;

public class NeighborJoining {

    private final static boolean VERBOSE_DEFAULT = false;
    private DistanceMatrix       _d;
    private DistanceMatrix       _m;
    private double[]             _r;
    private int                  _n;
    private PhylogenyNodeI[]     _external_nodes;
    private int[]                _mappings;
    private boolean              _verbose;

    private NeighborJoining() {
        init();
    }

    private void calculateDistancesFromNewNode( final int otu1, final int otu2, final double d ) {
        for( int i = 0; i < _n; ++i ) {
            if ( ( i == otu1 ) || ( i == otu2 ) ) {
                continue;
            }
            final double nd = ( getValueFromD( otu1, i ) + getValueFromD( i, otu2 ) - d ) / 2;
            setValueInD( nd, otu1, i );
        }
    }

    private double calculateM( final int i, final int j ) {
        return getValueFromD( i, j ) - ( _r[ i ] + _r[ j ] ) / ( _n - 2 );
    }

    private void calculateNetDivergences() {
        for( int i = 0; i < _n; ++i ) {
            double d = 0.0;
            for( int n = 0; n < _n; ++n ) {
                d += getValueFromD( i, n );
            }
            _r[ i ] = d;
        }
    }

    public Phylogeny execute( final DistanceMatrix distance ) {
        reset( distance );
        // TODO this should work with a "calculation bases" phylo factory!
        final Phylogeny phylogeny = new Phylogeny();
        while ( _n > 2 ) {
            updateM();
            final int[] s = findMinimalDistance();
            final int otu1 = s[ 0 ];
            final int otu2 = s[ 1 ];
            // It is a condition that otu1 < otu2.
            if ( otu1 > otu2 ) {
                throw new AssertionError( "NJ code is faulty: otu1 > otu2" );
            }
            final PhylogenyNodeI node = new PhylogenyNode();// TODO there should
            // be
            // a
            // node factory, instead of
            // new!
            final double d = getValueFromD( otu1, otu2 );
            final double d1 = ( d / 2 ) + ( ( _r[ otu1 ] - _r[ otu2 ] ) / ( 2 * ( _n - 2 ) ) );
            final double d2 = d - d1;
            getExternalPhylogenyNode( otu1 ).setDistanceToParent( d1 );
            getExternalPhylogenyNode( otu2 ).setDistanceToParent( d2 );
            node.addAsChild( getExternalPhylogenyNode( otu1 ) );
            node.addAsChild( getExternalPhylogenyNode( otu2 ) );
            if ( isVerbose() ) {
                printProgress( otu1, otu2 );
            }
            calculateDistancesFromNewNode( otu1, otu2, d );
            setExternalPhylogenyNode( node, otu1 );
            updateMappings( otu2 );
            --_n;
        }
        final double d = getValueFromD( 0, 1 ) / 2;
        getExternalPhylogenyNode( 0 ).setDistanceToParent( d );
        getExternalPhylogenyNode( 1 ).setDistanceToParent( d );
        // a
        // node factory, instead of
        // new!
        final PhylogenyNodeI root = new PhylogenyNode();
        root.addAsChild( getExternalPhylogenyNode( 0 ) );
        root.addAsChild( getExternalPhylogenyNode( 1 ) );
        if ( isVerbose() ) {
            printProgress( 0, 1 );
        }
        phylogeny.setRoot( ( PhylogenyNode ) root );
        phylogeny.setRooted( false );
        return phylogeny;
    }

    public List<Phylogeny> execute( final List<DistanceMatrix> distances_list ) {
        final List<Phylogeny> pl = new ArrayList<Phylogeny>();
        for( final DistanceMatrix distances : distances_list ) {
            pl.add( execute( distances ) );
        }
        return pl;
    }

    private int[] findMinimalDistance() {
        // if more than one minimal distances, always the first found is
        // returned
        // i could randomize this, so that any would be returned in a randomized
        // fashion...
        double minimum = Double.MAX_VALUE;
        int otu_1 = -1;
        int otu_2 = -1;
        for( int j = 1; j < _n; ++j ) {
            for( int i = 0; i < j; ++i ) {
                if ( _m.getValue( i, j ) < minimum ) {
                    minimum = _m.getValue( i, j );
                    otu_1 = i;
                    otu_2 = j;
                }
            }
        }
        return new int[] { otu_1, otu_2 };
    }

    private PhylogenyNodeI getExternalPhylogenyNode( final int i ) {
        return _external_nodes[ _mappings[ i ] ];
    }

    private double getValueFromD( final int otu1, final int otu2 ) {
        return _d.getValue( _mappings[ otu1 ], _mappings[ otu2 ] );
    }

    private void init() {
        setVerbose( VERBOSE_DEFAULT );
    }

    private void initExternalNodes() {
        _external_nodes = new PhylogenyNodeI[ _n ];
        for( int i = 0; i < _n; ++i ) {
            _external_nodes[ i ] = new PhylogenyNode();
            // TODO there should be a node factory, instead of new!
            final String id = _d.getIdentifier( i );
            if ( id != null ) {
                _external_nodes[ i ].setName( id );
            }
            else {
                _external_nodes[ i ].setName( "" + i );
            }
            _mappings[ i ] = i;
        }
    }

    private boolean isVerbose() {
        return _verbose;
    }

    private void printProgress( final int otu1, final int otu2 ) {
        final PhylogenyNodeI n1 = getExternalPhylogenyNode( otu1 );
        final PhylogenyNodeI n2 = getExternalPhylogenyNode( otu2 );
        System.out.println( "Node " + ( ForesterUtil.isEmpty( n1.getNodeName() ) ? n1.getNodeId() : n1.getNodeName() )
                + " joins " + ( ForesterUtil.isEmpty( n2.getNodeName() ) ? n2.getNodeId() : n2.getNodeName() ) );
    }

    // only the values in the lower triangle are used.
    // !matrix values will be changed!
    private void reset( final DistanceMatrix distances ) {
        _n = distances.getSize();
        _d = distances;
        _m = new BasicSymmetricalDistanceMatrix( _n );
        _r = new double[ _n ];
        _mappings = new int[ _n ];
        initExternalNodes();
    }

    private void setExternalPhylogenyNode( final PhylogenyNodeI node, final int i ) {
        _external_nodes[ _mappings[ i ] ] = node;
    }

    private void setValueInD( final double d, final int otu1, final int otu2 ) {
        _d.setValue( _mappings[ otu1 ], _mappings[ otu2 ], d );
    }

    public void setVerbose( final boolean verbose ) {
        _verbose = verbose;
    }

    private void updateM() {
        calculateNetDivergences();
        for( int j = 1; j < _n; ++j ) {
            for( int i = 0; i < j; ++i ) {
                _m.setValue( i, j, calculateM( i, j ) );
            }
        }
    }

    // otu2 will, in effect, be "deleted" from the matrix.
    private void updateMappings( final int otu2 ) {
        for( int i = otu2; i < _mappings.length - 1; ++i ) {
            _mappings[ i ] = _mappings[ i + 1 ];
        }
    }

    public static NeighborJoining createInstance() {
        return new NeighborJoining();
    }
}
