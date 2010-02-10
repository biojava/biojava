// $Id: DevelopmentTools.java,v 1.11 2009/10/28 19:11:22 cmzmasek Exp $
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

package org.forester.development;

import java.util.Random;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;

public final class DevelopmentTools {

    /**
     * Creates a completely unbalanced Phylogeny with i external nodes.
     * 
     * @return a newly created unbalanced Phylogeny
     */
    // public static Phylogeny createUnbalancedTree( int i ) {
    //
    // Phylogeny t1 = null;
    //
    // try {
    // PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
    // t1 = factory.create( ":S=", new SimpleNHXParser() );
    //            
    // t1.setRooted( true );
    //
    // for ( int j = 1; j < i; ++j ) {
    // t1.addNodeAndConnect( "", "" );
    // }
    // t1.setRoot( t1.getFirstExternalNode().getRoot() );
    // t1.calculateRealHeight();
    // }
    //
    // catch ( PhylogenyParserException e ) {
    // System.err
    // .println( "Unexpected exception during \"createUnbalancedTree\":" );
    // System.err.println( e.toString() );
    // System.exit( -1 );
    // }
    //
    // return t1;
    // }
    private DevelopmentTools() {
    }

    /**
     * Creates a completely balanced rooted phylogeny with a given number of levels and
     * children per node.
     * 
     * @param levels
     * @param children_per_node
     * @return a completely balanced rooted phylogeny
     */
    public static Phylogeny createBalancedPhylogeny( final int levels, final int children_per_node ) {
        final PhylogenyNode root = new PhylogenyNode();
        final Phylogeny p = new Phylogeny();
        p.setRoot( root );
        p.setRooted( true );
        DevelopmentTools.createBalancedPhylogenyRecursion( levels, children_per_node, root );
        return p;
    }

    private static void createBalancedPhylogenyRecursion( int current_level,
                                                          final int children_per_node,
                                                          final PhylogenyNode current_node ) {
        if ( current_level > 0 ) {
            --current_level;
            for( int i = 0; i < children_per_node; ++i ) {
                final PhylogenyNode new_node = new PhylogenyNode();
                current_node.addAsChild( new_node );
                DevelopmentTools.createBalancedPhylogenyRecursion( current_level, children_per_node, new_node );
            }
        }
    }

    /**
     * Sets the species name of the external Nodes of Phylogeny t to 1, 1+i, 2,
     * 2+i, 3, 3+i, .... Examples: i=2: 1, 3, 2, 4 i=4: 1, 5, 2, 6, 3, 7, 4, 8
     * i=8: 1, 9, 2, 10, 3, 11, 4, 12, ...
     */
    public static void intervalNumberSpecies( final Phylogeny t, final int i ) {
        if ( ( t == null ) || t.isEmpty() ) {
            return;
        }
        PhylogenyNode n = t.getFirstExternalNode();
        int j = 1;
        boolean odd = true;
        while ( n != null ) {
            if ( odd ) {
                PhylogenyMethods.setScientificName( n, j + "" );
            }
            else {
                PhylogenyMethods.setScientificName( n, ( j + i ) + "" );
                j++;
            }
            odd = !odd;
            n = n.getNextExternalNode();
        }
    }

    /**
     * Sets the species namea of the external Nodes of Phylogeny t to descending
     * integers, ending with 1.
     */
    public static void numberSpeciesInDescOrder( final Phylogeny t ) {
        if ( ( t == null ) || t.isEmpty() ) {
            return;
        }
        PhylogenyNode n = t.getFirstExternalNode();
        int j = t.getRoot().getNumberOfExternalNodes();
        while ( n != null ) {
            PhylogenyMethods.setTaxonomyCode( n, j + "" );
            j--;
            n = n.getNextExternalNode();
        }
    }

    /**
     * Sets the species namea of the external Nodes of Phylogeny t to ascending
     * integers, starting with 1.
     */
    public static void numberSpeciesInOrder( final Phylogeny t ) {
        if ( ( t == null ) || t.isEmpty() ) {
            return;
        }
        PhylogenyNode n = t.getFirstExternalNode();
        int j = 1;
        while ( n != null ) {
            PhylogenyMethods.setScientificName( n, j + "" );
            j++;
            n = n.getNextExternalNode();
        }
    }

    /**
     * Sets the species names of the external Nodes of Phylogeny t to a random
     * positive integer number between (and including) min and max.
     * 
     * @param t
     *            whose external species names are to be randomized
     * @param min
     *            minimal value for random numbers
     * @param max
     *            maximum value for random numbers
     */
    public static void randomizeSpecies( final int min, final int max, final Phylogeny t ) {
        if ( ( t == null ) || t.isEmpty() ) {
            return;
        }
        final int mi = Math.abs( min );
        final int ma = Math.abs( max );
        final Random r = new Random();
        PhylogenyNode n = t.getFirstExternalNode();
        while ( n != null ) {
            final String code = ( ( Math.abs( r.nextInt() ) % ( ma - mi + 1 ) ) + mi ) + "";
            PhylogenyMethods.setTaxonomyCode( n, code );
            n = n.getNextExternalNode();
        }
    }
}