// $Id: ORcount.java,v 1.20 2009/11/20 22:22:10 cmzmasek Exp $
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

package org.forester.sdi;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/*
 * Allows to <ul> <li> <li> <li> </ul>
 * 
 * @see SDIse
 * 
 * @see SDI
 * 
 * @author Christian M. Zmasek
 * 
 * @version 1.400 -- last modified: 10/29/2005
 */
public class ORcount {

    private static final String[]                     group_1              = { "ANOGA", "DROME", "CAEBR", "CAEEL" };
    private static final String[]                     group_2              = { "CIOIN", "FUGRU", "MOUSE", "RAT",
            "HUMAN"                                                       };
    private static final String[]                     all_species          = { "ANOGA", "DROME", "CAEBR", "CAEEL",
            "CIOIN", "FUGRU", "MOUSE", "RAT", "HUMAN"                     };
    private final Phylogeny[]                         _trees;
    private HashMap<String, HashMap<Object, Integer>> _species             = null;
    private ArrayList<String>                         _names               = null;
    private int                                       _group1_vs_2_counter = 0;

    /**
     * Default contructor which
     */
    public ORcount( final Phylogeny[] trees ) {
        _trees = trees;
    } // ORcount( final Phylogeny[] trees )

    private void count( final PhylogenyNode node ) {
        final List<PhylogenyNode> external_nodes = node.getAllExternalDescendants();
        for( int i = 1; i < external_nodes.size(); ++i ) {
            for( int j = 0; j < i; ++j ) {
                final PhylogenyNode node_i = external_nodes.get( i );
                final PhylogenyNode node_j = external_nodes.get( j );
                final String si = PhylogenyMethods.getSpecies( node_i );
                final String sj = PhylogenyMethods.getSpecies( node_j );
                count( si, sj, node_i.getNodeName(), node_j.getNodeName() );
            }
        }
    } // count( PhylogenyNode )

    private void count( final String a, final String b, final String seq_name_a, final String seq_name_b ) {
        HashMap<Object, Integer> h1 = _species.get( a );
        if ( h1 == null ) {
            throw new RuntimeException( "Unexpected error: Species \"" + a + "\" not present in species matrix." );
        }
        Object h2 = h1.get( b );
        String species_in_h1 = b;
        // We only look at the half matrix, and we do not know/care about the
        // order
        // of the keys (species).
        if ( h2 == null ) {
            h1 = _species.get( b );
            if ( h1 == null ) {
                throw new RuntimeException( "Unexpected error: Species \"" + b + "\" not present in species matrix." );
            }
            h2 = h1.get( a );
            species_in_h1 = a;
        }
        if ( h2 == null ) {
            throw new RuntimeException( "Unexpected error: Species \"" + a + "\" not present in species matrix." );
        }
        h1.put( species_in_h1, new Integer( ( ( Integer ) h2 ).intValue() + 1 ) );
        _names.add( a + "-" + seq_name_a + " = " + b + "-" + seq_name_b );
    } // count( String, String )

    public void countSharedAncestralClades( final Phylogeny tree,
                                            final int bootstrap_threshold,
                                            final String[] group_1,
                                            final String[] group_2 ) {
        if ( ( group_1 == null ) || ( group_2 == null ) ) {
            throw new IllegalArgumentException( "String[](s) in arguments to method \"ORcount.countSharedAncestralClades\" is (are) null." );
        }
        if ( !tree.isRooted() ) {
            throw new IllegalArgumentException( "Phylogeny must be rooted in order to count shared ancestral clades." );
        }
        final PhylogenyNodeIterator it = tree.iteratorPostorder();
        tree.setIndicatorsToZero();
        while ( it.hasNext() ) {
            final PhylogenyNode current_node = it.next();
            if ( current_node.getNumberOfDescendants() != 2 ) {
                throw new IllegalArgumentException( "Phylogeny can not contain multifurcations in order to count shared ancestral clades." );
            }
            if ( !current_node.isExternal() ) {
                final PhylogenyNode child1 = current_node.getChildNode1();
                final PhylogenyNode child2 = current_node.getChildNode2();
                if ( ( child1.getIndicator() == 1 ) || ( child2.getIndicator() == 1 ) ) {
                    current_node.setIndicator( ( byte ) 1 );
                }
                else {
                    final List<PhylogenyNode> external_nodes = current_node.getAllExternalDescendants();
                    final String[] external_species = new String[ external_nodes.size() ];
                    for( int i = 0; i < external_nodes.size(); ++i ) {
                        final PhylogenyNode n = external_nodes.get( i );
                        external_species[ i ] = PhylogenyMethods.getSpecies( n ).trim().toUpperCase();
                    }
                    if ( ForesterUtil.isIntersecting( external_species, group_1 )
                            && ForesterUtil.isIntersecting( external_species, group_2 ) ) {
                        current_node.setIndicator( ( byte ) 1 );
                        if ( ( group_1.length == 1 ) && ( group_2.length == 1 ) ) {
                            count( group_1[ 0 ], group_2[ 0 ], "name a", "name b" );
                        }
                        else {
                            increaseGroup1Vs2Counter();
                        }
                    }
                }
            }
        } // while
    } // countSharedAncestralClades( Phylogeny, int )

    public void countSharedAncestralClades( final Phylogeny[] trees, final int bootstrap_threshold ) {
        for( int i = 1; i < ORcount.all_species.length; ++i ) {
            for( int j = 0; j < i; ++j ) {
                final String all_i = ORcount.all_species[ i ].trim().toUpperCase();
                final String all_j = ORcount.all_species[ j ].trim().toUpperCase();
                final String[] a = { all_i };
                final String[] b = { all_j };
                for( int k = 0; k < trees.length; ++k ) {
                    countSharedAncestralClades( trees[ k ], bootstrap_threshold, a, b );
                }
            }
        }
        // print();
        if ( ( ORcount.group_1 != null ) && ( ORcount.group_2 != null ) && ( ORcount.group_1.length > 0 )
                && ( ORcount.group_2.length > 0 ) ) {
            setGroup1Vs2Counter( 0 );
            for( int k = 0; k < trees.length; ++k ) {
                countSharedAncestralClades( trees[ k ], bootstrap_threshold, ORcount.group_1, ORcount.group_2 );
            }
            System.out.println( "\nCount [(" + ForesterUtil.stringArrayToString( ORcount.group_1 ) + ") vs ("
                    + ForesterUtil.stringArrayToString( ORcount.group_2 ) + ")] = " + getGroup1Vs2Counter() );
        }
    }

    public void countSuperOrthologousRelations( final int bootstrap_threshold ) {
        reset();
        for( int i = 0; i < _trees.length; ++i ) {
            countSuperOrthologousRelations( _trees[ i ], bootstrap_threshold );
        }
    }

    private void countSuperOrthologousRelations( final Phylogeny tree, final int bootstrap_threshold ) {
        final PhylogenyNodeIterator it = tree.iteratorPostorder();
        if ( !tree.isRooted() ) {
            throw new IllegalArgumentException( "Phylogeny must be rooted in order to count 1:1 orthologous relationships." );
        }
        // The purpose of this is to find all substrees
        // which contain only speciation events on all their nodes.
        // All nodes in these subtrees are "painted" with 0's, wheres
        // the rest od the nodes in painted with 1's.
        tree.setIndicatorsToZero();
        it.reset();
        while ( it.hasNext() ) {
            final PhylogenyNode current_node = it.next();
            if ( current_node.getNumberOfDescendants() != 2 ) {
                throw new IllegalArgumentException( "Phylogeny can not contain multifurcations in order to count 1:1 orthologous relationships." );
            }
            if ( !current_node.isExternal() && !current_node.isHasAssignedEvent() ) {
                throw new IllegalArgumentException( "All nodes must have duplication or speciation assigned in order to count 1:1 orthologous relationships." );
            }
            if ( !current_node.isExternal()
                    && ( current_node.isDuplication() || ( current_node.getChildNode1().getIndicator() == 1 ) || ( current_node
                            .getChildNode2().getIndicator() == 1 ) ) ) {
                current_node.setIndicator( ( byte ) 1 );
            }
        }
        // These find the largest subtrees containing only speciations
        // and uses their largest nodes to count all possible species
        // combinations
        // in their extant external nodes.
        // ~~~ this could possibly be combined with the first iteration ~~
        // <<<<<<<<<<<~~~~~~~~~~~~~~~<<<<<<<<<<<<<<<
        it.reset();
        while ( it.hasNext() ) {
            final PhylogenyNode current_node = it.next();
            if ( !current_node.isExternal()
                    && ( current_node.getIndicator() == 0 )
                    && ( current_node.isRoot() || ( current_node.getParent().getIndicator() == 1 ) )
                    && ( ( bootstrap_threshold < 1 ) || ( ( PhylogenyMethods.getConfidenceValue( current_node ) >= bootstrap_threshold )
                            && ( PhylogenyMethods.getConfidenceValue( current_node.getChildNode1() ) >= bootstrap_threshold ) && ( PhylogenyMethods
                            .getConfidenceValue( current_node.getChildNode2() ) >= bootstrap_threshold ) ) ) ) {
                count( current_node );
            }
        }
    } // countOneToOneOrthologs( Phylogeny, int )

    // This puts all the species found in Phylogeny array _trees into
    // species HashMap.
    private void getAllSpecies() {
        if ( ( getTrees() == null ) || ( getTrees().length < 1 ) ) {
            throw new RuntimeException( "Phylogeny array in method \"getAllSpecies( HashMap hash )\" is null or empty." );
        }
        setSpecies( new HashMap<String, HashMap<Object, Integer>>() );
        for( int i = 0; i < getTrees().length; ++i ) {
            PhylogenyNode node = getTrees()[ i ].getFirstExternalNode();
            while ( node != null ) {
                getSpecies().put( PhylogenyMethods.getSpecies( node ), null );
                node = node.getNextExternalNode();
            }
        }
    } // void getAllSpecies( HashMap hash )

    private int getGroup1Vs2Counter() {
        return _group1_vs_2_counter;
    }

    private HashMap<String, HashMap<Object, Integer>> getSpecies() {
        return _species;
    }

    private Phylogeny[] getTrees() {
        return _trees;
    }

    private void increaseGroup1Vs2Counter() {
        _group1_vs_2_counter++;
    }

    private void printCount() {
        if ( ( _species == null ) || ( _species.size() < 2 ) ) {
            throw new RuntimeException( "Species HashMap in method \"setUpCountingMatrix()\" is null or contains less than two species." );
        }
        final Object[] species_array = _species.keySet().toArray();
        final int s = species_array.length;
        for( int i = 0; i < s - 1; ++i ) {
            final String species = ( String ) species_array[ i ];
            System.out.println();
            System.out.println( species + ":" );
            final HashMap<?, ?> h = _species.get( species );
            // Setting up HashMaps linked to by hash (=_species)
            // Diagonals are ignored, only half the matrix is needed.
            for( int j = 1 + i; j < s; ++j ) {
                final String sp = ( String ) species_array[ j ];
                final int c = ( ( Integer ) h.get( sp ) ).intValue();
                System.out.println( species + "-" + sp + ": " + c );
            }
        }
    }

    private void printNames() {
        for( int i = 0; i < _names.size(); ++i ) {
            System.out.println( i + ": " + _names.get( i ) );
        }
    }

    public void reset() {
        getAllSpecies();
        setUpCountingMatrix();
        setGroup1Vs2Counter( 0 );
        _names = new ArrayList<String>();
    }

    private void setGroup1Vs2Counter( final int group1_vs_2_counter ) {
        _group1_vs_2_counter = group1_vs_2_counter;
    }

    private void setSpecies( final HashMap<String, HashMap<Object, Integer>> species ) {
        _species = species;
    }

    private void setUpCountingMatrix() {
        if ( ( getSpecies() == null ) || ( getSpecies().size() < 2 ) ) {
            throw new RuntimeException( "Species HashMap in method \"setUpCountingMatrix()\" is null or contains less than two species." );
        }
        final Object[] species_array = getSpecies().keySet().toArray();
        final int s = species_array.length;
        for( int i = 0; i < s; ++i ) {
            final String species = ( String ) species_array[ i ];
            final HashMap<Object, Integer> h = new HashMap<Object, Integer>();
            // Setting up HashMaps linked to by hash (=_species)
            // Diagonals are ignored, only half the matrix is needed.
            for( int j = 1 + i; j < s; ++j ) {
                h.put( species_array[ j ], new Integer( 0 ) );
            }
            getSpecies().put( species, h );
        }
    }

    private static void errorInCommandLine() {
        System.out.println( "\nORcount: Error in command line.\n" );
        System.out.println( "Usage: \"\"" );
        System.out.println( "\nOptions:" );
        System.out.println( " -" );
        System.out.println( "" );
        System.exit( -1 );
    } // errorInCommandLine()

    /**
     * Main method for this class.
     * <p>
     * (Last modified: 11/26/03)
     * 
     * @param args[1or2]
     *            gene tree file name (in NHX format with species names in
     *            species name fields and sequence names in sequence name
     *            fields; unless -n option is used)
     */
    public static void main( final String args[] ) {
        if ( args.length == 0 ) {
            ORcount.errorInCommandLine();
        }
        final Phylogeny[] trees = new Phylogeny[ args.length ];
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        for( int i = 0; i < trees.length; ++i ) {
            try {
                System.out.println( "Reading tree #" + i + "  [" + args[ i ] + "]" );
                final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( new File( args[ i ] ), true );
                trees[ i ] = factory.create( new File( args[ i ] ), pp )[ 0 ];
            }
            catch ( final Exception e ) {
                System.out.println( "\nFailed to read \"" + args[ i ] + "\". Terminating.\n" );
                System.exit( -1 );
            }
        }
        System.out.println( "Finished reading in trees.\n\n" );
        final ORcount or_count = new ORcount( trees );
        try {
            System.out.println( "\n\n\n\"1:1 ORTHOLOGOUS GENE PAIRS\":\n" );
            System.out.println( "\n\n\n\"SUPER ORTHOLOGOUS GENE PAIRS\":\n" );
            or_count.countSuperOrthologousRelations( 0 );
            or_count.printNames();
            or_count.printCount();
            // System.out.println( "\n\n\n\"SHARED ANCESTRAL CLADES\":\n");
            // or_count.reset();
            // or_count.countSharedAncestralClades( trees, 0 );
        }
        catch ( final Exception e ) {
            System.out.println( "\nException. Terminating.\n" );
            System.out.println( "\nException is: " + e + "\n" );
            e.printStackTrace();
            System.exit( -1 );
        }
        System.out.println( "\nDone." );
        System.exit( 0 );
    } // main ( String )
} // End of class ORcount.
