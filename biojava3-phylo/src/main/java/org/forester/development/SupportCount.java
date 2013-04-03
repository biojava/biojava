// $Id: SupportCount.java,v 1.23 2009/10/26 23:29:40 cmzmasek Exp $
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/*
 * A simple class containing a static method to evaluate the topology of a given
 * phylogeny with a list of resampled phylogenies.
 * 
 * 
 * @author Christian M Zmasek
 */
public final class SupportCount {

    private SupportCount() {
    }

    public static double compare( final Phylogeny phylogeny,
                                  final Phylogeny evaluator_phylogeny,
                                  final boolean strip_evaluator_phylogeny,
                                  final boolean update_support_in_phylogeny,
                                  final boolean re_root ) {
        String[] seq_names_to_keep = null;
        if ( strip_evaluator_phylogeny ) {
            seq_names_to_keep = phylogeny.getAllExternalNodeNames();
            SupportCount.strip( seq_names_to_keep, evaluator_phylogeny );
        }
        if ( re_root ) {
            final String child0_name = phylogeny.getFirstExternalNode().getNodeName();
            phylogeny.reRoot( phylogeny.getNode( child0_name ) );
            evaluator_phylogeny.reRoot( evaluator_phylogeny.getNode( child0_name ) );
        }
        final Map<Integer, ArrayList<String>> phylogeny_external_names_per_node = SupportCount
                .extractExternalNamesPerNode( phylogeny );
        return ( SupportCount.compare( phylogeny,
                                       evaluator_phylogeny,
                                       phylogeny_external_names_per_node,
                                       update_support_in_phylogeny,
                                       -1 ) );
    }

    /**
     * 
     * Precondition: phylogeny and evaluator_phylogeny have to be rooted in the
     * same manner.
     * 
     * Returns a measure of the similarity ("average bootstrap similarity")
     * between the topologies of phylogeny and evaluator_phylogeny: (sum of
     * branches which divide phylogeny in a manner consitent with
     * evaluator_phylogeny)/sum of branches in phylogeny. Therefore, this
     * measure is 1.0 for indentical topologies and 0.0 for completely
     * incompatible topologies.
     * 
     * 
     * @param phylogeny
     * @param evaluator_phylogeny
     * @param external_names_per_node
     * @param update_support_in_phylogeny
     *            set to true to update support values in phylogeny, otherwise,
     *            just calculation of the "average bootstrap similarity"
     * @return a measure of the similarity ("average bootstrap similarity")
     *         between phylogeny and evaluator_phylogeny
     */
    private static double compare( final Phylogeny phylogeny,
                                   final Phylogeny evaluator_phylogeny,
                                   final Map<Integer, ArrayList<String>> phylogeny_external_names_per_node,
                                   final boolean update_support_in_phylogeny,
                                   final double similarity_threshold ) {
        int matching_branches = 0;
        int phylogeny_total_internal_branches = 0;
        for( final PhylogenyNodeIterator it = phylogeny.iteratorPostorder(); it.hasNext(); ) {
            if ( !it.next().isExternal() ) {
                ++phylogeny_total_internal_branches;
            }
        }
        final Map<PhylogenyNode, Double> support_values = new HashMap<PhylogenyNode, Double>();
        E: for( final PhylogenyNodeIterator evaluator_phylogeny_it = evaluator_phylogeny.iteratorPostorder(); evaluator_phylogeny_it
                .hasNext(); ) {
            final List<String> c1 = new ArrayList<String>();
            for( final Object element : evaluator_phylogeny_it.next().getAllExternalDescendants() ) {
                c1.add( ( ( PhylogenyNode ) element ).getNodeName() );
            }
            for( final Integer id : phylogeny_external_names_per_node.keySet() ) {
                final List<String> c2 = phylogeny_external_names_per_node.get( id );
                if ( ( c2.size() == c1.size() ) && c2.containsAll( c1 ) ) {
                    if ( c2.size() > 1 ) {
                        matching_branches++;
                    }
                    if ( update_support_in_phylogeny ) {
                        final PhylogenyNode node = phylogeny.getNode( id.intValue() );
                        double d = PhylogenyMethods.getConfidenceValue( node );
                        if ( d < 1.0 ) {
                            d = 1.0;
                        }
                        else {
                            ++d;
                        }
                        support_values.put( node, new Double( d ) );
                    }
                    continue E;
                }
            }
        }
        final double similarity = ( double ) matching_branches / phylogeny_total_internal_branches;
        if ( ( similarity_threshold < 0.0 ) || ( similarity >= similarity_threshold ) ) {
            for( final PhylogenyNode node : support_values.keySet() ) {
                double b = support_values.get( node ).doubleValue();
                if ( b < 0 ) {
                    b = 0.0;
                }
                PhylogenyMethods.setBootstrapConfidence( node, b );
            }
        }
        return similarity;
    }

    public static void count( final Phylogeny phylogeny,
                              final Phylogeny[] evaluator_phylogenies,
                              final boolean strip_evaluator_phylogenies,
                              final boolean verbose ) {
        SupportCount.count( phylogeny, evaluator_phylogenies, strip_evaluator_phylogenies, -1, verbose );
    }

    /**
     * This counts the support of topology phylogeny by the topologies in
     * phylogenies. If phylogenies contains topogies with names not present in
     * phylogeny, strip_phylogenies must be set to true. phylogeny must not
     * contain names not found in all phylogenies.
     * 
     * @param phylogeny
     *            the topology to be evaluated
     * @param evaluator_phylogenies
     *            the topologies used for evaluation
     * @param strip_evaluator_phylogenies
     *            set to true if phylogenies contains topologies with names not
     *            present in phylogeny
     */
    public static List<Phylogeny> count( final Phylogeny phylogeny,
                                         final Phylogeny[] evaluator_phylogenies,
                                         final boolean strip_evaluator_phylogenies,
                                         final double similarity_threshold,
                                         final boolean verbose ) {
        String[] seq_names_to_keep = null;
        final List<Phylogeny> evaluator_phylogenies_above_threshold = new ArrayList<Phylogeny>();
        if ( strip_evaluator_phylogenies ) {
            seq_names_to_keep = phylogeny.getAllExternalNodeNames();
        }
        final String child0_name = phylogeny.getFirstExternalNode().getNodeName();
        phylogeny.reRoot( phylogeny.getNode( child0_name ) );
        final Map<Integer, ArrayList<String>> phylogeny_external_names_per_node = SupportCount
                .extractExternalNamesPerNode( phylogeny );
        if ( verbose ) {
            System.out.println();
            System.out.println( "evaluator phylogeny #: similarity score (max is 1.0)" );
            System.out.println( "----------------------------------------------------" );
            System.out.println();
        }
        for( int i = 0; i < evaluator_phylogenies.length; ++i ) {
            final Phylogeny evaluator_phylogeny = evaluator_phylogenies[ i ];
            evaluator_phylogeny.reRoot( evaluator_phylogeny.getNode( child0_name ) );
            Phylogeny unstripped_evaluator_phylogeny = evaluator_phylogeny;
            if ( strip_evaluator_phylogenies ) {
                unstripped_evaluator_phylogeny = evaluator_phylogeny.copy();
                SupportCount.strip( seq_names_to_keep, evaluator_phylogeny );
                evaluator_phylogeny.orderAppearance( true ); // This is for
                // easer
                // comparison if
                // phylos are saved
                // to file.
            }
            final double s = SupportCount.compare( phylogeny,
                                                   evaluator_phylogenies[ i ],
                                                   phylogeny_external_names_per_node,
                                                   true,
                                                   similarity_threshold );
            if ( ( similarity_threshold < 0.0 ) || ( s >= similarity_threshold ) ) {
                unstripped_evaluator_phylogeny.orderAppearance( true );
                evaluator_phylogenies_above_threshold.add( unstripped_evaluator_phylogeny );
            }
            if ( verbose ) {
                if ( similarity_threshold < 0.0 ) {
                    System.out.println( i + ": " + s );
                }
                else if ( s >= similarity_threshold ) {
                    System.out.println( i + ": " + s + " <====" );
                }
                else {
                    System.out.println( i + ": " + s );
                }
            }
        }
        if ( verbose ) {
            System.out.println( "----------------------------------------------------" );
            System.out.println();
        }
        return evaluator_phylogenies_above_threshold;
    }

    private static Map<Integer, ArrayList<String>> extractExternalNamesPerNode( final Phylogeny phylogeny )
            throws NoSuchElementException {
        final HashMap<Integer, ArrayList<String>> phylogeny_external_names_per_node = new HashMap<Integer, ArrayList<String>>();
        for( final PhylogenyNodeIterator it = phylogeny.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final List<PhylogenyNode> l = n.getAllExternalDescendants();
            final ArrayList<String> c = new ArrayList<String>();
            phylogeny_external_names_per_node.put( new Integer( n.getNodeId() ), c );
            for( final PhylogenyNode phylogenyNode : l ) {
                c.add( phylogenyNode.getNodeName() );
            }
        }
        return phylogeny_external_names_per_node;
    }

    private static void strip( final String[] to_keep, final Phylogeny to_be_stripped ) {
        PhylogenyMethods.deleteExternalNodesPositiveSelection( to_keep, to_be_stripped );
    }
}
