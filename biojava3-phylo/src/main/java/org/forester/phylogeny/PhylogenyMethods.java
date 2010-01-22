// $Id: PhylogenyMethods.java,v 1.68 2009/10/28 19:11:22 cmzmasek Exp $
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

package org.forester.phylogeny;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class PhylogenyMethods {

    private static PhylogenyMethods  _instance      = null;
    private final Set<PhylogenyNode> _temp_hash_set = new HashSet<PhylogenyNode>();
    private PhylogenyNode            _farthest_1    = null;
    private PhylogenyNode            _farthest_2    = null;

    private PhylogenyMethods() {
        // Hidden constructor.
    }

    /**
     * Calculates the distance between PhylogenyNodes node1 and node2.
     * 
     * 
     * @param node1
     * @param node2
     * @return distance between node1 and node2
     */
    public double calculateDistance( final PhylogenyNode node1, final PhylogenyNode node2 ) {
        final PhylogenyNode lca = getLCA( node1, node2 );
        final PhylogenyNode n1 = node1;
        final PhylogenyNode n2 = node2;
        return ( PhylogenyMethods.getDistance( n1, lca ) + PhylogenyMethods.getDistance( n2, lca ) );
    }

    public double calculateFurthestDistance( final Phylogeny phylogeny ) {
        _farthest_1 = null;
        _farthest_2 = null;
        if ( phylogeny.getNumberOfExternalNodes() < 2 ) {
            return 0.0;
        }
        PhylogenyNode node_1 = null;
        PhylogenyNode node_2 = null;
        double farthest_d = -Double.MAX_VALUE;
        final PhylogenyMethods methods = PhylogenyMethods.getInstance();
        final List<PhylogenyNode> ext_nodes = phylogeny.getRoot().getAllExternalDescendants();
        for( int i = 1; i < ext_nodes.size(); ++i ) {
            for( int j = 0; j < i; ++j ) {
                final double d = methods.calculateDistance( ext_nodes.get( i ), ext_nodes.get( j ) );
                if ( d > farthest_d ) {
                    farthest_d = d;
                    node_1 = ext_nodes.get( i );
                    node_2 = ext_nodes.get( j );
                }
            }
        }
        _farthest_1 = node_1;
        _farthest_2 = node_2;
        return farthest_d;
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        throw new CloneNotSupportedException();
    }

    public PhylogenyNode getFarthestNode1() {
        return _farthest_1;
    }

    public PhylogenyNode getFarthestNode2() {
        return _farthest_2;
    }

    /**
     * Returns the LCA of PhylogenyNodes node1 and node2.
     * 
     * 
     * @param node1
     * @param node2
     * @return LCA of node1 and node2
     */
    public PhylogenyNode getLCA( final PhylogenyNode node1, final PhylogenyNode node2 ) {
        _temp_hash_set.clear();
        boolean done_with_1 = false;
        boolean done_with_2 = false;
        PhylogenyNode n1 = node1;
        PhylogenyNode n2 = node2;
        while ( !done_with_1 || !done_with_2 ) {
            if ( n1 == n2 ) {
                return n1;
            }
            if ( !done_with_1 ) {
                _temp_hash_set.add( n1 );
            }
            if ( !done_with_2 ) {
                _temp_hash_set.add( n2 );
            }
            if ( !done_with_1 ) {
                if ( !n1.isRoot() ) {
                    n1 = n1.getParent();
                    if ( _temp_hash_set.contains( n1 ) ) {
                        return n1;
                    }
                }
                else {
                    done_with_1 = true;
                }
            }
            if ( !done_with_2 ) {
                if ( !n2.isRoot() ) {
                    n2 = n2.getParent();
                    if ( _temp_hash_set.contains( n2 ) ) {
                        return n2;
                    }
                }
                else {
                    done_with_2 = true;
                }
            }
        }
        throw new IllegalArgumentException( "attempt to get LCA of two nodes which do not share a common root" );
    }

    /**
     * Returns all orthologs of the external PhylogenyNode n of this Phylogeny.
     * Orthologs are returned as List of node references.
     * <p>
     * PRECONDITION: This tree must be binary and rooted, and speciation -
     * duplication need to be assigned for each of its internal Nodes.
     * <p>
     * Returns null if this Phylogeny is empty or if n is internal.
     * @param n
     *            external PhylogenyNode whose orthologs are to be returned
     * @return Vector of references to all orthologous Nodes of PhylogenyNode n
     *         of this Phylogeny, null if this Phylogeny is empty or if n is
     *         internal
     */
    public List<PhylogenyNode> getOrthologousNodes( final Phylogeny phy, final PhylogenyNode node ) {
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        final PhylogenyNodeIterator it = phy.iteratorExternalForward();
        while ( it.hasNext() ) {
            final PhylogenyNode temp_node = it.next();
            if ( ( temp_node != node ) && isAreOrthologous( node, temp_node ) ) {
                nodes.add( temp_node );
            }
        }
        return nodes;
    }

    public boolean isAreOrthologous( final PhylogenyNode node1, final PhylogenyNode node2 ) {
        return !getLCA( node1, node2 ).isDuplication();
    }

    static double addPhylogenyDistances( final double a, final double b ) {
        if ( ( a >= 0.0 ) && ( b >= 0.0 ) ) {
            return a + b;
        }
        else if ( a >= 0.0 ) {
            return a;
        }
        else if ( b >= 0.0 ) {
            return b;
        }
        return PhylogenyNode.DISTANCE_DEFAULT;
    }

    // Helper for getUltraParalogousNodes( PhylogenyNode ).
    public static boolean areAllChildrenDuplications( final PhylogenyNode n ) {
        if ( n.isExternal() ) {
            return false;
        }
        else {
            if ( n.isDuplication() ) {
                //FIXME test me!
                for( final PhylogenyNode desc : n.getDescendants() ) {
                    if ( !areAllChildrenDuplications( desc ) ) {
                        return false;
                    }
                }
                return true;
            }
            else {
                return false;
            }
        }
    }

    public static int calculateDepth( final PhylogenyNode node ) {
        PhylogenyNode n = node;
        int steps = 0;
        while ( !n.isRoot() ) {
            steps++;
            n = n.getParent();
        }
        return steps;
    }

    public static double calculateDistanceToRoot( final PhylogenyNode node ) {
        PhylogenyNode n = node;
        double d = 0.0;
        while ( !n.isRoot() ) {
            if ( n.getDistanceToParent() > 0.0 ) {
                d += n.getDistanceToParent();
            }
            n = n.getParent();
        }
        return d;
    }

    public static short calculateMaxBranchesToLeaf( final PhylogenyNode node ) {
        if ( node.isExternal() ) {
            return 0;
        }
        short max = 0;
        for( PhylogenyNode d : node.getAllExternalDescendants() ) {
            short steps = 0;
            while ( d != node ) {
                if ( d.isCollapse() ) {
                    steps = 0;
                }
                else {
                    steps++;
                }
                d = d.getParent();
            }
            if ( max < steps ) {
                max = steps;
            }
        }
        return max;
    }

    public static int calculateMaxDepth( final Phylogeny phy ) {
        int max = 0;
        for( final PhylogenyNodeIterator iter = phy.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final int steps = calculateDepth( node );
            if ( steps > max ) {
                max = steps;
            }
        }
        return max;
    }

    public static double calculateMaxDistanceToRoot( final Phylogeny phy ) {
        double max = 0.0;
        for( final PhylogenyNodeIterator iter = phy.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final double d = calculateDistanceToRoot( node );
            if ( d > max ) {
                max = d;
            }
        }
        return max;
    }

    public static int calculateMaximumNumberOfDescendantsPerNode( final Phylogeny phy ) {
        int max = 0;
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.getNumberOfDescendants() > max ) {
                max = node.getNumberOfDescendants();
            }
        }
        return max;
    }

    /**
     * Returns the sum of different taxonomies of
     * all external nodes of node.
     * If at least one the external nodes has no taxonomy,
     * 0 is returned.
     * 
     */
    public static int calculateSumOfDistinctTaxonomies( final PhylogenyNode node ) {
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        final Set<Taxonomy> tax_set = new HashSet<Taxonomy>();
        for( final PhylogenyNode n : descs ) {
            if ( !n.getNodeData().isHasTaxonomy() || n.getNodeData().getTaxonomy().isEmpty() ) {
                return 0;
            }
            tax_set.add( n.getNodeData().getTaxonomy() );
        }
        return tax_set.size();
    }

    /**
     * Deep copies the phylogeny originating from this node.
     */
    static PhylogenyNode copySubTree( final PhylogenyNode source ) {
        if ( source == null ) {
            return null;
        }
        else {
            final PhylogenyNode newnode = source.copyNodeData();
            if ( !source.isExternal() ) {
                for( int i = 0; i < source.getNumberOfDescendants(); ++i ) {
                    newnode.setChildNode( i, PhylogenyMethods.copySubTree( source.getChildNode( i ) ) );
                }
            }
            return newnode;
        }
    }

    public static void deleteExternalNodesNegativeSelection( final Set<Integer> to_delete, final Phylogeny phy ) {
        phy.hashIDs();
        for( final Integer id : to_delete ) {
            phy.deleteSubtree( phy.getNode( id ), true );
        }
        phy.hashIDs();
    }

    public static void deleteExternalNodesNegativeSelection( final String[] node_names_to_delete, final Phylogeny p )
            throws IllegalArgumentException {
        for( int i = 0; i < node_names_to_delete.length; ++i ) {
            if ( ForesterUtil.isEmpty( node_names_to_delete[ i ] ) ) {
                continue;
            }
            List<PhylogenyNode> nodes = null;
            nodes = p.getNodes( node_names_to_delete[ i ] );
            final Iterator<PhylogenyNode> it = nodes.iterator();
            while ( it.hasNext() ) {
                final PhylogenyNode n = it.next();
                if ( !n.isExternal() ) {
                    throw new IllegalArgumentException( "attempt to delete non-external node \""
                            + node_names_to_delete[ i ] + "\"" );
                }
                p.deleteSubtree( n, true );
            }
        }
    }

    public static void deleteExternalNodesPositiveSelection( final Set<Taxonomy> species_to_keep, final Phylogeny phy ) {
        //   final Set<Integer> to_delete = new HashSet<Integer>();
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy() ) {
                if ( !species_to_keep.contains( n.getNodeData().getTaxonomy() ) ) {
                    //to_delete.add( n.getNodeId() );
                    phy.deleteSubtree( n, true );
                }
            }
            else {
                throw new IllegalArgumentException( "node " + n.getNodeId() + " has no taxonomic data" );
            }
        }
        phy.hashIDs();
        phy.externalNodesHaveChanged();
        //  deleteExternalNodesNegativeSelection( to_delete, phy );
    }

    public static void deleteExternalNodesPositiveSelection( final String[] node_names_to_keep, final Phylogeny p ) {
        final PhylogenyNodeIterator it = p.iteratorExternalForward();
        final String[] to_delete = new String[ p.getNumberOfExternalNodes() ];
        int i = 0;
        Arrays.sort( node_names_to_keep );
        while ( it.hasNext() ) {
            final String curent_name = it.next().getNodeName();
            if ( Arrays.binarySearch( node_names_to_keep, curent_name ) < 0 ) {
                to_delete[ i++ ] = curent_name;
            }
        }
        PhylogenyMethods.deleteExternalNodesNegativeSelection( to_delete, p );
    }

    public static Set<PhylogenyNode> getAllDescendants( final PhylogenyNode node ) {
        final Set<PhylogenyNode> descs = new HashSet<PhylogenyNode>();
        if ( !node.isExternal() ) {
            final List<PhylogenyNode> exts = node.getAllExternalDescendants();
            for( PhylogenyNode current : exts ) {
                descs.add( current );
                while ( current != node ) {
                    current = current.getParent();
                    if ( descs.contains( current ) ) {
                        continue;
                    }
                    descs.add( current );
                }
            }
        }
        return descs;
    }

    /**
     * 
     * Convenience method
     * 
     * @param node
     * @return
     */
    public static Color getBranchColorValue( final PhylogenyNode node ) {
        if ( node.getBranchData().getBranchColor() == null ) {
            return null;
        }
        return node.getBranchData().getBranchColor().getValue();
    }

    /**
     * Convenience method
     */
    public static double getBranchWidthValue( final PhylogenyNode node ) {
        if ( !node.getBranchData().isHasBranchWidth() ) {
            return BranchWidth.BRANCH_WIDTH_DEFAULT_VALUE;
        }
        return node.getBranchData().getBranchWidth().getValue();
    }

    /**
     * Convenience method
     */
    public static double getConfidenceValue( final PhylogenyNode node ) {
        if ( !node.getBranchData().isHasConfidences() ) {
            return Confidence.CONFIDENCE_DEFAULT_VALUE;
        }
        return node.getBranchData().getConfidence( 0 ).getValue();
    }

    /**
     * Convenience method
     */
    public static double[] getConfidenceValuesAsArray( final PhylogenyNode node ) {
        if ( !node.getBranchData().isHasConfidences() ) {
            return new double[ 0 ];
        }
        final double[] values = new double[ node.getBranchData().getConfidences().size() ];
        int i = 0;
        for( final Confidence c : node.getBranchData().getConfidences() ) {
            values[ i++ ] = c.getValue();
        }
        return values;
    }

    /**
     * Calculates the distance between PhylogenyNodes n1 and n2.
     * PRECONDITION: n1 is a descendant of n2.
     * 
     * @param n1
     *            a descendant of n2
     * @param n2
     * @return distance between n1 and n2
     */
    private static double getDistance( PhylogenyNode n1, final PhylogenyNode n2 ) {
        double d = 0.0;
        while ( n1 != n2 ) {
            if ( n1.getDistanceToParent() > 0.0 ) {
                d += n1.getDistanceToParent();
            }
            n1 = n1.getParent();
        }
        return d;
    }

    /**
     * Returns taxonomy t if all external descendants have 
     * the same taxonomy t, null otherwise.
     * 
     */
    public static Taxonomy getExternalDescendantsTaxonomy( final PhylogenyNode node ) {
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        Taxonomy tax = null;
        for( final PhylogenyNode n : descs ) {
            if ( !n.getNodeData().isHasTaxonomy() || n.getNodeData().getTaxonomy().isEmpty() ) {
                return null;
            }
            else if ( tax == null ) {
                tax = n.getNodeData().getTaxonomy();
            }
            else if ( n.getNodeData().getTaxonomy().isEmpty() || !tax.isEqual( n.getNodeData().getTaxonomy() ) ) {
                return null;
            }
        }
        return tax;
    }

    public static PhylogenyNode getFurthestDescendant( final PhylogenyNode node ) {
        final List<PhylogenyNode> children = node.getAllExternalDescendants();
        PhylogenyNode farthest = null;
        double longest = -Double.MAX_VALUE;
        for( final PhylogenyNode child : children ) {
            if ( PhylogenyMethods.getDistance( child, node ) > longest ) {
                farthest = child;
                longest = PhylogenyMethods.getDistance( child, node );
            }
        }
        return farthest;
    }

    public static PhylogenyMethods getInstance() {
        if ( PhylogenyMethods._instance == null ) {
            PhylogenyMethods._instance = new PhylogenyMethods();
        }
        return PhylogenyMethods._instance;
    }

    /**
     * Returns the largest confidence value found on phy.
     */
    static public double getMaximumConfidenceValue( final Phylogeny phy ) {
        double max = -Double.MAX_VALUE;
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final double s = PhylogenyMethods.getConfidenceValue( iter.next() );
            if ( ( s != Confidence.CONFIDENCE_DEFAULT_VALUE ) && ( s > max ) ) {
                max = s;
            }
        }
        return max;
    }

    /**
     * Convenience method for display purposes.
     * Not intended for algorithms.
     */
    public static String getSpecies( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasTaxonomy() ) {
            return "";
        }
        if ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
            return node.getNodeData().getTaxonomy().getTaxonomyCode();
        }
        else if ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) ) {
            return node.getNodeData().getTaxonomy().getScientificName();
        }
        else {
            return node.getNodeData().getTaxonomy().getCommonName();
        }
    }

    /**
     * Returns all Nodes which are connected to external PhylogenyNode n of this
     * Phylogeny by a path containing only speciation events. We call these
     * "super orthologs". Nodes are returned as Vector of references to Nodes.
     * <p>
     * PRECONDITION: This tree must be binary and rooted, and speciation -
     * duplication need to be assigned for each of its internal Nodes.
     * <p>
     * Returns null if this Phylogeny is empty or if n is internal.
     * @param n
     *            external PhylogenyNode whose strictly speciation related Nodes
     *            are to be returned
     * @return Vector of references to all strictly speciation related Nodes of
     *         PhylogenyNode n of this Phylogeny, null if this Phylogeny is
     *         empty or if n is internal
     */
    public static List<PhylogenyNode> getSuperOrthologousNodes( final PhylogenyNode n ) {
        // FIXME
        PhylogenyNode node = n, deepest = null;
        final List<PhylogenyNode> v = new ArrayList<PhylogenyNode>();
        if ( !node.isExternal() ) {
            return null;
        }
        while ( !node.isRoot() && !node.getParent().isDuplication() ) {
            node = node.getParent();
        }
        deepest = node;
        deepest.setIndicatorsToZero();
        do {
            if ( !node.isExternal() ) {
                if ( node.getIndicator() == 0 ) {
                    node.setIndicator( ( byte ) 1 );
                    if ( !node.isDuplication() ) {
                        node = node.getChildNode1();
                    }
                }
                if ( node.getIndicator() == 1 ) {
                    node.setIndicator( ( byte ) 2 );
                    if ( !node.isDuplication() ) {
                        node = node.getChildNode2();
                    }
                }
                if ( ( node != deepest ) && ( node.getIndicator() == 2 ) ) {
                    node = node.getParent();
                }
            }
            else {
                if ( node != n ) {
                    v.add( node );
                }
                if ( node != deepest ) {
                    node = node.getParent();
                }
                else {
                    node.setIndicator( ( byte ) 2 );
                }
            }
        } while ( ( node != deepest ) || ( deepest.getIndicator() != 2 ) );
        return v;
    }

    /**
     * Convenience method for display purposes.
     * Not intended for algorithms.
     */
    public static String getTaxonomyIdentifier( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasTaxonomy() || ( node.getNodeData().getTaxonomy().getIdentifier() == null ) ) {
            return "";
        }
        return node.getNodeData().getTaxonomy().getIdentifier().getValue();
    }

    /**
     * Returns all Nodes which are connected to external PhylogenyNode n of this
     * Phylogeny by a path containing, and leading to, only duplication events.
     * We call these "ultra paralogs". Nodes are returned as Vector of
     * references to Nodes.
     * <p>
     * PRECONDITION: This tree must be binary and rooted, and speciation -
     * duplication need to be assigned for each of its internal Nodes.
     * <p>
     * Returns null if this Phylogeny is empty or if n is internal.
     * <p>
     * (Last modified: 10/06/01)
     * 
     * @param n
     *            external PhylogenyNode whose ultra paralogs are to be returned
     * @return Vector of references to all ultra paralogs of PhylogenyNode n of
     *         this Phylogeny, null if this Phylogeny is empty or if n is
     *         internal
     */
    public static List<PhylogenyNode> getUltraParalogousNodes( final PhylogenyNode n ) {
        // FIXME test me
        PhylogenyNode node = n;
        if ( !node.isExternal() ) {
            return null;
        }
        while ( !node.isRoot() && node.getParent().isDuplication() && areAllChildrenDuplications( node.getParent() ) ) {
            node = node.getParent();
        }
        final List<PhylogenyNode> nodes = node.getAllExternalDescendants();
        nodes.remove( n );
        return nodes;
    }

    public static String inferCommonPartOfScientificNameOfDescendants( final PhylogenyNode node ) {
        final List<PhylogenyNode> descs = node.getDescendants();
        String sn = null;
        for( final PhylogenyNode n : descs ) {
            if ( !n.getNodeData().isHasTaxonomy()
                    || ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                return null;
            }
            else if ( sn == null ) {
                sn = n.getNodeData().getTaxonomy().getScientificName().trim();
            }
            else {
                String sn_current = n.getNodeData().getTaxonomy().getScientificName().trim();
                if ( !sn.equals( sn_current ) ) {
                    boolean overlap = false;
                    while ( ( sn.indexOf( ' ' ) >= 0 ) || ( sn_current.indexOf( ' ' ) >= 0 ) ) {
                        if ( ForesterUtil.countChars( sn, ' ' ) > ForesterUtil.countChars( sn_current, ' ' ) ) {
                            sn = sn.substring( 0, sn.lastIndexOf( ' ' ) ).trim();
                        }
                        else {
                            sn_current = sn_current.substring( 0, sn_current.lastIndexOf( ' ' ) ).trim();
                        }
                        if ( sn.equals( sn_current ) ) {
                            overlap = true;
                            break;
                        }
                    }
                    if ( !overlap ) {
                        return null;
                    }
                }
            }
        }
        return sn;
    }

    public static boolean isHasExternalDescendant( final PhylogenyNode node ) {
        for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            if ( node.getChildNode( i ).isExternal() ) {
                return true;
            }
        }
        return false;
    }

    private static boolean match( final String s,
                                  final String query,
                                  final boolean case_sensitive,
                                  final boolean partial ) {
        if ( ForesterUtil.isEmpty( s ) || ForesterUtil.isEmpty( query ) ) {
            return false;
        }
        String my_s = s.trim();
        String my_query = query.trim();
        if ( !case_sensitive ) {
            my_s = my_s.toLowerCase();
            my_query = my_query.toLowerCase();
        }
        if ( partial ) {
            return my_s.indexOf( my_query ) >= 0;
        }
        else {
            return my_s.equals( my_query );
        }
    }

    public static void midpointRoot( final Phylogeny phylogeny ) {
        if ( phylogeny.getNumberOfExternalNodes() < 2 ) {
            return;
        }
        final PhylogenyMethods methods = getInstance();
        final double farthest_d = methods.calculateFurthestDistance( phylogeny );
        final PhylogenyNode f1 = methods.getFarthestNode1();
        final PhylogenyNode f2 = methods.getFarthestNode2();
        if ( farthest_d <= 0.0 ) {
            return;
        }
        double x = farthest_d / 2.0;
        PhylogenyNode n = f1;
        if ( PhylogenyMethods.getDistance( f1, phylogeny.getRoot() ) < PhylogenyMethods.getDistance( f2, phylogeny
                .getRoot() ) ) {
            n = f2;
        }
        while ( x > n.getDistanceToParent() ) {
            x -= n.getDistanceToParent();
            n = n.getParent();
        }
        phylogeny.reRoot( n, x );
        phylogeny.recalculateNumberOfExternalDescendants( true );
        final PhylogenyNode a = getFurthestDescendant( phylogeny.getRoot().getChildNode1() );
        final PhylogenyNode b = getFurthestDescendant( phylogeny.getRoot().getChildNode2() );
        final double da = getDistance( a, phylogeny.getRoot() );
        final double db = getDistance( b, phylogeny.getRoot() );
        if ( Math.abs( da - db ) > 0.000001 ) {
            throw new IllegalStateException( " THIS SHOULD NOT HAVE HAPPENED \n midpoint rooting failed:  da=" + da
                    + ",  db=" + db + ",  diff=" + Math.abs( da - db ) );
        }
    }

    public static void normalizeBootstrapValues( final Phylogeny phylogeny,
                                                 final double max_bootstrap_value,
                                                 final double max_normalized_value ) {
        for( final PhylogenyNodeIterator iter = phylogeny.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isInternal() ) {
                final double confidence = getConfidenceValue( node );
                if ( confidence != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                    if ( confidence >= max_bootstrap_value ) {
                        setBootstrapConfidence( node, max_normalized_value );
                    }
                    else {
                        setBootstrapConfidence( node, ( confidence * max_normalized_value ) / max_bootstrap_value );
                    }
                }
            }
        }
    }

    public static Set<PhylogenyNode> obtainAllNodesAsSet( final Phylogeny phy ) {
        final Set<PhylogenyNode> nodes = new HashSet<PhylogenyNode>();
        if ( phy.isEmpty() ) {
            return nodes;
        }
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            nodes.add( iter.next() );
        }
        return nodes;
    }

    public static void postorderBranchColorAveragingExternalNodeBased( final Phylogeny p ) {
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            double red = 0.0;
            double green = 0.0;
            double blue = 0.0;
            int n = 0;
            if ( node.isInternal() ) {
                for( final PhylogenyNodeIterator iterator = node.iterateChildNodesForward(); iterator.hasNext(); ) {
                    final PhylogenyNode child_node = iterator.next();
                    final Color child_color = getBranchColorValue( child_node );
                    if ( child_color != null ) {
                        ++n;
                        red += child_color.getRed();
                        green += child_color.getGreen();
                        blue += child_color.getBlue();
                    }
                }
                setBranchColorValue( node, new Color( ForesterUtil.roundToInt( red / n ), ForesterUtil
                        .roundToInt( green / n ), ForesterUtil.roundToInt( blue / n ) ) );
            }
        }
    }

    public static void removeNode( final PhylogenyNode remove_me, final Phylogeny phylogeny ) {
        if ( remove_me.isRoot() ) {
            throw new IllegalArgumentException( "ill advised attempt to remove root node" );
        }
        if ( remove_me.isExternal() ) {
            phylogeny.deleteSubtree( remove_me, false );
        }
        else {
            final PhylogenyNode parent = remove_me.getParent();
            final List<PhylogenyNode> descs = remove_me.getDescendants();
            parent.removeChildNode( remove_me );
            for( final PhylogenyNode desc : descs ) {
                parent.addAsChild( desc );
                desc.setDistanceToParent( addPhylogenyDistances( remove_me.getDistanceToParent(), desc
                        .getDistanceToParent() ) );
            }
            remove_me.setParent( null );
            phylogeny.setIdHash( null );
            phylogeny.externalNodesHaveChanged();
        }
    }

    public static Set<PhylogenyNode> searchData( final String query,
                                                 final Phylogeny phy,
                                                 final boolean case_sensitive,
                                                 final boolean partial ) {
        final Set<PhylogenyNode> nodes = new HashSet<PhylogenyNode>();
        if ( phy.isEmpty() || ( query == null ) ) {
            return nodes;
        }
        if ( ForesterUtil.isEmpty( query ) ) {
            return nodes;
        }
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            boolean match = false;
            if ( match( node.getNodeName(), query, case_sensitive, partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasTaxonomy()
                    && match( node.getNodeData().getTaxonomy().getTaxonomyCode(), query, case_sensitive, partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasTaxonomy()
                    && match( node.getNodeData().getTaxonomy().getCommonName(), query, case_sensitive, partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasTaxonomy()
                    && match( node.getNodeData().getTaxonomy().getScientificName(), query, case_sensitive, partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasTaxonomy()
                    && ( node.getNodeData().getTaxonomy().getIdentifier() != null )
                    && match( node.getNodeData().getTaxonomy().getIdentifier().getValue(),
                              query,
                              case_sensitive,
                              partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasSequence()
                    && match( node.getNodeData().getSequence().getName(), query, case_sensitive, partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasSequence()
                    && match( node.getNodeData().getSequence().getSymbol(), query, case_sensitive, partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getAccession() != null )
                    && match( node.getNodeData().getSequence().getAccession().getValue(),
                              query,
                              case_sensitive,
                              partial ) ) {
                match = true;
            }
            else if ( node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
                final DomainArchitecture da = node.getNodeData().getSequence().getDomainArchitecture();
                I: for( int i = 0; i < da.getNumberOfDomains(); ++i ) {
                    if ( match( da.getDomain( i ).getName(), query, case_sensitive, partial ) ) {
                        match = true;
                        break I;
                    }
                }
            }
            if ( match ) {
                nodes.add( node );
            }
        }
        return nodes;
    }

    public static Set<PhylogenyNode> searchDataLogicalAnd( final String[] queries,
                                                           final Phylogeny phy,
                                                           final boolean case_sensitive,
                                                           final boolean partial ) {
        final Set<PhylogenyNode> nodes = new HashSet<PhylogenyNode>();
        if ( phy.isEmpty() || ( queries == null ) || ( queries.length < 1 ) ) {
            return nodes;
        }
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            boolean all_matched = true;
            for( final String query : queries ) {
                boolean match = false;
                if ( ForesterUtil.isEmpty( query ) ) {
                    continue;
                }
                if ( match( node.getNodeName(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasTaxonomy()
                        && match( node.getNodeData().getTaxonomy().getTaxonomyCode(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasTaxonomy()
                        && match( node.getNodeData().getTaxonomy().getCommonName(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasTaxonomy()
                        && match( node.getNodeData().getTaxonomy().getScientificName(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasTaxonomy()
                        && ( node.getNodeData().getTaxonomy().getIdentifier() != null )
                        && match( node.getNodeData().getTaxonomy().getIdentifier().getValue(),
                                  query,
                                  case_sensitive,
                                  partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasSequence()
                        && match( node.getNodeData().getSequence().getName(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasSequence()
                        && match( node.getNodeData().getSequence().getSymbol(), query, case_sensitive, partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getAccession() != null )
                        && match( node.getNodeData().getSequence().getAccession().getValue(),
                                  query,
                                  case_sensitive,
                                  partial ) ) {
                    match = true;
                }
                else if ( node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
                    final DomainArchitecture da = node.getNodeData().getSequence().getDomainArchitecture();
                    I: for( int i = 0; i < da.getNumberOfDomains(); ++i ) {
                        if ( match( da.getDomain( i ).getName(), query, case_sensitive, partial ) ) {
                            match = true;
                            break I;
                        }
                    }
                }
                if ( !match ) {
                    all_matched = false;
                    break;
                }
            }
            if ( all_matched ) {
                nodes.add( node );
            }
        }
        return nodes;
    }

    /**
     * Convenience method.
     * Sets value for the first confidence value (created if not present, values overwritten otherwise). 
     */
    public static void setBootstrapConfidence( final PhylogenyNode node, final double bootstrap_confidence_value ) {
        setConfidence( node, bootstrap_confidence_value, "bootstrap" );
    }

    public static void setBranchColorValue( final PhylogenyNode node, final Color color ) {
        if ( node.getBranchData().getBranchColor() == null ) {
            node.getBranchData().setBranchColor( new BranchColor() );
        }
        node.getBranchData().getBranchColor().setValue( color );
    }

    /**
     * Convenience method
     */
    public static void setBranchWidthValue( final PhylogenyNode node, final double branch_width_value ) {
        node.getBranchData().setBranchWidth( new BranchWidth( branch_width_value ) );
    }

    /**
     * Convenience method.
     * Sets value for the first confidence value (created if not present, values overwritten otherwise). 
     */
    public static void setConfidence( final PhylogenyNode node, final double confidence_value ) {
        setConfidence( node, confidence_value, "" );
    }

    /**
     * Convenience method.
     * Sets value for the first confidence value (created if not present, values overwritten otherwise). 
     */
    public static void setConfidence( final PhylogenyNode node, final double confidence_value, final String type ) {
        Confidence c = null;
        if ( node.getBranchData().getNumberOfConfidences() > 0 ) {
            c = node.getBranchData().getConfidence( 0 );
        }
        else {
            c = new Confidence();
            node.getBranchData().addConfidence( c );
        }
        c.setType( type );
        c.setValue( confidence_value );
    }

    public static void setScientificName( final PhylogenyNode node, final String scientific_name ) {
        if ( !node.getNodeData().isHasTaxonomy() ) {
            node.getNodeData().setTaxonomy( new Taxonomy() );
        }
        node.getNodeData().getTaxonomy().setScientificName( scientific_name );
    }

    /**
     * Convenience method to set the taxonomy code of a phylogeny node.
     * 
     * 
     * @param node
     * @param taxonomy_code
     */
    public static void setTaxonomyCode( final PhylogenyNode node, final String taxonomy_code ) {
        if ( !node.getNodeData().isHasTaxonomy() ) {
            node.getNodeData().setTaxonomy( new Taxonomy() );
        }
        node.getNodeData().getTaxonomy().setTaxonomyCode( taxonomy_code );
    }

    /**
     * Removes from Phylogeny to_be_stripped all external Nodes which are
     * associated with a species NOT found in Phylogeny reference.
     * 
     * @param reference
     *            a reference Phylogeny
     * @param to_be_stripped
     *            Phylogeny to be stripped
     * @return number of external nodes removed from to_be_stripped
     */
    public static int taxonomyBasedDeletionOfExternalNodes( final Phylogeny reference, final Phylogeny to_be_stripped ) {
        final Set<String> ref_ext_taxo = new HashSet<String>();
        final ArrayList<PhylogenyNode> nodes_to_delete = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator it = reference.iteratorExternalForward(); it.hasNext(); ) {
            ref_ext_taxo.add( getSpecies( it.next() ) );
        }
        for( final PhylogenyNodeIterator it = to_be_stripped.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !ref_ext_taxo.contains( getSpecies( n ) ) ) {
                nodes_to_delete.add( n );
            }
        }
        for( final PhylogenyNode phylogenyNode : nodes_to_delete ) {
            to_be_stripped.deleteSubtree( phylogenyNode, true );
        }
        return nodes_to_delete.size();
    }
}
