// $Id: SDIR.java,v 1.5 2009/06/19 16:27:30 cmzmasek Exp $
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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyBranch;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/*
 * Allows to infer duplications - speciations on a unrooted gene tree. It
 * reroots the gene trees on each of its branches and performs SDIse on each of
 * the resulting trees. Trees which minimize a certain criterion are returned as
 * the "correctly" rooted ones. The criterions are: <ul> <li>Sum of duplications
 * <li>Mapping cost L <li>Phylogeny height - which is the largest distance from
 * root to external node (minimizing of which is the same as "midpoint rooting")
 * </ul>
 * 
 * @see SDIse
 * 
 * @see SDI
 * 
 * @author Christian M. Zmasek
 */
public class SDIR {

    private final static double ZERO_DIFF = 1.0E-6; // Due to inaccurate
    // calculations on
    // Java's side, not
    // everything that should
    // be 0.0 is 0.0.
    private int                 _count;
    private int                 _min_dup;
    private int                 _min_cost;
    private double              _min_height;
    private double              _min_diff;
    private long                _time_sdi;

    /**
     * Default contructor which creates an "empty" object..
     */
    public SDIR() {
        init();
    }

    /**
     * Returns the number of differently rooted trees which minimize the
     * (rooting) "criterion" - as determined by method "infer".
     * 
     * @see #infer(Phylogeny,Phylogeny,boolean,boolean,boolean,boolean,int,boolean)
     * @return number of differently rooted trees which minimized the criterion
     */
    public int getCount() {
        return _count;
    }

    /**
     * Returns the (absolue value of the) minimal difference in tree heights of
     * the two subtrees at the root (of the (re)rooted gene tree) - as
     * determined by method "infer" - if minimize_height is set to true.
     * <p>
     * If a tree is midpoint rooted this number is zero.
     * <p>
     * <B>IMPORTANT </B>: If minimize_mapping_cost or minimize_sum_of_dup are
     * also set to true, then this returns the minimal difference in tree
     * heights of the trees which minimize the first criterion, and is therefore
     * not necessarily zero.
     * <p>
     * (Last modified: 01/22/00)
     * 
     * @see #infer(Phylogeny,Phylogeny,boolean,boolean,boolean,boolean,int,boolean)
     * @return the minimal difference in tree heights -- IF calculated by
     *         "infer"
     */
    public double getMinimalDiffInSubTreeHeights() {
        return _min_diff;
    }

    /**
     * Returns the minimal number of duplications - as determined by method
     * "infer".
     * <p>
     * <B>IMPORTANT </B>: If the tree is not rooted by minimizing the sum of
     * duplications or the mapping cost L, then this number is NOT NECESSARILY
     * the MINIMAL number of duplications.
     * 
     * @see #infer(Phylogeny,Phylogeny,boolean,boolean,boolean,boolean,int,boolean)
     * @return (minimal) number of duplications
     */
    public int getMinimalDuplications() {
        return _min_dup;
    }

    /**
     * Returns the minimal mapping cost L - as determined by method "infer" - if
     * minimize_mapping_cost is set to true.
     * <p>
     * (Last modified: 11/07/00)
     * 
     * @see #infer(Phylogeny,Phylogeny,boolean,boolean,boolean,boolean,int,boolean)
     * @return the minimal mapping cost "L" -- IF calculated by "infer"
     */
    public int getMinimalMappingCost() {
        return _min_cost;
    }

    /**
     * Returns the minimal tree height - as determined by method "infer" - if
     * minimize_height is set to true. <B>IMPORTANT </B>: If
     * minimize_mapping_cost or minimize_sum_of_dup are also set to true, then
     * this returns the minimal tree height of the trees which minimize the
     * first criterion.
     * <p>
     * (Last modified: 01/12/00)
     * 
     * @see #infer(Phylogeny,Phylogeny,boolean,boolean,boolean,boolean,int,boolean)
     * @return the minimal tree height -- IF calculated by "infer"
     */
    public double getMinimalTreeHeight() {
        return _min_height;
    }

    /**
     * Returns the sum of times (in ms) needed to run method infer of class SDI.
     * Final variable TIME needs to be set to true.
     * 
     * @return sum of times (in ms) needed to run method infer of class SDI
     */
    public long getTimeSumSDI() {
        return _time_sdi;
    }

    /**
     * Infers gene duplications on a possibly unrooted gene Phylogeny gene_tree.
     * The tree is rooted be minimizing either the sum of duplications, the
     * mapping cost L, or the tree height (or combinations thereof). If
     * return_trees is set to true, it returns an array of possibly more than
     * one differently rooted Trees. <br>
     * The maximal number of returned trees is set with max_trees_to_return.
     * <br>
     * Phylogeny species_tree is a species Phylogeny to which the gene Phylogeny
     * gene_tree is compared to. <br>
     * If both minimize_sum_of_dup and minimize_mapping_cost are true, the tree
     * is rooted by minimizing the mapping cost L.<br>
     * If minimize_sum_of_dup, minimize_mapping_cost, and minimize_height are
     * false tree gene_tree is assumed to be alreadty rooted and no attempts at
     * rooting are made, and only one tree is returned. <br>
     * <p>
     * Conditions:
     * </p>
     * <ul>
     * <li>Both Trees must be completely binary (except deepest node of gene
     * tree)
     * <li>The species Phylogeny must be rooted
     * <li>Both Trees must have species names in the species name fields of
     * their nodes
     * <li>Both Trees must not have any collapses nodes
     * </ul>
     * <p>
     * (Last modified: 10/01/01)
     * 
     * @param gene_tree
     *            a binary (except deepest node) gene Phylogeny
     * @param species_tree
     *            a rooted binary species Phylogeny
     * @param minimize_mapping_cost
     *            set to true to root by minimizing the mapping cost L (and also
     *            the sum of duplications)
     * @param minimize_sum_of_dup
     *            set to true to root by minimizing the sum of duplications
     * @param minimize_height
     *            set to true to root by minimizing the tree height - if
     *            minimize_mapping_cost is set to true or minimize_sum_of_dup is
     *            set to true, then out of the resulting trees with minimal
     *            mapping cost or minimal number of duplications the tree with
     *            the minimal height is chosen
     * @param return_trees
     *            set to true to return Array of Trees, otherwise null is
     *            returned
     * @param max_trees_to_return
     *            maximal number of Trees to return (=maximal size of returned
     *            Array) must be no lower than 1
     * @return array of rooted Trees with duplication vs. speciation assigned if
     *         return_trees is set to true, null otherwise
     */
    public Phylogeny[] infer( final Phylogeny gene_tree,
                              final Phylogeny species_tree,
                              final boolean minimize_mapping_cost,
                              boolean minimize_sum_of_dup,
                              final boolean minimize_height,
                              final boolean return_trees,
                              int max_trees_to_return ) {
        init();
        SDIse sdise = null;
        final ArrayList<Phylogeny> trees = new ArrayList<Phylogeny>();
        Phylogeny[] tree_array = null;
        List<PhylogenyBranch> branches = null;
        Phylogeny g = null;
        PhylogenyNode prev_root = null;
        PhylogenyNode prev_root_c1 = null;
        PhylogenyNode prev_root_c2 = null;
        int duplications = 0;
        int cost = 0;
        int counter = 0;
        int min_duplications = Integer.MAX_VALUE;
        int min_cost = Integer.MAX_VALUE;
        int j = 0;
        double height = 0.0;
        double diff = 0.0;
        double min_height = Double.MAX_VALUE;
        double min_diff = 0.0;
        double[] height__diff = new double[ 2 ];
        boolean smaller = false;
        boolean equal = false;
        boolean prev_root_was_dup = false;
        if ( max_trees_to_return < 1 ) {
            max_trees_to_return = 1;
        }
        if ( minimize_mapping_cost && minimize_sum_of_dup ) {
            minimize_sum_of_dup = false;
        }
        if ( !minimize_mapping_cost && !minimize_sum_of_dup && !minimize_height ) {
            throw new IllegalArgumentException( "parameter to minimize not given for rooting of phylogeny" );
        }
        g = gene_tree.copy();
        if ( g.getNumberOfExternalNodes() <= 1 ) {
            g.setRooted( true );
            setMinimalDuplications( 0 );
            setMinimalTreeHeight( 0.0 );
            tree_array = new Phylogeny[ 1 ];
            tree_array[ 0 ] = g;
            return tree_array;
        }
        for( final PhylogenyNodeIterator iter = g.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( n.isRoot() ) {
                if ( ( n.getNumberOfDescendants() != 2 ) && ( n.getNumberOfDescendants() != 3 ) ) {
                    throw new IllegalArgumentException( "attempt to run SDI on gene tree with "
                            + n.getNumberOfDescendants() + " child nodes at its root" );
                }
            }
            else if ( !n.isExternal() && ( n.getNumberOfDescendants() != 2 ) ) {
                throw new IllegalArgumentException( "attempt to run SDI on gene tree which is not completely binary [found node with "
                        + n.getNumberOfDescendants() + " child nodes]" );
            }
        }
        for( final PhylogenyNodeIterator iter = species_tree.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( !n.isExternal() && ( n.getNumberOfDescendants() != 2 ) ) {
                throw new IllegalArgumentException( "attempt to run SDI with a species tree which is not completely binary (after stripping) [found node with "
                        + n.getNumberOfDescendants() + " child nodes]" );
            }
        }
        g.reRoot( g.getFirstExternalNode() );
        branches = SDIR.getBranchesInPreorder( g );
        if ( minimize_mapping_cost || minimize_sum_of_dup ) {
            sdise = new SDIse( g, species_tree );
            duplications = sdise.getDuplicationsSum();
        }
        final Set<PhylogenyBranch> used_root_placements = new HashSet<PhylogenyBranch>();
        F: for( j = 0; j < branches.size(); ++j ) {
            prev_root = g.getRoot();
            prev_root_c1 = prev_root.getChildNode1();
            prev_root_c2 = prev_root.getChildNode2();
            prev_root_was_dup = prev_root.isDuplication();
            final PhylogenyBranch current_branch = branches.get( j );
            g.reRoot( current_branch );
            if ( minimize_mapping_cost || minimize_sum_of_dup ) {
                duplications = sdise.updateM( prev_root_was_dup, prev_root_c1, prev_root_c2 );
            }
            if ( !used_root_placements.contains( current_branch ) ) {
                if ( minimize_mapping_cost ) {
                    cost = sdise.computeMappingCostL();
                    if ( minimize_height && ( cost <= min_cost ) ) {
                        height__diff = SDIR.moveRootOnBranchToMinHeight( g );
                        height = height__diff[ 0 ];
                        diff = height__diff[ 1 ];
                    }
                    if ( cost == min_cost ) {
                        if ( minimize_height ) {
                            smaller = equal = false;
                            if ( height < min_height ) {
                                min_height = height;
                                counter = 1;
                                smaller = true;
                            }
                            else if ( height == min_height ) {
                                counter++;
                                equal = true;
                            }
                            if ( Math.abs( diff ) < min_diff ) {
                                min_diff = Math.abs( diff );
                            }
                        }
                        if ( return_trees ) {
                            if ( minimize_height ) {
                                if ( smaller ) {
                                    trees.clear();
                                    trees.add( g.copy() );
                                }
                                else if ( equal && ( trees.size() < max_trees_to_return ) ) {
                                    trees.add( g.copy() );
                                }
                            }
                            else {
                                counter++;
                                if ( trees.size() < max_trees_to_return ) {
                                    trees.add( g.copy() );
                                }
                            }
                        }
                        else if ( !minimize_height ) {
                            counter++;
                        }
                    }
                    else if ( cost < min_cost ) {
                        if ( minimize_height ) {
                            min_height = height;
                            min_diff = Math.abs( diff );
                        }
                        if ( return_trees ) {
                            trees.clear();
                            trees.add( g.copy() );
                        }
                        counter = 1;
                        min_cost = cost;
                    }
                    if ( duplications < min_duplications ) {
                        min_duplications = duplications;
                    }
                }
                else if ( minimize_sum_of_dup ) {
                    if ( minimize_height && ( duplications <= min_duplications ) ) {
                        height__diff = SDIR.moveRootOnBranchToMinHeight( g );
                        height = height__diff[ 0 ];
                        diff = height__diff[ 1 ];
                    }
                    if ( duplications == min_duplications ) {
                        if ( minimize_height ) {
                            smaller = equal = false;
                            if ( height < min_height ) {
                                min_height = height;
                                counter = 1;
                                smaller = true;
                            }
                            else if ( height == min_height ) {
                                counter++;
                                equal = true;
                            }
                            if ( Math.abs( diff ) < min_diff ) {
                                min_diff = Math.abs( diff );
                            }
                        }
                        if ( return_trees ) {
                            if ( minimize_height ) {
                                if ( smaller ) {
                                    trees.clear();
                                    trees.add( g.copy() );
                                }
                                else if ( equal && ( trees.size() < max_trees_to_return ) ) {
                                    trees.add( g.copy() );
                                }
                            }
                            else {
                                counter++;
                                if ( trees.size() < max_trees_to_return ) {
                                    trees.add( g.copy() );
                                }
                            }
                        }
                        else if ( !minimize_height ) {
                            counter++;
                        }
                    }
                    else if ( duplications < min_duplications ) {
                        if ( minimize_height ) {
                            min_height = height;
                            min_diff = Math.abs( diff );
                        }
                        if ( return_trees ) {
                            trees.clear();
                            trees.add( g.copy() );
                        }
                        counter = 1;
                        min_duplications = duplications;
                    }
                }
                else if ( minimize_height ) {
                    height__diff = SDIR.moveRootOnBranchToMinHeight( g );
                    height = height__diff[ 0 ];
                    diff = height__diff[ 1 ];
                    if ( Math.abs( diff ) < SDIR.ZERO_DIFF ) {
                        sdise = new SDIse( g, species_tree );
                        min_duplications = sdise.getDuplicationsSum();
                        min_height = height;
                        min_diff = Math.abs( diff );
                        counter = 1;
                        if ( return_trees ) {
                            trees.add( g.copy() );
                        }
                        break F;
                    }
                }
            } // if ( used_root_placements.containsKey( current_branch ) )
            used_root_placements.add( current_branch );
        } // End of huge for loop "F".
        if ( return_trees ) {
            trees.trimToSize();
            tree_array = new Phylogeny[ trees.size() ];
            for( int i = 0; i < trees.size(); ++i ) {
                tree_array[ i ] = trees.get( i );
                tree_array[ i ].recalculateNumberOfExternalDescendants( false );
            }
        }
        setCount( counter );
        setMinimalDuplications( min_duplications );
        setMinimalMappingCost( min_cost );
        setMinimalTreeHeight( min_height );
        setMinimalDiffInSubTreeHeights( Math.abs( min_diff ) );
        return tree_array;
    }

    private void init() {
        _count = -1;
        _min_dup = Integer.MAX_VALUE;
        _min_cost = Integer.MAX_VALUE;
        _min_height = Double.MAX_VALUE;
        _min_diff = Double.MAX_VALUE;
        _time_sdi = -1;
    }

    private void setCount( final int i ) {
        _count = i;
    }

    private void setMinimalDiffInSubTreeHeights( final double d ) {
        _min_diff = d;
    }

    private void setMinimalDuplications( final int i ) {
        _min_dup = i;
    }

    private void setMinimalMappingCost( final int i ) {
        _min_cost = i;
    }

    private void setMinimalTreeHeight( final double d ) {
        _min_height = d;
    }

    // This was totally changed on 2006/10/03.
    // Places references to all Branches of Phylogeny t into a List.
    // The order is preorder.
    // Trees are treated as if they were unrooted (i.e. child 1 and
    // child 2 of the root are treated as if they were connected
    // directly).
    // The resulting List allows to visit all branches without ever
    // traversing more than one node at a time.
    public static List<PhylogenyBranch> getBranchesInPreorder( final Phylogeny t ) {
        final ArrayList<PhylogenyBranch> branches = new ArrayList<PhylogenyBranch>();
        if ( t.isEmpty() || ( t.getNumberOfExternalNodes() <= 1 ) ) {
            return branches;
        }
        if ( t.getNumberOfExternalNodes() == 2 ) {
            branches.add( new PhylogenyBranch( t.getRoot().getChildNode1(), t.getRoot().getChildNode2() ) );
            return branches;
        }
        final Set<PhylogenyNode> one = new HashSet<PhylogenyNode>();
        final Set<PhylogenyNode> two = new HashSet<PhylogenyNode>();
        PhylogenyNode node = t.getRoot();
        while ( !node.isRoot() || !two.contains( node ) ) {
            if ( !node.isExternal() && !two.contains( node ) ) {
                if ( !one.contains( node ) && !two.contains( node ) ) {
                    one.add( node );
                    node = node.getChildNode1();
                }
                else {
                    two.add( node );
                    node = node.getChildNode2();
                }
                if ( !node.getParent().isRoot() ) {
                    branches.add( new PhylogenyBranch( node, node.getParent() ) );
                }
                else if ( !node.isExternal() ) {
                    branches.add( new PhylogenyBranch( t.getRoot().getChildNode1(), t.getRoot().getChildNode2() ) );
                }
            }
            else {
                if ( !node.getParent().isRoot() && !node.isExternal() ) {
                    branches.add( new PhylogenyBranch( node, node.getParent() ) );
                }
                node = node.getParent();
            }
        }
        return branches;
    }

    // This places the root of t on its branch in such a way that it
    // minimizes the tree height as good as possible.
    // Returns the height and the difference in heights of the resulting
    // modified Phylogeny t.
    private static double[] moveRootOnBranchToMinHeight( final Phylogeny t ) {
        final PhylogenyNode root = t.getRoot();
        if ( root.getNumberOfDescendants() != 2 ) {
            throw new IllegalArgumentException( "attempt to move root to minimize height on root where number of child nodes does not equal two" );
        }
        final PhylogenyNode child0 = root.getChildNode( 0 );
        final PhylogenyNode child1 = root.getChildNode( 1 );
        final double newdist = 0.5 * ( ( child0.getDistanceToParent() > 0 ? child0.getDistanceToParent() : 0 ) + ( child1
                .getDistanceToParent() > 0 ? child1.getDistanceToParent() : 0 ) );
        child0.setDistanceToParent( newdist );
        child1.setDistanceToParent( newdist );
        final double d = child0.getDistanceToParent();
        double diff = 0.0;
        double height = 0.0;
        final double[] height_diff = new double[ 2 ];
        final double l0 = t.calculateSubtreeHeight( t.getRoot().getChildNode( 0 ) );
        final double l1 = t.calculateSubtreeHeight( t.getRoot().getChildNode( 1 ) );
        diff = l0 - l1;
        height = t.getHeight();
        if ( d > 0.0 ) {
            if ( ( 2 * d ) > Math.abs( diff ) ) {
                child0.setDistanceToParent( d - ( diff / 2.0 ) );
                child1.setDistanceToParent( d + ( diff / 2.0 ) );
                height_diff[ 0 ] = height - Math.abs( diff / 2 );
                height_diff[ 1 ] = 0.0;
            }
            else {
                if ( diff > 0 ) {
                    child0.setDistanceToParent( 0.0 );
                    child1.setDistanceToParent( 2 * d );
                    height_diff[ 1 ] = diff - ( 2 * d );
                }
                else {
                    child0.setDistanceToParent( 2 * d );
                    child1.setDistanceToParent( 0.0 );
                    height_diff[ 1 ] = diff + ( 2 * d );
                }
                height_diff[ 0 ] = height - d;
            }
        }
        else {
            height_diff[ 0 ] = height;
            height_diff[ 1 ] = diff;
        }
        return height_diff;
    }
}
