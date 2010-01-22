// $Id: DistanceCalculator.java,v 1.13 2009/11/20 22:22:10 cmzmasek Exp $
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
import java.util.ListIterator;
import java.util.Vector;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.ForesterUtil;

/*
 * @author Christian M. Zmasek
 * 
 * @version 1.001 -- last modified: 12/04/00
 */
public class DistanceCalculator {

    public final static double       DEFAULT = -1.0;
    private Phylogeny                tree_;
    private ArrayList<PhylogenyNode> nodes_;
    private int                      n_;
    private double                   mean_, variance_, stand_dev_;
    private PhylogenyNode            lca_;                        // The LCA of the

    // Nodes in nodes_
    /**
     * Default constructor. (Last modified: 11/30/00)
     */
    public DistanceCalculator() {
        tree_ = null;
        nodes_ = null;
        n_ = 0;
        mean_ = DistanceCalculator.DEFAULT;
        variance_ = DistanceCalculator.DEFAULT;
        stand_dev_ = DistanceCalculator.DEFAULT;
        lca_ = null;
    }

    /**
     * Constructor. Sets the rooted Phylogeny t for which the mean distance to
     * the root and its variance and standard deviation are calculated. (Last
     * modified: 12/01/00)
     * 
     * @param t
     *            the rooted Phylogeny for which the mean distance to the root
     *            and its variance and standard deviation are calculated
     */
    public DistanceCalculator( final Phylogeny t ) {
        setTree( t );
    }

    /**
     * Constructor. Sets the rooted Phylogeny t and the external Nodes ext_nodes
     * for which the mean distance to their lowest common ancestor and its
     * variance and standard deviation are calculated. (Last modified: 12/01/00)
     * 
     * @param t
     *            the rooted Phylogeny containing Nodes in Vector ext_nodes
     * @param ext_nodes
     *            a Vector of Nodes of t, the mean distance to their lowest
     *            common ancestor and its variance and standard deviation are
     *            calculated
     */
    public DistanceCalculator( final Phylogeny t, final Vector<PhylogenyNode> ext_nodes ) {
        setTreeAndExtNodes( t, ext_nodes );
    }

    // (Last modified: 12/01/00)
    private PhylogenyNode calculateLCA( final ArrayList<PhylogenyNode> nodes ) {
        if ( ( nodes == null ) || nodes.isEmpty() ) {
            return null;
        }
        PhylogenyNode node = nodes.get( 0 );
        int c = node.getNumberOfExternalNodes();
        final int v = nodes.size();
        while ( !node.isRoot() && ( c < v ) ) {
            node = node.getParent();
            c = node.getNumberOfExternalNodes();
        }
        ArrayList<PhylogenyNode> current_nodes = new ArrayList<PhylogenyNode>( node.getAllExternalDescendants() );
        while ( !node.isRoot() && !current_nodes.containsAll( nodes ) ) {
            node = node.getParent();
            current_nodes = new ArrayList<PhylogenyNode>( node.getAllExternalDescendants() );
        }
        return node;
    }

    // (Last modified: 11/31/00)
    private void calculateMean() {
        if ( ( nodes_ == null ) || nodes_.isEmpty() || ( tree_ == null ) || tree_.isEmpty() ) {
            return;
        }
        double sum = 0.0;
        final ListIterator<PhylogenyNode> li = nodes_.listIterator();
        n_ = 0;
        try {
            while ( li.hasNext() ) {
                n_++;
                sum += getDistanceToNode( li.next(), lca_ );
            }
        }
        catch ( final Exception e ) {
            System.err.println( "calculateMean(): " + "Exception: " + e );
            System.exit( -1 );
        }
        setMean( sum / n_ );
    }

    // (Last modified: 11/30/00)
    private void calculateMeanDistToRoot() {
        if ( ( tree_ == null ) || tree_.isEmpty() ) {
            return;
        }
        double sum = 0.0;
        PhylogenyNode node = tree_.getFirstExternalNode();
        n_ = 0;
        while ( node != null ) {
            n_++;
            sum += getDistanceToRoot( node );
            node = node.getNextExternalNode();
        }
        setMean( sum / n_ );
    }

    // (Last modified: 11/31/00)
    private void calculateStandardDeviation() {
        if ( ( getVariance() == DistanceCalculator.DEFAULT ) || ( getVariance() < 0.0 ) ) {
            return;
        }
        setStandardDeviation( java.lang.Math.sqrt( getVariance() ) );
    }

    // (Last modified: 11/31/00)
    private void calculateVariance() {
        if ( ( getMean() == DistanceCalculator.DEFAULT ) || ( nodes_ == null ) || nodes_.isEmpty() || ( tree_ == null )
                || tree_.isEmpty() || ( n_ <= 1.0 ) ) {
            return;
        }
        double x = 0.0, sum = 0.0;
        final ListIterator<PhylogenyNode> li = nodes_.listIterator();
        try {
            while ( li.hasNext() ) {
                x = getDistanceToNode( li.next(), lca_ ) - getMean();
                sum += ( x * x );
            }
        }
        catch ( final Exception e ) {
            System.err.println( "calculateVariance(): " + "Exception: " + e );
            System.exit( -1 );
        }
        setVariance( sum / ( n_ - 1 ) );
    }

    // (Last modified: 11/31/00)
    private void calculateVarianceDistToRoot() {
        if ( ( getMean() == DistanceCalculator.DEFAULT ) || ( tree_ == null ) || tree_.isEmpty() || ( n_ <= 1.0 ) ) {
            return;
        }
        double x = 0.0, sum = 0.0;
        PhylogenyNode node = tree_.getFirstExternalNode();
        while ( node != null ) {
            x = getDistanceToRoot( node ) - getMean();
            sum += ( x * x );
            node = node.getNextExternalNode();
        }
        setVariance( sum / ( n_ - 1 ) );
    }

    /**
     * Calculates the distance of the PhylogenyNode with seq name seq_name to
     * the LCA of ext_nodes, which has been set either with constructor
     * DistanceCalculator(Phylogeny,Vector) or method
     * setTreeAndExtNodes(Phylogeny,Vector). Throws an exception if no
     * PhylogenyNode with seq name_seq name is found or if seq_name is not
     * unique. (Last modified: 12/03/00)
     * 
     * @param seq_name
     *            the seq name for the PhylogenyNode for which the distance to
     *            the LCA is to be calculated
     * @return distance of PhylogenyNode with seq name seq_name to the LCA of
     *         Nodes in ext_nodes
     * @see #DistanceCalculator(Phylogeny,Vector)
     * @see #setTreeAndExtNodes(Phylogeny,Vector)
     * @see #setTreeAndExtNodes(Phylogeny,ArrayList)
     */
    public double getDistanceToLCA( final String seq_name ) {
        if ( ( tree_ == null ) || tree_.isEmpty() || ( lca_ == null ) ) {
            return 0.0;
        }
        return getDistanceToNode( seq_name, lca_ );
    }

    /**
     * Calculates the distance of PhylogenyNode outer to PhylogenyNode inner.
     * PhylogenyNode inner must be closer to the root than PhylogenyNode outer
     * and on the same "path". (Last modified: 12/01/00)
     * 
     * @param outer
     *            a PhylogenyNode
     * @param inner
     *            a PhylogenyNode closer to the root than outer
     * @return distance of PhylogenyNode outer to PhylogenyNode inner
     */
    public double getDistanceToNode( PhylogenyNode outer, final PhylogenyNode inner ) {
        double d = 0.0, dist = 0.0;
        while ( ( inner != outer ) && !outer.isRoot() ) {
            d = outer.getDistanceToParent();
            if ( d > 0.0 ) {
                dist += d;
            }
            outer = outer.getParent();
        }
        if ( !inner.isRoot() && outer.isRoot() ) {
            throw new IllegalArgumentException( "getDistanceToNode(PhylogenyNode outer,PhylogenyNode inner): "
                    + "PhylogenyNode inner is not closer to the root than PhylogenyNode outer "
                    + "or is not on the same \"subtree\"" );
        }
        return dist;
    }

    /**
     * Calculates the distance of the PhylogenyNode with seq name seq_name to
     * PhylogenyNode inner. PhylogenyNode inner must be closer to the root than
     * the PhylogenyNode with seq name seq_name and on the same "path". Throws
     * an exception if no PhylogenyNode with seq name_seq name is found or if
     * seq_name is not unique. (Last modified: 12/01/00)
     * 
     * @param seq_name
     *            the seq name of a PhylogenyNode further from the root than
     *            PhylogenyNode inner
     * @param inner
     *            a PhylogenyNode
     * @return distance of PhylogenyNode with seq name seq_nam to PhylogenyNode
     *         inner
     */
    public double getDistanceToNode( final String seq_name, final PhylogenyNode inner ) {
        if ( ( tree_ == null ) || tree_.isEmpty() ) {
            return 0.0;
        }
        return getDistanceToNode( tree_.getNodeViaSequenceName( seq_name ), inner );
    }

    /**
     * Calculates the distance of PhylogenyNode n to the root of Phylogeny t
     * which has been set either with a constructor, setTree(Phylogeny), or
     * setTreeAndExtNodes(Phylogeny,Vector). (Last modified: 12/01/00)
     * 
     * @param n
     *            the PhylogenyNode for which the distance to the root is to be
     *            calculated
     * @return distance of PhylogenyNode n to the root
     * @see #DistanceCalculator(Phylogeny)
     * @see #DistanceCalculator(Phylogeny,Vector)
     * @see #setTree(Phylogeny)
     * @see #setTreeAndExtNodes(Phylogeny,Vector)
     */
    public double getDistanceToRoot( final PhylogenyNode n ) {
        if ( ( tree_ == null ) || tree_.isEmpty() ) {
            return 0.0;
        }
        double d = 0.0;
        try {
            d = getDistanceToNode( n, tree_.getRoot() );
        }
        catch ( final Exception e ) {
            System.err.println( "getDistanceToRoot(PhylogenyNode): Unexpected " + "exception: " + e );
            System.exit( -1 );
        }
        return d;
    }

    /**
     * Calculates the distance of the PhylogenyNode with seq name seq_name to
     * the root of Phylogeny t, which has been set either with a constructor,
     * setTree(Phylogeny), or setTreeAndExtNodes(Phylogeny,Vector). Throws an
     * exception if no PhylogenyNode with seq name_seq name is found or if
     * seq_name is not unique. (Last modified: 12/01/00)
     * 
     * @param seq_name
     *            the seq name for the PhylogenyNode for which the distance to
     *            the root is to be calculated
     * @return distance of PhylogenyNode with seq name seq_name to the root
     * @see #DistanceCalculator(Phylogeny)
     * @see #DistanceCalculator(Phylogeny,Vector)
     * @see #setTree(Phylogeny)
     * @see #setTreeAndExtNodes(Phylogeny,Vector)
     * @see #setTreeAndExtNodes(Phylogeny,ArrayList)
     */
    public double getDistanceToRoot( final String seq_name ) {
        if ( ( tree_ == null ) || tree_.isEmpty() ) {
            return 0.0;
        }
        return getDistanceToNode( seq_name, tree_.getRoot() );
    }

    /**
     * Returns the mean distance. If constructor DistanceCalculator(Phylogeny)
     * or method setTree(Phylogeny) have been used, it is the mean of the
     * distances from the root to all external Nodes. If constructor
     * DistanceCalculator(Phylogeny,Vector) or method
     * setTreeAndExtNodes(Phylogeny,Vector) have been used, it is the mean of
     * the distances from the external nodes ext_nodes to their lowest common
     * ancestor. (Last modified: 11/30/00)
     * 
     * @return mean distance
     * @see #DistanceCalculator(Phylogeny)
     * @see #DistanceCalculator(Phylogeny,Vector)
     * @see #setTree(Phylogeny)
     * @see #setTreeAndExtNodes(Phylogeny,Vector)
     * @see #setTreeAndExtNodes(Phylogeny,ArrayList)
     */
    public double getMean() {
        return mean_;
    }

    /**
     * Returns the sum of all Nodes used to calculate the mean. (Last modified:
     * 12/01/00)
     * 
     * @return n
     */
    public int getN() {
        return n_;
    }

    /**
     * Returns the standard deviation. If constructor
     * DistanceCalculator(Phylogeny) or method setTree(Phylogeny) have been
     * used, it is the standard deviation of the distances from the root to all
     * external Nodes. If constructor DistanceCalculator(Phylogeny,Vector) or
     * method setTreeAndExtNodes(Phylogeny,Vector) have been used, it is the
     * standard deviation of the distances from the external nodes ext_nodes to
     * their lowest common ancestor. (Last modified: 11/30/00)
     * 
     * @return standard deviation
     * @see #DistanceCalculator(Phylogeny)
     * @see #DistanceCalculator(Phylogeny,Vector)
     * @see #setTree(Phylogeny)
     * @see #setTreeAndExtNodes(Phylogeny,Vector)
     * @see #setTreeAndExtNodes(Phylogeny,ArrayList)
     */
    public double getStandardDeviation() {
        return stand_dev_;
    }

    /**
     * Returns the variance. ( 1/(N - 1) * Sum((x-mean)^2) ) If constructor
     * DistanceCalculator(Phylogeny) or method setTree(Phylogeny) have been
     * used, it is the variance of the distances from the root to all external
     * Nodes. If constructor DistanceCalculator(Phylogeny,Vector) or method
     * setTreeAndExtNodes(Phylogeny,Vector) have been used, it is the variance
     * of the distances from the external nodes ext_nodes to their lowest common
     * ancestor. (Last modified: 11/30/00)
     * 
     * @return variance
     * @see #DistanceCalculator(Phylogeny)
     * @see #DistanceCalculator(Phylogeny,Vector)
     * @see #setTree(Phylogeny)
     * @see #setTreeAndExtNodes(Phylogeny,Vector)
     * @see #setTreeAndExtNodes(Phylogeny,ArrayList)
     */
    public double getVariance() {
        return variance_;
    }

    // (Last modified: 11/30/00)
    private void setMean( final double d ) {
        mean_ = d;
    }

    // (Last modified: 11/30/00)
    private void setStandardDeviation( final double d ) {
        stand_dev_ = d;
    }

    /**
     * Sets the rooted Phylogeny t for which the mean distance to the root and
     * its variance and standard deviation are calculated. (Last modified:
     * 12/01/00)
     * 
     * @param t
     *            the rooted Phylogeny for which the mean distance to the root
     *            and its variance and standard deviation are calculated
     */
    public void setTree( final Phylogeny t ) {
        tree_ = t;
        nodes_ = null;
        n_ = 0;
        mean_ = DistanceCalculator.DEFAULT;
        variance_ = DistanceCalculator.DEFAULT;
        stand_dev_ = DistanceCalculator.DEFAULT;
        lca_ = null;
        calculateMeanDistToRoot();
        calculateVarianceDistToRoot();
        calculateStandardDeviation();
    }

    /**
     * Sets the rooted Phylogeny t and the external Nodes ext_nodes for which
     * the mean distance to their lowest common ancestor and its variance and
     * standard deviation are calculated. (Last modified: 12/03/00)
     * 
     * @param t
     *            the rooted Phylogeny containing Nodes in Vector ext_nodes
     * @param ext_nodes
     *            a ArrayList of Nodes of t, the mean distance to their lowest
     *            common ancestor and its variance and standard deviation are
     *            calculated
     */
    public void setTreeAndExtNodes( final Phylogeny t, final ArrayList<PhylogenyNode> ext_nodes ) {
        tree_ = t;
        nodes_ = ext_nodes;
        n_ = 0;
        mean_ = DistanceCalculator.DEFAULT;
        variance_ = DistanceCalculator.DEFAULT;
        stand_dev_ = DistanceCalculator.DEFAULT;
        lca_ = calculateLCA( nodes_ );
        calculateMean();
        calculateVariance();
        calculateStandardDeviation();
    }

    /**
     * Sets the rooted Phylogeny t and the external Nodes ext_nodes for which
     * the mean distance to their lowest common ancestor and its variance and
     * standard deviation are calculated. (Last modified: 12/03/00)
     * 
     * @param t
     *            the rooted Phylogeny containing Nodes in Vector ext_nodes
     * @param ext_nodes
     *            a Vector of Nodes of t, the mean distance to their lowest
     *            common ancestor and its variance and standard deviation are
     *            calculated
     */
    public void setTreeAndExtNodes( final Phylogeny t, final Vector<PhylogenyNode> ext_nodes ) {
        setTreeAndExtNodes( t, new ArrayList<PhylogenyNode>( ext_nodes ) );
    }

    // (Last modified: 11/30/00)
    private void setVariance( final double d ) {
        variance_ = d;
    }

    // Main for testing.
    public static void main( final String args[] ) {
        File tree_file = null;
        Phylogeny tree = null;
        DistanceCalculator dc = null;
        tree_file = new File( args[ 0 ] );
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( tree_file, true );
            tree = factory.create( tree_file, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( e.toString() );
            System.exit( -1 );
        }
        double time = System.currentTimeMillis();
        dc = new DistanceCalculator( tree );
        final double m = dc.getMean(), var = dc.getVariance(), sd = dc.getStandardDeviation();
        time = ( System.currentTimeMillis() - time );
        System.out.println( "\nn   = " + dc.getN() );
        System.out.println( "mea = " + m );
        System.out.println( "var = " + var );
        System.out.println( "sd  = " + sd + "\n" );
        System.out.println( "t=" + time + "\n" );
    }
}
