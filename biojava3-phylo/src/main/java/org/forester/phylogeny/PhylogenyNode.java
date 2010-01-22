// $Id: PhylogenyNode.java,v 1.127 2009/11/17 03:51:34 cmzmasek Exp $
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

import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.nhx.NHXFormatException;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.phylogeny.data.BranchData;
import org.forester.phylogeny.data.NodeData;
import org.forester.phylogeny.iterators.ChildNodeIteratorForward;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterUtil;

public class PhylogenyNode implements PhylogenyNodeI, Comparable<PhylogenyNode> {

    /** Value of -99.0 is used as default value. */
    public final static double       DISTANCE_DEFAULT = -1024.0;
    private static int               _node_count      = 0;
    private byte                     _indicator;
    private int                      _id;
    private int                      _sum_ext_nodes;
    private float                    _x;
    private float                    _y;
    private String                   _node_name;
    private double                   _distance_parent;
    private boolean                  _collapse;
    private PhylogenyNode            _parent;
    private PhylogenyNode            _link;
    private ArrayList<PhylogenyNode> _descendants;
    private NodeData                 _node_data;
    private BranchData               _branch_data;
    private float                    _x_secondary;
    private float                    _y_secondary;

    private boolean                  _pathToParent = false;

    /**
     * Default constructor for PhylogenyNode.
     */
    public PhylogenyNode() {
        init();
        setNodeId( PhylogenyNode.getNodeCount() );
        PhylogenyNode.increaseNodeCount();
        setSumExtNodes( 1 ); // For ext node, this number is 1 (not 0!!)
    }

    public PhylogenyNode( final String nhx ) throws NHXFormatException {
        this( nhx, TAXONOMY_EXTRACTION.NO );
    }

    public PhylogenyNode( final String nhx, final TAXONOMY_EXTRACTION taxonomy_extraction ) throws NHXFormatException {
        init();
        NHXParser.parseNHX( nhx, this, taxonomy_extraction, false );
        setNodeId( PhylogenyNode.getNodeCount() );
        PhylogenyNode.increaseNodeCount();
        setSumExtNodes( 1 ); // For ext node, this number is 1 (not 0!!)
    }

    /**
     * Constructor for PhylogenyNode.
     * <p>
     * 
     * @param s
     *            String representing one PhylogenyNode in New Hampshire (NH) or
     *            New Hampshire X (NHX) format.
     * @throws NHXFormatException
     * @throws PhylogenyParserException
     */
    public PhylogenyNode( final String nhx,
                          final TAXONOMY_EXTRACTION taxonomy_extraction,
                          final boolean replace_underscores ) throws NHXFormatException {
        init();
        NHXParser.parseNHX( nhx, this, taxonomy_extraction, replace_underscores );
        setNodeId( PhylogenyNode.getNodeCount() );
        PhylogenyNode.increaseNodeCount();
        setSumExtNodes( 1 ); // For ext node, this number is 1 (not 0!!)
    }

    public void setPathToParent(boolean value){
        _pathToParent = value;
    }

    public boolean getPathToParent(){
        return _pathToParent;
    }

    /**
     * Adds PhylogenyNode n to the list of child nodes and sets the _parent of n
     * to this.
     * 
     * @param n
     *            the PhylogenyNode to add
     */
    final public void addAsChild( final PhylogenyNodeI node ) {
        final PhylogenyNode n = ( PhylogenyNode ) node;
        // of
        // interface!!!!
        addChildNode( n );
        n.setParent( this );
    }

    /**
     * Adds PhylogenyNode n to the list of child nodes. But does NOT set the
     * _parent of n to this.
     * 
     * @see addAsChild( PhylogenyNode n )
     * @param n
     *            the PhylogenyNode to add
     */
    final private void addChildNode( final PhylogenyNode child ) {
        getDescendants().add( child );
    }

    final public int compareTo( final PhylogenyNode o ) {
        final PhylogenyNode n = o;
        if ( ( getNodeName() == null ) || ( n.getNodeName() == null ) ) {
            return 0;
        }
        return getNodeName().compareTo( n.getNodeName() );
    }

    // ---------------------------------------------------------
    // Copy and delete Nodes, copy subtress
    // ---------------------------------------------------------
    /**
     * Returns a new PhylogenyNode which has its data copied from this
     * PhylogenyNode. Links to the other Nodes in the same Phylogeny are NOT
     * copied (e.g. _link to _parent). Field "_link" IS copied.
     * 
     * @see #getLink() 
     */
    final public PhylogenyNode copyNodeData() {
        final PhylogenyNode node = new PhylogenyNode();
        PhylogenyNode.decreaseNodeCount();
        node._id = _id;
        node._sum_ext_nodes = _sum_ext_nodes;
        node._indicator = _indicator;
        node._x = _x;
        node._y = _y;
        if ( _node_name != null ) {
            node._node_name = new String( _node_name );
        }
        node._distance_parent = _distance_parent;
        node._collapse = _collapse;
        node._link = _link;
        if ( _node_data != null ) {
            node._node_data = ( NodeData ) _node_data.copy();
        }
        if ( _branch_data != null ) {
            node._branch_data = ( BranchData ) _branch_data.copy();
        }
        return node;
    } // copyNodeData()

    @Override
    /**
     * Based on node name, sequence, and taxonomy.
     * 
     * 
     */
    final public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            return false;
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            final PhylogenyNode other = ( PhylogenyNode ) o;
            if ( !getNodeName().equals( other.getNodeName() ) ) {
                return false;
            }
            final NodeData this_data = getNodeData();
            final NodeData other_data = other.getNodeData();
            if ( ( this_data.isHasSequence() && other_data.isHasSequence() )
                    && ( this_data.isHasTaxonomy() && other_data.isHasTaxonomy() ) ) {
                return ( this_data.getTaxonomy().isEqual( other_data.getTaxonomy() ) && this_data.getSequence()
                        .isEqual( other_data.getSequence() ) );
            }
            else if ( this_data.isHasSequence() && other_data.isHasSequence() ) {
                return ( this_data.getSequence().isEqual( other_data.getSequence() ) );
            }
            else if ( this_data.isHasTaxonomy() && other_data.isHasTaxonomy() ) {
                return ( this_data.getTaxonomy().isEqual( other_data.getTaxonomy() ) );
            }
            else if ( getNodeName().length() > 0 ) {
                // Node name is not empty, and equal.
                return true;
            }
            else {
                return false;
            }
        }
    }

    // ---------------------------------------------------------
    // Obtaining of Nodes
    // ---------------------------------------------------------
    /**
     * Returns a List containing references to all external children of this
     * PhylogenyNode.
     * 
     * @return List of references to external Nodes
     */
    final public List<PhylogenyNode> getAllExternalDescendants() {
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        if ( isExternal() ) {
            nodes.add( this );
            return nodes;
        }
        PhylogenyNode node1 = this;
        while ( !node1.isExternal() ) {
            node1 = node1.getFirstChildNode();
        }
        PhylogenyNode node2 = this;
        while ( !node2.isExternal() ) {
            node2 = node2.getLastChildNode();
        }
        //  do {
        //      nodes.add( node1 );
        //      node1 = node1.getNextExternalNode();
        //  } while ( ( node1 != node2 ) && ( node1 != null ) );
        //  if ( node1 == null ) {
        //      System.out.println( "node1 == null" );
        //  }
        while ( node1 != node2 ) {
            nodes.add( node1 );
            node1 = node1.getNextExternalNode();
        }
        nodes.add( node2 );
        return nodes;
    }

    /**
     * Returns a List containing references to all names of the external
     * children of this PhylogenyNode.
     * 
     * @return List of references to names of external Nodes
     */
    final public List<String> getAllExternalDescendantsNames() {
        final List<PhylogenyNode> c = getAllExternalDescendants();
        final List<String> n = new ArrayList<String>( c.size() );
        for( final PhylogenyNode phylogenyNode : c ) {
            n.add( phylogenyNode.getNodeName() );
        }
        return n;
    }

    final public BranchData getBranchData() {
        if ( _branch_data == null ) {
            _branch_data = new BranchData();
        }
        return _branch_data;
    }

    final BranchData getBranchDataDirectly() {
        return _branch_data;
    }

    /**
     * This return child node n of this node.
     * 
     * @param n
     *            the index of the child to get
     * @return the child node with index n
     * @throws IllegalArgumentException
     *             if n is out of bounds
     */
    final public PhylogenyNode getChildNode( final int i ) {
        if ( isExternal() ) {
            throw new UnsupportedOperationException( "attempt to get the child node of an external node." );
        }
        if ( ( i >= getNumberOfDescendants() ) || ( i < 0 ) ) {
            throw new IllegalArgumentException( "attempt to get child node " + i + " of a node with "
                    + getNumberOfDescendants() + " child nodes" );
        }
        return getDescendants().get( i );
    }

    /**
     * Convenience method. Returns the first child PhylogenyNode of this
     * PhylogenyNode.
     */
    final public PhylogenyNode getChildNode1() {
        return getChildNode( 0 );
    }

    /**
     * Convenience method. Returns the second child PhylogenyNode of this
     * PhylogenyNode.
     * <p>
     * [last modified May 18, 2005 by CMZ]
     */
    final public PhylogenyNode getChildNode2() {
        return getChildNode( 1 );
    }

    /**
     * This gets the child node index of this node.
     * <p>
     * 
     * @return the child node index of this node
     * @throws UnsupportedOperationException
     *             if this node is a root node
     */
    final public int getChildNodeIndex() {
        return getChildNodeIndex( getParent() );
    }

    /**
     * This gets the child node index of this node, given that parent is its
     * parent
     * <p>
     * [last modified Aug 14, 2006 by CMZ]
     * 
     * @return the child node index of this node
     * @throws UnsupportedOperationException
     *             if this node is a root node
     */
    final public int getChildNodeIndex( final PhylogenyNode parent ) {
        if ( isRoot() ) {
            throw new UnsupportedOperationException( "Cannot get the child index for a root node." );
        }
        for( int i = 0; i < parent.getNumberOfDescendants(); ++i ) {
            if ( parent.getChildNode( i ) == this ) {
                return i;
            }
        }
        throw new RuntimeException( "Unexpected exception: Could not determine the child index for node: " + this );
    }

    final public List<PhylogenyNode> getDescendants() {
        return _descendants;
    }

    /**
     * Returns the length of the branch leading to the _parent of this
     * PhylogenyNode (double).
     */
    final public double getDistanceToParent() {
        return _distance_parent;
    }

    /**
     * Convenience method. Returns the first child node of this node.
     * <p>
     * [last modified May 18, 2005 by CMZ]
     * 
     * @return the first child node of this node
     */
    public final PhylogenyNode getFirstChildNode() {
        return getChildNode( 0 );
    }

    /**
     * Returns the _indicator value of this PhylogenyNode.
     */
    public final byte getIndicator() {
        return _indicator;
    }

    /**
     * Convenience method. Returns the last child node of this node.
     * <p>
     * [last modified May 18, 2005 by CMZ]
     * 
     * @return the last child node of this node
     */
    public final PhylogenyNode getLastChildNode() {
        return getChildNode( getNumberOfDescendants() - 1 );
    }

    /**
     * Returns a refernce to the linked PhylogenyNode of this PhylogenyNode.
     * Currently, this method is only used for the speciation-_duplication
     * assignment algorithms.
     */
    public final PhylogenyNode getLink() {
        return _link;
    }

    /**
     * Returns a refernce to the next external PhylogenyNode of this
     * PhylogenyNode. TODO should be in Phylogeny. Returns null if no next
     * external node is available.
     */
    public final PhylogenyNode getNextExternalNode() {
        if ( isInternal() ) {
            throw new UnsupportedOperationException( "attempt to get next external node of an internal node" );
        }
        else if ( isLastExternalNode() ) {
            return null;
        }
        int index = getChildNodeIndex();
        PhylogenyNode previous_node = this;
        PhylogenyNode current_node = getParent();
        while ( !current_node.isRoot()
                && ( ( current_node.getNumberOfDescendants() == 1 ) || previous_node.isLastChildNode() ) ) {
            index = current_node.getChildNodeIndex();
            previous_node = current_node;
            current_node = current_node.getParent();
        }
        current_node = current_node.getChildNode( index + 1 );
        while ( current_node.isInternal() ) {
            current_node = current_node.getFirstChildNode();
        }
        return current_node;
    }

    public final NodeData getNodeData() {
        if ( _node_data == null ) {
            _node_data = new NodeData();
        }
        return _node_data;
    }

    final NodeData getNodeDataDirectly() {
        return _node_data;
    }

    // ---------------------------------------------------------
    // Set and get methods for Nodes
    // ---------------------------------------------------------
    /**
     * Returns the ID (int) of this PhylogenyNode.
     */
    final public int getNodeId() {
        return _id;
    }

    /**
     * Returns the name of this node.
     */
    final public String getNodeName() {
        return _node_name;
    }

    final public int getNumberOfDescendants() {
        return _descendants.size();
    }

    /**
     * Returns the total number of external Nodes originating from this
     * PhylogenyNode (int).
     */
    final public int getNumberOfExternalNodes() {
        return _sum_ext_nodes;
    }

    final public int getNumberOfParents() {
        return 1;
    }

    /**
     * Returns a refernce to the parent PhylogenyNode of this PhylogenyNode.
     */
    final public PhylogenyNode getParent() {
        return _parent;
    }

    /**
     * Returns a refernce to the next external PhylogenyNode of this
     * PhylogenyNode. TODO should be in Phylogeny. Returns null if no next
     * external node is available.
     */
    final public PhylogenyNode getPreviousExternalNode() {
        if ( isInternal() ) {
            throw new UnsupportedOperationException( "Cannot get the previous external node for an internal node." );
        }
        else if ( isRoot() /* TODO && tree is rooted */) {
            throw new UnsupportedOperationException( "Cannot get the previous external node for a root node." );
        }
        else if ( isFirstExternalNode() ) {
            throw new UnsupportedOperationException( "Attempt to get previous external node of the first external node." );
        }
        int index = getChildNodeIndex();
        PhylogenyNode previous_node = this;
        PhylogenyNode current_node = getParent();
        while ( !current_node.isRoot()
                && ( ( current_node.getNumberOfDescendants() == 1 ) || previous_node.isFirstChildNode() ) ) {
            index = current_node.getChildNodeIndex();
            previous_node = current_node;
            current_node = current_node.getParent();
        }
        current_node = current_node.getChildNode( index - 1 );
        while ( current_node.isInternal() ) {
            current_node = current_node.getLastChildNode();
        }
        return current_node;
    }

    /**
     * Used for drawing of Trees.
     */
    final public float getXcoord() {
        return _x;
    }

    final public float getXSecondary() {
        return _x_secondary;
    }

    /**
     * Used for drawing of Trees.
     */
    final public float getYcoord() {
        return _y;
    }

    final public float getYSecondary() {
        return _y_secondary;
    }

    @Override
    final public int hashCode() {
        final NodeData data = getNodeData();
        if ( ( getNodeName().length() < 1 ) && !data.isHasSequence() && !data.isHasTaxonomy() ) {
            return super.hashCode();
        }
        int result = getNodeName().hashCode();
        if ( data.isHasSequence() ) {
            result ^= data.getSequence().hashCode();
        }
        if ( data.isHasTaxonomy() ) {
            result ^= data.getTaxonomy().hashCode();
        }
        return result;
    }

    final private void init() {
        _descendants = new ArrayList<PhylogenyNode>();
        _parent = null;
        _id = 0;
        initializeData();
    }

    /**
     * Deletes data of this PhylogenyNode. Links to the other Nodes in the
     * Phylogeny, the ID and the sum of external nodes are NOT deleted. Field
     * "_link" (_link to Nodes in other Phylogeny) IS deleted.
     * 
     * @see #getLink() (Last modified: 12/20/03)
     */
    final public void initializeData() {
        _indicator = 0;
        _x = 0;
        _y = 0;
        _node_name = "";
        _distance_parent = PhylogenyNode.DISTANCE_DEFAULT;
        _collapse = false;
        _link = null;
        _branch_data = null;
        _node_data = null;
    }

    /**
     * Returns whether this PhylogenyNode should be drawn as collapsed.
     */
    final public boolean isCollapse() {
        return _collapse;
    }

    /**
     * Returns true if this PhylogenyNode represents a _duplication event, false
     * otherwise.
     */
    final public boolean isDuplication() {
        return getNodeData().isHasEvent() && getNodeData().getEvent().isDuplication();
    }

    /**
     * Checks whether this PhylogenyNode is external (tip).
     * 
     * @return true if this PhylogenyNode is external, false otherwise
     */
    final public boolean isExternal() {
        return ( getNumberOfDescendants() < 1 );
    }

    /**
     * DOCUMENT ME!
     * 
     * @return DOCUMENT ME!
     */
    final public boolean isFirstChildNode() {
        if ( isRoot() /* and tree is rooted TODO */) {
            throw new UnsupportedOperationException( "Cannot determine whether the root is the first child node of its _parent." );
        }
        return ( getChildNodeIndex() == 0 );
    }

    /**
     * DOCUMENT ME!
     * 
     * @return DOCUMENT ME!
     */
    final public boolean isFirstExternalNode() {
        if ( isInternal() ) {
            return false;
        }
        PhylogenyNode node = this;
        while ( !node.isRoot() ) {
            if ( !node.isFirstChildNode() ) {
                return false;
            }
            node = node.getParent();
        }
        return true;
    }

    /**
     * Returns whether a _duplication or speciation event has been assigned for
     * this PhylogenyNode.
     */
    final public boolean isHasAssignedEvent() {
        if ( !getNodeData().isHasEvent() ) {
            return false;
        }
        if ( ( getNodeData().getEvent() ).isUnassigned() ) {
            return false;
        }
        return true;
    }

    /**
     * Checks whether this PhylogenyNode is internal (tip).
     * 
     * @return true if this PhylogenyNode is external, false otherwise
     */
    final public boolean isInternal() {
        return ( !isExternal() );
    }

    /**
     * Returns true if this node is the last child node of its _parent.
     * <p>
     * [last modified June 01, 2005 by CMZ]
     * 
     * @return true if this node is the last child node of its _parent, false
     *         otherwise
     */
    final public boolean isLastChildNode() {
        if ( isRoot() /* and tree is rooted TODO */) {
            throw new UnsupportedOperationException( "Cannot determine whether the root is the last child node of its _parent." );
        }
        return ( getChildNodeIndex() == ( getParent().getNumberOfDescendants() - 1 ) );
    }

    /**
     * DOCUMENT ME!
     * 
     * @return DOCUMENT ME!
     */
    final public boolean isLastExternalNode() {
        if ( isInternal() ) {
            return false;
        }
        PhylogenyNode node = this;
        while ( !node.isRoot() ) {
            if ( !node.isLastChildNode() ) {
                return false;
            }
            node = node.getParent();
        }
        return true;
    }

    /**
     * Checks whether this PhylogenyNode is a root.
     * 
     * @return true if this PhylogenyNode is the root, false otherwise
     */
    final public boolean isRoot() {
        return _parent == null;
    }

    final public boolean isSpeciation() {
        return getNodeData().isHasEvent() && getNodeData().getEvent().isSpeciation();
    }

    // ---------------------------------------------------------
    // Iterator
    // ---------------------------------------------------------
    final public PhylogenyNodeIterator iterateChildNodesForward() {
        return new ChildNodeIteratorForward( this );
    }

    // ---------------------------------------------------------
    // Basic printing
    // ---------------------------------------------------------
    /**
     * Prints to the console the subtree originating from this PhylogenyNode in
     * preorder.
     */
    public void preorderPrint() {
        System.out.println( this + "\n" );
        if ( isInternal() ) {
            for( int i = 0; i < getNumberOfDescendants(); ++i ) {
                getChildNode( i ).preorderPrint();
            }
        }
    }

    final public void removeChildNode( final int i ) {
        if ( isExternal() ) {
            throw new UnsupportedOperationException( "cannot get the child node for a external node." );
        }
        if ( ( i >= getNumberOfDescendants() ) || ( i < 0 ) ) {
            throw new IllegalArgumentException( "attempt to get child node " + i + " of a node with "
                    + getNumberOfDescendants() + " child nodes." );
        }
        getDescendants().remove( i );
    }

    final public void removeChildNode( final PhylogenyNode remove_me ) {
        removeChildNode( remove_me.getChildNodeIndex() );
    }

    final public void setBranchData( final BranchData branch_data ) {
        _branch_data = branch_data;
    }

    /**
     * Sets the first child PhylogenyNode of this PhylogenyNode to n.
     */
    final public void setChild1( final PhylogenyNode n ) {
        setChildNode( 0, n );
    }

    /**
     * Sets the second child PhylogenyNode of this PhylogenyNode to n.
     */
    final public void setChild2( final PhylogenyNode n ) {
        setChildNode( 1, n );
    }

    /**
     * Inserts PhylogenyNode n at the specified position i into the list of
     * child nodes. This does not allow null slots in the list of child nodes:
     * If i is larger than the number of child nodes, n is just added to the
     * list, not place at index i.
     * 
     * @param i
     *            the index of position where to add the child
     * @param n
     *            the PhylogenyNode to add
     */
    final public void setChildNode( final int i, final PhylogenyNode node ) {
        node.setParent( this );
        if ( getNumberOfDescendants() <= i ) {
            addChildNode( node );
        }
        else {
            getDescendants().set( i, node );
        }
    }

    final void setChildNodeOnly( final int i, final PhylogenyNode node ) {
        if ( getNumberOfDescendants() <= i ) {
            addChildNode( node );
        }
        else {
            getDescendants().set( i, node );
        }
    }

    /**
     * Sets whether this PhylogenyNode should be drawn as collapsed.
     */
    final public void setCollapse( final boolean b ) {
        _collapse = b;
    }

    /**
     * Sets the length of the branch leading to the _parent of this
     * PhylogenyNode to double d.
     */
    final public void setDistanceToParent( final double d ) {
        _distance_parent = d;
    }

    /**
     * Sets the _indicator value of this PhylogenyNode to i.
     */
    final public void setIndicator( final byte i ) {
        _indicator = i;
    }

    // --------------------------------------------------------------------
    // Adjust methods (related to Phylogeny construction and
    // Phylogeny modification)
    // --------------------------------------------------------------------
    /**
     * Sets the indicators of all the children of this PhylogenyNode to zero.
     */
    final void setIndicatorsToZero() {
        for( final PreorderTreeIterator it = new PreorderTreeIterator( this ); it.hasNext(); ) {
            it.next().setIndicator( ( byte ) 0 );
        }
    }

    /**
     * Sets the linked PhylogenyNode of this PhylogenyNode to n. Currently, this
     * method is only used for the speciation-_duplication assignment
     * algorithms.
     */
    final public void setLink( final PhylogenyNode n ) {
        _link = n;
    }

    /**
     * Sets the name of this node.
     */
    final public void setName( final String node_name ) {
        _node_name = node_name;
    }

    /**
     * Sets the ID of this PhylogenyNode to i (int). In most cases, this number
     * should not be set to values lower than getNodeCount().
     */
    final protected void setNodeId( final int i ) {
        _id = i;
    }

    /**
     * Sets the _parent PhylogenyNode of this PhylogenyNode to n.
     */
    final public void setParent( final PhylogenyNode n ) {
        _parent = n;
    }

    /**
     * Sets the total number of external Nodes originating from this
     * PhylogenyNode to i (int).
     */
    final public void setSumExtNodes( final int i ) {
        _sum_ext_nodes = i;
    }

    /**
     * Used for drawing of Trees.
     */
    final public void setXcoord( final float x ) {
        _x = x;
    }

    final public void setXSecondary( final float x_secondary ) {
        _x_secondary = x_secondary;
    }

    // -----------
    /**
     * Used for drawing of Trees.
     */
    final public void setYcoord( final float y ) {
        _y = y;
    }

    final public void setYSecondary( final float y_secondary ) {
        _y_secondary = y_secondary;
    }

    // ---------------------------------------------------------
    // Writing of Nodes to Strings
    // ---------------------------------------------------------
    final public String toNewHampshire( final boolean simple_nh, final boolean write_distance_to_parent ) {
        final StringBuilder sb = new StringBuilder();
        String data = "";
        if ( !ForesterUtil.isEmpty( getNodeName() ) ) {
            data = getNodeName();
        }
        else if ( getNodeData().isHasTaxonomy() ) {
            if ( !ForesterUtil.isEmpty( getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
                data = getNodeData().getTaxonomy().getTaxonomyCode();
            }
            else if ( !ForesterUtil.isEmpty( getNodeData().getTaxonomy().getScientificName() ) ) {
                data = getNodeData().getTaxonomy().getScientificName();
            }
            else if ( !ForesterUtil.isEmpty( getNodeData().getTaxonomy().getCommonName() ) ) {
                data = getNodeData().getTaxonomy().getCommonName();
            }
            else if ( getNodeData().getTaxonomy().getTaxonomyCode() != null ) {
                data = getNodeData().getTaxonomy().getTaxonomyCode();
            }
        }
        else if ( getNodeData().isHasSequence() ) {
            if ( !ForesterUtil.isEmpty( getNodeData().getSequence().getName() ) ) {
                data = getNodeData().getSequence().getName();
            }
        }
        if ( data.length() > 0 ) {
            data = ForesterUtil.replaceIllegalNhCharacters( data );
            if ( simple_nh ) {
                if ( data.length() > 10 ) {
                    data = data.substring( 0, 11 );
                }
                if ( data.indexOf( '/' ) > 0 ) {
                    data = data.substring( 0, data.indexOf( '/' ) );
                }
                sb.append( data );
            }
            else {
                sb.append( data );
            }
        }
        if ( ( getDistanceToParent() != PhylogenyNode.DISTANCE_DEFAULT ) && write_distance_to_parent ) {
            sb.append( ":" );
            sb.append( getDistanceToParent() );
        }
        return sb.toString();
    }

    /**
     * Converts this PhylogenyNode to a New Hampshire X (NHX) String
     * representation.
     */
    final public String toNewHampshireX() {
        final StringBuffer s = new StringBuffer();
        final StringBuffer s_nhx = new StringBuffer();
        if ( !getNodeName().equals( "" ) ) {
            s.append( ForesterUtil.replaceIllegalNhxCharacters( getNodeName() ) );
        }
        if ( getDistanceToParent() != PhylogenyNode.DISTANCE_DEFAULT ) {
            s.append( ":" );
            s.append( getDistanceToParent() );
        }
        if ( getNodeDataDirectly() != null ) {
            s_nhx.append( getNodeDataDirectly().toNHX() );
        }
        if ( getBranchDataDirectly() != null ) {
            s_nhx.append( getBranchDataDirectly().toNHX() );
        }
        if ( s_nhx.length() > 0 ) {
            s.append( "[&&NHX" );
            s.append( s_nhx );
            s.append( "]" );
        }
        return s.toString();
    } // toNewHampshireX()

    @Override
    final public String toString() {
        final StringBuilder sb = new StringBuilder();
        if ( !ForesterUtil.isEmpty( getNodeName() ) ) {
            sb.append( getNodeName() );
            sb.append( " " );
        }
        sb.append( "[" );
        sb.append( getNodeId() );
        sb.append( "]" );
        return sb.toString();
    }

    /**
     * Decreases the total number of all Nodes created so far by one.
     */
    final static synchronized void decreaseNodeCount() {
        --PhylogenyNode._node_count;
    }

    /**
     * Returns the total number of all Nodes created so far.
     * 
     * @return total number of Nodes (int)
     */
    final public static int getNodeCount() {
        return PhylogenyNode._node_count;
    }

    /**
     * Increases the total number of all Nodes created so far by one.
     */
    final private static synchronized void increaseNodeCount() {
        ++PhylogenyNode._node_count;
    }

    /**
     * Sets the total number of all Nodes created so far to i (int).
     */
    final public static synchronized void setNodeCount( final int i ) {
        PhylogenyNode._node_count = i;
    }
}
