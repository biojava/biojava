// $Id: ChildNodeIteratorForward.java,v 1.9 2009/10/26 23:29:39 cmzmasek Exp $
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

package org.forester.phylogeny.iterators;

import java.util.NoSuchElementException;

import org.forester.phylogeny.PhylogenyNode;

/*
 * An iterator to forward iterate over child nodes of a PhylogenyNode. Created:
 * 10/23/2005 by Christian M. Zmasek. Last modified: 12/28/2006 by Christian M.
 * Zmasek.
 * 
 * @author Christian M. Zmasek
 * 
 * @version 1.000
 */
public class ChildNodeIteratorForward implements PhylogenyNodeIterator {

    // Instance variables
    // ------------------
    private int                 _i;
    final private PhylogenyNode _node;

    // Constructor
    // -----------
    /**
     * Creates a new ChildNodeIteratorForward.
     * 
     * @param node
     *            the parent of the PhylogenyNodes to iterate over.
     * @throws IllegalArgumentException
     *             if node has no child nodes
     */
    public ChildNodeIteratorForward( final PhylogenyNode node ) throws IllegalArgumentException {
        if ( node.getNumberOfDescendants() < 1 ) {
            throw new IllegalArgumentException( "Attempt to use ChildNodeIteratorForward on node with no child nodes." );
        }
        _node = node;
        reset();
    }

    // Private methods
    // ---------------
    /**
     * Returns the counter.
     */
    private int getI() {
        return _i;
    }

    /**
     * Returns the parent of the nodes to iterate over.
     * 
     * @return the parent of the nodes to iterate over.
     */
    private PhylogenyNode getNode() {
        return _node;
    }

    // Public methods
    // --------------
    /**
     * Returns true is this iterator has at least one more element, false
     * otherwise.
     * 
     * @return true is this iterator has at least one more element, false
     *         otherwise
     */
    public boolean hasNext() {
        return ( getI() < getNode().getNumberOfDescendants() );
    }

    /**
     * Increases the counter by one.
     */
    private void increaseI() {
        ++_i;
    }

    /**
     * Returns the next PhylogenyNode.
     * 
     * @return the next PhylogenyNode
     * @throws NoSuchElementException
     *             if iteration is complete
     */
    public PhylogenyNode next() throws NoSuchElementException {
        if ( !hasNext() ) {
            throw new NoSuchElementException( "Attempt to call \"next()\" on iterator which has no more next elements." );
        }
        final PhylogenyNode n = getNode().getChildNode( getI() );
        increaseI();
        return n;
    }

    /**
     * Not supported.
     * 
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Resets the iterator.
     */
    public void reset() {
        setI( 0 );
    }

    /**
     * Sets the counter.
     */
    private void setI( final int i ) {
        _i = i;
    }
} // end of class ChildNodeIteratorForward.
