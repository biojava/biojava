// $Id: LevelOrderTreeIterator.java,v 1.9 2009/10/26 23:29:39 cmzmasek Exp $
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

import org.forester.datastructures.Queue;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/*
 * An iterator to iterate a Phylogeny in level order.
 * 
 * Created: 10/23/2005 by Christian M. Zmasek. Last modified: 10/23/2005 by
 * Christian M. Zmasek.
 * 
 * @author Christian M. Zmasek
 * 
 * @version 1.000
 */
public class LevelOrderTreeIterator implements PhylogenyNodeIterator {

    // Instance variables
    // ------------------
    private final Queue         _queue;
    private final PhylogenyNode _root;

    // Constructors
    // ------------
    /**
     * Creates a new LevelOrderTreeIterator for iterating over all the nodes of
     * Phylogeny phylogeny
     * 
     * @param phylogeny
     *            the Phylogeny to iterate over
     * @throws IllegalArgumentException
     *             if phylogeny is empty
     */
    public LevelOrderTreeIterator( final Phylogeny phylogeny ) throws IllegalArgumentException {
        this( phylogeny.getRoot() );
        if ( phylogeny.isEmpty() ) {
            throw new IllegalArgumentException( "Attempt to use LevelOrderTreeIterator on an empty phylogeny." );
        }
    }

    /**
     * Creates a new LevelOrderTreeIterator for iterating over all the child
     * nodes of PhylogenyNode node (including node itself).
     * 
     * @param node
     *            the parent of the nodes to iterate over
     */
    public LevelOrderTreeIterator( final PhylogenyNode node ) {
        _queue = new Queue();
        _root = node;
        reset();
    }

    // Private methods
    // ---------------
    /**
     * Returns the queue upon which this iterator is based.
     * 
     */
    private Queue getQueue() {
        return _queue;
    }

    /**
     * Returns the root of the phylogeny this iterators parses over.
     * 
     * @return the root of the phylogeny this iterators parses over.
     */
    private PhylogenyNode getRoot() {
        return _root;
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
        return !getQueue().isEmpty();
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
        final PhylogenyNode node = ( PhylogenyNode ) getQueue().dequeue();
        for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            getQueue().enqueue( node.getChildNode( i ) );
        }
        return node;
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
        getQueue().clear();
        getQueue().enqueue( getRoot() );
    }
} // enod of class LevelOrderTreeIterator
