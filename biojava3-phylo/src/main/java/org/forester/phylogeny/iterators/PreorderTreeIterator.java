// $Id: PreorderTreeIterator.java,v 1.13 2009/10/26 23:29:39 cmzmasek Exp $
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
import java.util.Stack;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

// import java.util.Iterator; TODO should implement this, not some iterator of
// this package.
/*
 * @author Christian M. Zmasek
 * 
 * @version 1.020 -- last modified: 10/10/05
 */
public class PreorderTreeIterator implements PhylogenyNodeIterator {

    final private Phylogeny            _tree;
    final private Stack<PhylogenyNode> _stack;

    /**
     * @param tree
     *            Phylogeny for which a Iterator is to be constructed.
     */
    public PreorderTreeIterator( final Phylogeny tree ) throws IllegalArgumentException {
        if ( tree.isEmpty() ) {
            throw new IllegalArgumentException( "Attempt to use PreorderTreeIterator on empty tree." );
        }
        _stack = new Stack<PhylogenyNode>();
        _tree = tree;
        reset();
    }

    public PreorderTreeIterator( final PhylogenyNode node ) throws IllegalArgumentException {
        _stack = new Stack<PhylogenyNode>();
        _tree = null;
        reset( node );
    }

    private Stack<PhylogenyNode> getStack() {
        return _stack;
    }

    private Phylogeny getTree() {
        return _tree;
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.util.Iterator#hasNext()
     */
    public boolean hasNext() {
        return !getStack().isEmpty();
    }

    /**
     * Advances the Iterator by one.
     */
    public PhylogenyNode next() throws NoSuchElementException {
        if ( !hasNext() ) {
            throw new NoSuchElementException( "Attempt to call \"next()\" on iterator which has no more next elements." );
        }
        final PhylogenyNode node = getStack().pop();
        if ( !node.isExternal() ) {
            for( int i = node.getNumberOfDescendants() - 1; i >= 0; --i ) {
                getStack().push( node.getChildNode( i ) );
            }
        }
        return node;
    } // next()

    /**
     * Not supported.
     * 
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

    public void reset() {
        getStack().clear();
        getStack().push( getTree().getRoot() );
    }

    private void reset( final PhylogenyNode node ) {
        getStack().clear();
        getStack().push( node );
    }
} // End of class PreorderTreeIterator.
