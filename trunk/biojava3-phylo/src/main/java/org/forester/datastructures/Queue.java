// $Id: Queue.java,v 1.8 2009/10/26 23:29:41 cmzmasek Exp $
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

package org.forester.datastructures;

import java.util.LinkedList;
import java.util.NoSuchElementException;

/*
 * A simple Queue data structure. Created: 10/23/2005 by Christian M. Zmasek.
 * Last modified: 10/23/2005 by Christian M. Zmasek.
 * 
 * @author Christian M. Zmasek
 * 
 * @version 1.000
 */
public class Queue {

    // Instance variables
    // ------------------
    private final LinkedList<Object> _data;

    // Constructor
    // -----------
    /**
     * This created a new, empty Queue object.
     */
    public Queue() {
        _data = new LinkedList<Object>();
    }

    /**
     * Removes all elements from this queue.
     */
    public void clear() {
        getData().clear();
    }

    /**
     * Dequeues one element from this queue.
     * 
     * @return the dequeued object
     * @throws NoSuchElementException
     *             if this queue is empty
     */
    public Object dequeue() throws NoSuchElementException {
        if ( isEmpty() ) {
            throw new NoSuchElementException( "Attempt to dequeue from an empty Queue." );
        }
        return getData().removeFirst();
    }

    // Public methods
    // --------------
    /**
     * Adds Object element to thisqueue.
     * 
     * @param element
     *            the Object to be enqueued
     */
    public void enqueue( final Object element ) {
        getData().add( element );
    }

    // Private methods
    // ---------------
    /**
     * Returns the LinkedList upon which this queue is based.
     * 
     * @return the LinkedList upon which this queue is based
     */
    private LinkedList<Object> getData() {
        return _data;
    }

    /**
     * Returns whether or not this queue is empty.
     * 
     * @return true if this queue is empty, false otherwise
     */
    public boolean isEmpty() {
        return getData().isEmpty();
    }
} // end of class Queue.
