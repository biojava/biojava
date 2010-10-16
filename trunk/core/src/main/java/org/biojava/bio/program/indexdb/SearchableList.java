/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package org.biojava.bio.program.indexdb;

import java.util.Comparator;
import java.util.List;

import org.biojava.utils.Commitable;

/**
 * <code>SearchableList</code>s are ordered <code>List</code>s which
 * may be searched using an algorithm consistent with the sort
 * order. For example, using a binary search.
 *
 * @author Matthew Pocock
 * @author Keith James
 */
interface SearchableList
    extends
        List,
        Commitable
{
    /**
     * <code>search</code> returns an <code>Object</code> from the
     * list corresponding to the specified identifier.
     *
     * @param id a <code>String</code> ID.
     *
     * @return an <code>Object</code>.
     */
    public Object search(String id);
  
    /**
     * <code>searchAll</code> returns all <code>Object</code>s from
     * the list corresponding to the specified identifier.
     *
     * @param id a <code>String</code> ID.
     *
     * @return a <code>List</code> of <code>Object</code>s.
     */
    public List searchAll(String id);
  
    /**
     * <code>getComparator</code> returns the <code>Comparator</code>
     * responsible for maintaining list order.
     *
     * @return a <code>Comparator</code>.
     */
    public Comparator getComparator();
}
