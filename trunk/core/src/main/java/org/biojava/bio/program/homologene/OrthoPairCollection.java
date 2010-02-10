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

package org.biojava.bio.program.homologene;

import org.biojava.utils.ChangeVetoException;

/**
 * Interface for a Set of OrthoPairSets
 *
 * @author David Huen
 */
public interface OrthoPairCollection
{

    /**
     * Iterator for a OrthoPairCollection
     */
    public interface Iterator
    {
        /**
         * are more OrthoPairSets available?
         */
        public boolean hasNext();

        /**
         * returns the next OrthoPairSet
         */
        public OrthoPairSet nextSet();
    }


    public void add(OrthoPairSet group) throws ChangeVetoException;

    public boolean contains(OrthoPairSet group);

    public boolean isEmpty();

    public Iterator iterator();

    public OrthoPairCollection filter(OrthoPairSetFilter filters);
}

