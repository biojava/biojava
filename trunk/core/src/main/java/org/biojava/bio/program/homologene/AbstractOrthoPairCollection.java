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
 * An abstract implementation of the OrthoPairCollection
 * interface.  Its primary role is to implement the
 * filter() method.i
 *
 * @author David Huen
 */
public abstract class AbstractOrthoPairCollection
    implements OrthoPairCollection
{

    public abstract void add(OrthoPairSet group);

    public abstract boolean contains(OrthoPairSet group);

    public abstract boolean isEmpty();

    public abstract OrthoPairCollection.Iterator iterator();

    public OrthoPairCollection filter(OrthoPairSetFilter filters)
    {
        OrthoPairCollection results = new SimpleOrthoPairCollection();

        for (Iterator pairSetsI = iterator();
               pairSetsI.hasNext(); )
        {
            OrthoPairSet pairSet = pairSetsI.nextSet();

            if (filters.accept(pairSet)) {
                try {
                    results.add(pairSet);
                }
                catch (ChangeVetoException cve) {
                    // should be impossible as this group was created by me
                }
            }
        }
        return results;
    }
}

