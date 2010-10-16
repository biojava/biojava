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


import java.util.HashSet;
import java.util.Set;

public class SimpleOrthoPairCollection
    extends AbstractOrthoPairCollection
{
    protected Set groups;

    public class Iterator 
        implements OrthoPairCollection.Iterator
    {
        private java.util.Iterator groupsI;

        /**
         * constructor where the iterator is already created
         */
        private Iterator(java.util.Iterator groupsI) { this.groupsI = groupsI; }

        public boolean hasNext()
        {
            return groupsI.hasNext();
        }

        public OrthoPairSet nextSet()
        {
            return (OrthoPairSet) groupsI.next();
        }
    }

    public SimpleOrthoPairCollection()
    {
        groups = new HashSet();
    }

    SimpleOrthoPairCollection(Set groups)
    {
        this.groups = groups;
    }

    public void add(OrthoPairSet group)
    {
        groups.add(group);
    }

    public boolean contains(OrthoPairSet group)
    {
        return groups.contains(group);
    }

    public boolean isEmpty() { return groups.isEmpty(); }

    public OrthoPairCollection.Iterator iterator()
    {
        return new Iterator(groups.iterator());
    }
/*
    public OrthoPairCollection filter(OrthoPairSetFilter filters)
    {
        OrthoPairCollection results = new SimpleOrthoPairCollection();

        // this method uses its privileged access to groups
        for (java.util.Iterator groupsI = groups.iterator();
               groupsI.hasNext(); )
        {
            OrthoPairSet group = (OrthoPairSet) groupsI.next();

            if (filters.accept(group)) {
                try {
                    results.add(group);
                }
                catch (ChangeVetoException cve) {
                    // should be impossible as this group was created by me
                }
            }
        }
        return results;
    }
*/
}

