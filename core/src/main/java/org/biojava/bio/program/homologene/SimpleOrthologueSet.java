/*
 *                    BioJava development code
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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

public class SimpleOrthologueSet extends AbstractOrthologueSet
{

    public class Iterator implements OrthologueSet.Iterator
    {
        private java.util.Iterator orthoIterator;

        private Iterator(java.util.Iterator orthoIterator)
        {
            this.orthoIterator = orthoIterator;
        }

        public boolean hasNext()
        {
            return orthoIterator.hasNext();
        }

        public Orthologue nextOrthologue()
        {
            return (Orthologue) orthoIterator.next();
        }

    }

    // every Orthologue is stored in a Set
    private Set orthologueSet = new HashSet();
    private Map orthologueByHomologeneID = new HashMap();

    {
        generateChangeSupport();
    }

    public void addOrthologue(Orthologue ortho)
        throws ChangeVetoException
    {
        if (!hasListeners()) {
            orthologueSet.add(ortho);
        }
        else {
            // get the change support
            ChangeSupport cs = getChangeSupport(OrthologueSet.MODIFY);

            synchronized(cs) {
                ChangeEvent ce = new ChangeEvent(this, OrthologueSet.MODIFY);
                cs.firePreChangeEvent(ce);
                orthologueSet.add(ortho);
                cs.firePostChangeEvent(ce);
            }
        }
    }

    public void removeOrthologue(Orthologue ortho)
        throws ChangeVetoException
    {
        if (!hasListeners()) {
            orthologueSet.remove(ortho);
        }
        else {
            // get the change support
            ChangeSupport cs = getChangeSupport(OrthologueSet.MODIFY);

            synchronized(cs) {
                ChangeEvent ce = new ChangeEvent(this, OrthologueSet.MODIFY);
                cs.firePreChangeEvent(ce);
                orthologueSet.remove(ortho);
                cs.firePostChangeEvent(ce);
            }
        }
    }

    public Orthologue getOrthologue(String homologeneID)
    {
        return (Orthologue) orthologueByHomologeneID.get(homologeneID);
    }

    public OrthologueSet.Iterator iterator()
    {
        return new Iterator(orthologueSet.iterator());
    }
}

