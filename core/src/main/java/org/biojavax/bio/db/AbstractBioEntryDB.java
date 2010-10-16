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

package org.biojavax.bio.db;

import java.util.Iterator;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.BioEntry;
import org.biojavax.bio.BioEntryIterator;

/**
 * An abstract implementation of BioEntryDB that provides the getBioEntryIterator
 * method.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Richard Holland
 * @since 1.5
 */
public abstract class AbstractBioEntryDB extends AbstractChangeable implements BioEntryDB {
    
    public void addBioEntry(BioEntry seq) throws IllegalIDException, BioException, ChangeVetoException {
        throw new ChangeVetoException("Cannot add BioEntrys to a read-only database");
    }
    
    public void removeBioEntry(String id) throws IllegalIDException, BioException, ChangeVetoException {
        throw new ChangeVetoException("Cannot remove BioEntrys from a read-only database");
    }
           
    public BioEntryIterator getBioEntryIterator() {
        return new BioEntryIterator() {
            private Iterator pID = ids().iterator();
            
            public boolean hasNext() {
                return pID.hasNext();
            }
            
            public BioEntry nextBioEntry() throws BioException {
                return getBioEntry((String)pID.next());
            }
        };
    }
}
