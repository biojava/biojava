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


package org.biojava.bio.seq.db.biosql;

import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.IndexedChangeHub;

/**
 * Handles ChangeEvents for BioSQLEntry instances.
 *
 * @author Thomas Down (original implementation)
 * @author David Huen (refactoring)
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */
class BioSQLEntryChangeHub extends IndexedChangeHub
{
    BioSQLSequenceDB seqDB;

    public BioSQLEntryChangeHub(BioSQLSequenceDB seqDB)
    {
        super();
        this.seqDB = seqDB;
    }

    protected final boolean isMyChangeEvent(ChangeEvent cev, ListenerMemento lm)
    {
        ChangeType ct = cev.getType();
        return ct.isMatchingType(lm.type);
    }

    void addListener(int bioentry_id,
                     ChangeListener listener,
                     ChangeType ct)
    {
        Integer id = new Integer(bioentry_id);

        super.addListener(id, listener, ct);
    }

    void removeListener(int bioentry_id,
                            ChangeListener listener,
                            ChangeType ct)
    {
        Integer id = new Integer(bioentry_id);

        super.removeListener(id, listener, ct);
    }


    void firePreChange(ChangeEvent cev)
        throws ChangeVetoException
    {
        BioSQLSequenceI source = (BioSQLSequenceI) cev.getSource();
        Integer id = new Integer(source.getBioEntryID());

        super.firePreChange(id, cev);

        ChangeEvent pcev = new ChangeEvent(seqDB, SequenceDB.SEQUENCES, null, null, cev);
        seqDB.firePreChangeEvent(pcev);
    }

    void firePostChange(ChangeEvent cev)
    {
        BioSQLSequenceI source = (BioSQLSequenceI) cev.getSource();
        Integer id = new Integer(source.getBioEntryID());

        super.firePostChange(id, cev);

        ChangeEvent pcev = new ChangeEvent(seqDB, SequenceDB.SEQUENCES, null, null, cev);
        seqDB.firePostChangeEvent(pcev);
    }
}

