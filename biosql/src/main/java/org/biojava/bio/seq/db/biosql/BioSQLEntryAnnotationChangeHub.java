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

import org.biojava.bio.Annotatable;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.Sequence;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.IndexedChangeHub;

/**
 * Handles ChangeEvents for BioSQLEntryAnnotation instances.
 *
 * @author Thomas Down (original implementation)
 * @author David Huen (refactoring)
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */
class BioSQLEntryAnnotationChangeHub extends IndexedChangeHub
{
    BioSQLSequenceDB seqDB;
    BioSQLEntryChangeHub entryHub;

    public BioSQLEntryAnnotationChangeHub(BioSQLSequenceDB seqDB, BioSQLEntryChangeHub entryHub)
    {
        this.seqDB = seqDB;
        this.entryHub = entryHub;
    }

    protected final boolean isMyChangeEvent(ChangeEvent cev, ListenerMemento lm)
    {
        ChangeType ct = cev.getType();
        return ct.isMatchingType(lm.type);
    }

    void firePreChange(ChangeEvent cev)
        throws ChangeVetoException
    {
        BioSQLSequenceAnnotation source = (BioSQLSequenceAnnotation) cev.getSource();
        Integer bioentry_id = new Integer(source.getBioentryID());

        super.firePreChange(bioentry_id, cev);

        try {
            Sequence seq = seqDB.getSequence(null, bioentry_id.intValue());
            ChangeEvent pcev = new ChangeEvent(seq, Annotatable.ANNOTATION, null, null, cev);
            entryHub.firePreChange(pcev);
        } catch (BioException ex) {
            throw new BioRuntimeException("Sequence has gone missing");
        }
    }

    void firePostChange(ChangeEvent cev)
    {
        BioSQLSequenceAnnotation source = (BioSQLSequenceAnnotation) cev.getSource();
        Integer bioentry_id = new Integer(source.getBioentryID());

        super.firePostChange(bioentry_id, cev);

        // omitted in original BioSQLChangeHub implementation but appears obvious
        try {
            Sequence seq = seqDB.getSequence(null, bioentry_id.intValue());
            ChangeEvent pcev = new ChangeEvent(seq, Annotatable.ANNOTATION, null, null, cev);
            entryHub.firePostChange(pcev);
        } catch (BioException ex) {
            throw new BioRuntimeException("Sequence has gone missing");
        }
    }
}

