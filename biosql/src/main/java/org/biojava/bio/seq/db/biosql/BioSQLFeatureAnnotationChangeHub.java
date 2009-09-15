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
import org.biojava.bio.seq.Feature;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.IndexedChangeHub;

/**
 * Handles ChangeEvents for BioSQLFeatureAnnotation instances.
 *
 * @author Thomas Down (original implementation)
 * @author David Huen (refactoring)
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */

class BioSQLFeatureAnnotationChangeHub extends IndexedChangeHub
{
    BioSQLSequenceDB seqDB;
    BioSQLFeatureChangeHub featureHub;

    public BioSQLFeatureAnnotationChangeHub(BioSQLSequenceDB seqDB, BioSQLFeatureChangeHub featureHub)
    {
        this.seqDB = seqDB;
        this.featureHub = featureHub;
    }

    protected final boolean isMyChangeEvent(ChangeEvent cev, ListenerMemento lm)
    {
        ChangeType ct = cev.getType();
        return ct.isMatchingType(lm.type);
    }

    void firePreChange(ChangeEvent cev)
        throws ChangeVetoException
    {
        BioSQLFeatureAnnotation source = (BioSQLFeatureAnnotation) cev.getSource();
        Integer feature_id = new Integer(source.getFeatureID());

        super.firePreChange(feature_id, cev);

        Feature parent = seqDB.getFeatureByID(feature_id.intValue());
        ChangeEvent pcev = new ChangeEvent(parent, Annotatable.ANNOTATION, null, null, cev);
        featureHub.firePreChange(pcev);
    }

    void firePostChange(ChangeEvent cev)
    {
        BioSQLFeatureAnnotation source = (BioSQLFeatureAnnotation) cev.getSource();
        Integer feature_id = new Integer(source.getFeatureID());

        super.firePostChange(feature_id, cev);

        Feature parent = seqDB.getFeatureByID(feature_id.intValue());
        ChangeEvent pcev = new ChangeEvent(parent, Annotatable.ANNOTATION, null, null, cev);
        featureHub.firePostChange(pcev);
    }
}

