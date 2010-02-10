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

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.RealizingFeatureHolder;
import org.biojava.bio.seq.Sequence;

/**
 * Sequence keyed off a BioSQL biosequence.
 *
 * @author Thomas Down
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 * @since 1.3
 */

interface BioSQLSequenceI extends Sequence, RealizingFeatureHolder {
    void persistFeature(Feature f, int parent_id)
        throws BioException;
    
    BioSQLSequenceDB getSequenceDB();

    int getBioEntryID();
}
