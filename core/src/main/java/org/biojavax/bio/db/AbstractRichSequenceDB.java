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
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.MergeFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.BioEntry;
import org.biojavax.bio.BioEntryIterator;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 * An abstract implementation of RichSequenceDB that provides the getRichSequenceIterator
 * method.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Richard Holland
 * @since 1.5
 */
public abstract class AbstractRichSequenceDB extends AbstractBioEntryDB implements RichSequenceDB {
    
    public Sequence getSequence(String id) throws BioException, IllegalIDException {
        return this.getRichSequence(id);
    }
    
    public BioEntry getBioEntry(String id) throws BioException, IllegalIDException {
        return this.getRichSequence(id);
    }
    
    public BioEntryDB getBioEntrys(Set ids) throws BioException, IllegalIDException {
        return this.getRichSequences(ids);
    }
    
    public BioEntryDB getBioEntrys(Set ids, BioEntryDB db) throws BioException, IllegalIDException {
        RichSequenceDB db2 = this.getRichSequences(ids);
        for (BioEntryIterator i = db2.getBioEntryIterator(); i.hasNext(); ) {
            BioEntry be = i.nextBioEntry();
            try {
                db.addBioEntry(be);
            } catch (ChangeVetoException ce) {
                throw new BioException("Unexpectedly unable to add to a BioEntryDB", ce);
            }
        }
        return db;
    }
    
    public void addSequence(Sequence seq) throws IllegalIDException, BioException, ChangeVetoException {
        this.addRichSequence(RichSequence.Tools.enrich(seq));
    }
    
    public void removeSequence(String id) throws IllegalIDException, BioException, ChangeVetoException {
        this.removeRichSequence(id);
    }
    
    public void addBioEntry(BioEntry seq) throws IllegalIDException, BioException, ChangeVetoException {
        throw new ChangeVetoException("Cannot add BioEntrys to a RichSequence database - use addRichSequence");
    }
    
    public void removeBioEntry(String id) throws IllegalIDException, BioException, ChangeVetoException {
        throw new ChangeVetoException("Cannot remove BioEntrys from a RichSequence database - use addRichSequence");
    }
    
    public void addRichSequence(RichSequence seq) throws IllegalIDException, BioException, ChangeVetoException {
        throw new ChangeVetoException("Cannot add RichSequences to a read-only database");
    }
    
    public void removeRichSequence(String id) throws IllegalIDException, BioException, ChangeVetoException {
        throw new ChangeVetoException("Cannot remove RichSequences from a read-only database");
    }
    
    public SequenceIterator sequenceIterator() {
        return this.getRichSequenceIterator();
    }
    
    public BioEntryIterator getBioEntryIterator() {
        return this.getRichSequenceIterator();
    }
    
    public RichSequenceIterator getRichSequenceIterator() {
        return new RichSequenceIterator() {
            private Iterator pID = ids().iterator();
            
            public boolean hasNext() {
                return pID.hasNext();
            }
            
            public Sequence nextSequence() throws BioException {
                return nextRichSequence();
            }
            
            public BioEntry nextBioEntry() throws BioException {
                return nextRichSequence();
            }
            
            public RichSequence nextRichSequence() throws BioException {
                return getRichSequence((String)pID.next());
            }
        };
    }
    
    public FeatureHolder filter(FeatureFilter ff) {
        // Default behaviour is accept-all.
        MergeFeatureHolder results = new MergeFeatureHolder();
        try {
            for (RichSequenceIterator si = getRichSequenceIterator(); si.hasNext(); ) {
                RichSequence seq = si.nextRichSequence();
                results.addFeatureHolder(seq.filter(ff));
            }
        } catch (BioException ex) {
            throw new BioRuntimeException(ex);
        } catch (ChangeVetoException cve) {
            throw new BioError("Assertion failed: couldn't modify newly created MergeFeatureHolder",cve);
        }
        return results;
    }
}
