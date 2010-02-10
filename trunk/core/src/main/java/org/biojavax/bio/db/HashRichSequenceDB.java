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
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.IDMaker;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.seq.RichSequence;


/**
 * An implementation of RichSequenceDB that uses an underlying HashMap to store the
 * RichSequence objects.
 *
 * @author Matthew Pocock
 * @author Gerald Loeffler
 * @author Richard Holland
 * @since 1.5
 */
public class HashRichSequenceDB extends AbstractRichSequenceDB implements RichSequenceDB {
    
    /**
     * The sequence-by-id map.
     */
    final private Map sequenceByID;
    
    /**
     * An object to extract an ID for a sequence.
     */
    final private IDMaker idMaker;
    
    /**
     * The name of this sequence database.
     */
    private String name;
    
    /**
     * Generate a HashRichSequenceDB object that will use byName to generate ids for
     * sequences and have a null name.
     */
    public HashRichSequenceDB() {
        this(IDMaker.byName, null);
    }
    
    /**
     * Generate a HashRichSequenceDB object that will use idMaker to generate ids for
     * sequences and have a null name.
     *
     * @param idMaker the object that will work out the default id for a sequence
     */
    public HashRichSequenceDB(IDMaker idMaker) {
        this(idMaker, null);
    }
    
    /**
     * Generate a HashRichSequenceDB object that will use byName to generate ids and
     * will have the requested name.
     *
     * @param name the name for this database
     */
    public HashRichSequenceDB(String name) {
        this(IDMaker.byName, name);
    }
    
    /**
     * Generate a HashRichSequenceDB object that will use idMaker to generate ids for
     * sequences and have the requested name.
     *
     * @param idMaker the object that will work out the default id for a sequence
     * @param name the name for this database
     */
    public HashRichSequenceDB(IDMaker idMaker, String name) {
        this.idMaker = idMaker;
        this.name = name;
        this.sequenceByID = new HashMap();
    }
    
    public String getName() {
        return name;
    }
    
    /**
     * Retrieve the IDMaker associated with this database.
     *
     * @return the current IDMaker object
     */
    public IDMaker getIDMaker() {
        return idMaker;
    }
    
    public RichSequence getRichSequence(String id) throws BioException, IllegalIDException {
        RichSequence seq = (RichSequence)sequenceByID.get(id);
        if (seq == null) throw new IllegalIDException("Sequence with ID " + id + " could not be found");
        return seq;
    }
    
    public RichSequenceDB getRichSequences(Set ids) throws BioException, IllegalIDException {
        return this.getRichSequences(ids,null);
    }
    
    public RichSequenceDB getRichSequences(Set ids, RichSequenceDB db) throws BioException, IllegalIDException {
        if (db==null) db = new HashRichSequenceDB();
        for (Iterator i = ids.iterator(); i.hasNext(); ) {
            String id = (String)i.next();
            if (!sequenceByID.containsKey(id)) throw new IllegalIDException("Sequence with ID " + id + " could not be found");
            else {
                try {
                    db.addSequence((RichSequence)sequenceByID.get(id));
                } catch (ChangeVetoException ce) {
                    throw new BioException("Unexpectedly couldn't add to a HashRichSequenceDB", ce);
                }
            }
        }
        return db;
    }
    
    public Set ids() {
        return sequenceByID.keySet();
    }
    
    /**
     * Add a sequence.
     * @param seq the RichSequence to add
     * @throws ChangeVetoException if this addition was vetoed
     */
    public void addRichSequence(RichSequence seq) throws IllegalIDException, BioException, ChangeVetoException {
        this.addRichSequence(idMaker.calcID(seq), seq);
    }
    
    protected void addRichSequence(String id, RichSequence seq) throws IllegalIDException, BioException, ChangeVetoException {
        if(!hasListeners(RichSequenceDB.SEQUENCES)) {
            sequenceByID.put(id, seq);
        } else {
            ChangeSupport changeSupport = getChangeSupport(RichSequenceDB.SEQUENCES);
            synchronized(changeSupport) {
                ChangeEvent ce = new ChangeEvent(
                        this,
                        RichSequenceDB.SEQUENCES,
                        new Object[] { id, seq },
                        null);
                        changeSupport.firePreChangeEvent(ce);
                        sequenceByID.put(id, seq);
                        changeSupport.firePostChangeEvent(ce);
            }
        }
    }
    
    public void removeSequence(String id) throws IllegalIDException, BioException, ChangeVetoException {
        if (!sequenceByID.containsKey(id)) throw new IllegalIDException("Sequence with ID " + id + " could not be found");
        if(!hasListeners(RichSequenceDB.SEQUENCES)) {
            sequenceByID.remove(id);
        } else {
            ChangeSupport changeSupport = getChangeSupport(RichSequenceDB.SEQUENCES);
            synchronized(changeSupport) {
                ChangeEvent ce = new ChangeEvent(
                        this,
                        RichSequenceDB.SEQUENCES,
                        null,
                        id
                        );
                changeSupport.firePreChangeEvent(ce);
                sequenceByID.remove(id);
                changeSupport.firePostChangeEvent(ce);
            }
        }
    }
}
