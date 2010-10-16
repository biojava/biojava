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
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.BioEntry;


/**
 * An implementation of RichSequenceDB that uses an underlying HashMap to store the
 * RichSequence objects.
 *
 * @author Matthew Pocock
 * @author Gerald Loeffler
 * @author Richard Holland
 * @since 1.5
 */
public class HashBioEntryDB extends AbstractBioEntryDB implements BioEntryDB {
    
    /**
     * The sequence-by-id map.
     */
    final private Map sequenceByID;
        
    /**
     * The name of this sequence database.
     */
    private String name;
    
    /**
     * Generate a HashRichSequenceDB object that will use byName to generate ids for
     * sequences and have a null name.
     */
    public HashBioEntryDB() {
        this(null);
    }
    
    /**
     * Generate a HashRichSequenceDB object that will use byName to generate ids and
     * will have the requested name.
     *
     * @param name the name for this database
     */
    public HashBioEntryDB(String name) {
        this.name = name;
        this.sequenceByID = new HashMap();
    }
    
    public String getName() {
        return name;
    }
    
    public BioEntry getBioEntry(String id) throws BioException, IllegalIDException {
        BioEntry seq = (BioEntry)sequenceByID.get(id);
        if (seq == null) throw new IllegalIDException("BioEntry with ID " + id + " could not be found");
        return seq;
    }
    
    public BioEntryDB getBioEntrys(Set ids) throws BioException, IllegalIDException {
        return this.getBioEntrys(ids,null);
    }
    
    public BioEntryDB getBioEntrys(Set ids, BioEntryDB db) throws BioException, IllegalIDException {
        if (db==null) db = new HashBioEntryDB();
        for (Iterator i = ids.iterator(); i.hasNext(); ) {
            String id = (String)i.next();
            if (!sequenceByID.containsKey(id)) throw new IllegalIDException("BioEntry with ID " + id + " could not be found");
            else {
                try {
                    db.addBioEntry((BioEntry)sequenceByID.get(id));
                } catch (ChangeVetoException ce) {
                    throw new BioException("Unexpectedly couldn't add to a HashBioEntryDB", ce);
                }
            }
        }
        return db;
    }
    
    public Set ids() {
        return sequenceByID.keySet();
    }
    
    /**
     * Add a BioEntry, the name of the BioEntry will be used as the ID
     * @param seq the BioEntry to add
     * @throws ChangeVetoException if this addition was vetoed
     */
    public void addBioEntry(BioEntry seq) throws IllegalIDException, BioException, ChangeVetoException {
        this.addBioEntry(seq.getName(), seq);
    }
    
    protected void addBioEntry(String id, BioEntry seq) throws IllegalIDException, BioException, ChangeVetoException {
        if(!hasListeners(BioEntryDB.BIOENTRYS)) {
            sequenceByID.put(id, seq);
        } else {
            ChangeSupport changeSupport = getChangeSupport(BioEntryDB.BIOENTRYS);
            synchronized(changeSupport) {
                ChangeEvent ce = new ChangeEvent(
                        this,
                        BioEntryDB.BIOENTRYS,
                        new Object[] { id, seq },
                        null);
                        changeSupport.firePreChangeEvent(ce);
                        sequenceByID.put(id, seq);
                        changeSupport.firePostChangeEvent(ce);
            }
        }
    }
    
    public void removeBioEntry(String id) throws IllegalIDException, BioException, ChangeVetoException {
        if (!sequenceByID.containsKey(id)) throw new IllegalIDException("BioEntry with ID " + id + " could not be found");
        if(!hasListeners(BioEntryDB.BIOENTRYS)) {
            sequenceByID.remove(id);
        } else {
            ChangeSupport changeSupport = getChangeSupport(BioEntryDB.BIOENTRYS);
            synchronized(changeSupport) {
                ChangeEvent ce = new ChangeEvent(
                        this,
                        BioEntryDB.BIOENTRYS,
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
