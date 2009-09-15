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
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.BioEntry;

/**
 * A database of BioEntrys. This may have several implementations with
 * rich behaviour, but basically most of the time you will just use
 * the interface methods to do stuff. A BioEntry database contains a
 * finite number of BioEntrys stored under unique keys.
 *
 * @author Matthew Pocock
 * @author Gerald Loeffler
 * @author Thomas Down
 * @author Richard Holland
 * @since 1.5
 */
public interface BioEntryDBLite {
    /**
     * Signals that sequences are being added to or remove from the database.
     * The sequences being removed should be listed in the previous field by
     * id, either as a single String, an array or a Set. The sequences
     * being added should be listed in the change field as either an array
     * Object[] { id, seq}, or a Map of id->seq.
     */
    public static final ChangeType BIOENTRYS = new ChangeType(
            "BioEntrys have been added or removed from the database",
            "org.biojavax.bio.db.BioEntryDB",
            "BIOENTRYS"
            );
    
    /**
     * Get the name of this sequence database.
     *
     * @return the name of the sequence database, which may be null.
     */
    public String getName();
    
    /**
     * Retrieve a single BioEntry by its id.
     *
     * @param id the id to retrieve by
     * @return  the BioEntry with that id
     * @throws IllegalIDException if the database doesn't know about the id
     * @throws BioException if there was a failure in retrieving the BioEntry
     */
    public BioEntry getBioEntry(String id) throws IllegalIDException, BioException;
    
    /**
     * Retrieve multiple BioEntry by their ids.
     *
     * @param ids a set of ids to retrieve by
     * @return  the BioEntrys with those ids
     * @throws IllegalIDException if the database doesn't know about the id
     */
    public BioEntryDB getBioEntrys(Set ids) throws BioException,IllegalIDException;
    
    /**
     * Retrieve multiple BioEntry into a specific sequence database. If
     * that database is null, a new HashBioEntryDB is used.
     *
     * @param ids a set of ids to retrieve by
     * @param db a database to load the seqs into
     * @return  the BioEntrys with that id
     * @throws IllegalIDException if the database doesn't know about the id
     */
    public BioEntryDB getBioEntrys(Set ids, BioEntryDB db) throws BioException,IllegalIDException;
    
    /**
     * Adds a sequence to the database.
     *
     * @param seq the BioEntry to add
     * @throws IllegalIDException if a uniqe ID could not be generated for BioEntry
     * @throws BioException if something goes wrong with adding the BioEntry
     * @throws ChangeVetoException  if either the database does not allow
     *         BioEntrys to be added or the modification was vetoed
     */
    public void addBioEntry(BioEntry seq) throws IllegalIDException, BioException, ChangeVetoException;
    
    /**
     * Remove the BioEntry associated with an ID from the database.
     *
     * @param id  the ID of the BioEntry to remove
     * @throws  IllegalIDException if there is no BioEntry for the ID
     * @throws  BioException if something failed while removing the BioEntry for
     *          that ID
     * @throws  ChangeVetoException  if either the database does not allow
     *          BioEntrys to be removed or the modification was vetoed
     */
    public void removeBioEntry(String id) throws IllegalIDException, BioException, ChangeVetoException;
}
