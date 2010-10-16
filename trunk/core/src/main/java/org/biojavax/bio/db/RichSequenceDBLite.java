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
import org.biojava.bio.seq.db.SequenceDBLite;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.seq.RichSequence;

/**
 * A database of RichSequences. This may have several implementations with
 * rich behaviour, but basically most of the time you will just use
 * the interface methods to do stuff. A RichSequence database contains a
 * finite number of RichSequences stored under unique keys.
 *
 * @author Matthew Pocock
 * @author Gerald Loeffler
 * @author Thomas Down
 * @author Richard Holland
 * @since 1.5
 */
public interface RichSequenceDBLite extends BioEntryDBLite,SequenceDBLite {      
    /**
     * Retrieve a single RichSequence by its id.
     *
     * @param id the id to retrieve by
     * @return  the Sequence with that id
     * @throws IllegalIDException if the database doesn't know about the id
     */
    public RichSequence getRichSequence(String id) throws BioException,IllegalIDException;
    
    /**
     * Retrieve multiple RichSequence by its id.
     *
     * @param ids a set of ids to retrieve by
     * @return  the RichSequences with that id
     * @throws IllegalIDException if the database doesn't know about the id
     */
    public RichSequenceDB getRichSequences(Set ids) throws BioException,IllegalIDException;
    
    /**
     * Retrieve multiple RichSequence into a specific sequence database. If
     * that database is null, a new HashRichSequenceDB is used.
     *
     * @param ids a set of ids to retrieve by
     * @param db a database to load the seqs into
     * @return  the RichSequences with that id
     * @throws IllegalIDException if the database doesn't know about the id
     */
    public RichSequenceDB getRichSequences(Set ids, RichSequenceDB db) throws BioException,IllegalIDException;
    
    /**
     * Adds a sequence to the database.
     *
     * @param seq the RichSequence to add
     * @throws IllegalIDException if a uniqe ID could not be generated for RichSequence
     * @throws BioException if something goes wrong with adding the RichSequence
     * @throws ChangeVetoException  if either the database does not allow
     *         RichSequences to be added or the modification was vetoed
     */
    public void addRichSequence(RichSequence seq) throws IllegalIDException, BioException, ChangeVetoException;
    
    /**
     * Remove the RichSequence associated with an ID from the database.
     *
     * @param id  the ID of the RichSequence to remove
     * @throws  IllegalIDException if there is no RichSequence for the ID
     * @throws  BioException if something failed while removing the RichSequence for
     *          that ID
     * @throws  ChangeVetoException  if either the database does not allow
     *          RichSequences to be removed or the modification was vetoed
     */
    public void removeRichSequence(String id) throws IllegalIDException, BioException, ChangeVetoException;
}
