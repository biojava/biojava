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

import org.biojavax.bio.BioEntryIterator;

/**.
 * A database of RichSequences with accessible keys and iterators over all
 * sequences.
 * <p>
 * This may have several implementations with rich behaviour, but basically most
 * of the time you will just use the interface methods to do stuff. A sequence
 * database contains a finite number of sequences stored under unique keys.
 *
 * @author Matthew Pocock
 * @author Gerald Loeffler
 * @author Thomas Down
 * @author Richard Holland
 * @since 1.5
 */
public interface BioEntryDB extends BioEntryDBLite {
    /**
     * Get an immutable set of all of the IDs in the database. The ids are legal
     * arguments to getBioEntry.
     *
     * @return  a Set of ids - at the moment, strings
     */
    public Set ids();
    
    /**
     * Returns a BioEntryIterator over all BioEntrys in the database. The order
     * of retrieval is undefined.
     * @return a BioEntryIterator over all BioEntrys
     */
    public BioEntryIterator getBioEntryIterator();
}
