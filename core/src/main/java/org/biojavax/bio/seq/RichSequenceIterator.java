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

package org.biojavax.bio.seq;
import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.SequenceIterator;
import org.biojavax.bio.BioEntryIterator;


/**
 * Essentially the same as SequenceIterator. It provides a new
 * method that returns RichSequence objects without the need for
 * explicit casting. Implementations of this interface should <b>always</b>
 * return RichSequence objects for both the nextSequence() and 
 * nextRichSequence() methods.
 *
 * @author Mark Schreiber
 * @author Richard Holland
 * @see org.biojava.bio.seq.SequenceIterator
 * @since 1.5
 */
public interface RichSequenceIterator extends BioEntryIterator,SequenceIterator {
    public RichSequence nextRichSequence() throws NoSuchElementException, BioException;
}
