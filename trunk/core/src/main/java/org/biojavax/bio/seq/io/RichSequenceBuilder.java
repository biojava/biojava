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

package org.biojavax.bio.seq.io;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.SequenceBuilder;
import org.biojavax.bio.seq.RichSequence;

/**
 * An interface for objects that can build RichSequences.
 * @author Mark Schreiber
 * @since 1.5
 */
public interface RichSequenceBuilder extends RichSeqIOListener, SequenceBuilder {
    
    /**
     * {@inheritDoc}
     * Implementations of this for a RichSequenceBuilder should
     * delegate to makeRichSequence() and return only RichSequence objects.
     */
    public Sequence makeSequence() throws BioException;
    
    /**
     * Build a RichSequence.
     * @return a RichSequence
     * @throws BioException if it is not possible to build a RichSequence
     */
    public RichSequence makeRichSequence() throws BioException;
}
