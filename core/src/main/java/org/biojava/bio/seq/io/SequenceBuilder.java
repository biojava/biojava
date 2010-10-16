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

package org.biojava.bio.seq.io;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;

/**
 * Interface for objects which accumulate state via SeqIOListener,
 * then construct a Sequence object.
 *
 * <p>
 * It is possible to build `transducer' objects which implement this
 * interface and pass on filtered notifications to a second, underlying
 * SequenceBuilder.  In this case, they should provide a
 * <code>makeSequence</code> method which delegates to the underlying
 * SequenceBuilder.
 * </p>
 *
 * <b>Note:</b> These are one-shot objects that can be used just once to make
 * one sequence. After that, they should be discarded. The usual way to get a
 * supply of these is via a SequenceBuilderFactory.
 *
 * More functionality is offered by the 
 * {@link org.biojavax.bio.seq.io.RichSequenceBuilder RichSequenceBuilder}.
 * @author Thomas Down
 * @since 1.1 [newio proposal]
 * @see org.biojavax.bio.seq.io.RichSequenceBuilder
 */

public interface SequenceBuilder extends SeqIOListener {
    /**
     * Return the Sequence object which has been constructed
     * by this builder.  This method is only expected to succeed
     * after the endSequence() notifier has been called.
     */

    public Sequence makeSequence()
            throws BioException;
}
