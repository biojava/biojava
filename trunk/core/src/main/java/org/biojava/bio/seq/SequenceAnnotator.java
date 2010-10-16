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

package org.biojava.bio.seq;

import org.biojava.bio.BioException;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.utils.ChangeVetoException;

/**
 * An object which adds some additional information to a Sequence.
 *
 * <p>
 * There are two approaches which can be taken to adding features
 * to a sequence:
 *
 * <ol><li>Directly adding features to a Sequence which implements
 * MutableFeatureHolder</li>
 * <li>Creating a new Sequence object which acts as a view on an
 * underlying Sequence, and presents extra features.</li></ol>
 *
 * At present, this interface supports both these mechanisms.  It
 * is the responsibility of the implementor to document which approach
 * is taken.
 * </p>
 *
 * @author Thomas Down
 */

public interface SequenceAnnotator {
    /**
     * Return an annotated version of a sequence.
     *
     * @param seq The sequence to be annotated.
     * @return An annotated version of <code>seq</code> (may be the
     *          same object).
     * @throws IllegalAlphabetException If the sequence is over
     *                                  an inappropriate alphabet for
     *                                  the annotated method being
     *                                  encapsulated
     * @throws BioException if the sequence fails to be annotated
     * @throws ChangeVetoException if either the sequence doesn't allow
     *         annotation or if the change was vetoed
     */

    public Sequence annotate(Sequence seq) throws BioException,
                                              IllegalAlphabetException,
                                              ChangeVetoException;
    
}
