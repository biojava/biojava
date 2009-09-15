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

import java.util.NoSuchElementException;

import org.biojava.bio.BioException;

/**
 * An iterator over a bag of sequences.
 * <p>
 * java.util.Iterator was not appropriate here, as we need specific exceptions
 * to be thrown, and as much type-safety as possible. However, we have made it
 * as compliant with Iterator as we could so that there is a minimal learning
 * curve.
 * @see org.biojavax.bio.seq.RichSequenceIterator
 * @author Matthew Pocock
 */
public interface SequenceIterator {
  /**
   * Returns whether there are more sequences to iterate over.
   *
   * @return  true if there are more sequences to get and false otherwise
   */
  boolean hasNext();
  
  /**
   * Returns the next sequence in the iterator.
   *
   * @return the next Sequence
   * @throws NoSuchElementException if you call nextSequence when hasNext
   *                                returns false
   * @throws BioException if for any reason the sequence could not be retrieved
   */
  Sequence nextSequence() throws NoSuchElementException, BioException;
}
