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

package org.biojava.bio.seq.db;

import java.util.Set;

import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.SequenceIterator;

/**
 * A database of sequences with accessible keys and iterators over all
 * sequences.
 * <p>
 * This may have several implementations with rich behaviour, but basically most
 * of the time you will just use the interface methods to do stuff. A sequence
 * database contains a finite number of sequences stored under unique keys.
 *
 * @author Matthew Pocock
 * @author <A href="mailto:Gerald.Loeffler@vienna.at">Gerald Loeffler</A>
 * @author Thomas Down
 */
public interface SequenceDB extends SequenceDBLite {
  /**
   * Get an immutable set of all of the IDs in the database. The ids are legal
   * arguments to getSequence.
   *
   * @return  a Set of ids - at the moment, strings
   */
  Set ids();
  
  /**
   * Returns a SequenceIterator over all sequences in the database. The order
   * of retrieval is undefined.
   *
   * @return  a SequenceIterator over all sequences
   */
  SequenceIterator sequenceIterator();

  /**
   * Query features attached to all sequences in this database.
   * This is equivalent to applying <code>filter</code> to all
   * sequences then merging the results.
   *
   * @param filter a <code>FeatureFilter</code>.
   * @since 1.3
   */

   public FeatureHolder filter(FeatureFilter filter);
}
