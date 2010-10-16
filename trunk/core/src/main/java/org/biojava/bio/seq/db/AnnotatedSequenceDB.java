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

import java.io.Serializable;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceAnnotator;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.utils.ChangeVetoException;

/**
 * SequenceDB implementation which lazily applies a SequenceAnnotator
 * to sequences retrieved from a SequenceDB.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public class AnnotatedSequenceDB
extends AbstractSequenceDB
implements SequenceDB, Serializable {
  private final SequenceDB parent;
  private final SequenceAnnotator annotator;

  public AnnotatedSequenceDB(SequenceDB parent, SequenceAnnotator a) {
    this.parent = parent;
    this.annotator = a;
  }

   /**
    * Get the original sequenceDB from this annotated sequenceDB.
    */

  public SequenceDB getParent() {
    return this.parent;
  }

  public String getName() {
    return parent.getName() + " (" + annotator.toString() + ")";
  }

  public Sequence getSequence(String id)
  throws BioException {
    return doAnnotation(parent.getSequence(id));
  }

  public Set ids() {
    return parent.ids();
  }

  public SequenceIterator sequenceIterator() {
    return new SequenceIterator() {
      SequenceIterator pi = parent.sequenceIterator();

            public boolean hasNext() {
        return pi.hasNext();
            }

            public Sequence nextSequence() throws BioException {
        return doAnnotation(pi.nextSequence());
            }
    };
  }

   /**
    * Apply the annotation to a sequence.
    * @param seq the sequence to annotate.
    */

  protected Sequence doAnnotation(Sequence seq) throws BioException  {
    try {
      return annotator.annotate(seq);
    } catch (IllegalAlphabetException ex) {
      throw new BioException("Couldn't apply annotator " + annotator.toString() + " to " + seq.getURN(), ex);
    } catch (ChangeVetoException cve) {
      throw new BioException("Couldn't apply annotator " + annotator.toString() + " to " + seq.getURN(), cve);
    }
  }
}
