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

import java.util.Iterator;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.MergeFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeVetoException;

/**
 * An abstract implementation of SequenceDB that provides the sequenceIterator
 * method.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
public abstract class AbstractSequenceDB
  extends
    AbstractChangeable
  implements
    SequenceDB
{
  public SequenceIterator sequenceIterator() {
    return new SequenceIterator() {
      private Iterator pID = ids().iterator();

      public boolean hasNext() {
        return pID.hasNext();
      }

      public Sequence nextSequence() throws BioException {
        return getSequence((String) pID.next());
      }
    };
  }

  public FeatureHolder filter(FeatureFilter ff) {
      MergeFeatureHolder results = new MergeFeatureHolder();
      try {
          for (SequenceIterator si = sequenceIterator(); si.hasNext(); ) {
              Sequence seq = si.nextSequence();
              FeatureHolder fh = seq.filter(ff);
              if (fh != FeatureHolder.EMPTY_FEATURE_HOLDER) {
                  results.addFeatureHolder(fh);
              }
          }
      } catch (BioException ex) {
          throw new BioRuntimeException(ex);
      } catch (ChangeVetoException cve) {
          throw new BioError("Assertion failed: couldn't modify newly created MergeFeatureHolder",cve);
      }
      return results;
  }

  public void addSequence(Sequence seq)
  throws BioException, ChangeVetoException {
    throw new ChangeVetoException("AbstractSequenceDB is immutable");
  }

  public void removeSequence(String id)
  throws BioException, ChangeVetoException {
    throw new ChangeVetoException("AbstractSequenceDB is immutable");
  }
}
