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

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;

/**
 * An abstract implementation of SequenceDB that wraps up another database.
 *
 * @author Matthew Pocock
 */
public abstract class SequenceDBWrapper extends AbstractSequenceDB
implements java.io.Serializable{
  private final SequenceDB parent;
  private transient SequencesForwarder seqFor;
  
  protected ChangeSupport getChangeSupport(ChangeType ct) {
    ChangeSupport changeSupport = super.getChangeSupport(ct);
    
    if(ct.isMatchingType(SequenceDB.SEQUENCES)) {
      seqFor = new SequencesForwarder(this, changeSupport);
      parent.addChangeListener(seqFor, SequenceDB.SEQUENCES);
    }
    
    return changeSupport;
  }
  
  /**
   * Return the parent SequenceDB.
   *
   * @return the parent SequenceDB
   */
  public SequenceDB getParent() {
    return this.parent;
  }
  
  protected class SequencesForwarder extends ChangeForwarder {
    public SequencesForwarder(Object source, ChangeSupport cs) {
      super(source, cs);
    }
    
    public ChangeEvent generateEvent(ChangeEvent ce) {
      if(ce.getType() == SequenceDB.SEQUENCES) {
        Object previous = ce.getPrevious();
        if(previous != null) {
          return new ChangeEvent(
            getSource(),
            SequenceDB.SEQUENCES,
            null,
            previous,
            ce
          );
        }
      }
      return null;
    }
  }

  public SequenceDBWrapper(SequenceDB parent) {
    this.parent = parent;
  }
}
