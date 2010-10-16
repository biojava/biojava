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

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.impl.ViewSequence;

/**
 * SequenceDB implementation that returns new SequenceView instances
 * wrapping the sequences in an underlying database. One appropriate
 * use of this would be to wrap a DB in one of these and then wrap
 * this in an annotating db so that the annotation is added to views,
 * not the underlying sequences.
 * 
 * @author Matthew Pocock
 * @since 1.2
 */

public class ViewingSequenceDB extends SequenceDBWrapper {
  /**
   * Create a new ViewingSequenceDB that views the sequences in parent.
   *
   * @param parent the SequenceDB to view
   */
  public ViewingSequenceDB(SequenceDB parent) {
    super(parent);
  }
  
  public String getName() {
    return getParent().getName();
  }
  
  public Sequence getSequence(String id) throws BioException {
    Sequence seq = getParent().getSequence(id);
    return new ViewSequence(seq);
  }
  
  public Set ids() {
    return getParent().ids();
  }
}
