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

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;

/**
 * @author Matthew Pocock
 */
public class SubSequenceDB extends SequenceDBWrapper {
  private final Set ids;
  
  public SubSequenceDB(SequenceDB parent, Set ids)
  throws BioException {
    super(parent);
    this.ids = new HashSet(ids);
    Set pids = parent.ids();
    if(!pids.containsAll(ids)) {
      throw new BioException(
        "IDs must all be contained in the parent database " +
        parent.getName()
      );
    }
  }
  
  public String getName() {
    return getParent().getName() + " subset " + ids.toString();
  }
  
  public Sequence getSequence(String id) throws BioException {
    if(!ids.contains(id)) {
      throw new BioException(
        "No sequence for " + id + " found in database " + getName()
      );
    }
    return getParent().getSequence(id);
  }
  
  public Set ids() {
    return Collections.unmodifiableSet(ids);
  }
}
