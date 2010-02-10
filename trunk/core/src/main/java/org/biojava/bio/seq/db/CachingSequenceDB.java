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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.lang.ref.SoftReference;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;

/**
 * SequenceDB implementation that caches the results of another SequenceDB.
 * 
 * @author Matthew Pocock
 */

public class CachingSequenceDB extends SequenceDBWrapper {
  private transient Map cache;
  
  /**
   * Create a new CachingSequenceDB that caches the sequences in parent.
   *
   * @param parent the SequenceDB to cache
   */
  public CachingSequenceDB(SequenceDB parent) {
    super(parent);
    cache = new HashMap();
  }
  
  public String getName() {
    return getParent().getName();
  }
  
  public Sequence getSequence(String id) throws BioException {
    SoftReference ref = (SoftReference) cache.get(id);
    Sequence seq;
    if(ref == null) {
      seq = getParent().getSequence(id);
      cache.put(id, new SoftReference(seq));
    } else {
      seq = (Sequence) ref.get();
      if(seq == null) {
        seq = getParent().getSequence(id);
        cache.put(id, new SoftReference(seq));
      }
    }
    return seq;
  }
  
  public Set ids() {
    return getParent().ids();
  }
  
  private void readObject(ObjectInputStream in)
  throws IOException, ClassNotFoundException {
    in.defaultReadObject();
    this.cache = new HashMap();
  }
}
