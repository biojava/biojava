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


package org.biojava.bio.dp;

import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.GapDistribution;
import org.biojava.bio.symbol.Alphabet;

/**
 * Start/end state for HMMs.
 * <p>
 * All MagicalState objects emit over MAGICAL_ALPHABET, which only contains
 * MAGICAL_STATE.
 *
 * @author Matthew Pocock
 */
public final class MagicalState extends SimpleEmissionState {
  /**
   * A cache of magical state objects so that we avoid making the same
   * thing twice.
   */
  private static final Map stateCache;
  static {
    stateCache = new HashMap();
  }
  
  public static MagicalState getMagicalState(Alphabet alphabet, int heads) {
    AlphaHeads ah = new AlphaHeads(alphabet, heads);
    MagicalState ms = (MagicalState) stateCache.get(ah);
    if(ms == null) {
      ms = new MagicalState(alphabet, heads);
      stateCache.put(ah, ms);
    }
    return ms;
  }

  private MagicalState(Alphabet alpha, int heads) {
    super(
      "!-" + heads,
      Annotation.EMPTY_ANNOTATION,
      new int[heads],
      new GapDistribution(alpha)
    );
    int [] advance = getAdvance();
    for(int i = 0; i < heads; i++) {
      advance[i] = 1;
    }
  }

  private Object writeReplace() throws ObjectStreamException {
    return new PlaceHolder(getDistribution().getAlphabet(), getAdvance().length);
  }
  
  private static class AlphaHeads {
    public Alphabet alpha;
    public int heads;
    
    public AlphaHeads(Alphabet alpha, int heads) {
      this.alpha = alpha;
      this.heads = heads;
    }
    
    public boolean equals(Object o) {
      AlphaHeads ah = (AlphaHeads) o;
      return this.alpha == ah.alpha && this.heads == ah.heads;
    }
    
    public int hashCode() {
      return alpha.hashCode() ^ heads;
    }
  }
  
  private static class PlaceHolder implements Serializable {
    private Alphabet alpha;
    private int heads;
    
    public PlaceHolder(Alphabet alpha, int heads) {
      this.alpha = alpha;
      this.heads = heads;
    }
    
    public Object readResolve() throws ObjectStreamException {
      return MagicalState.getMagicalState(alpha, heads);
    }
  }
}

