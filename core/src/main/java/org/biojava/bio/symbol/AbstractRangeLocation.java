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
package org.biojava.bio.symbol;

import java.util.Collections;
import java.util.Iterator;

/**
 * Base class for simple contiguous Location implementations.
 *
 * Just implement <code>getMin</code> and <code>getMax</code>, and <code>translate</code>..
 *
 * @author Matthew Pocock
 * @author Keith James
 */
public abstract class AbstractRangeLocation extends AbstractLocation {
  public Iterator blockIterator() {
    return Collections.singleton(this).iterator();
  }
  
  public boolean isContiguous() {
    return true;
  }
  
  public SymbolList symbols(SymbolList seq) {
    return seq.subList(getMin(), getMax());
  }
  
  public boolean contains(int p) {
    return p >= getMin() && p <= getMax();
  }
}
