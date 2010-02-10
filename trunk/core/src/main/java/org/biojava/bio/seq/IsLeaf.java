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

/**
 * Accept features which themselves have no children.
 *
 * @author Matthew Pocock
 * @since 1.3
 */

final class IsLeaf implements OptimizableFilter {
  public boolean accept(Feature f) {
    return f.countFeatures() == 0;
  }
  
  public int hashCode() {
    return 41;
  }
  
  /**
  * All instances are equal (this should really be a singleton, but
  * that doesn't quite fit current </code>FeatureFilter</code>
  * patterns.
  */
  
  public boolean equals(Object o) {
    return (o instanceof IsLeaf);
  }
  
  public boolean isProperSubset(FeatureFilter ff) {
    return (ff instanceof IsLeaf) || (ff instanceof AcceptAllFilter);
  }
  
  public boolean isDisjoint(FeatureFilter ff) {
    return ff instanceof ByChild || ff instanceof ByDescendant || ff instanceof AcceptNoneFilter;
  }
}

