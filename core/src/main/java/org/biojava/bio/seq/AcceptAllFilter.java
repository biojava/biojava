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
 * The class that accepts all features.
 * <p>
 * Use the FeatureFilter.all member.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.0
 */
class AcceptAllFilter implements OptimizableFilter {
  protected AcceptAllFilter() {}
  
  public boolean accept(Feature f) { return true; }
  
  public boolean equals(Object o) {
    return o instanceof AcceptAllFilter;
  }
  
  public int hashCode() {
    return 0;
  }
  
  public boolean isProperSubset(FeatureFilter sup) {
    return sup instanceof AcceptAllFilter;
  }
  
  public boolean isDisjoint(FeatureFilter filt) {
    return filt instanceof AcceptNoneFilter;
  }
  
  public String toString() {
    return "All";
  }
}

