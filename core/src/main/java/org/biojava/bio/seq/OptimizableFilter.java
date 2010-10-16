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
 * The interface for filters that can potentialy optimize themselves, and
 * compare themselves with other filters.
 *
 * @author Matthew Pocock
 * @since 1.2
 */
public interface OptimizableFilter extends FeatureFilter {
  /**
   * Returns true if this filter is a proper subset of sup - that is, for every
   * feature that matches this, it also matches sup. The empty filter is a
   * proper subset of all filters. All filters are a proper subset of the all
   * filter. All filters are proper subsets of themselves.
   *
   * @param sup the potential super set
   * @return true if sup contains all features contained by this filter
   */
  public boolean isProperSubset(FeatureFilter sup);
  
  /**
   * Returns true if this filter is disjoint from filt - that is, there is no
   * Feature that is accepted by both filters. The empty filter is disjoint from
   * all other filters. The all filter is disjoint from none.
   */
    public boolean isDisjoint(FeatureFilter filt);
}
