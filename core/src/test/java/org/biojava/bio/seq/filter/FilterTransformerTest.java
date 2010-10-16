package org.biojava.bio.seq.filter;

import junit.framework.TestCase;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;

/**
 * Checks that FilterTransformer is not totaly nuts.
 *
 * @author Matthew Pocock
 */
public class FilterTransformerTest
        extends TestCase {
  public void testIdentity() {
    try {
      FeatureFilter filt = FilterUtils.and(
              FilterUtils.byType("pigs"),
              FilterUtils.or(
                      FilterUtils.bySource("iran"),
                      FilterUtils.not(
                              FilterUtils.onlyChildren(
                                      FilterUtils.hasAnnotation("beer")
                              ))));

      FeatureFilter filt2 = (FeatureFilter) FilterUtils.visitFilter(
              filt,
              new FilterTransformer());
      assertTrue("Non-moidfying transformer gives equal results:\n\t" + filt + "\n\t" + filt2,
                 FilterUtils.areEqual(filt, filt2));
    } catch (BioException be) {
      throw (AssertionError) new AssertionError("Couldn't make walker").initCause(be);
    }
  }

  public void testReplace() {
    try {
      int trans = 2000;
      Location loc = new RangeLocation(20, 50);
      Location loc2 = loc.translate(trans);
      FeatureFilter filt = FilterUtils.containedByLocation(loc);
      FeatureFilter filt2 = FilterUtils.containedByLocation(loc2);

      FeatureFilter filt3 = (FeatureFilter) FilterUtils.visitFilter(
              filt,
              new Translater(trans));
      assertTrue("Non-moidfying transformer gives equal results:\n\t" + filt2 + "\n\t" + filt3,
                 FilterUtils.areEqual(filt2, filt3));
    } catch (BioException be) {
      throw (AssertionError) new AssertionError("Couldn't make walker").initCause(be);
    }
  }

  public class Translater
  extends FilterTransformer {
    private int trans;

    public Translater(int trans) {
      this.trans = trans;
    }

    public FeatureFilter containedByLocation(FeatureFilter.ContainedByLocation filt) {
      return FilterUtils.containedByLocation(filt.getLocation().translate(trans));
    }
  }
}
