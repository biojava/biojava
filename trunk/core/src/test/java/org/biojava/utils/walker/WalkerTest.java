package org.biojava.utils.walker;

import junit.framework.TestCase;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.RangeLocation;

/**
 * Test some walkers for some viewers, ensuring that they get the events we
 * would expect.
 *
 * @author Matthew Pocock
 */
public class WalkerTest
extends TestCase {
  private FeatureFilter booring1;
  private FeatureFilter booring2;
  private FeatureFilter booring4;
  private FeatureFilter and;
  private FeatureFilter andOr;

  protected void setUp() {
    booring1 = new FeatureFilter.OverlapsLocation(new RangeLocation(20, 50));
    booring2 = new FeatureFilter.ByClass(StrandedFeature.class);
    booring4 = new FeatureFilter.ByClass(ComponentFeature.class);
    and = new FeatureFilter.And(booring1, booring2);
    andOr = new FeatureFilter.And(booring1,
                                  new FeatureFilter.Or(booring2, booring4));
  }

  public void testCountAll() {
    try {
      CountAll ca = new CountAll();
      Walker walker = WalkerFactory.getInstance().getWalker(ca);

      walker.walk(booring1, ca);
      assertEquals("One filter: " + booring1, 1, ca.count);

      ca.count = 0;
      walker.walk(and, ca);
      assertEquals("Three filters: " + and, 3, ca.count);

      ca.count = 0;
      walker.walk(andOr, ca);
      assertEquals("Five filters: " + andOr, 5, ca.count);
    } catch (BioException be) {
      throw (AssertionError) new AssertionError(
              "Could not instantiate visitor").initCause(be);
    }
  }

  public void testCountSome() {
    try {
      CountByClass cbc = new CountByClass();
      Walker walker = WalkerFactory.getInstance().getWalker(cbc);

      walker.walk(booring1, cbc);
      assertEquals("One filter, none interesting", 0, cbc.count);

      cbc.count = 0;
      walker.walk(booring2, cbc);
      assertEquals("One filter, one interesting", 1, cbc.count);

      cbc.count = 0;
      walker.walk(andOr, cbc);
      assertEquals("Five filter, two interesting", 2, cbc.count);
    } catch (BioException be) {
      throw (AssertionError) new AssertionError(
              "Could not instantiate visitor").initCause(be);
    }
  }

  public void testCountWithFallback() {
    try {
      CountWithFallback cwf = new CountWithFallback();
      Walker walker = WalkerFactory.getInstance().getWalker(cwf);

      walker.walk(booring1, cwf);
      assertEquals("One filter, one booring", 0, cwf.others);
      assertEquals("One filter, zero interesting", 1, cwf.byLoc);

      cwf.others = cwf.byLoc = 0;
      walker.walk(booring2, cwf);
      assertEquals("One filter, one booring", 1, cwf.others);
      assertEquals("One filter, zero interesting", 0, cwf.byLoc);

      cwf.others = cwf.byLoc = 0;
      walker.walk(andOr, cwf);
      assertEquals("One filter, four booring", 4, cwf.others);
      assertEquals("One filter, one interesting", 1, cwf.byLoc);
    } catch (BioException be) {
      throw (AssertionError) new AssertionError(
              "Could not instantiate visitor").initCause(be);
    }
  }

  public void testReturnAll() {
    try {
      ReturnOne ra = new ReturnOne();
      Walker walker = WalkerFactory.getInstance().getWalker(ra);

      walker.walk(booring1, ra);
      assertEquals("One filter", new Integer(1), walker.getValue());

      walker.walk(andOr, ra);
      assertEquals("Five filters", new Integer(1), walker.getValue());
    } catch (BioException be) {
      throw (AssertionError) new AssertionError(
              "Could not instantiate visitor").initCause(be);
    }
  }

  public class CountAll implements Visitor {
    int count = 0;

    public void featureFilter(FeatureFilter filter) {
      System.err.println("Increasing counter: " + filter);
      count++;
    }
  }

  public class CountByClass
  implements Visitor {
    int count = 0;

    public void byClass(FeatureFilter.ByClass byClass) {
      count++;
    }
  }

  public class CountWithFallback
  implements Visitor {
    int byLoc = 0;
    int others = 0;

    public void featureFilter(FeatureFilter filter) {
      System.err.println("Feature: " + filter);
      others++;
    }

    public void overlapsLocation(FeatureFilter.OverlapsLocation overlaps) {
      System.err.println("OverlapsLocation: " + overlaps);
      byLoc++;
    }
  }

  public class ReturnOne
  implements Visitor {
    public Integer featureFilter(FeatureFilter filter) {
      return new Integer(1);
    }
  }
}
