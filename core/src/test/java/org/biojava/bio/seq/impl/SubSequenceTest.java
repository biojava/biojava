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

package org.biojava.bio.seq.impl;

import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.RemoteFeature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;

/**
 * Tests for SimpleAssembly.  By dependancy, this also
 * tests ProjectedFeatureHolder and SimpleAssembly.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.3
 */

public class SubSequenceTest extends TestCase {
  protected Sequence seq;
  protected Sequence subseq;

  public SubSequenceTest(String name) {
    super(name);
  }

  protected void setUp() throws Exception {
    seq = new SimpleSequence(DNATools.createDNA("aacgtaggttccatgc"),
                             "fragment1",
                             "fragment1",
                             Annotation.EMPTY_ANNOTATION);

    Feature.Template sft = new Feature.Template();
    sft.source = "test1";
    sft.annotation = Annotation.EMPTY_ANNOTATION;
    sft.location = new RangeLocation(1, 3);
    seq.createFeature(sft);

    sft.type = "test2";
    sft.location = new RangeLocation(10, 12);
    seq.createFeature(sft);

    sft.type = "test3";
    sft.location = new RangeLocation(5, 13);
    Feature choppedFeature = seq.createFeature(sft);

    sft.type = "test4";
    sft.location = new RangeLocation(5, 6);
    choppedFeature.createFeature(sft);

    sft.type = "test5";
    sft.location = new RangeLocation(9, 10);
    choppedFeature.createFeature(sft);

    subseq = new SubSequence(seq, 8, 14);
  }

  public void testSymbols()
          throws Exception {
    assertTrue(compareSymbolList(subseq,
                                 DNATools.createDNA("gttccat")));
  }

  public void testFeatureClipping()
          throws Exception {
    assertEquals("subseq.countFeatures()", 2, subseq.countFeatures());
  }

  public void testFeatureProjection()
          throws Exception {
    Feature f = (Feature) subseq.filter(new FeatureFilter.Not(new FeatureFilter.ByClass(RemoteFeature.class)), false).features().next();
    Location fl = f.getLocation();
    assertEquals("projection.getMin()", 3, fl.getMin());
    assertEquals("projection.getMax()", 5, fl.getMax());
  }

  public void testRemoteFeature()
          throws Exception {
    FeatureHolder remotes = subseq.filter(new FeatureFilter.ByClass(RemoteFeature.class), false);
    assertEquals("One remote feature: ", 1, remotes.countFeatures());

    Feature f = (Feature) remotes.features().next();
    assertTrue("Feature implements RemoteFeature", f instanceof RemoteFeature);

    RemoteFeature rf = (RemoteFeature) f;
    Location fl = rf.getLocation();
    assertEquals("remote feature getMin()", 1, fl.getMin());
    assertEquals("remote feature getMax()", 6, fl.getMax());
    assertEquals("remote feature name", seq.getName(), rf.getRemoteFeature().getSequence().getName());
  }

  public void testRemoteChildFeature()
          throws Exception {
    FeatureHolder remotes = subseq.filter(new FeatureFilter.ByClass(RemoteFeature.class), false);
    assertEquals("One remote feature: ", 1, remotes.countFeatures());

    Feature f = (Feature) remotes.features().next();
    assertTrue("Feature implements RemoteFeature", f instanceof RemoteFeature);

    RemoteFeature rf = (RemoteFeature) f;
    assertEquals("child remote features should be pruned - counting them", 1, rf.countFeatures());

    Feature cf = (Feature) rf.features().next();
    Location cfl = cf.getLocation();
    assertEquals("cf.getMin()", 2, cfl.getMin());
    assertEquals("cf.getMax()", 3, cfl.getMax());
  }

  public void testCreateOnSubsequence()
          throws Exception {
    Feature.Template templ = new Feature.Template();
    templ.type = "create_on_subsequence";
    templ.source = "test";
    templ.location = new RangeLocation(2, 3);
    templ.annotation = Annotation.EMPTY_ANNOTATION;
    synchronized (subseq) {
    	subseq.createFeature(templ);	
	}
    

    Feature f = (Feature) seq.filter(new FeatureFilter.ByType("create_on_subsequence"), false).features().next();
    Location fl = f.getLocation();
    assertEquals("fl.getMin()", 9, fl.getMin());
    assertEquals("fl.getMax()", 10, fl.getMax());
  }

  public void testCreateOnSubsequenceFeature()
          throws Exception {
    Feature.Template templ = new Feature.Template();
    templ.type = "create_on_subsequence_feature";
    templ.source = "test";
    templ.location = new RangeLocation(3, 4);
    templ.annotation = Annotation.EMPTY_ANNOTATION;

    Feature subf = (Feature) subseq.filter(new FeatureFilter.Not(new FeatureFilter.ByClass(RemoteFeature.class)), false).features().next();
    synchronized(subf){
    	subf.createFeature(templ);
    }

    Feature f = (Feature) seq.filter(new FeatureFilter.ByType("create_on_subsequence_feature"), true).features().next();
    Location fl = f.getLocation();
    assertEquals(fl.getMin(), 10);
    assertEquals(fl.getMax(), 11);
  }

  public void testRemoveFeatureFromSubsequence()
          throws Exception {
    FeatureHolder fh = subseq.filter(new FeatureFilter.Not(new FeatureFilter.ByClass(RemoteFeature.class)), false);
    assertEquals(fh.countFeatures(), 1);
    Feature f = (Feature) fh.features().next();
    subseq.removeFeature(f);

    fh = subseq.filter(new FeatureFilter.Not(new FeatureFilter.ByClass(RemoteFeature.class)), false);
    assertEquals(fh.countFeatures(), 0);
  }

  private boolean compareSymbolList(SymbolList sl1, SymbolList sl2) {
    if (sl1.length() != sl2.length()) {
      return false;
    }

    Iterator si1 = sl1.iterator();
    Iterator si2 = sl2.iterator();
    while (si1.hasNext()) {
      if (!(si1.next() == si2.next())) {
        return false;
      }
    }

    return true;
  }
}
