/*
 *                BioJava development code
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

import java.util.Set;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;

/**
 * JUnit test for FeatureHolderUtils
 * 
 * @author Markus Brosch
 */
public class FeatureHolderUtilsTest extends TestCase {

  protected Feature f1;
  protected Feature f2;
  protected Feature f3;

  protected SimpleFeatureHolder fh1; //has f1, f2
  protected SimpleFeatureHolder fh2; //has f2, f3

  protected void setUp() throws Exception {
    Sequence seq = DNATools.createDNASequence("attagagg", "seq");

    f1 = generateFeatures(seq, "f1");
    f2 = generateFeatures(seq, "f2");
    f3 = generateFeatures(seq, "f3");

    fh1 = new SimpleFeatureHolder();
    fh1.addFeature(f1);
    fh1.addFeature(f2); //also in fh2

    fh2 = new SimpleFeatureHolder();
    fh2.addFeature(f2); //also in fh1
    fh2.addFeature(f3);
  }

  public void testUnion() throws ChangeVetoException {
    FeatureHolder fh = FeatureHolderUtils.union(fh1, fh2);
    assertTrue(fh.containsFeature(f1));
    assertTrue(fh.containsFeature(f2));
    assertTrue(fh.containsFeature(f3));
    assertEquals(3, fh.countFeatures());
  }

  public void testIntersect() throws ChangeVetoException {
    FeatureHolder fh = FeatureHolderUtils.intersect(fh1, fh2);
    assertTrue(fh.containsFeature(f2));
    assertFalse(fh.containsFeature(f1));
    assertFalse(fh.containsFeature(f3));
    assertEquals(1, fh.countFeatures());
  }

  public void testNot() throws ChangeVetoException {
    FeatureHolder fh = FeatureHolderUtils.not(fh1, fh2);
    assertTrue(fh.containsFeature(f1));
    assertFalse(fh.containsFeature(f2));
    assertFalse(fh.containsFeature(f3));
    assertEquals(1, fh.countFeatures());
  }

  public void testFeatureHolderAsSet() {
    Set set = FeatureHolderUtils.featureHolderAsSet(fh1);
    assertTrue(set.contains(f1));
    assertTrue(set.contains(f2));
    assertFalse(set.contains(f3));

    set.add(f3);
    assertTrue(set.contains(f3));
    set.remove(f3);
    assertFalse(set.contains(f3));
  }

  private Feature generateFeatures(FeatureHolder fh, String type) throws BioException,
      ChangeVetoException {
    StrandedFeature.Template templ = new StrandedFeature.Template();
    templ.type = type;
    templ.location = new RangeLocation(1, 5);
    templ.strand = StrandedFeature.POSITIVE;
    templ.source = "foo";
    templ.annotation = Annotation.EMPTY_ANNOTATION;
    return fh.createFeature(templ);
  }
}