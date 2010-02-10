package org.biojava.bio.dp;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Test emission state events.
 *
 * @author Matthew Pocock
 */
public class EmissionStateEventTest
extends TestCase {
  public void testDistribution_setWeight() {
    try {
      Distribution dist = DistributionFactory.DEFAULT.createDistribution(
              DNATools.getDNA());
      EmissionState es = new SimpleEmissionState("test",
                                                 Annotation.EMPTY_ANNOTATION,
                                                 new int[] { 1 },
                                                 dist);
      EventCounter everything = new EventCounter("Everything");
      EventCounter distC = new EventCounter("Distribution counter");
      EventCounter advC = new EventCounter("Advance counter");

      es.addChangeListener(everything, ChangeType.UNKNOWN);
      es.addChangeListener(distC, EmissionState.DISTRIBUTION);
      es.addChangeListener(advC, EmissionState.ADVANCE);

      // test setWeight - should cause a DISTRIBUTION event, and no Advance event
      dist.setWeight(DNATools.a(), 0.3);
      dist.setWeight(DNATools.a(), 0.7);
      assertEquals("No distribution events vetoed: " + distC + "\n\t" + everything, distC.getPreCounts(), distC.getPostCounts());
      assertEquals("No advance events vetoed: " + advC + "\n\t" + everything, advC.getPreCounts(), advC.getPostCounts());
      assertEquals("Two distribution events: " + distC + "\n\t" + everything, 2, distC.getPostCounts());
      assertEquals("No advance events: " + advC + "\n\t" + everything, 0, advC.getPostCounts());
    } catch (IllegalAlphabetException iae) {
      throw (AssertionError) new AssertionError("Can not create distribution.").initCause(iae);
    } catch (IllegalSymbolException ise) {
      throw (AssertionError) new AssertionError("Can not set weight.").initCause(ise);
    } catch (ChangeVetoException cve) {
      throw (AssertionError) new AssertionError("Prevented from setting weight.").initCause(cve);
    }
  }

  public void testSetDistribution() {
    try {
      Distribution dist = DistributionFactory.DEFAULT.createDistribution(
              DNATools.getDNA());
      Distribution dist2 = DistributionFactory.DEFAULT.createDistribution(
              DNATools.getDNA());
      EmissionState es = new SimpleEmissionState("test",
                                                 Annotation.EMPTY_ANNOTATION,
                                                 new int[] { 1 },
                                                 dist);
      EventCounter everything = new EventCounter("Everything");
      EventCounter distC = new EventCounter("Distribution counter");
      EventCounter advC = new EventCounter("Advance counter");

      es.addChangeListener(everything, ChangeType.UNKNOWN);
      es.addChangeListener(distC, EmissionState.DISTRIBUTION);
      es.addChangeListener(advC, EmissionState.ADVANCE);

      // test setDistribution - should cause a DISTRIBUTION event, and no Advance event
      es.setDistribution(dist2);
      assertEquals("No distribution events vetoed: " + distC + "\n\t" + everything, distC.getPreCounts(), distC.getPostCounts());
      assertEquals("No advance events vetoed: " + advC + "\n\t" + everything, advC.getPreCounts(), advC.getPostCounts());
      assertEquals("One distribution event: " + distC + "\n\t" + everything, 1, distC.getPostCounts());
      assertEquals("No advance events: " + advC + "\n\t" + everything, 0, advC.getPostCounts());
    } catch (IllegalAlphabetException iae) {
      throw (AssertionError) new AssertionError("Can not create distribution.").initCause(iae);
    } catch (ChangeVetoException cve) {
      throw (AssertionError) new AssertionError("Prevented from setting weight.").initCause(cve);
    }
  }
}
