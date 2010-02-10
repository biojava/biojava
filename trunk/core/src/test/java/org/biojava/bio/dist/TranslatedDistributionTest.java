package org.biojava.bio.dist;

import junit.framework.AssertionFailedError;
import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * Tests that check whether translated distributions behave properly.
 *
 * @author Thomas Down
 * @since 1.4
 */
public class TranslatedDistributionTest extends TestCase {
  private double a = 0.1;
  private double g = 0.3;
  private double c = 0.25;
  private double t = 0.35;
  private double delta = 0.000001;
    
  public TranslatedDistributionTest(String name) {
    super(name);
  }
  
  public void testWeights() 
    throws Exception
  {
    try {
        Distribution d1 = DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
        d1.setWeight(DNATools.a(), a);
        d1.setWeight(DNATools.c(), c);
        d1.setWeight(DNATools.g(), g);
        d1.setWeight(DNATools.t(), t);
        Distribution d2 = new TranslatedDistribution(
            DNATools.complementTable(),
            d1,
            DistributionFactory.DEFAULT
        );
        assertEquals(a, d2.getWeight(DNATools.t()), delta);
        assertEquals(c, d2.getWeight(DNATools.g()), delta);
        assertEquals(g, d2.getWeight(DNATools.c()), delta);
        assertEquals(t, d2.getWeight(DNATools.a()), delta);
        
        d1.setWeight(DNATools.a(), t);
        d1.setWeight(DNATools.c(), g);
        d1.setWeight(DNATools.g(), c);
        d1.setWeight(DNATools.t(), a);

        assertEquals(a, d2.getWeight(DNATools.a()), delta);
        assertEquals(c, d2.getWeight(DNATools.c()), delta);
        assertEquals(g, d2.getWeight(DNATools.g()), delta);
        assertEquals(t, d2.getWeight(DNATools.t()), delta);
    } catch (IllegalSymbolException ise) {
      throw new AssertionFailedError("Can't retrieve weight: "
      + ise.getMessage());
    } 
  }
}
