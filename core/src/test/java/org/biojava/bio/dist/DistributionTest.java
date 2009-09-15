package org.biojava.bio.dist;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import junit.framework.AssertionFailedError;
import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;

/**
 * Tests that simple distributions work as advertised.
 *
 * @author Matthew Pocock
 * @since 1.2
 */
public class DistributionTest extends TestCase {
  private double delta = 0.00001;
  private double a = 0.1;
  private double g = 0.3;
  private double c = 0.25;
  private double t = 0.35;

  private Distribution dist;
  private Distribution gap;
  private Distribution gapDist;
  
  public DistributionTest(String name) {
    super(name);
  }
  
  protected void setUp() {
    // create distributions
    try {
      dist = DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
      gap = new GapDistribution(DNATools.getDNA());
      gapDist = new PairDistribution(gap, dist);
    } catch (IllegalAlphabetException iae) {
      throw new AssertionFailedError("Can't initialize test distributions " + iae.getMessage());
    }

    // set the weights in dist
    try {
      dist.setWeight(DNATools.a(), a);
      dist.setWeight(DNATools.g(), g);
      dist.setWeight(DNATools.c(), c);
      dist.setWeight(DNATools.t(), t);
    } catch (IllegalSymbolException ise) {
      throw new AssertionFailedError("Unable to set weights: "
      + ise.getMessage());
    } catch (ChangeVetoException cve) {
      throw new AssertionFailedError("Unable to set weights: "
      + cve.getMessage());
    }
  }
  
  // make sure that we get out the weights out that we put in
  public void testSimpleDistribution() {
    try {
      assertEquals("gap", 0.0, dist.getWeight(dist.getAlphabet().getGapSymbol()), delta);
      assertEquals(DNATools.a().toString(), a, dist.getWeight(DNATools.a()), delta);
      assertEquals(DNATools.g().toString(), g, dist.getWeight(DNATools.g()), delta);
      assertEquals(DNATools.c().toString(), c, dist.getWeight(DNATools.c()), delta);
      assertEquals(DNATools.t().toString(), t, dist.getWeight(DNATools.t()), delta);
    } catch (IllegalSymbolException ise) {
      throw new AssertionFailedError("Can't retrieve weight: "
      + ise.getMessage());
    }
  }
  
  // check that the ambiguity stuff works
  public void testAmbiguities() {
    try {
      // implicit from sets
      Set syms = new HashSet();
      double tot = 0.0;
      assertEquals(syms.toString(), tot, dist.getWeight(dist.getAlphabet().getAmbiguity(syms)), delta);
      
      syms.add(DNATools.a());
      tot += a;
      assertEquals(syms.toString(), tot, dist.getWeight(dist.getAlphabet().getAmbiguity(syms)), delta);
      
      syms.add(DNATools.g());
      tot += g;
      assertEquals(syms.toString(), tot, dist.getWeight(dist.getAlphabet().getAmbiguity(syms)), delta);
      
      syms.add(DNATools.c());
      tot += c;
      assertEquals(syms.toString(), tot, dist.getWeight(dist.getAlphabet().getAmbiguity(syms)), delta);
      
      syms.add(DNATools.t());
      tot += t;
      assertEquals(syms.toString(), tot, dist.getWeight(dist.getAlphabet().getAmbiguity(syms)), delta);
    } catch (IllegalSymbolException ise) {
      throw new AssertionFailedError("Can't retrieve weight: "
      + ise.getMessage());
    }
  }
  
  // check gap distribution
  public void testGapDistribution() {
    try {
      assertEquals("real-gap", 1.0
      , gap.getWeight(AlphabetManager.getGapSymbol()), delta);
      assertEquals("gap", 1.0, gap.getWeight(gap.getAlphabet().getGapSymbol()), delta);
      assertEquals(DNATools.a().toString(), 1.0, gap.getWeight(DNATools.a()), delta);
      assertEquals(DNATools.g().toString(), 1.0, gap.getWeight(DNATools.g()), delta);
      assertEquals(DNATools.c().toString(), 1.0, gap.getWeight(DNATools.c()), delta);
      assertEquals(DNATools.t().toString(), 1.0, gap.getWeight(DNATools.t()), delta);
    } catch (IllegalSymbolException ise) {
      throw new AssertionFailedError("Can't retrieve weight: "
      + ise.getMessage());
    }
  }
  
  public void testGapDist() {
    try {
      Symbol[] syms = new Symbol[2];
      syms[0] = gap.getAlphabet().getGapSymbol();
      List symL = Arrays.asList(syms);
      
      for(Iterator si = ((FiniteAlphabet) dist.getAlphabet()).iterator(); si.hasNext(); ) {
        Symbol s = (Symbol) si.next();
        syms[1] = s;
        Symbol sym = gapDist.getAlphabet().getSymbol(symL);
        
        assertEquals(sym.toString(), gapDist.getWeight(sym), dist.getWeight(s), delta);
      }
      
    } catch (IllegalSymbolException ise) {
      throw new AssertionFailedError("Can't retrieve weight: "
      + ise.getMessage());
    }
  }
}
