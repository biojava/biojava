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
package org.biojavax.ga.util;

import junit.framework.TestCase;

import org.biojava.bio.dist.OrderNDistribution;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;


/**
 * @author Mark Schreiber
 */
public class GAToolsTest
    extends TestCase {
  Symbol a = DNATools.a();
  Symbol c = DNATools.c();
  Symbol g = DNATools.g();
  Symbol t = DNATools.t();

  public GAToolsTest(String name) {
    super(name);
  }

  protected void setUp() throws Exception {
    super.setUp();
    /**don't really need to do this but hey, why not*/
   // gATools = new GATools();
  }

  protected void tearDown() throws Exception {
   // gATools = null;
    super.tearDown();
  }

  public void testStandardMutationDistribution() throws Exception{
    OrderNDistribution d =
        GATools.standardMutationDistribution(DNATools.getDNA());

    assertTrue(d.getDistribution(a).getWeight(a) == 0.0);
    assertTrue(d.getDistribution(a).getWeight(c) == 1.0/3.0);
    assertTrue(d.getDistribution(a).getWeight(g) == 1.0/3.0);
    assertTrue(d.getDistribution(a).getWeight(t) == 1.0/3.0);

    assertTrue(d.getDistribution(c).getWeight(c) == 0.0);
    assertTrue(d.getDistribution(c).getWeight(a) == 1.0 / 3.0);
    assertTrue(d.getDistribution(c).getWeight(g) == 1.0 / 3.0);
    assertTrue(d.getDistribution(c).getWeight(t) == 1.0/3.0);

    assertTrue(d.getDistribution(g).getWeight(g) == 0.0);
    assertTrue(d.getDistribution(g).getWeight(a) == 1.0 / 3.0);
    assertTrue(d.getDistribution(g).getWeight(c) == 1.0 / 3.0);
    assertTrue(d.getDistribution(g).getWeight(t) == 1.0/3.0);

    assertTrue(d.getDistribution(t).getWeight(t) == 0.0);
    assertTrue(d.getDistribution(t).getWeight(a) == 1.0 / 3.0);
    assertTrue(d.getDistribution(t).getWeight(g) == 1.0 / 3.0);
    assertTrue(d.getDistribution(t).getWeight(c) == 1.0/3.0);

  }

  public void testUniformMutationDistribution() throws Exception {
    OrderNDistribution d =
        GATools.uniformMutationDistribution(DNATools.getDNA());

    assertTrue(d.getDistribution(a).getWeight(a) == 0.25);
    assertTrue(d.getDistribution(a).getWeight(c) == 0.25);
    assertTrue(d.getDistribution(a).getWeight(g) == 0.25);
    assertTrue(d.getDistribution(a).getWeight(t) == 0.25);

    assertTrue(d.getDistribution(c).getWeight(a) == 0.25);
    assertTrue(d.getDistribution(c).getWeight(c) == 0.25);
    assertTrue(d.getDistribution(c).getWeight(g) == 0.25);
    assertTrue(d.getDistribution(c).getWeight(t) == 0.25);

    assertTrue(d.getDistribution(g).getWeight(a) == 0.25);
    assertTrue(d.getDistribution(g).getWeight(c) == 0.25);
    assertTrue(d.getDistribution(g).getWeight(g) == 0.25);
    assertTrue(d.getDistribution(g).getWeight(t) == 0.25);

    assertTrue(d.getDistribution(t).getWeight(a) == 0.25);
    assertTrue(d.getDistribution(t).getWeight(c) == 0.25);
    assertTrue(d.getDistribution(t).getWeight(g) == 0.25);
    assertTrue(d.getDistribution(t).getWeight(t) == 0.25);

  }

  public void testBinaryAlphabet(){
    FiniteAlphabet bin = GATools.getBinaryAlphabet();
    assertEquals(bin.getName(), "GA_Binary");
    FiniteAlphabet bin2 =
        (FiniteAlphabet)AlphabetManager.alphabetForName("GA_Binary");

    assertNotNull(bin);
    assertNotNull(bin2);
    assertTrue(bin == bin2);
    assertTrue(bin2 == bin);
    assertEquals(bin, bin2);
    assertEquals(bin2, bin);

    //System.out.println(bin.size());
    assertTrue(bin.size() == 2);
    assertTrue(bin.contains(GATools.one()));
    assertTrue(bin.contains(GATools.one()));

    Symbol one = GATools.one();
    Symbol x = GATools.one();
    assertTrue(one == x);
    assertTrue(x == one);
    assertEquals(one, x);
    assertEquals(x, one);

    Symbol zero = GATools.zero();
              x = GATools.zero();
    assertTrue(zero == x);
    assertTrue(x == zero);
    assertEquals(zero, x);
    assertEquals(x, zero);
  }

  public void testCreateBinary(){
    String bin = "00101000101010101010100101010010101010";
    try {
      SymbolList symL = GATools.createBinary(bin);
      assertNotNull(symL);
      assertTrue(bin.length() == 38);

      String out = symL.seqString();
      assertNotNull(out);
      assertEquals(bin, out);
    }
    catch (IllegalSymbolException ex) {
      fail(ex.getMessage());
    }

  }
}
