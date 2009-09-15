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

package org.biojavax.ga.functions;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.SymbolList;


/**
 * @author Mark Schreiber
 */
public class SimpleCrossOverFunctionTest
    extends TestCase {
  private SimpleCrossOverFunction simpleCrossOverFunction = null;

  public SimpleCrossOverFunctionTest(String name) {
    super(name);
  }

  protected void setUp() throws Exception {
    super.setUp();
    simpleCrossOverFunction = new SimpleCrossOverFunction();
    simpleCrossOverFunction.setMaxCrossOvers(1);
    simpleCrossOverFunction.setCrossOverProbs(new double[]{0.0, 0.0, 1.0, 0.0});
  }

  protected void tearDown() throws Exception {
    simpleCrossOverFunction = null;
    super.tearDown();
  }


  public void testPerformCrossOver() throws Exception {
    SymbolList chromA = DNATools.createDNA("aaaaaaaa");
    SymbolList chromB = DNATools.createDNA("gggggggg");
    GACrossResult cross = simpleCrossOverFunction.performCrossOver(
        chromA, chromB);
    assertTrue(cross.getCrossOverPositions().length == 1);
    assertEquals(cross.getCrossOverPositions()[0], new PointLocation(3));

    assertEquals(cross.getChromosomes()[0] , DNATools.createDNA("aagggggg"));
    assertEquals(chromA , DNATools.createDNA("aagggggg"));
    assertEquals(cross.getChromosomes()[1] , DNATools.createDNA("ggaaaaaa"));
    assertEquals(chromB , DNATools.createDNA("ggaaaaaa"));
  }
}
