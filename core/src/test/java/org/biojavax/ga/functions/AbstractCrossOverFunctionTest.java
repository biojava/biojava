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

import org.biojava.utils.ChangeVetoException;


/**
 * @author Mark Schreiber
 */
public class AbstractCrossOverFunctionTest
    extends TestCase {
  private AbstractCrossOverFunction abstractCrossOverFunction = null;

  public AbstractCrossOverFunctionTest(String name) {
    super(name);
  }

  protected void setUp() throws Exception {
    super.setUp();

    abstractCrossOverFunction = new SimpleCrossOverFunction();
  }

  protected void tearDown() throws Exception {
    abstractCrossOverFunction = null;
    super.tearDown();
  }

  public void testAbstractCrossOverFunction() {
    abstractCrossOverFunction = new SimpleCrossOverFunction();
    assertTrue(abstractCrossOverFunction.getCrossOverProbs() == CrossOverFunction.DEFAULT_CROSS_PROB);
    assertTrue(abstractCrossOverFunction.getMaxCrossOvers() == CrossOverFunction.DEFAULT_MAX_CROSS);
  }

  public void testGetCrossOverProbs() throws Exception{
    double[] expectedReturn = new double[]{0.1, 0.4, 0.34};
    abstractCrossOverFunction.setCrossOverProbs(expectedReturn);
    double[] actualReturn = abstractCrossOverFunction.getCrossOverProbs();
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testGetMaxCrossOvers() throws Exception{
    int expectedReturn = 12;
    abstractCrossOverFunction.setMaxCrossOvers(expectedReturn);
    int actualReturn = abstractCrossOverFunction.getMaxCrossOvers();
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testSetCrossOverProbs() throws ChangeVetoException {
    double[] crossOverProbs = new double[]{0.5, 0.9};
    abstractCrossOverFunction.setCrossOverProbs(crossOverProbs);
  }

  public void testSetMaxCrossOvers() throws ChangeVetoException {
    int maxCrossOvers = 0;
    abstractCrossOverFunction.setMaxCrossOvers(maxCrossOvers);
  }
}
