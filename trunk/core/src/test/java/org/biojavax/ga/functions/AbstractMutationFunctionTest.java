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
public class AbstractMutationFunctionTest extends TestCase{

  private AbstractMutationFunction func = null;

  public AbstractMutationFunctionTest( String name ){
    super(name);
  }

  protected void setUp() throws Exception {
    super.setUp();

    func = new SimpleMutationFunction();
  }

  protected void tearDown() throws Exception {
    func = null;
    super.tearDown();
  }

  public void TestAbstractCrossOverFunction(){
    assertTrue(
        func.getMutationProbs() == MutationFunction.DEFAULT_MUTATION_PROBS);
  }

  public void testSetAndGetMutationProbs() {
    try {
      double[] probs = new double[]{0.1, 0.9, 0.5};
      func.setMutationProbs(probs);
      assertEquals( func.getMutationProbs(), probs );
    }
    catch (ChangeVetoException ex) {
      fail("Some how the Function has been locked, message is: "+ex.getMessage());
    }
  }
}
