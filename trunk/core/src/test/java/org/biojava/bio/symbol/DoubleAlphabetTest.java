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

package org.biojava.bio.symbol;

import junit.framework.TestCase;

/**
 * Tests IntegerAlphabet to make sure we can get it, that it is cannonical, that
 * symbols it returns are cannonical, that finite ranges are cannonical and that
 * symbols in finite ranges are cannonical.
 *
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a>
 * @since 1.2
 */
public class DoubleAlphabetTest extends TestCase
{
  public DoubleAlphabetTest(String name) {
    super(name);
  }

  public void testCanonicalSymbols() {
    Symbol[] syms = new Symbol[10];
    for(int i = 0; i < syms.length; i++) {
      syms[i] = DoubleAlphabet.getInstance().getSymbol(0.5 + i);
    }

    for(int i = 0; i < syms.length; i++) {
      assertEquals(
        syms[i],
        DoubleAlphabet.getInstance().getSymbol(0.5 + i)
      );
    }
  }

  public void testCanonicalSubAlphabets(){
    assertTrue(DoubleAlphabet.getSubAlphabet(1.0,2.5) ==
               DoubleAlphabet.getSubAlphabet(1.0,2.5));

    assertTrue(DoubleAlphabet.getSubAlphabet(1.0,2.5) !=
               DoubleAlphabet.getSubAlphabet(1.00001,2.5));
  }

  public void testSubAlphabetSymbolsCanonical() throws Exception{
    DoubleAlphabet d = DoubleAlphabet.getInstance();
    DoubleAlphabet.SubDoubleAlphabet sd = DoubleAlphabet.getSubAlphabet(2.0,99.0);

    assertTrue(d.getSymbol(3.0) == sd.getSymbol(3.0));
  }
}
