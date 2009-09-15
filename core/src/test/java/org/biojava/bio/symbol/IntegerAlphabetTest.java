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
 * Tests IntegerAlphabet to make sure we can get it, that it is
 * canonical, that symbols it returns are canonical, that finite
 * ranges are canonical and that symbols in finite ranges are
 * canonical.
 *
 * @since 1.2
 */
public class IntegerAlphabetTest extends TestCase
{
  public IntegerAlphabetTest(String name) {
    super(name);
  }

  public void testCanonicalAlphabet() {
    assertEquals(
      AlphabetManager.alphabetForName("INTEGER"),
      IntegerAlphabet.getInstance()
    );
    
    assertEquals(
      IntegerAlphabet.getInstance(),
      IntegerAlphabet.getInstance()
    );
  }

  public void testCanonicalSymbols() {
    Symbol[] syms = new Symbol[10];
    for(int i = 0; i < syms.length; i++) {
      syms[i] = IntegerAlphabet.getInstance().getSymbol(i);
    }

    for(int i = 0; i < syms.length; i++) {
      assertEquals(
        syms[i],
        IntegerAlphabet.getInstance().getSymbol(i)
      );
    }
  }

  public void testCanonicalSubAlphabet() {
    int min = 200;
    int max = 300;
 
    Alphabet a1 = IntegerAlphabet.getSubAlphabet(min, max);
    Alphabet a2 = IntegerAlphabet.getSubAlphabet(min, max);
 
    assertEquals(a1, a2);
  }

  public void testSubAlphabetSymbolsCanonical()
  throws IllegalSymbolException {
    int min = 400;
    int max = 500;

    IntegerAlphabet alpha = IntegerAlphabet.getInstance();

    IntegerAlphabet.SubIntegerAlphabet a1 = IntegerAlphabet.getSubAlphabet(min, max);
    for(int i = min; i <= max; i++) {
      assertEquals(
        a1.getSymbol(i),
        a1.getSymbol(i)
      );

      assertEquals(
        a1.getSymbol(i),
        alpha.getSymbol(i)
      );
    }
  }
}
