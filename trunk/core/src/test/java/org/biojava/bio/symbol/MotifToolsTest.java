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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.AssertionFailure;

public class MotifToolsTest
    extends TestCase {

    protected String n;

    protected void setUp() {
      try {
          StringBuffer sb = new StringBuffer();
          sb.append("[");
          SymbolTokenization sTok = DNATools.getDNA().getTokenization("token");
          FiniteAlphabet na = (FiniteAlphabet) DNATools.n().getMatches();

          Set rawSyms = AlphabetManager.getAllSymbols(na);
          List gapSyms = new ArrayList();

          for (Iterator si = rawSyms.iterator(); si.hasNext();) {
              Symbol rawSym = (Symbol) si.next();
              // Crude check for gap symbol
              if (((FiniteAlphabet) rawSym.getMatches()).size() == 0) {
                  gapSyms.add(rawSym);
              }
          }

          rawSyms.removeAll(gapSyms);

          // getAllSymbols returns a Set (i.e. unordered) so
          // we convert to char array so we can sort tokens
          Symbol [] nSyms = (Symbol []) rawSyms.toArray(new Symbol [0]);
          char [] nChars = new char [nSyms.length];

          for (int i = 0; i < nSyms.length; i++) {
              nChars[i] = sTok.tokenizeSymbol(nSyms[i]).charAt(0);
          }

          Arrays.sort(nChars);
          sb.append(nChars);
          sb.append("]");
          n = sb.toString();
      } catch (Exception e) {
          throw new AssertionFailure("Couldn't initialize motif tools test", e);
      }
    }

    public MotifToolsTest(String name) {
        super(name);
    }

    public void testPlain() {
        doTest("atcg", "atcg");
    }

    public void testTwoStart() {
        doTest("aatcg", "a{2}tcg");
    }

    public void testThreeStart() {
        doTest("aaatcg", "a{3}tcg");
    }

    public void testTwoInternal() {
        doTest("attcg", "at{2}cg");
    }

    public void testThreeInternal() {
        doTest("atttcg", "at{3}cg");
    }

    public void testTwoEnd() {
        doTest("atcgg", "atcg{2}");
    }

    public void testThreeEnd() {
        doTest("atcggg", "atcg{3}");
    }

    public void testTwoOnly() {
        doTest("aa", "a{2}");
    }

    public void testThreeOnly() {
        doTest("aaa", "a{3}");
    }

    public void testAmbStart() {
        doTest("ngct", n + "gct");
    }

    public void testAmbMiddle() {
        doTest("anct", "a" + n + "ct");
    }

    public void testAmbEnd() {
        doTest("agcn", "agc" + n);
    }

    public void testTwoAmbOnly() {
        doTest("nn", n + "{2}");
    }

    void doTest(String pattern, String target) {
        try {
            assertEquals(target, MotifTools.createRegex(DNATools.createDNA(pattern)));
        } catch (IllegalSymbolException ise) {
            throw new AssertionFailure(ise);
        }
    }
}
