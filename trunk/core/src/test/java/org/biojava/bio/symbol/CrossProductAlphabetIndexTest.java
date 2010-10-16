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
 */


package org.biojava.bio.symbol;

import java.util.Arrays;
import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;

/**
 * <p>Title: CrossProductSybmolIndexTest</p>
 * <p>Description: Tests indexing of cross product alphabets</p>
 * @author Matthew Pocock
 */

public class CrossProductAlphabetIndexTest extends TestCase {
    IntegerAlphabet.SubIntegerAlphabet subint;
    FiniteAlphabet dna;
    FiniteAlphabet prot;
    FiniteAlphabet product;
    AlphabetIndex index;

    public CrossProductAlphabetIndexTest(String name){
        super(name);
    }

    protected void setUp() throws java.lang.Exception {
        super.setUp();

        subint = IntegerAlphabet.getSubAlphabet(20,99);
        dna = DNATools.getDNA();
        prot = ProteinTools.getAlphabet();

        product = (FiniteAlphabet) AlphabetManager.getCrossProductAlphabet(
          Arrays.asList(new FiniteAlphabet[]
            { subint, dna, prot }
          )
        );

        index = new CrossProductAlphabetIndex(product);
    }

    protected void tearDown() throws java.lang.Exception {
        super.tearDown();
    }

    public void testIndex2Symbol() throws Exception {
      for(int i = 0; i < product.size(); i++) {
        product.validate(index.symbolForIndex(i));
      }
    }

    public void testSymbol2Index() throws Exception {
      for(Iterator i = product.iterator(); i.hasNext(); ) {
        Symbol s = (Symbol) i.next();
        int j = index.indexForSymbol(s);
        assertTrue(j >= 0 && j < product.size());
      }
    }

    public void testISI() throws Exception {
      for(int i = 0; i < product.size(); i++) {
        assertTrue(
          "index = index->symbol->index\t" +
          i + " -> " + index.symbolForIndex(i) + " -> " +
          index.indexForSymbol(index.symbolForIndex(i)),
          i == index.indexForSymbol(index.symbolForIndex(i))
        );
      }
    }

    public void testSIS() throws Exception {
      for(Iterator i = product.iterator(); i.hasNext(); ) {
        Symbol s = (Symbol) i.next();
        assertTrue(
          "symbol = symbol -> index -> symbol\t" +
          s + " -> " +
          index.indexForSymbol(s) + " -> " +
          index.symbolForIndex(index.indexForSymbol(s)),
          s == index.symbolForIndex(index.indexForSymbol(s))
        );
      }
    }

    public static void main(String[] args)
    throws Throwable {
      CrossProductAlphabetIndexTest test = new CrossProductAlphabetIndexTest("test");
      test.setUp();
      test.testIndex2Symbol();
      test.testSymbol2Index();
      test.testISI();
      test.testSIS();
      test.tearDown();
    }
}
