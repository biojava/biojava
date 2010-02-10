/**
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

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.seq.io.ChunkedSymbolListFactory;
import org.biojava.utils.ChangeVetoException;

/**
 * JUnit test for SymbolList objects
 * @author David Huen
 * @since 1.3
 */
public class SymbolListTest extends TestCase
{
    // SymbolList lengths to run tests at.
    int testLengths[] = {100, 16384, 32000, 100000};

    // number of times to repeat each test to deal with chance
    // matches in last symbol.
    int noRepeats = 8;

    /**
     * an interface for internal classes
     */
    private interface SymListFactory {
        public SymbolList createSymbolList(Symbol [] array, FiniteAlphabet alpha, int length);
    }

    public SymbolListTest(String string)
    {
        super(string);
    }

    /**
     * creates a random SymbolList
     *
     * @param the Alphabet from which Symbols are to be drawn.  Can include ambiguity symbols.
     */
    protected Symbol [] createRandomSymbolArray(FiniteAlphabet alpha, int length)
        throws Exception
    {
        int alfaSize = alpha.size();
        AlphabetIndex indx = AlphabetManager.getAlphabetIndex(alpha);
        Random rand = new Random();

        Symbol [] array = new Symbol [length];

        for (int i=0; i < length; i++) {
            array[i] = indx.symbolForIndex(rand.nextInt(alfaSize));
        }

        return array;
    }

    /**
     * compares a SymbolList against a Symbol array
     */
    protected boolean compareSymbolLists(SymbolList list, Symbol [] array)
    {
        // array must be at least as long as SymbolList
        int length = list.length();
        if (length > array.length) return false;

        // compare symbol lists across length
        for (int i =1; i <= length; i++) {
            if (list.symbolAt(i) != array[i-1]) return false;
        }

        return true;
    }

    /**
     * Note that arrayAlpha <b>MUST NOT</b> have Symbols that are incompatible with symListAlpha!!!
     * e.g. arrayAlpha may have ambiguity symbols in it that are only implicitly defined by symListAlpha.
     * @param arrayAlpha the Alphabet from which the Symbols in the Symbol[] are to be drawn.
     * @param symListAlpha the Alphabet on which the symbolList is to be defined.
     */
    protected boolean runSymbolListTest(FiniteAlphabet arrayAlpha, FiniteAlphabet symListAlpha, int length, SymListFactory factory)
        throws Exception
    {
            // create a Symbol array of the kind required
            Symbol [] array = createRandomSymbolArray(arrayAlpha, length);
            assertNotNull(array);

            // create the required SymbolList
            SymbolList symList = factory.createSymbolList(array, symListAlpha, length);
//            if (symList == null) return false;
            assertNotNull(symList);

            // verify and return result.
            return compareSymbolLists(symList, array);
    }


    /**
     * runs repeated tests for the constructor
     * that takes a SymbolList argument
     */
    private boolean runRepeatedSymbolListTests(FiniteAlphabet arrayAlpha, FiniteAlphabet symListAlpha, SymListFactory factory)
        throws Exception
    {
        for (int i=0; i < testLengths.length; i++) {

            // setup test for specified length
            int length = testLengths[i];

            for (int j=0; j < noRepeats; j++ ) {
                assertTrue(runSymbolListTest(arrayAlpha, symListAlpha, length, factory));
            }
        }

        return true;
    }

    /**
     * set that generates a DNA alphabet including ambiguity symbols.
     */
    private FiniteAlphabet generateAmbiguousDNA()
    {
        FiniteAlphabet dna = DNATools.getDNA();

        FiniteAlphabet ambiguous = new SimpleAlphabet();

        try {
        ambiguous.addSymbol(DNATools.a());
        ambiguous.addSymbol(DNATools.c());
        ambiguous.addSymbol(DNATools.g());
        ambiguous.addSymbol(DNATools.t());

        Set chars = new HashSet();
        chars.add(DNATools.a());
        chars.add(DNATools.c());
        ambiguous.addSymbol(dna.getAmbiguity(chars));

        chars = new HashSet();
        chars.add(DNATools.a());
        chars.add(DNATools.g());
        ambiguous.addSymbol(dna.getAmbiguity(chars));

        chars = new HashSet();
        chars.add(DNATools.a());
        chars.add(DNATools.t());
        ambiguous.addSymbol(dna.getAmbiguity(chars));

        chars = new HashSet();
        chars.add(DNATools.c());
        chars.add(DNATools.g());
        ambiguous.addSymbol(dna.getAmbiguity(chars));

        chars = new HashSet();
        chars.add(DNATools.c());
        chars.add(DNATools.t());
        ambiguous.addSymbol(dna.getAmbiguity(chars));

        chars = new HashSet();
        chars.add(DNATools.g());
        chars.add(DNATools.t());
        ambiguous.addSymbol(dna.getAmbiguity(chars));

        chars = new HashSet();
        chars.add(DNATools.a());
        chars.add(DNATools.c());
        chars.add(DNATools.g());
        ambiguous.addSymbol(dna.getAmbiguity(chars));

        chars = new HashSet();
        chars.add(DNATools.a());
        chars.add(DNATools.c());
        chars.add(DNATools.t());
        ambiguous.addSymbol(dna.getAmbiguity(chars));

        chars = new HashSet();
        chars.add(DNATools.a());
        chars.add(DNATools.g());
        chars.add(DNATools.t());
        ambiguous.addSymbol(dna.getAmbiguity(chars));

        chars = new HashSet();
        chars.add(DNATools.c());
        chars.add(DNATools.g());
        chars.add(DNATools.t());
        ambiguous.addSymbol(dna.getAmbiguity(chars));

        chars.add(DNATools.n());

        return ambiguous;
        }
        catch (IllegalSymbolException ise) {
            return null;
        }
        catch (ChangeVetoException cve) {
            return null;
        }
    }

    /**
     * test for PackedSymbolList that implements ambiguity symbols.
     */
    public void testPackedSymbolListWithAmbiguitySymbols()
        throws Exception
    {
        // create an alphabet with ambiguity symbols
        FiniteAlphabet symListAlpha = (FiniteAlphabet) DNATools.getDNA();
        FiniteAlphabet arrayAlpha = generateAmbiguousDNA();
        assertNotNull(arrayAlpha);
        assertNotNull(symListAlpha);

        // create a PackedSymbolList that supports ambiguity symbols
        SymListFactory factory = new SymListFactory () {
                public SymbolList createSymbolList(Symbol [] array, FiniteAlphabet alpha, int length)
                {
                    try {
                        // get a packing
                        Packing pack = PackingFactory.getPacking(DNATools.getDNA(), true);
                        assertNotNull(pack);
                        assertNotNull(array);
                        assertNotNull(alpha);
                        assertTrue(length > 0);

                        return new PackedSymbolList(pack, array, length, alpha);
                    }
                    catch (IllegalAlphabetException iae) {
                        iae.printStackTrace();
                        return null;
                    }
                }
            };

        // exercise the PackedSymbolList implementation
        assertTrue(runRepeatedSymbolListTests(arrayAlpha, symListAlpha, factory));
    }

    /**
     * test for SimpleSymbolList that implements ambiguity symbols.
     */
    public void testSimpleSymbolListWithAmbiguitySymbols()
        throws Exception
    {
        // create an alphabet with ambiguity symbols
        FiniteAlphabet symListAlpha = (FiniteAlphabet) DNATools.getDNA();
        FiniteAlphabet arrayAlpha = generateAmbiguousDNA();
        assertNotNull(arrayAlpha);
        assertNotNull(symListAlpha);

        // create a SimpleSymbolList that supports ambiguity symbols
        SymListFactory factory = new SymListFactory () {
                public SymbolList createSymbolList(Symbol [] array, FiniteAlphabet alpha, int length)
                {
                    assertNotNull(array);
                    assertNotNull(alpha);
                    assertTrue(length > 0);

                    return new SimpleSymbolList(array, length, alpha);
                }
            };

        // exercise the PackedSymbolList implementation
        assertTrue(runRepeatedSymbolListTests(arrayAlpha, symListAlpha, factory));
    }

    /**
     * test for PackedSymbolList that does not implement ambiguity symbols.
     */
    public void testPackedSymbolListWithoutAmbiguitySymbols()
        throws Exception
    {
        // create an alphabet without ambiguity symbols
        FiniteAlphabet symListAlpha = (FiniteAlphabet) DNATools.getDNA();
        assertNotNull(symListAlpha);

        // create a PackedSymbolList that supports ambiguity symbols
        SymListFactory factory = new SymListFactory () {
                public SymbolList createSymbolList(Symbol [] array, FiniteAlphabet alpha, int length)
                {
                    try {
                        // get a packing
                        Packing pack = PackingFactory.getPacking(DNATools.getDNA(), false);
                        assertNotNull(pack);
                        assertNotNull(array);
                        assertNotNull(alpha);
                        assertTrue(length > 0);

                        return new PackedSymbolList(pack, array, length, alpha);
                    }
                    catch (IllegalAlphabetException iae) {
                        iae.printStackTrace();
                        return null;
                    }
                }
            };

        // exercise the PackedSymbolList implementation
        assertTrue(runRepeatedSymbolListTests(symListAlpha, symListAlpha, factory));
    }
    
    public void testSymbolListsEqual() throws Exception {
        SymbolList dna = DNATools.createDNA("");
        SymbolList rna = RNATools.createRNA("");
        assertFalse(dna.equals(rna));
    }

    /**
     * test for SimpleSymbolList that implements ambiguity symbols.
     */
    public void testChunkedPackedSymbolListWithAmbiguitySymbols()
        throws Exception
    {
        // create an alphabet with ambiguity symbols
        FiniteAlphabet symListAlpha = (FiniteAlphabet) DNATools.getDNA();
        FiniteAlphabet arrayAlpha = generateAmbiguousDNA();
        assertNotNull(arrayAlpha);
        assertNotNull(symListAlpha);

        // create a SimpleSymbolList that supports ambiguity symbols
        SymListFactory factory = new SymListFactory () {
                public SymbolList createSymbolList(Symbol [] array, FiniteAlphabet alpha, int length)
                {
                    assertNotNull(array);
                    assertNotNull(alpha);
                    assertTrue(length > 0);

                    try {
                        // create the factory that will make our sequence
                        ChunkedSymbolListFactory chunker = new ChunkedSymbolListFactory(new PackedSymbolListFactory(true));

                        // add symbols to it.
                        chunker.addSymbols(alpha, array, 0, length);
                        return chunker.makeSymbolList();
                    }
                        catch (IllegalAlphabetException iae) {
                        return null;
                    }
                }
            };

        // exercise the PackedSymbolList implementation
        assertTrue(runRepeatedSymbolListTests(arrayAlpha, symListAlpha, factory));
    }

    // creates a suite
    public static Test suite()
    {
        TestSuite suite = new TestSuite(SymbolListTest.class);

//        TestSuite suite = new TestSuite();
//        suite.add(new SymbolListTest("PackedSymbolListWithAmbiguitySymbols"));

        return suite;
    }

    // harness for tests
    public static void main(String [] args)
    {
        junit.textui.TestRunner.run(suite());
    }
}
