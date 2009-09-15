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

package org.biojava.bio.seq;

import java.util.Collections;
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.SymbolList;

/**
 * <code>DNAToolsTest</code> tests are to ensure that the class can be
 * instantiated and that the results are internally consistent. The
 * functionality of the classes to which <code>DNATools</code>
 * delegates e.g. <code>AlphabetManager</code> and
 * <code>SymbolTokenization</code> should be tested in their own unit
 * tests.
 *
 * @author Keith James
 */
public class DNAToolsTest extends TestCase
{
    protected SymbolTokenization dnaTokens;

    public DNAToolsTest(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        dnaTokens =
            AlphabetManager.alphabetForName("DNA").getTokenization("token");
    }

    public void testSymbols() throws Exception
    {
        assertEquals("a", dnaTokens.tokenizeSymbol(DNATools.a()));
        assertEquals("g", dnaTokens.tokenizeSymbol(DNATools.g()));
        assertEquals("c", dnaTokens.tokenizeSymbol(DNATools.c()));
        assertEquals("t", dnaTokens.tokenizeSymbol(DNATools.t()));
        assertEquals("n", dnaTokens.tokenizeSymbol(DNATools.n()));
    }

    public void testGetDNA()
    {
        assertEquals(AlphabetManager.alphabetForName("DNA"),
                     DNATools.getDNA());
    }

    public void testGetCodon(){
      List l = Collections.nCopies(3, DNATools.getDNA());
      assertEquals(DNATools.getCodonAlphabet(), AlphabetManager.getCrossProductAlphabet(l));
    }

    public void testCreateDNA() throws Exception
    {
        SymbolList dnaSyms =  DNATools.createDNA("agctn");

        assertEquals(5, dnaSyms.length());
        assertEquals("a", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(1)));
        assertEquals("g", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(2)));
        assertEquals("c", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(3)));
        assertEquals("t", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(4)));
        assertEquals("n", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(5)));
    }

    public void testIndex() throws Exception
    {
        assertEquals(0, DNATools.index(DNATools.a()));
        assertEquals(1, DNATools.index(DNATools.g()));
        assertEquals(2, DNATools.index(DNATools.c()));
        assertEquals(3, DNATools.index(DNATools.t()));
    }

    public void testForIndex()
    {
        assertEquals(DNATools.a(), DNATools.forIndex(0));
        assertEquals(DNATools.g(), DNATools.forIndex(1));
        assertEquals(DNATools.c(), DNATools.forIndex(2));
        assertEquals(DNATools.t(), DNATools.forIndex(3));
    }

    public void testComplement() throws Exception
    {
        assertEquals(DNATools.a(), DNATools.complement(DNATools.t()));
        assertEquals(DNATools.g(), DNATools.complement(DNATools.c()));
        assertEquals(DNATools.c(), DNATools.complement(DNATools.g()));
        assertEquals(DNATools.t(), DNATools.complement(DNATools.a()));
        assertEquals(DNATools.n(), DNATools.complement(DNATools.n()));

        SymbolList dnaSyms =
            DNATools.complement(DNATools.createDNA("agctn"));

        assertEquals(5, dnaSyms.length());
        assertEquals("t", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(1)));
        assertEquals("c", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(2)));
        assertEquals("g", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(3)));
        assertEquals("a", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(4)));
        assertEquals("n", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(5)));
    }

    public void testForSymbol() throws Exception
    {
        assertEquals("a", dnaTokens.tokenizeSymbol(DNATools.forSymbol('a')));
        assertEquals("g", dnaTokens.tokenizeSymbol(DNATools.forSymbol('g')));
        assertEquals("c", dnaTokens.tokenizeSymbol(DNATools.forSymbol('c')));
        assertEquals("t", dnaTokens.tokenizeSymbol(DNATools.forSymbol('t')));
    }

    public void testReverseComplement() throws Exception
    {
        SymbolList dnaSyms =
            DNATools.reverseComplement(DNATools.createDNA("agctn"));

        assertEquals(5, dnaSyms.length());
        assertEquals("n", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(1)));
        assertEquals("a", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(2)));
        assertEquals("g", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(3)));
        assertEquals("c", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(4)));
        assertEquals("t", dnaTokens.tokenizeSymbol(dnaSyms.symbolAt(5)));
    }

    public void testDNAToken() throws Exception
    {
        assertEquals('a', DNATools.dnaToken(DNATools.a()));
        assertEquals('g', DNATools.dnaToken(DNATools.g()));
        assertEquals('c', DNATools.dnaToken(DNATools.c()));
        assertEquals('t', DNATools.dnaToken(DNATools.t()));
        assertEquals('n', DNATools.dnaToken(DNATools.n()));
    }
}
