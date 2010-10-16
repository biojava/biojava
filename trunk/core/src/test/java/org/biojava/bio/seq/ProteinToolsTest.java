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

import junit.framework.TestCase;

import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SymbolList;

/**
 * <code>ProteinToolsTest</code> tests are to ensure that the class can be
 * instantiated and that the results are internally consistent. The
 * functionality of the classes to which <code>ProteinTools</code>
 * delegates e.g. <code>AlphabetManager</code> and
 * <code>SymbolTokenization</code> should be tested in their own unit
 * tests.
 *
 * @author Mark Schreiber
 * @author gwaldon (pyrrolysine)
 */
public class ProteinToolsTest extends TestCase
{
    protected SymbolTokenization tokens;
    protected FiniteAlphabet proteinAlpha;
    protected FiniteAlphabet proteinTAlpha;

    public ProteinToolsTest(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        tokens =
            AlphabetManager.alphabetForName("PROTEIN-TERM").getTokenization("token");
        proteinAlpha = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN");
        proteinTAlpha = (FiniteAlphabet) AlphabetManager.alphabetForName(
            "PROTEIN-TERM");
    }


    public void testGetAlpha()
    {
        assertEquals(proteinAlpha,
                     ProteinTools.getAlphabet());
        assertEquals(proteinTAlpha,
                     ProteinTools.getTAlphabet());

        assertTrue(proteinAlpha == ProteinTools.getAlphabet());
        assertTrue(proteinTAlpha == ProteinTools.getTAlphabet());
    }

    public void testCreateProtein() throws Exception
    {
        SymbolList aa = ProteinTools.createProtein(
          "arndcqeghilkmfpstwyvuo*x");

        assertTrue(aa.symbolAt(1) == ProteinTools.a());
        assertTrue(aa.symbolAt(1) == ProteinTools.ala());

        assertTrue(aa.symbolAt(2) == ProteinTools.r());
        assertTrue(aa.symbolAt(2) == ProteinTools.arg());

        assertTrue(aa.symbolAt(3) == ProteinTools.n());
        assertTrue(aa.symbolAt(3) == ProteinTools.asn());

        assertTrue(aa.symbolAt(4) == ProteinTools.d());
        assertTrue(aa.symbolAt(4) == ProteinTools.asp());

        assertTrue(aa.symbolAt(5) == ProteinTools.c());
        assertTrue(aa.symbolAt(5) == ProteinTools.cys());

        assertTrue(aa.symbolAt(6) == ProteinTools.q());
        assertTrue(aa.symbolAt(6) == ProteinTools.gln());

        assertTrue(aa.symbolAt(7) == ProteinTools.e());
        assertTrue(aa.symbolAt(7) == ProteinTools.glu());

        assertTrue(aa.symbolAt(8) == ProteinTools.g());
        assertTrue(aa.symbolAt(8) == ProteinTools.gly());

        assertTrue(aa.symbolAt(9) == ProteinTools.h());
        assertTrue(aa.symbolAt(9) == ProteinTools.his());

        assertTrue(aa.symbolAt(10) == ProteinTools.i());
        assertTrue(aa.symbolAt(10) == ProteinTools.ile());

        assertTrue(aa.symbolAt(11) == ProteinTools.l());
        assertTrue(aa.symbolAt(11) == ProteinTools.leu());

        assertTrue(aa.symbolAt(12) == ProteinTools.k());
        assertTrue(aa.symbolAt(12) == ProteinTools.lys());

        assertTrue(aa.symbolAt(13) == ProteinTools.m());
        assertTrue(aa.symbolAt(13) == ProteinTools.met());

        assertTrue(aa.symbolAt(14) == ProteinTools.f());
        assertTrue(aa.symbolAt(14) == ProteinTools.phe());

        assertTrue(aa.symbolAt(15) == ProteinTools.p());
        assertTrue(aa.symbolAt(15) == ProteinTools.pro());

        assertTrue(aa.symbolAt(16) == ProteinTools.s());
        assertTrue(aa.symbolAt(16) == ProteinTools.ser());

        assertTrue(aa.symbolAt(17) == ProteinTools.t());
        assertTrue(aa.symbolAt(17) == ProteinTools.thr());

        assertTrue(aa.symbolAt(18) == ProteinTools.w());
        assertTrue(aa.symbolAt(18) == ProteinTools.trp());

        assertTrue(aa.symbolAt(19) == ProteinTools.y());
        assertTrue(aa.symbolAt(19) == ProteinTools.tyr());

        assertTrue(aa.symbolAt(20) == ProteinTools.v());
        assertTrue(aa.symbolAt(20) == ProteinTools.val());

        assertTrue(aa.symbolAt(21) == ProteinTools.u());
        assertTrue(aa.symbolAt(21) == ProteinTools.sec());

        assertTrue(aa.symbolAt(22) == ProteinTools.o());
        assertTrue(aa.symbolAt(22) == ProteinTools.pyl());
        
        assertTrue(aa.symbolAt(23) == ProteinTools.ter());
    }
}

