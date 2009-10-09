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

package org.biojava.bio.molbio;

import junit.framework.TestCase;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;

/**
 * <code>RestrictionEnzymeTest</code> tests enzyme functionality.
 *
 * @author Keith James
 * @author George Waldon - fix upstream cutsites
 */
public class RestrictionEnzymeTest extends TestCase
{
    public RestrictionEnzymeTest(String name)
    {
        super(name);
    }

    public void testGetName() throws BioException
    {
        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        assertEquals("EcoRI", ecoRI.getName());
    }

    public void testGetRecognitionSite() throws BioException
    {
        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        SymbolList site = ecoRI.getRecognitionSite();

        SymbolList test = null;
        try
        {
            test = DNATools.createDNA("GAATTC");
        }
        catch (IllegalSymbolException ise)
        {
            throw new BioError(ise, "Internal error in test");
        }

        for (int i = 1; i <= site.length(); i++)
        {
            assertEquals(test.symbolAt(i), site.symbolAt(i));
        }
    }

    public void testIsPalindromic()
    {
        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        assertTrue(ecoRI.isPalindromic());

        RestrictionEnzyme bsp24I = RestrictionEnzymeManager.getEnzyme("Bsp24I");
        assertTrue(! bsp24I.isPalindromic());
    }

    public void testGetCutType()
    {
        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        assertEquals(RestrictionEnzyme.CUT_SIMPLE, ecoRI.getCutType());

        RestrictionEnzyme bsp24I = RestrictionEnzymeManager.getEnzyme("Bsp24I");
        assertEquals(RestrictionEnzyme.CUT_COMPOUND, bsp24I.getCutType());
    }

    public void testGetDownstreamCut() throws BioException
    {
        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        int [] ds = ecoRI.getDownstreamCut();

        assertEquals(1, ds[0]);
        assertEquals(5, ds[1]);
    }

    public void testGetUpstreamCut() throws BioException
    {
        RestrictionEnzyme bsp24I = RestrictionEnzymeManager.getEnzyme("Bsp24I");
        int [] us = bsp24I.getUpstreamCut();

        assertEquals(-8,  us[0]);
        assertEquals(-13, us[1]);

        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        try
        {
            us = ecoRI.getUpstreamCut();
        }
        catch (BioException be)
        {
            return;
        }

        fail("Expected BioException");
    }

    public void testGetDownstreamEndType() throws BioException
    {
        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        assertEquals(RestrictionEnzyme.OVERHANG_5PRIME,
                     ecoRI.getDownstreamEndType());

        RestrictionEnzyme apaI = RestrictionEnzymeManager.getEnzyme("ApaI");
        assertEquals(RestrictionEnzyme.OVERHANG_3PRIME,
                     apaI.getDownstreamEndType());

        RestrictionEnzyme smaI = RestrictionEnzymeManager.getEnzyme("SmaI");
        assertEquals(RestrictionEnzyme.BLUNT,
                     smaI.getDownstreamEndType());
    }

    public void testGetUpstreamEndType() throws BioException
    {
        RestrictionEnzyme bsp24I = RestrictionEnzymeManager.getEnzyme("Bsp24I");
        assertEquals(RestrictionEnzyme.OVERHANG_3PRIME,
                     bsp24I.getUpstreamEndType());

        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        try
        {
            ecoRI.getUpstreamEndType();
        }
        catch (BioException be)
        {
            return;
        }

        fail("Expected BioException");
    }
}
