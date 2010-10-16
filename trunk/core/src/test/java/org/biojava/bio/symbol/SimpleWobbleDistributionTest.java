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

import java.util.Iterator;
import java.util.Set;

import junit.framework.TestCase;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.RNATools;

/**
 * tests the SimpleWobbleDistributionTest object
 *
 * @author David Huen
 * @since 1.3
 */
public class SimpleWobbleDistributionTest extends TestCase
{
    CodonPref testPref;
    FiniteAlphabet aaAlfa = ProteinTools.getTAlphabet();
    WobbleDistribution wobbleDist;
    Symbol serine;

    public SimpleWobbleDistributionTest(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        testPref = CodonPrefTools.getCodonPreference(CodonPrefTools.JUNIT);
        assertNotNull(testPref);

        // check that I can get back a WobbleDistribution
        // all subsequent tests done on the serine distribution in the JUNIT codon preference
        SymbolList serineL = ProteinTools.createProtein("S");
        assertNotNull(serineL);

        serine = serineL.symbolAt(1);
        wobbleDist = testPref.getWobbleDistributionForSynonyms(serine);
        assertNotNull(wobbleDist);
    }

    public void testWobbleDistribution()
    {
        try {
            // check that the residue association is correct
            assertEquals(wobbleDist.getResidue().getName(), "SER");

            // just test that the stats for the serine distro
            // look sensible
            Distribution nonWobbleDist = wobbleDist.getFrequencyOfNonWobbleBases();
            assertNotNull(nonWobbleDist);

            // ensure it sums to one
            Set nonWobbleBases = wobbleDist.getNonWobbleBases();
            assertNotNull(nonWobbleBases);

            double sum = 0.0;

            for (Iterator nonWobbleBasesI = nonWobbleBases.iterator(); nonWobbleBasesI.hasNext(); ) {
                AtomicSymbol dinucl = (AtomicSymbol) nonWobbleBasesI.next();
                sum += nonWobbleDist.getWeight(dinucl);
            }

            assertTrue(Math.abs(sum - 1.0) < 0.0001);

            // look up probability of a specific non-wobble dinucleotide
            Symbol testNonWobble = CodonPrefTools.getDinucleotideAlphabet().getSymbol(RNATools.createRNA("ag").toList());
            assertTrue(
                Math.abs(
                    nonWobbleDist.getWeight(testNonWobble) - 0.333333
                )
                < 0.0001
            );

            // now check on one wobble base frequency
            Distribution wobbleD = wobbleDist.getWobbleFrequency(testNonWobble);
            assertTrue(
                Math.abs(
                    wobbleD.getWeight(RNATools.c()) - 0.75
                )
                < 0.0001
            );
        }
        catch (IllegalSymbolException ise)
        {
            ise.printStackTrace();
            fail("unexpected IllegalSymbolException");
        }
    }
}
