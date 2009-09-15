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
 * Tests the CodonPrefTools class and
 * the CodonPref functionality.
 * @author David Huen
 * @since 1.3
 */

public class SimpleCodonPrefTest extends TestCase
{
    CodonPref testPref;
    FiniteAlphabet aaAlfa = ProteinTools.getTAlphabet();

    public SimpleCodonPrefTest(String name)
    {
	super(name);
    }

    protected void setUp() throws Exception
    {
        testPref = CodonPrefTools.getCodonPreference(CodonPrefTools.JUNIT);
        assertNotNull(testPref);
    }

    public void testGetGeneticCode()
    {
        // get a copy of the genetic code and confirm
        // correct one retrieved
        ManyToOneTranslationTable prefGCode = testPref.getGeneticCode();

        // retrieve the UNIVERSAL code and check
        // that we have the same object
        ManyToOneTranslationTable correctGCode = RNATools.getGeneticCode(TranslationTable.UNIVERSAL);

        assertSame(prefGCode, correctGCode);
    }

    public void testGetFrequency()
    {

        // check that the sum is correct
        try {
            Distribution codonUse = testPref.getFrequency();
            assertNotNull(codonUse);

            // sum frequencies over all codons and check that it comes to one
            FiniteAlphabet codons = RNATools.getCodonAlphabet();

            double sum = 0.0;

            for (Iterator codonI = codons.iterator(); codonI.hasNext(); ) {
                sum += codonUse.getWeight((Symbol) codonI.next());
            }

            assertTrue( Math.abs(sum - 1.0) < 0.00001 );

            // pick up a specific known frequency and check
            SymbolList agc = RNATools.createRNA("agc");
            Symbol testCodon = codons.getSymbol(agc.toList());

            assertTrue(Math.abs(codonUse.getWeight(testCodon) - 0.0234375) < 0.0001);
        }
        catch (IllegalSymbolException ise) {
            fail("IllegalSymbolException occurred on codon frequency lookup");
        }
    }

    public void testGetFrequencyForSynonyms()
    {
        try {
            // get synonyms for this residue
            ManyToOneTranslationTable gCode = testPref.getGeneticCode();

            for (Iterator residueI = ProteinTools.getTAlphabet().iterator(); residueI.hasNext(); ) {
                Symbol residue = (Symbol) residueI.next();

                // filter out selenocysteine!
                if (residue.getName().equals("SEC")) continue;
                // filter out pyrrolysine!
                if (residue.getName().equals("PYL")) continue;
                
                Distribution synonymUse = testPref.getFrequencyForSynonyms(residue);
                assertNotNull(synonymUse);

                // get set of synonyms
                Set synonyms = gCode.untranslate(residue);

                double sum = 0.0;

                for (Iterator synonymsI = synonyms.iterator(); synonymsI.hasNext(); ) {
                    Symbol synonym = (Symbol) synonymsI.next();

                    sum += synonymUse.getWeight(synonym);
                }
                assertTrue(Math.abs(sum - 1.0) < 0.0001);
            }

            // check specific value: serine and agc
            SymbolList serine = ProteinTools.createProtein("S");
            assertNotNull(serine);

            Distribution serineDist = testPref.getFrequencyForSynonyms(serine.symbolAt(1));

            SymbolList agc = RNATools.createRNA("agc");
            Symbol testCodon = RNATools.getCodonAlphabet().getSymbol(agc.toList());

            assertTrue(Math.abs(serineDist.getWeight(testCodon) - 0.25) < 0.0001);

        }
        catch (IllegalSymbolException ise) {
            fail("IllegalSymbolException occurred on codon frequency lookup");
        }
    }

    public void testGetWobbleDistributionForSynonyms()
    {
        try {
            // check that I can get back a WobbleDistribution
            SymbolList serine = ProteinTools.createProtein("S");
            assertNotNull(serine);

            WobbleDistribution wobbleDist = testPref.getWobbleDistributionForSynonyms(serine.symbolAt(1));
            assertNotNull(wobbleDist);
        }
        catch (IllegalSymbolException ise) {
            fail("IllegalSymbolException occurred on wobble frequency lookup");
        }
    }
}
