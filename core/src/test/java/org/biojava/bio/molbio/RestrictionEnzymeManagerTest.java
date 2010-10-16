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

import java.io.InputStream;
import java.util.Iterator;
import java.util.Set;
import java.util.regex.Pattern;

import junit.framework.TestCase;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;

/**
 * <code>RestrictionEnzymeManagerTest</code> tests manager
 * functionality.
 *
 * @author Keith James
 * @author G. Waldon
 */
public class RestrictionEnzymeManagerTest extends TestCase
{
    public RestrictionEnzymeManagerTest(String name)
    {
        super(name);
    }

    public void testGetAllEnzymes()
    {
        Set allEnz = RestrictionEnzymeManager.getAllEnzymes();
        assertEquals(40, allEnz.size());
    }

    public void testGetEnzyme()
    {
        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        assertEquals("EcoRI", ecoRI.getName());

        try
        {
            RestrictionEnzyme invalid = RestrictionEnzymeManager.getEnzyme("xxxx");
            invalid=invalid==null?null:invalid;//trick
        }
        catch (IllegalArgumentException iae)
        {
            return;
        }

        fail("Expected IllegalArgumentException");
    }

    public void testGetIsoschizomers()
    {
        Set isoAvaI = RestrictionEnzymeManager.getIsoschizomers("AvaI");
        assertEquals(1, isoAvaI.size());

        Set isoAvrI = RestrictionEnzymeManager.getIsoschizomers("AvrI");
        assertEquals(1, isoAvrI.size());

        RestrictionEnzyme avaI = RestrictionEnzymeManager.getEnzyme("AvaI");
        RestrictionEnzyme avrI = RestrictionEnzymeManager.getEnzyme("AvrI");

        assertTrue(isoAvaI.contains(avrI));
        assertTrue(isoAvrI.contains(avaI));

        try
        {
            Set invalid = RestrictionEnzymeManager.getIsoschizomers("xxxx");
            invalid=invalid==null?null:invalid;//trick
        }
        catch (IllegalArgumentException iae)
        {
            return;
        }

        fail("Expected IllegalArgumentException");
    }

    public void testGetNCutters()
    {
        Set all6Cutters = RestrictionEnzymeManager.getNCutters(6);
        assertEquals(29, all6Cutters.size());

        for (Iterator ei = all6Cutters.iterator(); ei.hasNext();)
        {
            RestrictionEnzyme e = (RestrictionEnzyme) ei.next();
            assertEquals(6, e.getRecognitionSite().length());
        }
    }

    public void testGetPatterns()
    {
        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        Pattern [] pat = RestrictionEnzymeManager.getPatterns(ecoRI);

        assertEquals("ga{2}t{2}c", pat[0].pattern());
        assertEquals("ga{2}t{2}c", pat[1].pattern());

        SymbolList site = null;
        try
        {
            site = DNATools.createDNA("a");
        }
        catch (IllegalSymbolException ise)
        {
            throw new BioError(ise, "Internal error in test");
        }

        RestrictionEnzyme custom = null;
        try
        {
            custom = new RestrictionEnzyme("custom", site, 1, 1);
        }
        catch (IllegalAlphabetException iae)
        {
            throw new BioError(iae, "Internal error in test");
        }

        try
        {
            pat = RestrictionEnzymeManager.getPatterns(custom);
        }
        catch (IllegalArgumentException iae)
        {
            return;
        }

        fail("Expected IllegalArgumentException");
    }
    
    /**
     * @author G. Waldon
     * @since 1.5
     */
    public void testGetRecognitionSequence() {
        String recognition = RestrictionEnzymeManager.getRecognitionSequence("EcoRI");
        assertEquals("G^AATTC",recognition);
    }
    
    /** Note: suppliers vary between REBASE releases.
     *
     * @author G. Waldon
     * @since 1.5
     */
    public void testgetSuppliers() {
        RestrictionEnzyme ecoRI = RestrictionEnzymeManager.getEnzyme("EcoRI");
        String suppliers = RestrictionEnzymeManager.getSuppliers(ecoRI);
        assertEquals("ABCEFGHIJKLMNOQRSTUVX",suppliers);
    }
    
    /**
     * @author G. Waldon
     * @since 1.5
     */
    public void testloadEnzymeFile() {
        String rebaseDataFileName = "org/biojava/bio/molbio/rebase.dat";
        InputStream is = getClass().getClassLoader().getResourceAsStream(rebaseDataFileName);
        RestrictionEnzymeManager.loadEnzymeFile(is,true);
        Set re = RestrictionEnzymeManager.getAllEnzymes();
        assertEquals(re.size(),1);
        try
        {
            RestrictionEnzyme hind3 = RestrictionEnzymeManager.getEnzyme("HindIII");
            assertEquals("HindIII", hind3.getName());
        }
        catch (IllegalArgumentException iae)
        {
            fail("IllegalArgumentException not expected");
        }
        is = getClass().getClassLoader().getResourceAsStream(rebaseDataFileName);
        RestrictionEnzymeManager.loadEnzymeFile(is,false);
        re = RestrictionEnzymeManager.getAllEnzymes();
        assertEquals(re.size(),2);
        try
        {
            RestrictionEnzyme hind3 = RestrictionEnzymeManager.getEnzyme("HindIII");
            assertEquals("HindIII", hind3.getName());
            RestrictionEnzyme ecor1 = RestrictionEnzymeManager.getEnzyme("EcoRI");
            assertEquals("EcoRI", ecor1.getName());
        }
        catch (IllegalArgumentException iae)
        {
            fail("IllegalArgumentException not expected");
        }
        rebaseDataFileName = "/org/biojava/bio/molbio/rebase_common.dat";
        is = RestrictionEnzymeManager.class.getResourceAsStream(rebaseDataFileName);
        RestrictionEnzymeManager.loadEnzymeFile(is,false);
    }
}
