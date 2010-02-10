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

package org.biojava.bio.seq.db;

import java.io.File;
import java.net.URI;
import java.net.URL;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import junit.framework.TestCase;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.FastaDescriptionLineParser;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.seq.io.SimpleSequenceBuilder;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;

/**
 * <code>EmblCDROMIndexStoreTest</code> contains unit tests for EMBL
 * CD-ROM format binary index <code>IndexStore</code>.
 *
 * @author Keith James
 * @since 1.2
 */
public class EmblCDROMIndexStoreTest extends TestCase
{
    protected SequenceFormat         format;
    protected Alphabet               alpha;
    protected SymbolTokenization     parser;
    protected SequenceBuilderFactory factory;
    protected SequenceDB             sequenceDB;

    protected IndexStore emblCDIndexStore;

    public EmblCDROMIndexStoreTest(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        URL divURL = EmblCDROMIndexStoreTest.class.getResource("/org/biojava/bio/seq/db/emblcd/division.lkp");
        URL entURL = EmblCDROMIndexStoreTest.class.getResource("/org/biojava/bio/seq/db/emblcd/entrynam.idx");

        File divisionLkp = new File(divURL.toURI());
        File entryNamIdx = new File(entURL.toURI());

        format  = new FastaFormat();
        alpha   = ProteinTools.getAlphabet();
        parser  = alpha.getTokenization("token");
        factory =
            new FastaDescriptionLineParser.Factory(SimpleSequenceBuilder.FACTORY);

        emblCDIndexStore = new EmblCDROMIndexStore(divisionLkp,
                                                   entryNamIdx,
                                                   format,
                                                   factory,
                                                   parser);

        ((EmblCDROMIndexStore) emblCDIndexStore).setPathPrefix(entryNamIdx.getParentFile().getAbsoluteFile());

        sequenceDB = new IndexedSequenceDB(emblCDIndexStore);
    }

    protected void tearDown() throws Exception
    {
        ((EmblCDROMIndexStore) emblCDIndexStore).close();
    }

    public void testCommit()
    {
        try
        {
            emblCDIndexStore.commit();
        }
        catch (BioException be)
        {
            return;
        }

        fail("Expected BioException");
    }

    public void testFetch() throws Exception
    {
        // Fetch from file 1
        Index i1 = emblCDIndexStore.fetch("NMA0007");
        // Fetch from file 2
        Index i2 = emblCDIndexStore.fetch("NMA0020");
        // Fetch from file 3
        Index i3 = emblCDIndexStore.fetch("NMA0030");

        // Test index fetches
        assertEquals(1811, i1.getStart());
        assertEquals("protDB1.aa", i1.getFile().getName());

        assertEquals(2379, i2.getStart());
        assertEquals("protDB2.aa", i2.getFile().getName());

        assertEquals(2510, i3.getStart());
        assertEquals("protDB3.aa", i3.getFile().getName());

        // Test actual sequence fetches
        Sequence seq = sequenceDB.getSequence("NMA0007");
        assertEquals("NMA0007", seq.getName());
        assertEquals(235, seq.length());

        seq = sequenceDB.getSequence("NMA0020");
        assertEquals("NMA0020", seq.getName());
        assertEquals(494, seq.length());

        seq = sequenceDB.getSequence("NMA0030");
        assertEquals("NMA0030", seq.getName());
        assertEquals(245, seq.length());
    }

    public void testIllegalIDFetch() throws BioException
    {
        try
        {
            Index i = emblCDIndexStore.fetch("xxxxx");
           i = i==null?null:i; //trick
        }
        catch (IllegalIDException iie)
        {
            return;
        }

        fail("Expected IllegalIDException");
    }

    public void testGetFiles()
    {
        Set files = emblCDIndexStore.getFiles();

        Set names = new HashSet();
        names.add("protDB1.aa");
        names.add("protDB2.aa");
        names.add("protDB3.aa");

        for (Iterator i = files.iterator(); i.hasNext();)
        {
            String name = (String) i.next();

            assertTrue(names.contains(name));
        }
    }

    public void testGetFormat()
    {
        assertEquals(format, emblCDIndexStore.getFormat());
    }

    public void testGetIDs()
    {
        Set ids = emblCDIndexStore.getIDs();
        assertTrue(30 == ids.size());
    }

    public void testGetName()
    {
        assertEquals("protDB", emblCDIndexStore.getName());
    }

    public void testGetSBFactory()
    {
        assertEquals(factory, emblCDIndexStore.getSBFactory());
    }

    public void testGetSymbolParser()
    {
        assertEquals(parser, emblCDIndexStore.getSymbolParser());
    }

    public void testRollback()
    {
        emblCDIndexStore.rollback();
    }

    public void testStore()
    {
        try
        {
            emblCDIndexStore.store(new SimpleIndex(new File("dummy"),
                                                   0,
                                                   -1,
                                                   "dummyID"));
        }
        catch (BioException be)
        {
            return;
        }

        fail("Expected BioException");
    }
}
