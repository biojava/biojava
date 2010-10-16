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

package org.biojava.bio.seq.db.emblcd;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URI;
import java.net.URL;

import junit.framework.TestCase;

/**
 * <code>EntryNamIdxReaderTest</code> is a unit test of reading the binary
 * entrynam.idx file type.
 *
 * @author Keith James
 * @since 1.2
 */
public class EntryNamIdxReaderTest extends TestCase
{
    protected EntryNamIdxReader ent;

    protected String [] seqID;
    protected   long [] rPos;
    protected   long [] sPos;
    protected    int [] fileNum;

    public EntryNamIdxReaderTest(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        URL url = EntryNamIdxReaderTest.class.getResource("entrynam.idx");

        BufferedInputStream bis = new BufferedInputStream(new
            FileInputStream(new File(url.toURI())));

        ent = new EntryNamIdxReader(bis);

        // First 10 sequence IDs.
        seqID   = new String [] {"NMA0001", "NMA0003", "NMA0004", "NMA0007", "NMA0011",
                                 "NMA0012", "NMA0013", "NMA0020", "NMA0021", "NMA0022" };

        // First 10 sequence byte offsets. There is more than one 0
        // byte offset because the database consists of more than one
        // sequence file.
        rPos    = new long [] { 0, 917, 1097, 1811, 0, 156, 270, 2379, 0, 283};

        // First 10 sequence file numbers.
        fileNum = new int [] { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3 };
    }

    protected void tearDown() throws Exception
    {
        ent.close();
        ent = null;
    }

    public void testReadFileLength()
    {
        assertTrue(1200 == ent.readFileLength());
    }

    public void testReadRecordCount()
    {
        assertTrue(30 == ent.readRecordCount());
    }

    public void testReadRecordLength()
    {
        assertTrue(30 == ent.readRecordLength());
    }

    public void testReadDBName()
    {
        assertEquals("protDB", ent.readDBName());
    }

    public void testReadDBRelease()
    {
        assertEquals("0.1", ent.readDBRelease());
    }

    public void testReadDBDate()
    {
        assertEquals("0:0:0", ent.readDBDate());
    }

    public void testReadRecord() throws IOException
    {
        for (int i = 0; i < 10; i++)
        {
            Object [] rec = ent.readRecord();

            assertEquals(seqID[i],     (String)  rec[0]);
            assertEquals(rPos[i],     ((Long)    rec[1]).longValue());
            assertEquals(0,           ((Long)    rec[2]).longValue());
            assertEquals(fileNum[i],  ((Integer) rec[3]).intValue());
        }
    }
}
