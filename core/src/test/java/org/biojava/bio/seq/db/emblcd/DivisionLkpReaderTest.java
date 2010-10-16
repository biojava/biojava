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
 * <code>DivisionLkpReaderTest</code> is a unit test of reading the
 * binary division.lkp file type.
 *
 * @author Keith James
 * @since 1.2
 */
public class DivisionLkpReaderTest extends TestCase
{
    protected DivisionLkpReader div;

    public DivisionLkpReaderTest(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        URL url = DivisionLkpReaderTest.class.getResource("division.lkp");
        
        BufferedInputStream bis = new BufferedInputStream(new
            FileInputStream(new File(url.toURI())));

        div = new DivisionLkpReader(bis);
    }

    protected void tearDown() throws Exception
    {
        div.close();
    }

    public void testReadFileLength()
    {
        assertTrue(366 == div.readFileLength());
    }

    public void testReadRecordCount()
    {
        assertTrue(3 == div.readRecordCount());
    }

    public void testReadRecordLength()
    {
        assertTrue(22 == div.readRecordLength());
    }

    public void testReadDBName()
    {
        assertEquals("protDB", div.readDBName());
    }

    public void testReadDBRelease()
    {
        assertEquals("0.1", div.readDBRelease());
    }

    public void testReadDBDate()
    {
        assertEquals("0:0:0", div.readDBDate());
    }

    public void testReadRecord() throws IOException
    {
        for (int i = 1; i <= 3; i++)
        {
            Object [] rec = div.readRecord();

            assertEquals(i, ((Integer) rec[0]).intValue());
            assertEquals("protDB" + i + ".aa", (String) rec[1]);
        }
    }
}
