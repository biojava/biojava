/*    
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *    AAindexStreamReaderTest.java
 *    Copyright (C) 2005 BioWeka.org
 */
package org.biojava.bio.proteomics.aaindex;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.Map;

import junit.framework.TestCase;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.SymbolList;

/**
 * Test class for {@link org.biojava.bio.proteomics.aaindex.AAindexStreamReader}. 
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class AAindexStreamReaderTest extends TestCase {

    /**
     * Main entry point.
     * @param args command line arguments
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(AAindexStreamReaderTest.class);
    }
    
    /**
     * Test for reading AAindex files.
     * @throws IOException if aaindex1 file could not be loaded.
     * @throws BioException if the format of the file is invalid.
     */
    public void testReading() throws IOException, BioException {
        InputStream input = getClass().getClassLoader().getResourceAsStream(
                "org/biojava/bio/proteomics/aaindex/aaindex1");
        AAindexStreamReader reader = new AAindexStreamReader(
                new InputStreamReader(input));
        while (reader.hasNext()) {
            AAindex aaindex = (AAindex) reader.nextTable();
            SymbolList symbols = ProteinTools.createProtein(
                    "ARNDCEQGHILKMFPSTWYV");
            for (int i = 1; i <= symbols.length(); i++) {
                aaindex.getDoubleValue(symbols.symbolAt(i));
            }
            assertNotNull(aaindex.accessionNumber());
            assertEquals(ProteinTools.getAlphabet(), aaindex.getAlphabet());
            assertNotNull(aaindex.getArticleAuthors());
            assertNotNull(aaindex.getArticleTitle());            
            assertNotNull(aaindex.getDescription());
            assertNotNull(aaindex.getJournalReference());
            assertTrue(aaindex.getLITDBEntryNumbers().length > 0);
            assertNotNull(aaindex.getName());
            Map entries = aaindex.similarEntries();
            Iterator keys = entries.keySet().iterator();
            while (keys.hasNext()) {
                String key = (String) keys.next();
                assertTrue(Double.class.isInstance(entries.get(key)));
            }
        }
        assertTrue(reader.eof());
    }

}
