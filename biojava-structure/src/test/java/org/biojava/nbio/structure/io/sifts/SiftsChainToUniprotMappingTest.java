/**
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
 * Created on Dec 7, 2013
 * Created by Douglas Myers-Turnbull
 *
 */
package org.biojava.nbio.structure.io.sifts;

import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import static org.junit.Assert.assertEquals;

/**
 * Tests {@link SiftsChainToUniprotMapping}.
 * @author dmyersturnbull
 * @author Jose Duarte
 * @since 3.0.7
 * @since 5.0.0 uses mock SIFTS data to avoid download of full file
 */
public class SiftsChainToUniprotMappingTest {

	@Test
	public void test() throws IOException {
		SiftsChainToUniprotMapping.DEFAULT_FILE = File.createTempFile("biojavaSiftsTest-", "");
		SiftsChainToUniprotMapping.DEFAULT_FILE.deleteOnExit();
		
		BufferedReader br = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("mock_sifts.tsv")));
		PrintWriter pw = new PrintWriter(SiftsChainToUniprotMapping.DEFAULT_FILE);
		String line;
		while ( (line = br.readLine())!=null)
			pw.println(line);
		pw.close();
		br.close();
		
		SiftsChainToUniprotMapping sifts = SiftsChainToUniprotMapping.build();
		SiftsChainEntry entry = sifts.getByChainId("1hiv", "A");
		assertEquals("1hiv", entry.getPdbId());
		assertEquals("A", entry.getChainId());
		assertEquals("P04585", entry.getUniProtId());
	}

}
