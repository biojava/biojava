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
package org.biojava.nbio.structure.chem;

import org.biojava.nbio.core.util.FlatFileCache;
import org.biojava.nbio.structure.chem.ChemComp;
import org.biojava.nbio.structure.chem.DownloadChemCompProvider;
import org.biojava.nbio.structure.io.LocalPDBDirectory;
import org.junit.Test;
import static org.junit.Assert.*;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.zip.GZIPOutputStream;

public class TestDownloadChemCompProvider {

	@Test
	public void testProtectedIDs(){

		String id = "CON";

		DownloadChemCompProvider prov = new DownloadChemCompProvider();
		ChemComp cc = prov.getChemComp(id);
		assertNotNull(cc);

		assertEquals(cc.getId(), id);
	}

	@Test
	public void testWeDontCacheGarbage() {
		// see #703

		File file = new File(DownloadChemCompProvider.getLocalFileName("HEM"));

		file.delete();

		DownloadChemCompProvider prov = new DownloadChemCompProvider();

		// a fake URL that should give a 404
		DownloadChemCompProvider.serverBaseUrl = "https://www.rcsb.org/non-existent-ligand-url/";

		ChemComp cc = prov.getChemComp("HEM");

		// we got a 404 back from server so we shouldn't have cached a file
		assertTrue(!file.exists());

		file.delete();

		// very important: we have a memory cache of files, we need to reset it not to pollute the cache for later tests
		FlatFileCache.clear();

		// we couldn't parse, thus there should be no content
		assertNull(cc.getName());

		// reset to default URL or otherwise we could affect other tests
		DownloadChemCompProvider.serverBaseUrl = DownloadChemCompProvider.DEFAULT_SERVER_URL;


	}

	@Test
	public void testIfWeCachedGarbageWeCanDetectIt() throws IOException {
		// see #703
		// TODO this test for the moment only asserts that we get an empty chemcomp, since we can't detect bad cached files yet

		// very important: we have a memory cache of files, we need to reset it
		FlatFileCache.clear();

		File file = new File(DownloadChemCompProvider.getLocalFileName("HEM"));

		PrintWriter pw = new PrintWriter(new GZIPOutputStream(new FileOutputStream(file)));
		pw.println("This must produce a compressed file of at least LocalPDBDirectory.MIN_PDB_FILE_SIZE bytes to avoid deletion.");
		pw.close();

		DownloadChemCompProvider prov = new DownloadChemCompProvider();

		ChemComp cc = prov.getChemComp("HEM");

		assertTrue(file.exists());

		file.delete();

		// very important: we have a memory cache of files, we need to reset it not to pollute the cache for later tests
		// we've got to reset here before asserting, in case the assertion fails
		FlatFileCache.clear();

		assertNull(cc.getName());
	}

	@Test
	public void testPathUrlTemplateRegex() {
		String[] shouldMatch = {"{ccd_id}", "{ccd_id:1_2}", "{ccd_id:1}", "{ccd_id:-1}", "abcde{ccd_id}abcde", "abcde{ccd_id:1_2}abcde", "abcde{ccd_id:-1}abcde"};
		String[] expectedCaptures = {null, "1_2", "1", "-1", null, "1_2", "-1"};
		for (int i=0; i<shouldMatch.length; i++) {
			Matcher m = DownloadChemCompProvider.CCD_ID_TEMPLATE_REGEX.matcher(shouldMatch[i]);
			assertTrue("String '"+shouldMatch[i]+"' should match the regex",m.find());
			assertEquals(expectedCaptures[i], m.group(1));
		}
		String[] shouldntMatch = {"{ccd_id:}", "{ccd_id:-1_2}", "{ccd_id:x1}", "{ccd_id:1_-2}"};
		for (String testStr : shouldntMatch) {
			Matcher m = DownloadChemCompProvider.CCD_ID_TEMPLATE_REGEX.matcher(testStr);
			assertFalse("String '"+testStr+"' should not match the regex",m.find());
		}
	}

	@Test
	public void testPathUrlTemplateExpansion() {
		String templateStr = "/my/path/{ccd_id:1_2}/hello/{ccd_id:-1}/dir/abcdef/{ccd_id:2}/12345/{ccd_id}.cif";

		String e1 = "/my/path/T/hello/P/dir/abcdef/AT/12345/ATP.cif";
		String r1 = DownloadChemCompProvider.expandPathUrlTemplate(templateStr,"ATP");
		assertEquals(e1, r1);

		String e2 = "/my/path/D/hello/D/dir/abcdef/D/12345/D.cif";
		String r2 = DownloadChemCompProvider.expandPathUrlTemplate(templateStr,"D");
		assertEquals(e2, r2);

		String e3 = "/my/path//hello//dir/abcdef//12345/.cif";
		String r3 = DownloadChemCompProvider.expandPathUrlTemplate(templateStr,"");
		assertEquals(e3, r3);

		templateStr = "/my/fixed/dir";
		assertEquals(templateStr, DownloadChemCompProvider.expandPathUrlTemplate(templateStr, "" ));
	}

}
