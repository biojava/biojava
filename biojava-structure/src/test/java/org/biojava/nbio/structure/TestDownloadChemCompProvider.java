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
package org.biojava.nbio.structure;

import org.biojava.nbio.core.util.FlatFileCache;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.junit.Test;
import static org.junit.Assert.*;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
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
	public void testRedirectWorks() {
		// since August 2017, RCSB is redirecting:
		// http://rcsb.org/pdb/files/ligand/HEM.cif  ----> http://files.org/ligands/HEM.cif
		// see #703
		
		File file = new File(DownloadChemCompProvider.getLocalFileName("HEM"));
		file.delete();
		
		DownloadChemCompProvider prov = new DownloadChemCompProvider();
		
		DownloadChemCompProvider.serverBaseUrl = "http://www.rcsb.org/pdb/files/ligand/";
		
		ChemComp cc = prov.getChemComp("HEM");		
		
		//System.out.println(file.toString());
		
		assertTrue(file.exists());

		// just in case the we did get garbage, let's clean up
		file.delete();
		
		// very important: we have a memory cache of files, we need to reset it not to pollute the cache for later tests
		FlatFileCache.clear();
		
		assertNotNull(cc);
		
		assertNotNull(cc.getName());
		
		// reset to default URL or otherwise we could affect other tests
		DownloadChemCompProvider.serverBaseUrl = DownloadChemCompProvider.DEFAULT_SERVER_URL;
	}
	
	@Test
	public void testWeDontCacheGarbage() {
		// see #703
		
		File file = new File(DownloadChemCompProvider.getLocalFileName("HEM"));
		
		file.delete();
		
		DownloadChemCompProvider prov = new DownloadChemCompProvider();
		
		// a fake URL that should give a 404
		DownloadChemCompProvider.serverBaseUrl = "http://www.rcsb.org/non-existent-ligand-url/";
		
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
		pw.println("A lot of garbage"); 
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

}
