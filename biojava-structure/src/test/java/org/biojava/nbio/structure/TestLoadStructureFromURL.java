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

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import static org.junit.Assert.assertNotNull;


public class TestLoadStructureFromURL {

	public static final String lineSplit = System.getProperty("file.separator");


	@Test
	public void testLoadStructureFromURL() throws IOException, StructureException{

		// we use the cache path because there's no guarantee that the PDB dir is writable
		String path = new UserConfiguration().getCacheFilePath();

		File f = new File(path, "TEST DIR");
		f.deleteOnExit();
		if ( ! f.exists()) {
			System.out.println("making dir with space:" + f);
			f.mkdir();
		}
		AtomCache c = new AtomCache(f.toString(), f.toString());
		c.setUseMmCif(false);
		c.setUseMmtf(false);
		// fetch a random small structure

		c.getStructure("1znf");

		//and now create a URL for this file
		File subdir = f;
		for(String dir :PDBFileReader.PDB_SPLIT_DIR) {
			subdir = new File(subdir,dir);
			subdir.deleteOnExit();
		}
		subdir = new File(subdir,"zn");
		File newFile = new File(subdir, "pdb1znf.ent.gz");

		subdir.deleteOnExit();
		newFile.deleteOnExit();

		URL u = newFile.toURI().toURL();

		Structure s = c.getStructure(u.toString()+"?args=test");

		assertNotNull(s);

	}
}
