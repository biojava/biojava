package org.biojava.bio.structure;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.io.PDBFileReader;
import org.junit.Test;

import static org.junit.Assert.*;


public class TestLoadStructureFromURL {
	
	public static final String lineSplit = System.getProperty("file.separator");


	@Test
	public void testLoadStructureFromURL() throws IOException, StructureException{

		// we use the cache path because there's no guarantee that the PDB dir is writable
		String path = new UserConfiguration().getCacheFilePath();

		File f = new File(path, "SUB DIR"); 
		f.deleteOnExit();
		if ( ! f.exists()) {
			System.out.println("making dir with space:" + f);
			f.mkdir();
		}
		AtomCache c = new AtomCache(f.toString(), f.toString());
		c.setUseMmCif(false);

		// fetch a random small structure

		c.getStructure("1znf");

		//and now create a URL for this file
		File subdir = f;
		for(String dir :PDBFileReader.PDB_SPLIT_DIR) {
			subdir = new File(subdir,dir);
		}
		subdir = new File(subdir,"zn");
		File newFile = new File(subdir, "pdb1znf.ent.gz");

		subdir.deleteOnExit();
		newFile.deleteOnExit();
		
		URL u = newFile.toURI().toURL();

		System.out.println(u);

		Structure s = c.getStructure(u.toString()+"?args=test");

		System.out.println(s);

		assertNotNull(s);
		

	}
}
