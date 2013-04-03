package org.biojava.bio.structure;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import org.biojava.bio.structure.align.util.AtomCache;

import junit.framework.TestCase;


public class TestLoadStructureFromURL extends TestCase{
	
	private static final String FILE_SEPARATOR = System.getProperty("file.separator");
	public void testLoadStructureFromURL() throws IOException, StructureException{
		AtomCache cache = new AtomCache();
		String path = cache.getPath();
		
		path += "SUB DIR" + FILE_SEPARATOR;
		
		File f = new File(path);
		if ( ! f.exists()) {
			System.out.println("making dir with space:" + f);
			f.mkdir();
		}
		AtomCache c = new AtomCache(path, true);
		// fetch a random small structure
		
		c.getStructure("1znf");
		
		//and now create a URL for this file
		
		File newFile = new File(path+FILE_SEPARATOR+"zn"+ FILE_SEPARATOR + "pdb1znf.ent.gz");
		
		URL u = newFile.toURI().toURL();
		
		System.out.println(u);
		
		Structure s = c.getStructure(u.toString()+"?args=test");
		
		System.out.println(s);
		
		
	}
}
