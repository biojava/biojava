package org.biojava.bio.structure.scop;

import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;

public class ScopFactory {

	private static ScopDatabase scop ;
	
	static {
		
		String path = System.getProperty(AbstractUserArgumentProcessor.PDB_DIR);
		if ( path == null){
			String property = "java.io.tmpdir";
			String tempdir = System.getProperty(property);
			path = tempdir;
		}
		scop = new ScopInstallation(path);
	}
	
	public static ScopDatabase getSCOP(){
		return scop;
	}
	
	public static void setScopDatabase(ScopDatabase s){
		scop = s;
	}
		
		
		
}
