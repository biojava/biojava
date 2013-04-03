package org.biojava.bio.structure.scop;

public class ScopFactory {

	private static ScopDatabase scop = new ScopInstallation();
	
	public static ScopDatabase getSCOP(){
		return scop;
	}
	
	public static void setScopDatabase(ScopDatabase s){
		scop = s;
	}
		
		
		
}
