package org.biojava.bio.structure.scop;

public class ScopFactory {

	//private static ScopDatabase scop = new ScopInstallation();
	
	// by default we now request SCOP via web services -> much less memory consumption
	private static ScopDatabase scop = new RemoteScopInstallation();
	
	public static ScopDatabase getSCOP(){
		return scop;
	}
	
	public static void setScopDatabase(ScopDatabase s){
		scop = s;
	}
		
		
		
}
