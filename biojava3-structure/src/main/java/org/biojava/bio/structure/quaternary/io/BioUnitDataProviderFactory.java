package org.biojava.bio.structure.quaternary.io;

public class BioUnitDataProviderFactory {

	public static final String mmcifProviderClassName = "org.biojava.bio.structure.quaternary.io.MmCifBiolAssemblyProvider";
	
	public static final String remoteProviderClassName = "org.biojava.bio.structure.quaternary.io.RemoteBioUnitDataProvider";
	
	public static final String pdbProviderClassName 	= "org.biojava.bio.structure.quaternary.io.PDBBioUnitDataProvider";
	
	public static final String fileBasedProviderClassName = "org.biojava.bio.structure.quaternary.io.FileBasedPDBBioUnitDataProvider";
	
	public static String DEFAULT_PROVIDER_CLASSNAME =  pdbProviderClassName;
	
	private static String providerClassName = DEFAULT_PROVIDER_CLASSNAME;
	
	private BioUnitDataProviderFactory(){
		
	}
	
	public static BioUnitDataProvider getBioUnitDataProvider() {
		
		// use reflection to return a new instance...
		
		try {
			Class cls = Class.forName(providerClassName);
			//System.out.println("Using BioUnitProvider: " + providerClassName);
			return (BioUnitDataProvider) cls.newInstance();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return null;
		
	}

	public static void setBioUnitDataProvider(String className) {
		
		try {
			Class cls = Class.forName(providerClassName);
			Class interfaceClass = Class.forName("org.biojava.bio.structure.quaternary.io.BioUnitDataProvider");
			Class[] ifs = cls.getInterfaces();
			boolean found = false;
			for ( Class c : ifs){
				if ( c.equals(interfaceClass)){
					found = true;
					break;
				}
				
			}
			if ( ! found){
				System.err.println("The provided class " + className + " does not implement the correct interface!");
				return ;
			}
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			
			// we do not set classes that don't exist!
			return;
		
		}
		providerClassName = className;
	}
	
	
}
