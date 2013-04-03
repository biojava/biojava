package org.biojava3.alignment.aaindex;


/** Factory class to get Providers for substitution matrices the are provided by the AAINDEX database.
 *  
 * @author Andreas Prlic
 *
 */
public class AAindexFactory {

	
	private static AAIndexProvider provider = null;
	
	public static AAIndexProvider getAAIndexProvider() {
		if ( provider == null)
			provider = new DefaultAAIndexProvider();
		return provider;
	}

	public static void setAAIndexProvider(AAIndexProvider provider) {
		AAindexFactory.provider = provider;
	}
	
	
	
	
}
