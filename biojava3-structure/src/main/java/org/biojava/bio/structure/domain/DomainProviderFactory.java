package org.biojava.bio.structure.domain;


/** A simple factory object that returns the system wide default DomainProvider
 * 
 * @author andreas
 *
 */
public class DomainProviderFactory {

	private DomainProviderFactory(){
		
	}
	
	static DomainProvider domainProvider ;
	
	
	
	public static void setDomainProvider(DomainProvider provider){
		domainProvider = provider;
		
	}
	
	public static DomainProvider getDomainProvider(){
		if ( domainProvider == null)
			domainProvider = new RemoteDomainProvider(true);
		
		return domainProvider;
	}
}
