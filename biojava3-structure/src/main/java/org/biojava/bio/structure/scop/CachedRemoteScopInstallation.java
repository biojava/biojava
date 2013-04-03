/**
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
 * Created on Oct 12, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.scop;

import java.io.InputStream;
import java.net.URL;
import java.util.List;

import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.domain.SerializableCache;
import org.biojava.bio.structure.scop.server.ScopDomains;


/** An extension of the RemoteScopInstallation that caches some of the data locally.
 * 
 * @author Andreas Prlic
 *
 */
public class CachedRemoteScopInstallation extends SerializableCache<String,ScopDomain> implements ScopDatabase {

	private static String CACHE_FILE_NAME = "remotescopinstallation.ser";
	
	RemoteScopInstallation proxy ;
	
	SerializableCache<Integer,ScopDescription> scopDescriptionCache ;
	
	public CachedRemoteScopInstallation() {
		this(true);
		
		

	}
	
	public CachedRemoteScopInstallation(boolean useCache) {

		super(CACHE_FILE_NAME);

		proxy = new RemoteScopInstallation();
		
		scopDescriptionCache = new SerializableCache<Integer,ScopDescription>("scopDescriptionCache.ser");
		//scopDescriptionCache.setDebug(true);
		
		if ( ! useCache) {
			System.err.println("CachedRemoteScopInstallation disableing cache");
			disableCache();
			scopDescriptionCache.disableCache();
		} else {
			
			if ( serializedCache.size() < 8000){
				loadRepresentativeDomains();
			}
		}

	}


	/** get the ranges of representative domains from the centralized server
	 * 
	 */
	private void loadRepresentativeDomains() {
		
			ScopDomains results = null;
			try {
				URL u = new URL(RemoteScopInstallation.DEFAULT_SERVER + "getRepresentativeScopDomains");
				System.out.println(u);
				InputStream response = HTTPConnectionTools.getInputStream(u);
				String xml = JFatCatClient.convertStreamToString(response);			
				//System.out.println(xml);
				results  = ScopDomains.fromXML(xml);

				
				System.out.println("got " + results.getScopDomain().size() + " domain ranges for Scop domains from server.");
				for (ScopDomain dom : results.getScopDomain()){
					String scopId = dom.getScopId();
					serializedCache.put(scopId, dom);
				}
				
				
			} catch (Exception e){
				e.printStackTrace();
			}
			return ;
			
		

		
	}
	
	
	
	public List<ScopDescription> getByCategory(ScopCategory category) {
		return proxy.getByCategory(category);
	}

	
	public List<ScopDescription> filterByClassificationId(String query) {
		return proxy.filterByClassificationId(query);
	}

	
	public List<ScopNode> getTree(ScopDomain domain) {
		return proxy.getTree(domain);
	}

	
	public List<ScopDomain> filterByDomainName(String query) {
		return proxy.filterByDomainName(query);
	}

	
	public List<ScopDescription> filterByDescription(String query) {
		return proxy.filterByClassificationId(query);
	}

	
	public ScopDescription getScopDescriptionBySunid(int sunid) {
		
		ScopDescription desc = scopDescriptionCache.get(sunid);
		if ( desc != null)
			return desc;
		
		
		desc =  proxy.getScopDescriptionBySunid(sunid);
		if ( desc != null)
			scopDescriptionCache.cache(sunid,desc);
		return desc;
	}

	
	public List<ScopDomain> getDomainsForPDB(String pdbId) {
		
		return proxy.getDomainsForPDB(pdbId);
	}

	
	public ScopDomain getDomainByScopID(String scopId) {
		ScopDomain dom;
		
		if ( serializedCache != null){			
			if ( serializedCache.containsKey(scopId)) {
				dom = serializedCache.get(scopId);
				if ( dom != null) {
					return dom;
				}
			}			
		}
		
		 dom = proxy.getDomainByScopID(scopId);
		
		if ( dom != null)
			cache(scopId, dom);
		
		
		return dom;
	}


	public ScopNode getScopNode(int sunid) {
		return proxy.getScopNode(sunid);
	}


	public String getScopVersion() {
		return proxy.getScopVersion();
	}


	public List<ScopDomain> getScopDomainsBySunid(Integer sunid) {
		return proxy.getScopDomainsBySunid(sunid);
	}

	@Override
	public void flushCache() {
		System.out.println("flushing CachedRemoteScopInstallation");
		super.flushCache();
		scopDescriptionCache.flushCache();
	}

	
}
