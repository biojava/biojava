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

import java.util.List;


import org.biojava.bio.structure.domain.SerializableCache;


/** An extension of the RemoteScopInstallation that caches some of the data locally.
 * 
 * @author Andreas Prlic
 *
 */
public class CachedRemoteScopInstallation extends SerializableCache<String,ScopDomain> implements ScopDatabase {

	private static String CACHE_FILE_NAME = "remotescopinstallation.ser";
	
	RemoteScopInstallation proxy ;
	
	SerializableCache<Integer,ScopDescription> scopDescriptionCache = new SerializableCache<Integer,ScopDescription>("scopDescriptionCache.ser");
	
	public CachedRemoteScopInstallation() {
		this(true);

	}
	
	public CachedRemoteScopInstallation(boolean useCache) {

		super(CACHE_FILE_NAME);

		proxy = new RemoteScopInstallation();
		
		if ( ! useCache) {
			System.err.println("CachedRemoteScopInstallation disableing cache");
			disableCache();
		}

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
		
		scopDescriptionCache.cache(sunid,desc);
		return desc;
	}

	
	public List<ScopDomain> getDomainsForPDB(String pdbId) {
		
		return proxy.getDomainsForPDB(pdbId);
	}

	
	public ScopDomain getDomainByScopID(String scopId) {
		
		if ( serializedCache != null){			
			if ( serializedCache.containsKey(scopId)) {
				return serializedCache.get(scopId);
			}			
		}
		
		ScopDomain dom = proxy.getDomainByScopID(scopId);
		
		
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
	}

	
}
