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
package org.biojava.nbio.structure.scop;

import org.biojava.nbio.structure.align.client.JFatCatClient;
import org.biojava.nbio.structure.align.util.URLConnectionTools;
import org.biojava.nbio.structure.domain.SerializableCache;
import org.biojava.nbio.structure.scop.server.ScopDomains;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;


/** An extension of the RemoteScopInstallation that caches some of the data locally.
 *
 * @author Andreas Prlic
 *
 */
public class CachedRemoteScopInstallation extends SerializableCache<String,ScopDomain> implements ScopDatabase {

	private static final Logger logger = LoggerFactory.getLogger(CachedRemoteScopInstallation.class);

	private static final String CACHE_FILE_NAME = "remotescopinstallation.ser";

	RemoteScopInstallation proxy ;

	SerializableCache<Integer,ScopDescription> scopDescriptionCache ;

	public CachedRemoteScopInstallation() throws IOException {
		this(true);
	}

	public CachedRemoteScopInstallation(boolean useCache) throws IOException {

		super(CACHE_FILE_NAME);

		proxy = new RemoteScopInstallation();

		scopDescriptionCache = new SerializableCache<Integer,ScopDescription>("scopDescriptionCache.ser");

		if ( ! useCache) {
			logger.warn(getClass().getSimpleName() + " disabling cache");
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
	private void loadRepresentativeDomains() throws IOException {

		URL u = null;
		try {
			u = new URL(RemoteScopInstallation.DEFAULT_SERVER + "getRepresentativeScopDomains");
		} catch (MalformedURLException e) {
			throw new IOException("URL " + RemoteScopInstallation.DEFAULT_SERVER + "getRepresentativeScopDomains" + " is wrong", e);
		}
		logger.info("Using " + u + " to download representative domains");
		InputStream response = URLConnectionTools.getInputStream(u);
		String xml = JFatCatClient.convertStreamToString(response);
		ScopDomains results  = ScopDomains.fromXML(xml);

		logger.info("got " + results.getScopDomain().size() + " domain ranges for Scop domains from server.");
		for (ScopDomain dom : results.getScopDomain()){
			String scopId = dom.getScopId();
			serializedCache.put(scopId, dom);
		}

	}



	@Override
	public List<ScopDescription> getByCategory(ScopCategory category) {
		return proxy.getByCategory(category);
	}


	@Override
	public List<ScopDescription> filterByClassificationId(String query) {
		return proxy.filterByClassificationId(query);
	}


	@Override
	public List<ScopNode> getTree(ScopDomain domain) {
		return proxy.getTree(domain);
	}


	@Override
	public List<ScopDomain> filterByDomainName(String query) {
		return proxy.filterByDomainName(query);
	}


	@Override
	public List<ScopDescription> filterByDescription(String query) {
		return proxy.filterByClassificationId(query);
	}


	@Override
	public ScopDescription getScopDescriptionBySunid(int sunid) {

		ScopDescription desc = scopDescriptionCache.get(sunid);
		if ( desc != null)
			return desc;


		desc =  proxy.getScopDescriptionBySunid(sunid);
		if ( desc != null)
			scopDescriptionCache.cache(sunid,desc);
		return desc;
	}


	@Override
	public List<ScopDomain> getDomainsForPDB(String pdbId) {

		return proxy.getDomainsForPDB(pdbId);
	}


	@Override
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


	@Override
	public ScopNode getScopNode(int sunid) {
		return proxy.getScopNode(sunid);
	}


	@Override
	public String getScopVersion() {
		return proxy.getScopVersion();
	}

	@Override
	public void setScopVersion(String version) {
		proxy.setScopVersion(version);
	}


	@Override
	public List<ScopDomain> getScopDomainsBySunid(Integer sunid) {
		return proxy.getScopDomainsBySunid(sunid);
	}

	@Override
	public void flushCache() {
		logger.info("flushing " + getClass().getSimpleName());
		super.flushCache();
		scopDescriptionCache.flushCache();
	}

	@Override
	public List<String> getComments(int sunid) {
		return new ArrayList<String>(1);
	}


}
