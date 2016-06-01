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
 * Created on Aug 31, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure.domain;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.biojava.nbio.structure.ResidueRange;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.SubstructureIdentifier;
import org.biojava.nbio.structure.align.client.JFatCatClient;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.URLConnectionTools;
import org.biojava.nbio.structure.scop.server.XMLUtil;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/** A class that provided PDP assignments that are loaded from a remote web server
 *
 * @author Andreas Prlic
 *
 */
public class RemotePDPProvider extends SerializableCache<String,SortedSet<String>>  implements PDPProvider{

	private static final Logger logger = LoggerFactory.getLogger(RemotePDPProvider.class);

	public static final String DEFAULT_SERVER = "http://source.rcsb.org/jfatcatserver/domains/";

	String server = DEFAULT_SERVER;

	private static String CACHE_FILE_NAME = "remotepdpdomaindefs.ser";


	public static void main(String[] args) throws IOException, StructureException{
		RemotePDPProvider me = new RemotePDPProvider(true);

		//System.out.println(scop.getByCategory(ScopCategory.Superfamily));
		SortedSet<String> pdpdomains = me.getPDPDomainNamesForPDB("4HHB");
		System.out.println(pdpdomains);

		AtomCache cache = new AtomCache();
		Structure s = me.getDomain(pdpdomains.first(), cache);
		System.out.println(s);

		me.flushCache();

	}


	public RemotePDPProvider(){
		// equivalent to this(false) but without IOException
		super(CACHE_FILE_NAME);
		disableCache();
	}


	/**
	 *
	 * @param useCache
	 * @throws IOException
	 */
	public RemotePDPProvider(boolean useCache) throws IOException {

		super(CACHE_FILE_NAME);

		if ( ! useCache) {
			disableCache();
			//else if ( serializedCache.keySet().size() < 10000){
		} else {
			// make sure we always have the latest assignments...
			loadRepresentativeDomains();
		}

	}



	/** get the ranges of representative domains from the centralized server
	 * @throws IOException if the server cannot be reached
	 */
	private void loadRepresentativeDomains() throws IOException {

		AssignmentXMLSerializer results = null;
		try {
			URL u = new URL(server + "getRepresentativePDPDomains");
			logger.info("Fetching {}",u);
			InputStream response = URLConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);
			results  = AssignmentXMLSerializer.fromXML(xml);

			Map<String,String> data = results.getAssignments();
			logger.info("got {} domain ranges for PDP domains from server.",data.size());
			for (String key: data.keySet()){
				String range = data.get(key);

				// work around list in results;

				String[] spl = range.split(",");
				SortedSet<String> value = new TreeSet<String>();

				for (String s : spl){
					value.add(s);

				}
				serializedCache.put(key, value);
			}

		} catch (MalformedURLException e){
			logger.error("Malformed PDP server: "+server,e);
			throw new IllegalArgumentException("Invalid Server: "+server, e);
		}
	}


	public String getServer() {
		return server;
	}

	public void setServer(String server) {
		this.server = server;
	}

	/**
	 * Get the structure for a particular PDP domain
	 * @param pdpDomainName PDP identifier, e.g. "PDP:4HHBAa"
	 * @param cache AtomCache, responsible for fetching and storing the coordinates
	 * @return Structure representing the PDP domain
	 * @throws IOException if the server cannot be reached
	 * @throws StructureException For errors parsing the structure
	 */
	@Override
	public Structure getDomain(String pdpDomainName, AtomCache cache) throws IOException, StructureException {
		return cache.getStructure(getPDPDomain(pdpDomainName));
	}

	/**
	 * Get a StructureIdentifier representing the specified PDP domain.
	 *
	 * @param pdpDomainName PDP domain name
	 * @return a PDPDomain representing this domain name
	 * @throws IOException if the server cannot be reached
	 */
	@Override
	public PDPDomain getPDPDomain(String pdpDomainName) throws IOException{
		SortedSet<String> domainRanges = null;
		if ( serializedCache != null){
			if ( serializedCache.containsKey(pdpDomainName)){
				domainRanges= serializedCache.get(pdpDomainName);

			}
		}


		boolean shouldRequestDomainRanges = checkDomainRanges(domainRanges);

		try {
			if (shouldRequestDomainRanges){
				URL u = new URL(server + "getPDPDomain?pdpId="+pdpDomainName);
				logger.info("Fetching {}",u);
				InputStream response = URLConnectionTools.getInputStream(u);
				String xml = JFatCatClient.convertStreamToString(response);
				domainRanges = XMLUtil.getDomainRangesFromXML(xml);
				if ( domainRanges != null)
					cache(pdpDomainName,domainRanges);
			}
		} catch (MalformedURLException e){
			logger.error("Problem generating PDP request URL for "+pdpDomainName,e);
			throw new IllegalArgumentException("Invalid PDP name: "+pdpDomainName, e);
		}

		String pdbId = null;
		List<ResidueRange> ranges = new ArrayList<ResidueRange>();
		for(String domainRange : domainRanges) {
			SubstructureIdentifier strucId = new SubstructureIdentifier(domainRange);
			if(pdbId == null) {
				pdbId = strucId.getPdbId();
			} else if(!pdbId.equals(strucId.getPdbId())) {
				// should never happen with correct server implementation
				throw new RuntimeException("Don't know how to take the union of domains from multiple PDB IDs.");
			}

			ranges.addAll(strucId.getResidueRanges());
		}
		return new PDPDomain(pdpDomainName,ranges);
	}

	/** returns true if client should fetch domain definitions from server
	 *
	 * @param domainRanges
	 * @return
	 */
	private boolean checkDomainRanges(SortedSet<String> domainRanges) {

		if ( (domainRanges == null) || (domainRanges.size() == 0)){
			return true;
		}

		for ( String d : domainRanges){
			//System.out.println("domainRange: >" + d +"< " + d.length());
			if ( (d != null) && (d.length() >0)){
				return false;
			}
		}

		return true;
	}

	/**
	 * Get a list of all PDP domains for a given PDB entry
	 * @param pdbId PDB ID
	 * @return Set of domain names, e.g. "PDP:4HHBAa"
	 * @throws IOException if the server cannot be reached
	 */
	@Override
	public SortedSet<String> getPDPDomainNamesForPDB(String pdbId) throws IOException{
		SortedSet<String> results = null;
		try {
			URL u = new URL(server + "getPDPDomainNamesForPDB?pdbId="+pdbId);
			logger.info("Fetching {}",u);
			InputStream response = URLConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);
			results  = XMLUtil.getDomainRangesFromXML(xml);

		} catch (MalformedURLException e){
			logger.error("Problem generating PDP request URL for "+pdbId,e);
			throw new IllegalArgumentException("Invalid PDB name: "+pdbId, e);
		}
		return results;
	}
}
