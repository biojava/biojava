/*
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
 */
package org.biojava.nbio.structure.domain;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.client.JFatCatClient;
import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.align.util.URLConnectionTools;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.scop.server.XMLUtil;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * A DomainProvider that uses a mixture of SCOP and PDP domains.
 *
 * SCOP domains are preferred, with PDP providing a backup for structures where
 * SCOP has not been assigned.
 *
 * As of 2015, this class is equivalent to the method used by RCSB to define
 * representatives for structural similarity comparisons.
 */
public class RemoteDomainProvider extends SerializableCache<String,SortedSet<String>> implements DomainProvider{
	private static final Logger logger = LoggerFactory.getLogger(RemoteDomainProvider.class);

	public String url = RemotePDPProvider.DEFAULT_SERVER;

	ScopDatabase scop;
	PDPProvider pdp;

	private static String CACHE_FILE_NAME = "remotedomaincache.ser";


	public RemoteDomainProvider(){
		// equivalent to this(false) but without IOException
		super(CACHE_FILE_NAME);
		disableCache();
		scop = ScopFactory.getSCOP();
		pdp = new RemotePDPProvider();
	}

	/** initialize this provider with caching enabled
	 *
	 * @param cache
	 * @throws IOException
	 */
	public RemoteDomainProvider(boolean cache) throws IOException{
		super(CACHE_FILE_NAME);

		if( ! cache) {
			disableCache();
			//} else if ( serializedCache.keySet().size() < 20000){
		} else {
			// always load the representative assignments from server...
			// this makes sure we always have the latest assignments
			loadRepresentativeDomainAssignments();
		}

		scop = ScopFactory.getSCOP();
		pdp = new RemotePDPProvider(cache);
	}

	/** Requests the domain assignments for the current PDB IDs from the PDB.
	 * @throws IOException if the server cannot be reached
	 *
	 */
	private void loadRepresentativeDomainAssignments() throws IOException {
		AssignmentXMLSerializer results = null;
		try {
			URL u = new URL(url + "getRepresentativeDomains");
			logger.info("Fetching {}",u);
			InputStream response = URLConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);
			results  = AssignmentXMLSerializer.fromXML(xml);

			Map<String,String> data = results.getAssignments();
			logger.info("got {} ranges from server.",data.size());
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
			logger.error("Malformed Domain server: "+url,e);
			throw new IllegalArgumentException("Invalid Server: "+url, e);
		}
	}

	@Override
	public SortedSet<String> getDomainNames(String name) throws IOException, StructureException {


		if ( name.length() < 4)
			throw new IllegalArgumentException("Can't interpret IDs that are shorter than 4 residues!");

		if ( serializedCache != null){
			if ( serializedCache.containsKey(name)){
				return serializedCache.get(name);
			}
		}

		StructureName n = new StructureName(name);

		List<ScopDomain>scopDomains = scop.getDomainsForPDB(n.getPdbId());

		String chainID = n.getChainId();

		if ( scopDomains == null || scopDomains.size() == 0){
			SortedSet<String> data= getPDPDomains(n);
			cache(name,data);
			return data;
		} else {
			SortedSet<String> r = new TreeSet<String>();
			for ( ScopDomain d: scopDomains){
				StructureName s = new StructureName(d.getScopId());

				if( chainID == null){
					r.add(s.getIdentifier());

				} else if( s.getChainId().equalsIgnoreCase(n.getChainId())) {
					// SCOP IDS are case insensitive...
					r.add(s.getIdentifier());
				}
			}
			cache(name,r);
			return r;
		}



	}




	private SortedSet<String> getPDPDomains(StructureName n) throws IOException, StructureException {
		SortedSet<String> pdpDomains = pdp.getPDPDomainNamesForPDB(n.getPdbId());

		SortedSet<String> r = new TreeSet<String>();
		String chainID = n.getChainId();
		for ( String s : pdpDomains){
			StructureName d = new StructureName(s);
			if ( chainID == null)
				r.add(s);
			else if ( d.getChainId().equals(n.getChainId())){
				r.add(s);
			}
		}
		logger.info(n + " got PDP domains: "+ r);
		return r;
	}

	public static void main(String[] args) throws IOException, StructureException{
		String name ="3KIH.A";
		RemoteDomainProvider me = new RemoteDomainProvider(true);
		System.out.println(me.getDomainNames(name));
		StructureName n = new StructureName(name);
		System.out.println(n);
		//System.out.println(new  AtomCache().getStructure(name));
		me.flushCache();
	}

	@Override
	public void flushCache() {
		super.flushCache();
		if ( pdp instanceof RemotePDPProvider){
			RemotePDPProvider remotePDP = (RemotePDPProvider)pdp;
			remotePDP.flushCache();
		}
	}

	@Override
	public SortedSet<String> getRepresentativeDomains() throws IOException {

		String url = "http://source.rcsb.org/jfatcatserver/domains/getRepresentativeDomainNames";
		SortedSet<String> domainRanges = null;
		try {
			URL u = new URL(url);
			logger.info("Fetching {}",url);
			InputStream response = URLConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);
			//System.out.println(xml);
			domainRanges = XMLUtil.getDomainRangesFromXML(xml);
		} catch (MalformedURLException e){
			logger.error("Malformed Domain server: "+url,e);
			throw new IllegalArgumentException("Invalid Server: "+url, e);
		}
		return domainRanges;
	}




}
