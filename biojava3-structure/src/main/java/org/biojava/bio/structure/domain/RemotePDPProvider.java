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
package org.biojava.bio.structure.domain;

import java.io.InputStream;

import java.net.URL;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;


import org.biojava.bio.structure.Structure;

import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.scop.server.XMLUtil;


/** A class that provided PDP assignments that are loaded from a remote web server
 * 
 * @author Andreas Prlic
 *
 */
public class RemotePDPProvider extends SerializableCache<String,SortedSet<String>>  implements PDPProvider{

	public static final String DEFAULT_SERVER = "http://source.rcsb.org/jfatcatserver/domains/";

	String server = DEFAULT_SERVER;

	private static String CACHE_FILE_NAME = "remotepdpdomaindefs.ser";


	public static void main(String[] args){

		System.setProperty(AbstractUserArgumentProcessor.CACHE_DIR,"/Users/ap3/WORK/PDB");	

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
		this(false);
	}


	public RemotePDPProvider(boolean useCache) {

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
	 * 
	 */
	private void loadRepresentativeDomains() {

		AssignmentXMLSerializer results = null;
		try {
			URL u = new URL(server + "getRepresentativePDPDomains");
			System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);			
			//System.out.println(xml);
			results  = AssignmentXMLSerializer.fromXML(xml);

			Map<String,String> data = results.getAssignments();
			System.out.println("got " + data.size() + " domain ranges for PDP domains from server.");
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

		} catch (Exception e){
			e.printStackTrace();
		}
		return ;

	}


	public String getServer() {
		return server;
	}

	public void setServer(String server) {
		this.server = server;
	}

	public Structure getDomain(String pdpDomainName, AtomCache cache){

		SortedSet<String> domainRanges = null;
		if ( serializedCache != null){
			if ( serializedCache.containsKey(pdpDomainName)){
				domainRanges= serializedCache.get(pdpDomainName);

			}
		}



		Structure s = null;
		try {

			

			boolean shouldRequestDomainRanges = checkDomainRanges(domainRanges);
			
			if (shouldRequestDomainRanges){
				URL u = new URL(server + "getPDPDomain?pdpId="+pdpDomainName);
				System.out.println(u);
				InputStream response = HTTPConnectionTools.getInputStream(u);
				String xml = JFatCatClient.convertStreamToString(response);
				//System.out.println(xml);
				domainRanges = XMLUtil.getDomainRangesFromXML(xml);
				if ( domainRanges != null)
					cache(pdpDomainName,domainRanges);
			}

			
			int i =0 ;
			StringBuffer r = new StringBuffer();
			for (String domainRange : domainRanges){
				if ( ! domainRange.contains("."))
					r.append(domainRange);
				else {
					String[] spl = domainRange.split("\\.");
					if ( spl.length>1)
						r.append(spl[1]);
					else {
						System.out.println("not sure what to do with " + domainRange);
					}
				}
				i++;
				if ( i < domainRanges.size()) {
					r.append(",");
				}
			}

			String ranges = r.toString();

			StructureName sname = new  StructureName(pdpDomainName);
			Structure tmp = cache.getStructure(sname.getPdbId());
			//System.out.println(ranges);
			s = StructureTools.getSubRanges(tmp, ranges);


			s.setName(pdpDomainName);

		} catch (Exception e){
			e.printStackTrace();
		}

		return s;
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


	public SortedSet<String> getPDPDomainNamesForPDB(String pdbId){
		SortedSet<String> results = null;
		try {
			URL u = new URL(server + "getPDPDomainNamesForPDB?pdbId="+pdbId);
			System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);
			//System.out.println(xml);
			results  = XMLUtil.getDomainRangesFromXML(xml);

		} catch (Exception e){
			e.printStackTrace();
		}
		return results;
	}
}
