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
import java.util.SortedSet;


import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.scop.server.XMLUtil;


/** A class that provided PDP assignments that are loaded from a remote web server
 * 
 * @author Andreas Prlic
 *
 */
public class RemotePDPProvider extends SerializableCache<String,SortedSet<String>> {

	public static final String DEFAULT_SERVER = "http://source.rcsb.org/jfatcatserver/domains/";

	String server = DEFAULT_SERVER;

	private static String CACHE_FILE_NAME = "remotepdpdomaindefs.ser";


	public static void main(String[] args){

		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,"/Users/andreas/WORK/PDB");	

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

		if ( ! useCache)
			disableCache();

	}




	public String getServer() {
		return server;
	}

	public void setServer(String server) {
		this.server = server;
	}

	public Structure getDomain(String pdbDomainName, AtomCache cache){

		SortedSet<String> domainRanges = null;
		if ( serializedCache != null){
			if ( serializedCache.containsKey(pdbDomainName)){
				domainRanges= serializedCache.get(pdbDomainName);
				//System.out.println("got ranges: " + domainRanges);
			}
		}



		Structure s = null;
		try {

			if ( domainRanges == null || domainRanges.size() == 0){
				URL u = new URL(server + "getPDPDomain?pdpId="+pdbDomainName);
				System.out.println(u);
				InputStream response = HTTPConnectionTools.getInputStream(u);
				String xml = JFatCatClient.convertStreamToString(response);
				//System.out.println(xml);
				domainRanges = XMLUtil.getDomainRangesFromXML(xml);
				if ( domainRanges != null)
					cache(pdbDomainName,domainRanges);
			}
			

			if ( domainRanges.size() >0 ){
				String domainRange = domainRanges.first();
				s= cache.getStructure(domainRange);
				s.setName(pdbDomainName);								
			}

		} catch (Exception e){
			e.printStackTrace();
		}

		return s;
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
