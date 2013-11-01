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
 * Created on Aug 30, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.scop;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.List;

import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.scop.server.ScopDescriptions;
import org.biojava.bio.structure.scop.server.ScopDomains;
import org.biojava.bio.structure.scop.server.ScopNodes;
import org.biojava.bio.structure.scop.server.XMLUtil;


/** A class that fetches information about SCOP from a remote data-source. It requires port 80 to open for HTTP connection.
 *
 * @author Andreas Prlic
 *
 */
public class RemoteScopInstallation implements ScopDatabase {

	public static final String DEFAULT_SERVER = "http://source.rcsb.org/jfatcatserver/domains/";

	String server = DEFAULT_SERVER;

	private String version = null;

	public static void main(String[] args){

		ScopDatabase scop = new RemoteScopInstallation();
		ScopFactory.setScopDatabase(scop);

		//System.out.println(scop.getByCategory(ScopCategory.Superfamily));

		System.out.println(scop.getDomainsForPDB("4HHB"));
	}


	public String getServer() {
		return server;
	}

	public void setServer(String server) {
		this.server = server;
	}

	@Override
	public List<ScopDescription> getByCategory(ScopCategory category) {
		List<ScopDescription> results = null;
		try {
			URL u = new URL(server + "getByCategory?category="+category.toString()+"&version="+getScopVersion());
			//System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);

			if(! xml.trim().isEmpty()) {
				ScopDescriptions container = ScopDescriptions.fromXML(xml);
				results = container.getScopDescription();
			}
		} catch (Exception e) {
			System.err.println("Unable to reach "+ server + "getByCategory?category="+category.toString()+"&version="+getScopVersion());
		}
		return results;
	}

	@Override
	public List<ScopDescription> filterByClassificationId(String query) {
		List<ScopDescription> results = null;
		try {
			URL u = new URL(server + "filterByClassificationId?query="+query+"&version="+getScopVersion());
			//System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);

			if(! xml.trim().isEmpty()) {
				ScopDescriptions container = ScopDescriptions.fromXML(xml);
				results = container.getScopDescription();
			}
		} catch (Exception e){
			System.err.println("Unable to reach "+ server + "filterByClassificationId?query="+query+"&version="+getScopVersion());
		}
		return results;
	}

	@Override
	public List<ScopNode> getTree(ScopDomain domain) {
		List<ScopNode> results = null;
		try {
			URL u = new URL(server + "getTree?scopId="+domain.getScopId()+"&version="+getScopVersion());
			//System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);

			if(! xml.trim().isEmpty()) {
				ScopNodes container = ScopNodes.fromXML(xml);
				results = container.getScopNode();
			}
		} catch (Exception e){
			System.err.println("Unable to reach "+ server + "getTree?scopId="+domain.getScopId()+"&version="+getScopVersion());
		}
		return results;
	}

	@Override
	public List<ScopDomain> filterByDomainName(String query) {
		List<ScopDomain> results = null;
		try {
			URL u = new URL(server + "filterByDomainName?query="+query+"&version="+getScopVersion());
			//System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);

			if(! xml.trim().isEmpty()) {
				ScopDomains container = ScopDomains.fromXML(xml);
				results = container.getScopDomain();
			}
		} catch (Exception e){
			System.err.println("Unable to reach "+ server + "filterByDomainName?query="+query+"&version="+getScopVersion());
		}
		return results;
	}

	@Override
	public List<ScopDescription> filterByDescription(String query) {
		List<ScopDescription> results = null;
		try {
			URL u = new URL(server + "filterByDescription?query="+query+"&version="+getScopVersion());
			//System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);

			if(! xml.trim().isEmpty()) {
				ScopDescriptions container = ScopDescriptions.fromXML(xml);
				results = container.getScopDescription();
			}
		} catch (Exception e){
			System.err.println("Unable to reach "+ server + "filterByDescription?query="+query+"&version="+getScopVersion());
		}
		return results;
	}

	@Override
	public ScopDescription getScopDescriptionBySunid(int sunid) {

		ScopDescription desc = null;


		try {

			URL u = new URL(server + "getScopDescriptionBySunid?sunid="+sunid+"&version="+getScopVersion());
			//System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);

			if(! xml.trim().isEmpty()) {
				desc = XMLUtil.getScopDescriptionFromXML(xml);
			}

		} catch (Exception e){
			System.err.println("Unable to reach "+ server + "getScopDescriptionBySunid?sunid="+sunid+"&version="+getScopVersion());
		}
		return desc;
	}

	@Override
	public List<ScopDomain> getDomainsForPDB(String pdbId) {
		List<ScopDomain> results = null;
		try {
			URL u = new URL(server + "getDomainsForPDB?pdbId="+pdbId+"&version="+getScopVersion());
			//System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);

			if( !xml.trim().isEmpty()) {
				ScopDomains container = ScopDomains.fromXML(xml);
				results = container.getScopDomain();
			}
		} catch (Exception e){
			System.err.println("Unable to reach "+ server + "getDomainsForPDB?pdbId="+pdbId+"&version="+getScopVersion());
		}
		return results;
	}

	private ScopDomain requestRemoteDomainByScopID(String scopId)
	throws IOException{
		URL u = new URL(server + "getDomainByScopID?scopId="+scopId+"&version="+getScopVersion());
		//System.out.println(u);
		InputStream response = HTTPConnectionTools.getInputStream(u);
		String xml = JFatCatClient.convertStreamToString(response);

        if( !xml.trim().isEmpty()) {
        	return XMLUtil.getScopDomainFromXML(xml);
        }
        return null;
	}

	@Override
	public ScopDomain getDomainByScopID(String scopId) {
		ScopDomain desc = null;
		int i = 0;
		while ( desc == null && i < 3){
			i++;
			try {
				desc = requestRemoteDomainByScopID(scopId);
				i = 100;
				break;
			} catch (Exception e){
				System.err.println("Unable to reach "+ server + "getDomainByScopID?scopId="+scopId+"&version="+getScopVersion());
				// sleep 3 seconds and try again
				try {
					Thread.sleep(3000);
				} catch (InterruptedException e1) {
					e1.printStackTrace();
				}
			}

		}
		return desc;
	}

	@Override
	public ScopNode getScopNode(int sunid) {
		ScopNode desc = null;
		try {
			URL u = new URL(server + "getScopNode?sunid="+sunid+"&version="+getScopVersion());
			//System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);

			if( !xml.trim().isEmpty()) {
				desc = XMLUtil.getScopNodeFromXML(xml);
			}
		} catch (Exception e){
			System.err.println("Unable to reach "+ server + "getScopNode?sunid="+sunid+"&version="+getScopVersion());
		}
		return desc;
	}

	@Override
	public String getScopVersion() {
		// If no version is set, request the default version from the website
		if( version == null) {
			try {
				URL u = new URL(server + "getScopVersion");
				//System.out.println(u);
				InputStream response = HTTPConnectionTools.getInputStream(u);
				version = JFatCatClient.convertStreamToString(response);
				if( version != null)
					version = version.trim();

			} catch (Exception e){
				System.err.println("Unable to reach "+ server + "getScopVersion");
			}
		}
		return version;
	}

	public void setScopVersion(String version) {
		this.version = version;
	}

	@Override
	public List<ScopDomain> getScopDomainsBySunid(Integer sunid) {
		List<ScopDomain> results = null;
		try {
			URL u = new URL(server + "getScopDomainsBySunid?sunid="+sunid+"&version="+getScopVersion());
			//System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);

			if( !xml.trim().isEmpty()) {
				ScopDomains container = ScopDomains.fromXML(xml);
				results = container.getScopDomain();
			}
		} catch (Exception e){
			System.err.println("Unable to reach "+ server + "getScopDomainsBySunid?sunid="+sunid+"&version="+getScopVersion());
		}
		return results;
	}


	@Override
	public List<String> getComments(int sunid) {
		List<String> results = null;
		try {
			URL u = new URL(server + "getComments?sunid="+sunid+"&version="+getScopVersion());
			//System.out.println(u);
			InputStream response = HTTPConnectionTools.getInputStream(u);
			String xml = JFatCatClient.convertStreamToString(response);

			if( !xml.trim().isEmpty()) {
				results = XMLUtil.getCommentsFromXML(xml);
			}
		} catch (Exception e){
			System.err.println("Unable to reach "+ server + "getComments?sunid="+sunid+"&version="+getScopVersion());
		}
		return results;
	}

}
