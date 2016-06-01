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

import org.biojava.nbio.structure.align.util.URLConnectionTools;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.SortedSet;
import java.util.TreeSet;


/**
 * Class to fetch domains through the RCSB's REST API.
 *
 * @author Spencer Bliven
 *
 */
public class PDBDomainProvider implements DomainProvider{
	public static final String DEFAULT_PDB_HOST = "http://www.rcsb.org";
	public static final String DEFAULT_PDB_API_URL = DEFAULT_PDB_HOST + "/pdb/rest/";

	private String base;
	private int cutoff;

	/**
	 */
	public PDBDomainProvider() {
		this(DEFAULT_PDB_API_URL,40);
	}
	/**
	 * @param base
	 * @param cutoff
	 */
	public PDBDomainProvider(String base, int cutoff) {
		this.base = base;
		this.cutoff = cutoff;
	}


	/**
	 * Gets a list of domain representatives for a given PDB ID.
	 */
	@Override
	public SortedSet<String> getDomainNames(String name) {
		if ( name.length() < 4)
			throw new IllegalArgumentException("Can't interpret IDs that are shorter than 4 residues!");

		String url = String.format("%srepresentativeDomains?cluster=%s&structureId=%s",
				base, cutoff, name);
		return requestRepresentativeDomains(url);
	}
	/**
	 * Gets a list of all domain representatives
	 */
	@Override
	public SortedSet<String> getRepresentativeDomains() {
		String url = base + "representativeDomains?cluster="+ cutoff;
		return requestRepresentativeDomains(url);
	}

	/**
	 * Handles fetching and parsing XML from representativeDomains requests
	 * @param url Eg "http://www.rcsb.org/pdb/rest/representativeDomains"
	 * @return The names of all domain representatives
	 */
	private SortedSet<String> requestRepresentativeDomains(String url) {
		try {

			//System.out.println(url);

			final SortedSet<String> results = new TreeSet<String>();
			DefaultHandler handler = new DefaultHandler() {
				@Override
				public void startElement(String uri, String localName,String qName,
						Attributes attributes) throws SAXException {

					//System.out.println("Start Element :" + qName);

					if (qName.equalsIgnoreCase("representative")) {
						String name = attributes.getValue("name");
						results.add(name);
					}
				}
			};
			handleRestRequest(url,handler);
			return results;
		} catch (MalformedURLException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (SAXException e) {
			e.printStackTrace();
		} catch (ParserConfigurationException e) {
			e.printStackTrace();
		}
		return null;
	}
	/**
	 * Handles fetching and processing REST requests. The actual XML parsing is handled
	 * by the handler, which is also in charge of storing interesting data.
	 * @param url REST request
	 * @param handler SAX XML parser
	 * @throws SAXException
	 * @throws IOException
	 * @throws ParserConfigurationException
	 */
	private static void handleRestRequest(String url, DefaultHandler handler) throws SAXException, IOException, ParserConfigurationException {
		// Fetch XML stream
		URL u = new URL(url);
		InputStream response = URLConnectionTools.getInputStream(u);
		InputSource xml = new InputSource(response);

		// Parse XML
		SAXParserFactory factory = SAXParserFactory.newInstance();
		SAXParser saxParser = factory.newSAXParser();
		saxParser.parse(xml, handler);

	}


	//TODO Add methods to access http://www.rcsb.org/pdb/rest/representatives

	public static void main(String[] args){
		PDBDomainProvider dom = new PDBDomainProvider();
		String name;
		name = "2CDG";

		SortedSet<String> domains = dom.getDomainNames(name);

		System.out.println("Domains for "+name+":");
		for(String s : domains) {
			System.out.println(s);
		}

		SortedSet<String> reprs = dom.getRepresentativeDomains();
		System.out.format("%nFound %d clusters.%n",reprs.size());

		try {
			File outfile  = new File("/Users/blivens/Downloads/representativeDomainsJava.xml");
			Writer out = new BufferedWriter(new FileWriter(outfile));

			for(String repr : reprs) {
				out.write(String.format("  <representative name=\"%s\"/>%n", repr));
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}


}
