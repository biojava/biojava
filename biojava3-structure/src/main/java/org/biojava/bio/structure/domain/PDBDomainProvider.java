package org.biojava.bio.structure.domain;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.Writer;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;


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
		InputStream response = HTTPConnectionTools.getInputStream(u);
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
