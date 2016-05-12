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
package org.biojava.nbio.structure.quaternary.io;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.client.JFatCatClient;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.HTTPConnectionTools;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

public class RemoteBioUnitDataProvider implements BioUnitDataProvider {


	private static final Logger logger = LoggerFactory.getLogger(RemoteBioUnitDataProvider.class);

	public static final String DEFAULT_SERVERNAME = "http://pepper.rcsb.org:8080/pdb/rest/biolassembly/";

	public static final String NR_BIOL_APPEND = "nrBiolAssemblies?structureId=%s";

	public static final String BIO_ASSEMBLY = "bioAssembly?structureId=%s&nr=%s";

	private static final int DEFAULT_TIMEOUT = 5000;

	private String serverName;

	private int timeout;

	public RemoteBioUnitDataProvider(){
		serverName = DEFAULT_SERVERNAME;
		timeout = DEFAULT_TIMEOUT;
	}

	@Override
	public List<BiologicalAssemblyTransformation> getBioUnitTransformationList(
			String pdbId, int biolAssemblyNr) {

		String serverURL = serverName + BIO_ASSEMBLY;

		String u = String.format(serverURL,pdbId, biolAssemblyNr) ;
		List<BiologicalAssemblyTransformation> transformations = new ArrayList<BiologicalAssemblyTransformation>();

		try {
			URL url = new URL(u);
			System.out.println("requesting biol assemblies from server..."  + url);
			// have a short timeout for this...
			// 5 sec
			InputStream stream = HTTPConnectionTools.getInputStream(url,timeout);

			String xml = null;


			xml = JFatCatClient.convertStreamToString(stream);

			if ( stream != null) {

				System.out.println(xml);
				transformations = BiologicalAssemblyTransformation.fromMultiXML(xml);
			}
		} catch (IOException e){
			logger.error("Exception caught while getting biological assemblies",e);
		} catch (SAXException e) {
			logger.error("Exception caught while getting biological assemblies",e);
		} catch (ParserConfigurationException e) {
			logger.error("Exception caught while getting biological assemblies",e);
		}

		return transformations;

	}

	@Override
	public int getNrBiolAssemblies(String pdbId) {
		String serverURL = serverName + NR_BIOL_APPEND;
		int nrBiolAssemblies = -1;
		try {
			String u = String.format(serverURL,pdbId) ;


			URL url = new URL(u);
			System.out.println("requesting nr biol assemblies from server..."  + url);
			// have a short timeout for this...
			// 5 sec
			InputStream stream = HTTPConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = JFatCatClient.convertStreamToString(stream);
				System.out.println("got XML from server: " + xml);

				nrBiolAssemblies = extractNrBiolAssemblies(xml);
			}
		} catch (IOException e){
			logger.error("Exception caught while getting number of biological assemblies",e);
		}
		return nrBiolAssemblies;
	}

	@Override
	public boolean hasBiolAssembly(String pdbId) {
		int nrBiolAssemblies =  getNrBiolAssemblies(pdbId);
		if ( nrBiolAssemblies > 0)
			return true;
		return false;
	}

	private static int extractNrBiolAssemblies(String xml) {
		int nrBiolAssemblies = -1;

		try {
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder db = factory.newDocumentBuilder();
			InputSource inStream = new InputSource();
			inStream.setCharacterStream(new StringReader(xml));
			Document doc = db.parse(inStream);

			// normalize text representation
			doc.getDocumentElement().normalize();


			//Element rootElement = doc.getDocumentElement();

			NodeList listOfPairs = doc.getElementsByTagName("nrBiolAssemblies");
			//int numArrays = listOfArrays.getLength();

			// go over the blocks
			for(int i=0; i<listOfPairs.getLength() ; i++)
			{
				Node pair       = listOfPairs.item(i);
				//NodeList valList = pair.getChildNodes();
				//int numChildren  = valList.getLength();

				NamedNodeMap map = pair.getAttributes();

				String count =  map.getNamedItem("count").getTextContent();
				nrBiolAssemblies = Integer.parseInt(count);
			}

		} catch (IOException e){
			logger.error("Exception caught while getting number of biological assemblies",e);
		} catch (SAXException e) {
			logger.error("Exception caught while getting number of biological assemblies",e);
		} catch (ParserConfigurationException e) {
			logger.error("Exception caught while getting number of biological assemblies",e);
		}
		return nrBiolAssemblies;
	}

	@Override
	public Structure getAsymUnit(String pdbId) {
		logger.error("RemoteBioUnitDataProvider getAsymUnit Not implemented yet!");
		return null;


	}

	@Override
	public void setAsymUnit(Structure s){
		// nothing to be done here so far..
	}

	@Override
	public void setAtomCache(AtomCache cache) {
		// TODO Auto-generated method stub

	}

	@Override
	public AtomCache getAtomCache() {
		// TODO Auto-generated method stub
		return null;
	}

}
