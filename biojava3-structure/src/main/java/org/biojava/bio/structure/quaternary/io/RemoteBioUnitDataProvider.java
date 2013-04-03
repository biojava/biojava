package org.biojava.bio.structure.quaternary.io;

import java.io.InputStream;
import java.io.StringReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

public class RemoteBioUnitDataProvider implements BioUnitDataProvider {

	public static String DEFAULT_SERVERNAME = "http://pepper.rcsb.org:8080/pdb/rest/biolassembly/";

	public static String NR_BIOL_APPEND = "nrBiolAssemblies?structureId=%s";

	public static String BIO_ASSEMBLY = "bioAssembly?structureId=%s&nr=%s";

	String serverName;
	private static final int DEFAULT_TIMEOUT = 5000;

	int timeout;

	public RemoteBioUnitDataProvider(){
		serverName = DEFAULT_SERVERNAME;
		timeout = DEFAULT_TIMEOUT;
	}

	@Override
	public List<ModelTransformationMatrix> getBioUnitTransformationList(
			String pdbId, int biolAssemblyNr) {

		String serverURL = serverName + BIO_ASSEMBLY;

		String u = String.format(serverURL,pdbId, biolAssemblyNr) ;
		List<ModelTransformationMatrix> transformations = new ArrayList<ModelTransformationMatrix>();

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
				transformations = ModelTransformationMatrix.fromMultiXML(xml);
			}
		} catch (Exception e){
			e.printStackTrace();
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
		} catch (Exception e){
			e.printStackTrace();
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

		} catch (Exception e){
			e.printStackTrace();
		}
		return nrBiolAssemblies;
	}

	@Override
	public Structure getAsymUnit(String pdbId) {
		System.err.println("RemoteBioUnitDataProvider getAsymUnit Not implemented yet!");
		return null;
	
	
	}
	
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
