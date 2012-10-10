package org.biojava.bio.structure.quaternary.io;

import java.io.InputStream;
import java.io.StringReader;
import java.net.URL;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssembly;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGen;

import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyGenXMLContainer;

import org.biojava.bio.structure.io.mmcif.model.PdbxStructAssemblyXMLContainer;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructOperList;
import org.biojava.bio.structure.io.mmcif.model.PdbxStructOperListXMLContainer;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

/** A BioUnitDataProvider that fetches the symmetry operations via remote calls to servers from RCSB PDB
 * 
 * @author Andreas Prlic
 *
 */
public class RemoteRawBioUnitDataProvider implements RawBioUnitDataProvider {

	String pdbId;
	
	public static String DEFAULT_SERVERNAME = "http://pepper.rcsb.org:8080/pdb/rest/biolassembly/";
	
	public static String NR_BIOL_APPEND = "nrBiolAssemblies?structureId=%s";
	
	public static String GET_ASSEMBLY =  "pdbxStructAssemblies?structureId=%s";
	
	public static String GET_ASSEMBLY_GENS =  "pdbxStructAssemblyGens?structureId=%s";
	
	public static String GET_STRUCT_OPER = "pdbxStructOperList?structureId=%s";
	
	String serverName;
	private static final int DEFAULT_TIMEOUT = 5000;
	
	int timeout;
	public RemoteRawBioUnitDataProvider(){
		serverName = DEFAULT_SERVERNAME;
		timeout = DEFAULT_TIMEOUT;
	}
	
	public static void main(String[] args){
		
		RemoteRawBioUnitDataProvider me = new RemoteRawBioUnitDataProvider();
		
		me.setPdbId("4hhb");
		
		System.out.println("Nr biol assemblies: " + me.getNrBiolAssemblies());
		System.out.println("has biol assembly:" + me.hasBiolAssembly());
		System.out.println("assemblies: " + me.getPdbxStructAssemblies());
		System.out.println("assemblygens:" + me.getPdbxStructAssemblyGens());
		System.out.println("operations:" + me.getPdbxStructOperList());
		
		
	}
	
	
	public int getTimeout() {
		return timeout;
	}




	public void setTimeout(int timeout) {
		this.timeout = timeout;
	}




	@Override
	public void setPdbId(String pdbId) {
		this.pdbId = pdbId;

	}

	

	@Override
	public int getNrBiolAssemblies() {
		
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
	public boolean hasBiolAssembly() {
		int nrBiolAssemblies =  getNrBiolAssemblies();
		if ( nrBiolAssemblies > 0)
			return true;
		return false;
	}

	@Override
	public PdbxStructAssembly getPdbxStructAssembly(int biolAssemblyNr) {
		PdbxStructAssembly pdbxStructAssembly = null;
		
		return pdbxStructAssembly;
	}

	@Override
	public List<PdbxStructAssemblyGen> getPdbxStructAssemblyGen(int biolAssemblyNr) {
		
		List<PdbxStructAssemblyGen> pdbxStructAssemblyGen = null;
		
		return pdbxStructAssemblyGen;
	}
	@Override
	public List<PdbxStructAssembly> getPdbxStructAssemblies() {
		String serverURL = serverName + GET_ASSEMBLY;
		List<PdbxStructAssembly> assemblies = null;
		try {
			String u = String.format(serverURL,pdbId) ;

			
			URL url = new URL(u);
			System.out.println("requesting biol assemblies from server..."  + url);
			// have a short timeout for this...
			// 5 sec
			InputStream stream = HTTPConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = JFatCatClient.convertStreamToString(stream);
				System.out.println("got XML from server: " + xml);
				PdbxStructAssemblyXMLContainer tmp = PdbxStructAssemblyXMLContainer.fromXML(xml);
				assemblies = tmp.getPdbxStructAssemblies();
				//pdbxStructAssembly = extractNrBiolAssemblies(xml);
			}
		} catch (Exception e){
			e.printStackTrace();
		}
		return assemblies;
	}

	@Override
	public List<PdbxStructAssemblyGen> getPdbxStructAssemblyGens() {
		String serverURL = serverName + GET_ASSEMBLY_GENS;
		List<PdbxStructAssemblyGen> assemblies = null;
		try {
			String u = String.format(serverURL,pdbId) ;

			
			URL url = new URL(u);
			System.out.println("requesting  biol assembly gens from server..."  + url);
			// have a short timeout for this...
			// 5 sec
			InputStream stream = HTTPConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = JFatCatClient.convertStreamToString(stream);
				System.out.println("got XML from server: " + xml);
				PdbxStructAssemblyGenXMLContainer tmp = PdbxStructAssemblyGenXMLContainer.fromXML(xml);
				assemblies = tmp.getPdbxStructAssemblyGens();
				//pdbxStructAssembly = extractNrBiolAssemblies(xml);
			}
		} catch (Exception e){
			e.printStackTrace();
		}
		
		return assemblies;
		
	}

	@Override
	public List<PdbxStructOperList> getPdbxStructOperList() {
		String serverURL = serverName + GET_STRUCT_OPER;
		List<PdbxStructOperList> oper = null;
		try {
			String u = String.format(serverURL,pdbId) ;

			
			URL url = new URL(u);
			System.out.println("requesting operators from server..."  + url);
			// have a short timeout for this...
			// 5 sec
			InputStream stream = HTTPConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = JFatCatClient.convertStreamToString(stream);
				System.out.println("got XML from server: " + xml);
				PdbxStructOperListXMLContainer container =  PdbxStructOperListXMLContainer.fromXML(xml);
				oper = container.getPdbxStructOperLists();
			}
		} catch (Exception e){
			e.printStackTrace();
		}
		
		return oper;
	}

}
