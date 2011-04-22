/**
 * 
 */
package org.biojava.bio.structure;

import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Methods for getting the status of a PDB file (current, obsolete, etc)
 * and for accessing different versions of the structure.
 * 
 * All methods query the PDB website.
 * 
 * <p><b>TODO</b> Keep a small cache of queries around, to reduce server load
 * 
 * @author Spencer Bliven <sbliven@ucsd.edu>
 *
 */
public class PDBStatus {
	public static final String DEFAULT_PDB_SERVER = "www.rcsb.org";
	public static final String PDB_SERVER_PROPERTY = "PDB.SERVER";
	
	/**
	 * Represents the status of PDB IDs. 'OBSOLETE' and 'CURRENT' are the most
	 * common.
	 * @author Spencer Bliven <sbliven@ucsd.edu>
	 *
	 */
	public enum Status {
		OBSOLETE,
		CURRENT,
		AUTH,
		HOLD,
		HPUB,
		POLC,
		PROC,
		REFI,
		REPL,
		WAIT,
		WDRN,
		UNKNOWN;
		
		
		public static Status fromString(String statusStr) {
			Status status;
			if(statusStr.equalsIgnoreCase("OBSOLETE"))
				status = Status.OBSOLETE;
			else if(statusStr.equalsIgnoreCase("CURRENT"))
				status = Status.CURRENT;
			else if(statusStr.equalsIgnoreCase("AUTH"))
				status = Status.AUTH;
			else if(statusStr.equalsIgnoreCase("HOLD"))
				status = Status.HOLD;
			else if(statusStr.equalsIgnoreCase("HPUB"))
				status = Status.HPUB;
			else if(statusStr.equalsIgnoreCase("POLC"))
				status = Status.POLC;
			else if(statusStr.equalsIgnoreCase("PROC"))
				status = Status.PROC;
			else if(statusStr.equalsIgnoreCase("REFI"))
				status = Status.REFI;
			else if(statusStr.equalsIgnoreCase("REPL"))
				status = Status.REPL;
			else if(statusStr.equalsIgnoreCase("WAIT"))
				status = Status.WAIT;
			else if(statusStr.equalsIgnoreCase("WDRN"))
				status = Status.WDRN;
			else {
				status = null;
			}
			return status;
		}
	}
	
	/**
	 * Get the status of the PDB in question.
	 * 
	 * <p>Possible return values are:

	 * @param pdbId
	 * @return The status, or null if an error occurred.
	 */
	public static Status getStatus(String pdbId) {
		List<Attributes> attrList = getStatusIdRecords(new String[] {pdbId});
		//Expect a single record
		if(attrList == null || attrList.size() != 1) {
			System.err.println("Error getting Status for "+pdbId+" from the PDB website.");
			return null;
		}
		
		Attributes attrs = attrList.get(0);
		
		//Check that the record matches pdbId
		String id = attrs.getValue("structureId");
		if(id == null || !id.equals(pdbId)) {
			System.err.println("Error: Results returned from the query don't match "+pdbId);
			return null;
		}
		
		//Check that the status is given
		String statusStr = attrs.getValue("status");
		if(statusStr == null ) {
			System.err.println("Error: No status returned for "+pdbId);
			return null;
		}
		
		Status status = Status.fromString(statusStr);
		if(status == null) {
			System.err.println("Error: Unknown status '"+statusStr+"'");
			return null;
		}
	
		return status;
	}
	
	/**
	 * Gets the current version of a PDB ID. This is equivalent to calling
	 * {@link #getReplacement(String,boolean) getReplacement(oldPdbId,true)}
	 * 
	 * @param oldPdbId
	 * @return 
	 */
	public static String getCurrent(String oldPdbId) {
		return getReplacement(oldPdbId,true);
	}
	
	/**
	 * Gets the PDB which superseded oldPdbId. For CURRENT ids, this will
	 * be itself.
	 * 
	 * @param oldPdbId A pdb ID
	 * @param recurse If true, return the most current version of oldPdbId.
	 * 		Otherwise, just go one step newer than oldPdbId.
	 * @return The PDB which replaced oldPdbId. This may be oldPdbId itself, for
	 * 		current records. A return value of null indicates that the ID has
	 * 		been removed from the PDB.
	 */
	public static String getReplacement(String oldPdbId, boolean recurse) {
		List<Attributes> attrList = getStatusIdRecords(new String[] {oldPdbId});
		//Expect a single record
		if(attrList == null || attrList.size() != 1) {
			System.err.println("Error getting Status for "+oldPdbId+" from the PDB website.");
			return null;
		}
		
		Attributes attrs = attrList.get(0);
		
		//Check that the record matches pdbId
		String id = attrs.getValue("structureId");
		if(id == null || !id.equals(oldPdbId)) {
			System.err.println("Error: Results returned from the query don't match "+oldPdbId);
			return null;
		}
		
		//Check that the status is given
		String statusStr = attrs.getValue("status");
		if(statusStr == null ) {
			System.err.println("Error: No status returned for "+oldPdbId);
			return null;
		}
		
		Status status = Status.fromString(statusStr);
		if(status == null ) {
			System.err.println("Error: Unknown status '"+statusStr+"'");
			return null;
		}
		
		// If we're current, just return
		switch(status) {
			case CURRENT:
				return oldPdbId;
			case OBSOLETE: {
				String replacement = attrs.getValue("replacedBy");
				if(replacement == null) {
					System.err.format("Error: %s is OBSOLETE but lacks a replacedBy attribute.\n",oldPdbId);
					return null;
				}
				// Some PDBs are not replaced.
				if(replacement.equals("NONE")) {
					return null;
				}
				
				// Return the replacement.
				if(recurse) {
					return PDBStatus.getReplacement(replacement, recurse);
				}
				else {
					return replacement;
				}
			}
			case UNKNOWN:
				return null;
			default: { //TODO handle other cases explicitly. They might have other syntax than "replacedBy"
				String replacement = attrs.getValue("replacedBy");
				if(replacement == null) {
					// If no "replacedBy" attribute, assume at the root.
					// TODO is this correct?
					return oldPdbId;
				}
				// Some PDBs are not replaced.
				if(replacement.equals("NONE")) {
					return null;
				}
				
				return replacement;
			}
		}
	}
	/**
	 * Get the ID of the protein which was made obsolete by newPdbId.
	 * Equivalent to {@link #getReplaces(String,boolean) getReplaces(newPdbId, false)}
	 * 
	 * @param newPdbId PDB ID of the newer structure
	 * @return The ID of the direct ancestor of newPdbId, or <tt>null</tt> if 
	 * 		<tt>newPdbId</tt> is the original structure.
	 */
	public static String getReplaces(String newPdbId) {
		return getReplaces(newPdbId, false);
	}
	/**
	 * Get the ID of the protein which was made obsolete by newPdbId.
	 * 
	 * @param newPdbId PDB ID of the newer structure
	 * @param recurse If true, return the most oldest version of newPdbId.
	 * 		Otherwise, just go one step newer than oldPdbId.
	 * @return The ID of the direct ancestor of newPdbId, or <tt>null</tt> if 
	 * 		<tt>newPdbId</tt> is the original structure.
	 */
	public static String getReplaces(String newPdbId, boolean recurse) {
		List<Attributes> attrList = getStatusIdRecords(new String[] {newPdbId});
		//Expect a single record
		if(attrList == null || attrList.size() != 1) {
			//TODO is it possible to have multiple results? If so, should return String[] rather than null
			System.err.println("Error getting Status for "+newPdbId+" from the PDB website.");
			return null;
		}
		
		Attributes attrs = attrList.get(0);
		
		//Check that the record matches pdbId
		String id = attrs.getValue("structureId");
		if(id == null || !id.equals(newPdbId)) {
			System.err.println("Error: Results returned from the query don't match "+newPdbId);
			return null;
		}
		
		
		String replaced = attrs.getValue("replaces");
		if(replaced == null) {
			// no replaces value; assume root
			return null;
		}

		// Not the root! Return the replaced PDB.
		if(recurse) {
			String root = PDBStatus.getReplaces(replaced, recurse);
			if(root == null) {
				//replaced was already root
				return replaced;
			} else {
				return root;
			}
		} else {
			return replaced;
		}
	}
	
	
	/**
	 * Fetches the status of one or more pdbIDs from the server.
	 * 
	 * <p>Returns the results as a list of Attributes.
	 * Each attribute should contain "structureId" and "status" attributes, and
	 * possibly more.
	 * 
	 * <p>Example:</br>
	 * <tt>http://www.rcsb.org/pdb/rest/idStatus?structureID=1HHB,4HHB</tt></br>
	 *<pre>&lt;idStatus&gt;
	 *  &lt;record structureId="1HHB" status="OBSOLETE" replacedBy="4HHB"/&gt;
	 *  &lt;record structureId="4HHB" status="CURRENT" replaces="1HHB"/&gt;
	 *&lt;/idStatus&gt;
	 * </pre>
	 * 
	 * @param pdbIDs
	 * @return
	 */
	private static List<Attributes> getStatusIdRecords(String[] pdbIDs) {
		String serverName = System.getProperty(PDB_SERVER_PROPERTY);

		if ( serverName == null)
			serverName = DEFAULT_PDB_SERVER;
		else 
			System.out.format("Got System property %s=%s\n",PDB_SERVER_PROPERTY,serverName);

		// Build REST query URL
		if(pdbIDs.length < 1) {
			throw new IllegalArgumentException("No pdbIDs specified");
		}
		String urlStr = String.format("http://%s/pdb/rest/idStatus?structureId=%s",serverName,pdbIDs[0]);
		for(int i=1;i<pdbIDs.length;i++) {
			urlStr += "," + pdbIDs[i];
		}
		
		//System.out.println("Fetching " + urlStr);

		try {
			URL url = new URL(urlStr);

			InputStream uStream = url.openStream();
			
			/* // Print file directly
			BufferedReader r = new BufferedReader(new InputStreamReader(uStream));
			String line = r.readLine();
			while(line != null) {
				System.out.println(line);
				line = r.readLine();
			}
			r.close();
			
			uStream = url.openStream();
			*/
			
			
			InputSource source = new InputSource(uStream);
			SAXParserFactory parserFactory = SAXParserFactory.newInstance();
			SAXParser parser = parserFactory.newSAXParser();
			XMLReader reader = parser.getXMLReader();
			
			PDBStatusXMLHandler handler = new PDBStatusXMLHandler();
			
			reader.setContentHandler(handler);
			reader.parse(source);
			
			return handler.getRecords();
		} catch (Exception e){
			System.err.println("Problem getting status for " + pdbIDs.toString() + " from PDB server." );
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * Handles idStatus xml by storing attributes for all record elements.
	 * 
	 * @author Spencer Bliven <sbliven@ucsd.edu>
	 *
	 */
	private static class PDBStatusXMLHandler extends DefaultHandler {
		private List<Attributes> records;
		
		public PDBStatusXMLHandler() {
			records = new ArrayList<Attributes>();
		}
		
		/**
		 * @param uri
		 * @param localName
		 * @param qName
		 * @param attributes
		 * @throws SAXException
		 * @see org.xml.sax.helpers.DefaultHandler#startElement(java.lang.String, java.lang.String, java.lang.String, org.xml.sax.Attributes)
		 */
		@Override
		public void startElement(String uri, String localName, String qName,
				Attributes attributes) throws SAXException {
			//System.out.format("Starting element: uri='%s' localName='%s' qName='%s'\n", uri, localName, qName);
			if(qName.equals("record")) {
				records.add(attributes);
			}
		}


		/**
		 * @param e
		 * @throws SAXException
		 * @see org.xml.sax.helpers.DefaultHandler#error(org.xml.sax.SAXParseException)
		 */
		@Override
		public void error(SAXParseException e) throws SAXException {
			System.err.println(e.getMessage());
			super.error(e);
		}

		
		public List<Attributes> getRecords() {
			return records;
		}
	}

}