
/*
 *                  BioJava development code
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
 * Created on Jan 16, 2006
 *
 */
package org.biojava.dasobert.das2.io;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.program.das.dasalignment.DASException;
import org.biojava.dasobert.das2.Das2Capability;
import org.biojava.dasobert.das2.Das2CapabilityImpl;
import org.biojava.dasobert.das2.Das2Source;
import org.biojava.dasobert.das2.Das2SourceImpl;
import org.biojava.dasobert.das2.DasSourceConverter;
import org.biojava.dasobert.dasregistry.Das1Source;
import org.biojava.dasobert.dasregistry.DasCoordinateSystem;
import org.biojava.dasobert.dasregistry.DasSource;
import org.xml.sax.Attributes;
import org.xml.sax.helpers.DefaultHandler;

/** a parser for the DAS2 sources response
 * 
 * @author Jonathan Warren
 * @since 
 * @version %I% %G%
 */
public class Das1RegistryCoordinatesHandler  extends DefaultHandler{

	
	List coordinates;
	DasCoordinateSystem dcs;
	StringBuffer characterdata=new StringBuffer();
	//Das2Source currentSource;
	
	public  Das1RegistryCoordinatesHandler(){
		super();
		
	}

	private void startCoordinates (String uri, String name, String qName, Attributes atts){



	}

	private DasCoordinateSystem getCoordinateSystem(String uri, String name, String qname, Attributes atts){
		// e.g. uri="http://das.sanger.ac.uk/dasregistry/coordsys/CS_LOCAL6" 
		// source="Protein Sequence" authority="UniProt" test_range="P06213" />
		characterdata=new StringBuffer();
		 dcs = new DasCoordinateSystem();
		String id = atts.getValue("uri");
		dcs.setUniqueId(id);

		String source = atts.getValue("source");
		dcs.setCategory(source);

		String authority = atts.getValue("authority");
		dcs.setName(authority);

		String test_range = atts.getValue("test_range");
		dcs.setTestCode(test_range);
		
//		String organism=atts.getValue("organism");
//		dcs.setOrganismName(organism);

		try {
			String taxidstr = atts.getValue("taxid");
			int taxid = Integer.parseInt(taxidstr);
			dcs.setNCBITaxId(taxid);
		} catch (Exception e){}

		String version = atts.getValue("version");
		if ( version != null)
			dcs.setVersion(version);

		return dcs;
	}

	public void startElement (String uri, String name, String qName, Attributes atts){
		//System.out.println("new element "+qName);

		if (qName.equals("DASCOORDINATESYSTEM")) {
			//System.out.println("started root elemtent ");
			
			
			
			

		} else
		if ( qName.equals("COORDINATES")){
			//System.out.println("coordinate qname found");
			dcs = getCoordinateSystem(uri,name,qName,atts);
			

		}       
	}

	

	public void startDocument(){
		System.out.println("starting document");
		coordinates = new ArrayList();
		
	}

	public void endElement(String uri, String name, String qName) {
		//System.out.println("end qname="+qName);
		if ( qName.equals("DASCOORDINATESYSTEM")) {
			//System.out.println("ending DASCOORDINATESYSTEM");
			//currentSource.setCoordinateSystem((DasCoordinateSystem[])coordinates.toArray(new DasCoordinateSystem[coordinates.size()]));

		}
		if ( qName.equals("COORDINATES")) {
			//System.out.println("ending COORDINATES");
			String elementContent=characterdata.toString();
			String []coordinateComponents=elementContent.split(",");
			if(coordinateComponents.length>2){
				//System.out.println("setting organism="+coordinateComponents[2]);
				dcs.setOrganismName(coordinateComponents[2]);
			}else{
			dcs.setOrganismName("");
			}
			coordinates.add(dcs);
			
		}
	}

	public DasCoordinateSystem[] getRegistryCoordinates(){    
		//System.out.println("Das2SourceHandler: source size: " + sources.size());
		return (DasCoordinateSystem[])coordinates.toArray(new DasCoordinateSystem[coordinates.size()]);
	}

	
public void characters (char ch[], int start, int length){
		
		for (int i = start; i < start + length; i++) {

			characterdata.append(ch[i]);
		}

	}


}
