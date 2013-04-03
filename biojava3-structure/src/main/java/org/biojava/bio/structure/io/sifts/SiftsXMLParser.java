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
 * Created on Feb 22, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.io.sifts;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;


import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class SiftsXMLParser {

	Document dom;
	List<SiftsEntity> entities;

	public SiftsXMLParser(){
		entities = new ArrayList<SiftsEntity>();
	}

	public List<SiftsEntity> getEntities(){
		return entities;
	}


	public void parseXmlFile(InputStream is){
		entities = new ArrayList<SiftsEntity>();

		//get the factory
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

		try {

			//Using factory get an instance of document builder
			DocumentBuilder db = dbf.newDocumentBuilder();

			//parse using builder to get DOM representation of the XML file
			dom = db.parse(is);

			parseDocument();

		}catch(ParserConfigurationException pce) {
			pce.printStackTrace();
		}catch(SAXException se) {
			se.printStackTrace();
		}catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}

	

		private void parseDocument(){
			//get the root element
			Element docEle = dom.getDocumentElement();

			//get a nodelist of  entities

			NodeList nl = docEle.getElementsByTagName("entity");
			if(nl != null && nl.getLength() > 0) {
				for(int i = 0 ; i < nl.getLength();i++) {

					//get the entity element
					Element el = (Element)nl.item(i);
					//get the Employee object
					SiftsEntity e = getSiftsEntity(el);

					//add it to list
					entities.add(e);
				}
			}
		}

		/**
		 * <entity type="protein" entityId="A">
		 */
		private SiftsEntity getSiftsEntity(Element empEl) {

			//for each <employee> element get text or int values of
			//name ,id, age and name

			String type = empEl.getAttribute("type");
			String entityId = empEl.getAttribute("entityId");

			//Create a new Employee with the value read from the xml nodes
			SiftsEntity entity = new SiftsEntity(type,entityId);

			// get nodelist of segments...
			NodeList nl = empEl.getElementsByTagName("segment");
			if(nl != null && nl.getLength() > 0) {
				for(int i = 0 ; i < nl.getLength();i++) {

					//get the entity element
					Element el = (Element)nl.item(i);

					SiftsSegment s = getSiftsSegment(el);

					entity.addSegment(s);

				}
			}

			return entity;
		}

		/** segId="4hhb_A_1_140" start="1" end="140"
		 * 
		 * @param el
		 * @return
		 */
		private SiftsSegment getSiftsSegment(Element el) {

			String segId = el.getAttribute("segId");
			String start = el.getAttribute("start");
			String end = el.getAttribute("end");
			SiftsSegment seg = new SiftsSegment(segId,start,end);

			// get nodelist of segments...
			NodeList nl = el.getElementsByTagName("listResidue");
			if(nl != null && nl.getLength() > 0) {
				for(int i = 0 ; i < nl.getLength();i++) {
					//get the entity element
					Element listResidueEl = (Element)nl.item(i);

					NodeList residueNodes = listResidueEl.getElementsByTagName("residue");
					if(residueNodes != null && residueNodes.getLength() > 0) {
						for(int j = 0 ; j < residueNodes.getLength();j++) {
							Element residue = (Element) residueNodes.item(j);

							SiftsResidue pos = getResidue(residue);
							seg.addResidue(pos);
						}
					}

				}
			}


			return seg;
		}

		/**
		 *  <residue dbResNum="1" dbResName="THR">
          <crossRefDb dbSource="PDB" dbVersion="20101103"
          dbCoordSys="PDBresnum" dbAccessionId="1a4w" dbResNum="1H"
          dbResName="THR" dbChainId="L"></crossRefDb>
          <crossRefDb dbSource="UniProt" dbVersion="157-2"
          dbCoordSys="UniProt" dbAccessionId="P00734"
          dbResNum="328" dbResName="T"></crossRefDb>
          <crossRefDb dbSource="SCOP" dbVersion="1.75"
          dbCoordSys="PDBresnum" dbAccessionId="26083"
          dbResNum="1H" dbResName="THR" dbChainId="L"></crossRefDb>
          <residueDetail dbSource="MSD" property="Annotation">
          Not_Observed</residueDetail>
        </residue>

		 */
		private SiftsResidue getResidue(Element residue) {

			SiftsResidue res = new SiftsResidue();
			
			String dbResNumS = residue.getAttribute("dbResNum");
			res.setNaturalPos(Integer.parseInt(dbResNumS));
			
			String seqResName = residue.getAttribute("dbResName");
			res.setSeqResName(seqResName);
			
			boolean observed = true;
			
			String detail = getTextValue(residue, "residueDetail");
			if ( detail != null)
			//System.out.println(">"+detail+"<");
			if ( detail != null && detail.trim().equalsIgnoreCase("Not_Observed")){
				observed = false;
			}
			res.setNotObserved(! observed);
			//else if ( detail != null && detail.trim().equalsIgnoreCase("Conflict")){
				//
			//}
			
			NodeList nl = residue.getElementsByTagName("crossRefDb");
			if(nl != null && nl.getLength() > 0) {
				for(int i = 0 ; i < nl.getLength();i++) {
					//get the entity element
					Element crossRefEl = (Element)nl.item(i);

					String dbSource = crossRefEl.getAttribute("dbSource");
					String dbCoordSys = crossRefEl.getAttribute("dbCoordSys");
					String dbAccessionId = crossRefEl.getAttribute("dbAccessionId");
					String dbResNum = crossRefEl.getAttribute("dbResNum");
					String dbResName = crossRefEl.getAttribute("dbResName");
					String dbChainId = crossRefEl.getAttribute("dbChainId");

				//	System.out.println(dbSource + " " + dbCoordSys + " " + dbAccessionId + " " + dbResNum + " " + dbResName + " " + dbChainId);
					
					if ( dbSource.equals("PDB") && ( dbCoordSys.equals("PDBresnum"))){
						res.setPdbResNum(dbResNum);
						res.setPdbResName(dbResName);
						res.setChainId(dbChainId);
						res.setPdbId(dbAccessionId);
					} else if ( dbCoordSys.equals("UniProt")){
						res.setUniProtPos(Integer.parseInt(dbResNum));
						res.setUniProtResName(dbResName);
						res.setUniProtAccessionId(dbAccessionId);
					}
				}
			}
			return res;
		}



		/**
		 * I take a xml element and the tag name, look for the tag and get
		 * the text content
		 * i.e for <employee><name>John</name></employee> xml snippet if
		 * the Element points to employee node and tagName is 'name' I will return John
		 */
		private String getTextValue(Element ele, String tagName) {
			String textVal = null;
			NodeList nl = ele.getElementsByTagName(tagName);
			if(nl != null && nl.getLength() > 0) {
				Element el = (Element)nl.item(0);
				textVal = el.getFirstChild().getNodeValue();
			}

			return textVal;
		}


		

		


	}
