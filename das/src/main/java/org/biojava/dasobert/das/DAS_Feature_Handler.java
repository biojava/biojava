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
 * Created on 19.03.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.dasobert.das;

import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.Attributes;

import java.util.ArrayList ;
import java.util.HashMap ;
import java.util.List;
import java.util.Map;

/**
 * a class to parse the response of a DAS - Feature request
 * @author Andreas Prlic
 *
 */
public class DAS_Feature_Handler  extends DefaultHandler{

	/**
	 * 
	 */
	List<Map<String,String>> features ;
	boolean first_flag ;
	Map<String,String> feature ;
	String featurefield ;
	StringBuffer characterdata ;
	String dasCommand ;

	int comeBackLater ;

	int maxFeatures ;

	String segmentId ;
	String version;
	String type_id;
	String type_category;
	
	
	public DAS_Feature_Handler() {
		super();

		features= new ArrayList<Map<String, String>>() ;
		first_flag = true ;
		featurefield = "" ;
		characterdata = new StringBuffer();
		dasCommand = "" ;
		comeBackLater = -1; 
		maxFeatures = -1;
		segmentId = "";
		version   = "";
		type_id = "";
		type_category="";
	}

	
	
	/** get the id information specified int the SEGMENT field of the DAS response
	 * 
	 * @return the segmentId or an emtpy string if not available
	 */
	public String getSegmentId() {
		return segmentId;
	}

	/** get the version informationspecified in the SEGMENT field of the DAS response
	 * 
	 * @return the version information of an empty string if not available 
	 */
	public String getVersion() {
		return version;
	}
	
	public boolean isMD5Checksum(){
		
		if ((version != null) && (version.length() == 32))
			return true;
		return false;
	}


	/** specifies a maximum number of features to be downloaded. if a
	server returns more, they will be ignored.  default is to load
	all features 
    @param max the maximium number of features to be downloaded
	 */

	public void setMaxFeatures(int max) {
		maxFeatures = max;
	}

	public int getMaxFeatures() {
		return maxFeatures;
	}

	public void setDASCommand(String cmd) { dasCommand = cmd ;}
	public String getDASCommand() { return dasCommand; }

	public List<Map<String,String>> get_features() {
		return features ;
	}

	public int getComBackLater(){
		return comeBackLater;
	}

	void start_feature(String uri, String name, String qName, Attributes atts) {

		if (( maxFeatures > 0 ) && ( features.size() > maxFeatures ) ) {
			characterdata = new StringBuffer();
			return;
		}
		feature = new HashMap<String,String>() ;
		String id 	= atts.getValue("id");
		feature.put("id",id);
		feature.put("dassource",dasCommand);
		characterdata = new StringBuffer();
	}

	void add_featuredata(String uri, String name, String qName) {
		//System.out.println("featurefield "+featurefield+ " data "+characterdata);
		// NOTE can have multiple lines ..

		if (( maxFeatures > 0 ) && ( features.size() > maxFeatures ) ) {
			return;
		}


		String data = (String)feature.get(featurefield);
		String featureText = characterdata.toString();
		if (data != null){
			featureText = data + " " + featureText;
		}

		if ( qName.equals("TYPE")){
			if ( featureText.length() < 1)
				featureText = type_id;
		
			feature.put("TYPE_ID",type_id);
			feature.put("TYPE_CATEGORY", type_category);
			type_id = "";
			type_category="";
		}
		
		
		feature.put(featurefield,featureText);
		featurefield = "";
		characterdata = new StringBuffer();
	}

	private void addLink(String uri, String name, String qName, Attributes atts) {
		String href = atts.getValue("href");
		feature.put("LINK",href);
		characterdata= new StringBuffer();
		featurefield = "LINK-TEXT";

	}

	private void addGroup(String uri, String name, String qName, Attributes atts) {
		String id = atts.getValue("id");
		feature.put("GROUP",id);
		characterdata= new StringBuffer();
		featurefield = "GROUP";
	}
	
	public void startElement (String uri, String name, String qName, Attributes atts){
		//System.out.println("new element "+qName);

		if (qName.equals("FEATURE")) 
			start_feature(uri,  name,  qName,  atts);
		else if ( qName.equals("LINK"))
			addLink(uri,name,qName, atts);
		else if ( qName.equals("GROUP"))
			addGroup(uri,name,qName, atts);
		else if ( qName.equals("METHOD") || 
				qName.equals("TYPE") ||
				qName.equals("START") ||
				qName.equals("END") ||
				qName.equals("NOTE") ||                
				qName.equals("SCORE") ||
				qName.equals("ORIENTATION") 				
		){
			characterdata = new StringBuffer();
			featurefield = qName ;
		} else if (qName.equals("SEGMENT")){
			String id = atts.getValue("id");
			if (id != null)
				segmentId = id;
			
			String v = atts.getValue("version");
			if ( v != null)
				version = v;
			
			
		}
		if ( qName.equals("TYPE")){
			type_id      = atts.getValue("id");
			type_category= atts.getValue("category");
		}

	}



	public void endElement(String uri, String name, String qName) {

		if ( qName.equals("METHOD") || 
				qName.equals("TYPE") ||
				qName.equals("START") ||
				qName.equals("END") ||
				qName.equals("NOTE") || 
				qName.equals("LINK") || 
				qName.equals("SCORE") ||
				qName.equals("ORIENTATION") ||
				qName.equals("GROUP")
		) {
			add_featuredata(uri,name,qName);
		}
		else if ( qName.equals("FEATURE")) {

			if ( maxFeatures > 0 ) {
				if ( features.size() < maxFeatures ) {
					features.add(feature);
				} 
			} else {
				// no restriction
				features.add(feature);
			}
		}
	}

	public void characters (char ch[], int start, int length){
		if ( maxFeatures > 0)
			if ( features.size() > maxFeatures )
				return;

		for (int i = start; i < start + length; i++) {

			characterdata.append(ch[i]);
		}

	}

}
