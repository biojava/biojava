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
 * Created on 15.04.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.dasobert.dasregistry;

import java.util.ArrayList;
import java.util.Date ;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;


import org.biojava.dasobert.das.Capabilities;
import org.biojava.dasobert.das2.io.DasSourceWriter;
import org.biojava.dasobert.das2.io.DasSourceWriterImpl;
import org.biojava.utils.xml.PrettyXMLWriter;


/** a simple Bean class to be returned via SOAP
 * @author Andreas Prlic
 */

public class Das1Source implements DasSource {
	String url;
	protected String nickname;
	String adminemail;
	String description ;
	List<DasCoordinateSystem> coordinateSystem;
	List<Capabilities> capabilities;
	List<String> labels;
	String helperurl;    
	Date   registerDate; 
	Date   leaseDate;
	String id;
	boolean local;
	Map<String,String> properties;
	boolean alertAdmin;
	private List<String> validCapabilities;

	public static String EMPTY_ID = "UNK:-1" ;

	public Das1Source () {
		id               = EMPTY_ID;
		url              = "empty";
		adminemail       = "" ;
		description      = "" ;
		//String empty     = "" ;
		nickname         = "" ;
		coordinateSystem =new  ArrayList<DasCoordinateSystem>();

		capabilities     = new ArrayList<Capabilities>();
		labels 	         = new ArrayList<String>();

		registerDate     = new Date() ;
		leaseDate        = new Date() ;
		helperurl        = "";	
		local=true;
		properties = new HashMap<String, String>();
		
		
	}


	/** do a quick comparison if URL and nickname are the same
	 * 
	 */
	public boolean equals(DasSource other){
		//System.out.println("Das1Source equals, comparing with other DasSource");
		if (! (other instanceof Das1Source))
			return false;

		Das1Source ods = (Das1Source) other;

		if ( ods.getUrl().equals(url))
			return true;

		if ( ods.getNickname().equals(nickname))
			return true;

		return false;
	}

	/** makes a precise comparison between two das sources
	 * also compares description, coordsys, caps
	 * @param other
	 * @return if the two are exactly the same
	 */
	public boolean equalsExact(DasSource other){
		if ( ! this.equals(other))
			return false;

		if ( ! description.equals(other.getDescription()))
			return false;

		// test coordinate systems
		List<DasCoordinateSystem> newCoords = this.getCoordinateSystem();
		List<DasCoordinateSystem> oldCoords = other.getCoordinateSystem();

		if ( ! (newCoords.size() == oldCoords.size()))
			return false;

		//System.out.println("testing coords");
		for ( int i=0 ; i< newCoords.size() ; i++) {

			DasCoordinateSystem ncs = newCoords.get(i);

			boolean found = false ;

			for ( int j = 0 ; j < oldCoords.size(); j++){
				DasCoordinateSystem ocs = oldCoords.get(j);

				if ( ncs.getName().equals(ocs.getName())) {
					if ( ncs.getCategory().equals(ocs.getCategory())) {
						found = true;
						break;
					}
				}
				if (! found )
					return false;
			}
		}

		//System.out.println("testing capabs");

		// test capabilities
		List<Capabilities> otherCaps = other.getCapabilities();
		for (int i=0 ; i < capabilities.size();i++){
			Capabilities cap = capabilities.get(i);

			boolean found = false;
			for (int j=0; j < otherCaps.size();j++){
				Capabilities oCap = otherCaps.get(j);
				if ( oCap.equals(cap)) {
					found = true;
					break;
				}
			}
			if (! found)
				return false;
		}
		return true;
	}

	public int hashCode() {
		int h = 7;

		h = 31 * h + ( null == nickname ? 0 : nickname.hashCode());
		h = 31 * h + ( null == url ? 0 : url.hashCode());
System.out.println("hashcode="+h);
		return h;
	}

	/** the DAS2 string representation of this DAS source
	 * 
	 */
	public String toString() {

		StringWriter writer = new StringWriter();

		PrintWriter pw = new PrintWriter(writer);
		PrettyXMLWriter xw = new PrettyXMLWriter(pw);

		DasSourceWriter dswriter = new DasSourceWriterImpl();
		try {
			dswriter.writeDasSource(xw,this);
		} catch (IOException e){
			e.printStackTrace();
		}

		return writer.toString();

	}
	public void setLocal(boolean flag){ local = flag;}
	public boolean isLocal(){return local;}

	public void setId(String i) { id = i; }

	/** get a the Id of the DasSource. The Id is a unique db
	 * identifier. The public DAS-Registry has Auto_Ids that look like
	 * DASSOURCE:12345; public look like XYZ:12345, where the XYZ
	 * prefix can be configured in the config file.
	 */
	public String getId() { return id;}

	public void setNickname(String name) {
		nickname = name ;
	}
	public String getNickname(){
		return nickname;
	}
	public void setUrl(String u) {
		char lastChar = u.charAt(u.length()-1);
		if ( lastChar  != '/')
			u += "/";

		url = u ;
	}

	public void setAdminemail (String u) {
		adminemail = u ;
	}

	public void setDescription (String u) {
		description = u;
	}

	public void setCoordinateSystem (List<DasCoordinateSystem> u){
		coordinateSystem=u ;
	}

	public void setCapabilities (List<Capabilities> u){
		capabilities = u ;
	}

    /** test if a this source has a particular capability
     * 
     * @param testCapability
     * @return <code>true</code> if the server has this capability.
     */
    public boolean hasCapability(String testCapability){
        for (int i=0 ; i< capabilities.size();i++){
            Capabilities cap = capabilities.get(i);
            if ( cap.equals(testCapability))
                return true;
        }
        return false;
    }
    
	public String getUrl(){return url;}
    
	public String getAdminemail(){return adminemail;}
	
    public String getDescription(){return description;}
	
    public List<Capabilities> getCapabilities(){return capabilities;}
	
    public List<DasCoordinateSystem> getCoordinateSystem(){return coordinateSystem;}

	public void setRegisterDate(Date d) {
		registerDate = d;
	}
	public Date getRegisterDate() {
		return registerDate ;
	}
	public void setLeaseDate(Date d) {
		leaseDate =d ;
	}
	public Date getLeaseDate() {
		return leaseDate ;
	}

	public void setLabels(List<String> ls) {
		labels = ls ;
	}

	public List<String> getLabels() {
		return labels;
	}

	public void setHelperurl(String url) {
		helperurl = url;
	}

	public String getHelperurl() {
		return helperurl;
	}

	public void setAlertAdmin(boolean flag) {
		alertAdmin = flag;
	}

	public boolean getAlertAdmin() {
		return alertAdmin;
	}
	public Map<String,String> getProperties() {
		return properties;
	}
	public void setProperties(Map<String,String> properties) {
		this.properties = properties;
	}


	
	public List <String> getValidCapabilities(){
		
//			String []standin={"features"};
//			return  standin;
		
		return validCapabilities;
	}


	public void setValidCapabilities(List validCapabilities) {
		this.validCapabilities=validCapabilities;
		
	}
	

}
