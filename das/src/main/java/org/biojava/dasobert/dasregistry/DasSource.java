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
 * Created on Feb 8, 2006
 *
 */
package org.biojava.dasobert.dasregistry;

import java.util.Date;
import java.util.List;
import java.util.Map;

import org.biojava.dasobert.das.Capabilities;

public interface DasSource {

    public  void setLocal(boolean flag);

    public  boolean isLocal();

    /** compare if two das sources are equal
     * 
     * @param ds
     * @return returns true if two DAS sources are equivalent
     */
    public boolean equals(DasSource ds);
    
    /** classes that implement equals, should also implement hashKey
     * 
     * @return the hash code of a das source
     */
    public int hashCode();
    
    
    public  void setId(String i);

    /** get a the Id of the DasSource. The Id is a unique db
     * identifier. The public DAS-Registry has Auto_Ids that look like
     * DASSOURCE:12345; public look like XYZ:12345, where the XYZ
     * prefix can be configured in the config file.
     * @return String the ID of a Das Source
     */
    public  String getId();

    public  void setNickname(String name);

    public  String getNickname();

    public  void setUrl(String u);

    public  void setAdminemail(String u);

    public  void setDescription(String u);

    public  void setCoordinateSystem(List<DasCoordinateSystem> u);

    public  void setCapabilities(List<Capabilities> u);

    /** test if a this source has a particular capability
     * 
     * @param testCapability
     * @return <code>true</code> if the server has this capability.
     */
    public boolean hasCapability(String testCapability);
    
    public  String getUrl();

    public  String getAdminemail();

    public  String getDescription();

    public  List<Capabilities> getCapabilities();

    public  List<DasCoordinateSystem> getCoordinateSystem();

    public  void setRegisterDate(Date d);

    public  Date getRegisterDate();

    public  void setLeaseDate(Date d);

    public  Date getLeaseDate();

    public  void setLabels(List<String> ls);

    public  List<String> getLabels();

    public  void setHelperurl(String url);

    public  String getHelperurl();

    // TestCode is now part of the coordinate system!
    //public  void setTestCode(String code);
    //public  String getTestCode();

    public  void setAlertAdmin(boolean flag);

    public  boolean getAlertAdmin();
    
    /** set Properties for this DAS source, e.g. project name
     * 
     * @param properties
     */
    public void setProperties(Map<String,String> properties);
    
    /** get Properties for this DAS source
     * 
     * @return Properties
     */
    public Map<String,String> getProperties();

	public void setValidCapabilities(String[] validCapabilities);
	
	public String[] getValidCapabilities();

}