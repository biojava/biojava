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
 * Created on Oct 28, 2005
 *
 */
package org.biojava.dasobert.eventmodel;

import java.util.List;
import java.util.Map;

import org.biojava.dasobert.dasregistry.Das1Source;
public class FeatureEvent 
extends AbstractDasEvent {
	
    List<Map<String, String>> features;
   
    int comeBackLater;
    String version;
    
    public FeatureEvent(List<Map<String, String>> feats,Das1Source dasSource,String version) {
        super();
        this.features =feats;
        this.dasSource = dasSource;
        comeBackLater = -1;
        this.version = version;
    }
    
    public int getComeBackLater(){
        return comeBackLater;
    }
    
    public void setComeBackLater(int comeBackLater){
        this.comeBackLater = comeBackLater;
    }
    
    
    /** get the features that have been found.
     * 
     * do something like
     * Map[] features = event.getFeatures();
     * <pre>
     * for (int i = 0 ; i< features;i++) {            
     *      Map f = features[i];
     *      String type = (String) f.get("TYPE") ;
     *      System.out.println(type);
     * }      
     * </pre>
     * @return a Map containng the features
     */
    public List<Map<String,String>> getFeatures(){
        return features;
    }

    /** Get the version of the reference object that has been annotated.
     * Compare the version string with the version string obtained from the reference server.
     * If they don;t match there is a version problem between the annotation and the reference!
     * @return the version string (e.g. an MD5 digest of the reference sequence)
     */
	public String getVersion() {
		return version;
	}

	public void setVersion(String version) {
		this.version = version;
	}
    
   
    
}

