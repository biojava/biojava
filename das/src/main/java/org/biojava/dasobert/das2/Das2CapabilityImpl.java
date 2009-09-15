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
 * Created on Feb 9, 2006
 *
 */
package org.biojava.dasobert.das2;

public class Das2CapabilityImpl 
implements Das2Capability{

    String capability;
    String[] formats;
    String queryId;
    
    public static String DAS1_CAPABILITY_PREFIX = "das1:";
    
    public Das2CapabilityImpl() {
        super();
        capability = "undef";
        queryId = "";
        formats = new String[0];

    }
    
    public boolean isDas1Style(){
        
        if ( capability == null)
            return false;
        if ( capability.length() < DAS1_CAPABILITY_PREFIX.length())
            return false;
        if ( capability.substring(0,DAS1_CAPABILITY_PREFIX.length()).equals(DAS1_CAPABILITY_PREFIX))
            return true;
        return false;
    
    }
    
    public boolean equals(Das2Capability other){
        
        boolean status = true;
        
        if (!  capability.equals(other.getCapability()))
            status = false;
        if ( ! queryId.equals(other.getQueryUri()))
            status = false;
        
        return status;
    }
    
    public int hashCode(){
        int h = 7;
        h = 31 * h + ( null == capability ? 0 : capability.hashCode()) ;
        h = 31 * h + ( null == queryId    ? 0 : queryId.hashCode()) ;
        
        return h;
    }
    
    public String toString(){
        String txt ="capability " + capability + " queryId " + queryId;
        return txt;
    }

    public String getCapability() {
        
        return capability;
    }

    public String[] getFormats() {      
        return formats;
    }

    public String getQueryUri() {      
        return queryId;
    }

    public void setCapability(String type) {
       capability = type;
        
    }

    public void setFormats(String[] formats) {
       
        this.formats = formats; 
    }

    public void setQueryUri(String id) {
       queryId = id;
        
    }
    
    

}
