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

/**
 * @author ap3
 *
 */
public interface Das2Capability {

    /**
     * @param other
     * @return a flag if the two capabilities match
     */
    public boolean equals(Das2Capability other);
    public int hashCode();
    
    /**
     * @param type
     */
    public void setCapability(String type);
   
    /**
     * @return the capability
     */
    public String getCapability();
    
    /**
     * @param id
     */
    public void setQueryUri(String id);
    
    /**
     * @return the URI of the query
     */
    public String getQueryUri();
    
    /**
     * @param formats
     */
    public void setFormats(String[] formats);
    
    /**
     * @return an array of String[] that shows the supported file formats
     */
    public String[] getFormats();
    
    /** checks if this capability is actually of das1 style
     * 
     * @return boolean true if the capability is in DAS1 style
     */
    public boolean isDas1Style();
    
}
