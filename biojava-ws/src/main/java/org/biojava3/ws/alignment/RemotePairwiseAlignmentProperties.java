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
 */
package org.biojava3.ws.alignment;

import java.io.Serializable;
import java.util.Set;

/**
 * RemotePairwiseAlignmentProperties is a interface that contains the barest of 
 * methods for setting and getting Alignment properties.
 * 
 * Ideally, one would extend this class if creating a service by creating 
 * wrapper methods that actually call either getAlignementOption or setAlignementOption
 * with specific values for parameter names and checking values for options.
 * 
 * For an example, go see NCBIQBlastProperties
 * 
 * @author Sylvain Foisy, Diploide BioIT
 * @since Biojava 3
 *
 */
public interface RemotePairwiseAlignmentProperties extends Serializable{

	public static final long serialVersionUID = 1L;
	
    /**
     * Method that returns the value associated with the key given in parameter.
     * 
     * @param key :a String with the required key for this map.
     * @return a String with the value associated with this key
     * @throws Exception if key is not in the map of output options.
     */
	public String getAlignmentOption(String key) throws Exception;
	
	/**
	 * Method to set the value for a specific alignment parameter using a key to store in a map.
	 * 
	 * @param key :the key use to designate the value to be stored
	 * @param val :the actual value matched to key
	 */
	public void setAlignementOption(String key,String val);
	
	/**
	 * Method to get all keys to the information stored in this object.
	 * 
	 * @return a <code>Set</code> with all keys held in this instance of the object 
	 */
	public Set<String> getAlignmentOptions();
}