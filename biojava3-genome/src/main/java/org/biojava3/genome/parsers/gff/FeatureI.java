package org.biojava3.genome.parsers.gff;

import java.util.*;

/**
* A feature on a sequence (for example, an exon or a gene), defined by a location
* and a set of attributes encoded as key/value pairs.
*
* @author Hanno Hinsch
*/
public interface FeatureI
{
	
	/**
	 * Get the location of the feature.
	 *
	 * @return The location.
	 */
	public Location location();
	
	/**
	 * Get the group id of the feature. The group id is defined in the GFF1
	 * format; features that share a common group id are presumed to be part
	 * of some logical group. For example, a gene might be represented by several
	 * exon features that all share the same group id.
	 *<br><br> 
	 * The exact meaning of a feature's group id, or even its existence, is not guaranteed
	 * by this interface. An understanding of a particular file's data format is necessary to properly
	 * interpret the group id.
	 *
	 * @return The group id. This may be an empty string.
	 */
	public String group();
	
	/**
	 * Get the feature type, for example, "exon", "CDS", etc.
	 *
	 * @return The type.
	 */
	public String type();
	
    /**
     * Get the sequence name.  
     *
     * @return Sequence name.
     */
    public String seqname();
    
    /**
	 * Get the attribute value for this key.
	 *
	 * @param key The key.
	 * @return The corresponding value. Null if the key has no value defined .
	 */ 
	public String getAttribute( String key );
	
	
    /**
     * Check if the feature has a value defined for the specified key.
     *
     * @param key The key.
     * @return True if a value is defined for this key.
     */ 
    public boolean hasAttribute( String key );
	
    /**
     * Check if the feature attributes include the specified key/value pair.
     *
     * @param key The key.
     * @param value The value.
     * @return True if the feature's value for this key matches the specified value.
     */ 
    public boolean hasAttribute( String key, String value );
	
	/**
	 * A string representation of the feature.
	 * 
	 * @return The string.
	 */
	public String toString();
	
    /**
     * Get HashMap of user data.
     *
     * @return The user HashMap.
     */
    public  HashMap<String, String> userData();

	public HashMap<String, String> getAttributes();
    
}