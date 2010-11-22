package org.biojava.bio.structure.align.util;


public class ConfigurationException extends Exception {
    
    /**
	 * 
	 */
	private static final long serialVersionUID = -8047100079715000276L;

	/**
     * Constructs a ConfigurationException object.
     *
     * @param s  a String ...
     */
    public ConfigurationException(String s) {
	super(s);
    }
    
    /**
     * Constructs a ConfigurationException object.
     *
     * @param t  a Throwable object
     */
    public ConfigurationException (Throwable t) {
	super(t);
    }
}

