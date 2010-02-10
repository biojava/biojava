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
 * Created on Oct 26, 2006
 *
 */
package org.biojava.dasobert.util;


public class DasobertDefaults {

   

    public static final String VERSION = "0.3";
    
    public static final String PROJECTNAME = "dasobert - Java ";
    
    public static final String LOGGERNAME = "org.biojava.spice";
        
    public static final String USERAGENT;
    
    static {
        //e.g. "Mozilla/5.0 (Windows; U; Win98; en-US; rv:1.7.2) Gecko/20040803"
        String os_name    = java.lang.System.getProperty("os.name");
        String os_version = java.lang.System.getProperty("os.version");
        String os_arch    = java.lang.System.getProperty("os.arch");    
        
        USERAGENT = PROJECTNAME+ " " + VERSION + "("+os_name+"; "+os_arch + " ; "+ os_version+")";
    }
    
    
}
