/*
 *                    PDB web development code
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
 *
 * Created on Jul 21, 2009
 * Created by ap3
 *
 */

package org.biojava.bio.structure.align.webstart;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.logging.Level;
import java.util.logging.Logger;


public class BrowserOpener {
    public static Logger logger =  Logger.getLogger("org.biojava.spice");
   
    /** open a URL in the browser that was used to launch SPICE
     * 
     * @param url URL to be opened
     * @return true if this was successfull
     */
    public static boolean showDocument(URL url) 
    {
        if ( url != null ){
            boolean success = JNLPProxy.showDocument(url); 
            if ( ! success)
                logger.info("could not open URL "+url+" in browser. check your config or browser version.");
        return success;
        
        }
        else
            return false;
    }
    

    /** open a URL in the browser that was used to launch SPICE
     * 
     * @param urlstring string represntation of URL to be opened
     * @return true if this was successfull
     */
    public static boolean showDocument(String urlstring){
        try{
            URL url = new URL(urlstring);
            
            return showDocument(url);
        } catch (MalformedURLException e){
            logger.log(Level.WARNING,"malformed URL "+urlstring, e.getCause());
            return false;
        }
    }

}

