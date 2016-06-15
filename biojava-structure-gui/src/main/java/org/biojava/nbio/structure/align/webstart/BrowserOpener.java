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

package org.biojava.nbio.structure.align.webstart;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.net.MalformedURLException;
import java.net.URL;


public class BrowserOpener {

	private static final Logger logger = LoggerFactory.getLogger(BrowserOpener.class);

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
			logger.warn("malformed URL {}", urlstring, e);
			return false;
		}
	}

}

