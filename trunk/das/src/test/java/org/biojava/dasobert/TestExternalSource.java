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
 * Created on Sep 19, 2007
 * 
 */

package org.biojava.dasobert;

import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.List;

import org.biojava.dasobert.das2.Das2Source;
import org.biojava.dasobert.das2.io.DasSourceReaderImpl;
import org.biojava.dasobert.dasregistry.DasCoordinateSystem;
import org.biojava.dasobert.dasregistry.DasSource;

import junit.framework.TestCase;


/** a test case for parsing the SOURCES response
 * 
 * @author Andreas Prlic
 *
 */
public class TestExternalSource extends TestCase {
	
	String sourcesCmd = "http://vega.sanger.ac.uk/das/sources";
	
	public void testExternalSource() throws Exception{
		URL url = new URL(sourcesCmd);
		URLConnection conn = url.openConnection();
		InputStream inStream = conn.getInputStream();
		assertNotNull(inStream);

		DasSource[] sources = null;

		boolean errorHappend = false;
		try {

			DasSourceReaderImpl reader = new DasSourceReaderImpl();
			
			sources = reader.readDasSource(inStream);


		} catch (Exception ex) {
			ex.printStackTrace();
			errorHappend = true;
		}

		assertNotNull(sources);
		
		assertTrue(sources.length > 30);
		
		for (int i=0; i < sources.length;i++){
			DasSource s = sources[i];
			boolean isDAS2Source = false;
			if (  s instanceof Das2Source){
				isDAS2Source = true;
			}
			
			assertFalse("Found a DAS/2 source while parsing a DAS/1 sources listing",  isDAS2Source);
			
			
			List<DasCoordinateSystem> dcsses = s.getCoordinateSystem();
			assertTrue(dcsses != null);
			assertTrue(dcsses.size() > 0);
		}
		
		
		
		assertTrue(errorHappend != true);
		
	}
	
	
}
