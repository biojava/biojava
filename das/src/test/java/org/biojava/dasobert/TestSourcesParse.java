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
 * Created on Jun 13, 2007
 * 
 */

package org.biojava.dasobert;

import java.io.InputStream;

import org.biojava.dasobert.das2.Das2Source;
import org.biojava.dasobert.das2.io.DasSourceReaderImpl;
import org.biojava.dasobert.dasregistry.DasSource;

import junit.framework.TestCase;

public class TestSourcesParse extends TestCase {

	public void testParseSourcesResponse() {
		
		InputStream inStream = this.getClass().getResourceAsStream("/sources.xml");

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
		
		assertEquals("did not find the expected number of DAS sources.", 267,sources.length);
		
		for (int i=0; i < sources.length;i++){
			DasSource s = sources[i];
			boolean isDAS2Source = false;
			if (  s instanceof Das2Source){
				isDAS2Source = true;
			}
			
			assertFalse("Found a DAS/2 source while parsing a DAS/1 sources listing",  isDAS2Source);
		}
		
		
		
		assertTrue(errorHappend != true);
	}
}
