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
 * Created on Mar 19, 2014
 * Author: andreas
 *
 */

package org.biojava.nbio.ontology;

import junit.framework.TestCase;
import org.biojava.nbio.ontology.io.OboParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Set;

public class TestOboFileParsing extends TestCase{

	private static final Logger logger = LoggerFactory.getLogger(TestOboFileParsing.class);

	public void testParsingBioSapiensOBO(){
		OboParser parser = new OboParser();
		InputStream inStream = parser.getClass().getResourceAsStream("/ontology/biosapiens.obo");

		assertNotNull(inStream);

		BufferedReader oboFile = new BufferedReader ( new InputStreamReader ( inStream ) );

		Ontology ontology;
		try {
			ontology = parser.parseOBO(oboFile, "BioSapiens", "the BioSapiens ontology");
			Set<Term> keys = ontology.getTerms();

			assertTrue(keys.size() >4000);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			logger.error("Exception: ", e);
			fail(e.getMessage());
		}




	}

}
