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

package org.biojava3.ontology;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.util.Set;

import junit.framework.TestCase;

import org.biojava3.ontology.io.OboParser;

public class TestOboFileParsing extends TestCase{


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
			e.printStackTrace();
			fail(e.getMessage());
		}

		
		
		
	}

}
