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
 * Created on Jun 30, 2008
 * 
 */

package org.biojava.ontology;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Set;

import org.biojava.ontology.io.OboParser;

import junit.framework.TestCase;

public class ParseOBOFileTest extends TestCase{

	public void testsOBOFileParsing(){


		OboParser parser = new OboParser();
		InputStream inStream = this.getClass().getResourceAsStream("/org/biojava/bio/ontology/biosapiens.obo");

		assertNotNull("could not find biosapiens.obo file",inStream);
		BufferedReader oboFile = new BufferedReader ( new InputStreamReader ( inStream ) );
		try {
			Ontology ontology = parser.parseOBO(oboFile, "BioSapiens", "the BioSapiens ontology");

			Set keys = ontology.getTerms();
			assertTrue("did not find the expected nr of Ontology Terms!" , keys.size() > 50 );
			String soID = "SO:0001081";
			assertTrue("does not contain Term ", ontology.containsTerm(soID));
			Term hth = ontology.getTerm(soID);
			String soDesc = "helix_turn_helix";
			assertTrue("description does not match expectation ("+soDesc + ")", hth.getDescription().equals(soDesc));
		}
		catch (Exception e){
			fail (e.getMessage());
		}

	}
}
