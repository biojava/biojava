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
 * created at 20 Mar 2014
 * Author: ap3
 */

package org.biojava.nbio.ontology;

import org.biojava.nbio.ontology.io.OboParser;
import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Set;

public class TestParseSBO {

	private static final Logger logger = LoggerFactory.getLogger(TestParseSBO.class);

	@Test
	public void testSBO(){
		OboParser parser = new OboParser();


		try {

			InputStream inStream = parser.getClass().getResourceAsStream("/ontology/SBO_OBO.obo");

			Assert.assertNotNull(inStream);

			BufferedReader oboFile = new BufferedReader ( new InputStreamReader ( inStream ) );

			Ontology ontology;

			ontology = parser.parseOBO(oboFile, "SBO", "systems biology representation");
			Set<Term> keys = ontology.getTerms();

			Assert.assertTrue(keys.size() > 1200);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			logger.error("Exception: ", e);
			Assert.fail(e.getMessage());
		}
	}
}
