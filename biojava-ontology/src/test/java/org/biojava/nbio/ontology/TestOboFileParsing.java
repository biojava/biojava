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

import org.biojava.nbio.ontology.io.OboParser;
import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Set;
import java.util.Iterator;
import org.biojava.nbio.ontology.utils.Annotation;

public class TestOboFileParsing {

	private static final Logger logger = LoggerFactory.getLogger(TestOboFileParsing.class);

	@Test
	public void testParsingBioSapiensOBO() throws Exception {
		OboParser parser = new OboParser();
		InputStream inStream = parser.getClass().getResourceAsStream("/ontology/biosapiens.obo");

		Assert.assertNotNull(inStream);

		BufferedReader oboFile = new BufferedReader ( new InputStreamReader ( inStream ) );

		Ontology ontology;

		ontology = parser.parseOBO(oboFile, "BioSapiens", "the BioSapiens ontology");
		Set<Term> keys = ontology.getTerms();

		Assert.assertTrue(keys.size() > 4000);
	}

	@Test
	public void testParsingHPOOBO() throws Exception {
		OboParser parser = new OboParser();
		InputStream inStream = parser.getClass().getResourceAsStream("/ontology/hp.obo");

		Assert.assertNotNull(inStream);

		BufferedReader oboFile = new BufferedReader ( new InputStreamReader ( inStream ) );

		Ontology ontology;

		ontology = parser.parseOBO(oboFile, "Human_phenotype", "the Human Phenotype ontology");
		Set<Term> keys = ontology.getTerms();
		Iterator iter = keys.iterator();

		while (iter.hasNext()){
			Term term = (Term) iter.next();
			if(term.getName().equals("HP:0000057")) {
				Annotation anno = term.getAnnotation();
				Assert.assertTrue(anno.containsProperty("replaced_by"));
				Assert.assertEquals("HP:0008665", anno.getProperty("replaced_by"));
			}
		}
	}

}
