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
 * Created at 08 Aug 2017
 */

package org.biojava.nbio.ontology;

import org.biojava.nbio.ontology.io.OboParser;
import org.biojava.nbio.ontology.utils.Annotation;
import org.junit.Before;
import org.junit.Test;

import java.io.*;
import java.text.ParseException;
import java.util.List;
import java.util.Set;

import static org.biojava.nbio.ontology.obo.OboFileHandler.NAMESPACE;
import static org.biojava.nbio.ontology.obo.OboFileHandler.ALT_ID;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

public class TestParseOBO {

	private OboParser parser;

	final String testTermEntry = "\n[Term]\n" + "id: SO:0000691\n" + "name: cleaved_initiator_methionine \n"
			+ "namespace: sequence\n" + "alt_id: BS:00067\n"
			+ "def: \"The initiator methionine that has been cleaved from a mature polypeptide sequence.\" [EBIBS:GAR]\n"
			+ "subset: biosapiens\n" + "synonym: \"cleaved initiator methionine\" EXACT []\n"
			+ "synonym: \"init_met\" RELATED [uniprot:feature_type]\n"
			+ "synonym: \"initiator methionine\" RELATED []\n" + "is_a: SO:0100011 ! cleaved_peptide_region\n\n";
	private String replace;

	public Ontology readObo(String input) throws ParseException, IOException {
		parser = new OboParser();
		InputStream inStream = new ByteArrayInputStream(input.getBytes());
		assertNotNull(inStream);
		BufferedReader oboFile = new  BufferedReader(new InputStreamReader(inStream));
		return parser.parseOBO(oboFile, "so-xp/subsets/biosapiens",
				"snippet from biosapiens protein feature ontology");
	}

	@Test
	public void testNamespace() throws IOException, ParseException {
		Ontology ontology = readObo(testTermEntry);

		Set<Term> keys = ontology.getTerms();

		assertTrue(keys.size() > 1);
		assertTrue(getAnnotationForTerm(ontology).containsProperty(NAMESPACE));
		assertEquals("sequence", getAnnotationForTerm(ontology).getProperty(NAMESPACE));
		assertTrue(getAnnotationForTerm(ontology).getProperty(ALT_ID) instanceof List);
	}

	@Test
	public void testMultipleAltIds() throws IOException, ParseException {
		
		String replace = testTermEntry.replace("BS:00067", "BS:00067\nalt_id: BS:00068");
		Ontology ontology = readObo(replace);
		List<String> altIds = (List<String>) getAnnotationForTerm(ontology).getProperty(ALT_ID);
		assertEquals(2, altIds.size());
		assertEquals("BS:00067", altIds.get(0));
		assertEquals("BS:00068", altIds.get(1));
	}

	private Annotation getAnnotationForTerm(Ontology ontology) {
		return ontology.getTerm("SO:0000691").getAnnotation();

	}
}
