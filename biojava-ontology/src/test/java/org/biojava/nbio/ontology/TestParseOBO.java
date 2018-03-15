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
import org.junit.Assert;
import org.junit.Test;

import java.io.*;
import java.text.ParseException;
import java.util.Set;

import static org.biojava.nbio.ontology.obo.OboFileHandler.NAMESPACE;

public class TestParseOBO {

	@Test
	public void testNamespace() throws IOException, ParseException {

		String testTermEntry = "\n[Term]\n" +
                "id: SO:0000691\n" +
                "name: cleaved_initiator_methionine \n" +
                "namespace: sequence\n" +
                "alt_id: BS:00067\n" +
                "def: \"The initiator methionine that has been cleaved from a mature polypeptide sequence.\" [EBIBS:GAR]\n" +
                "subset: biosapiens\n" +
                "synonym: \"cleaved initiator methionine\" EXACT []\n" +
                "synonym: \"init_met\" RELATED [uniprot:feature_type]\n" +
                "synonym: \"initiator methionine\" RELATED []\n" +
                "is_a: SO:0100011 ! cleaved_peptide_region\n\n";

        OboParser parser = new OboParser();
        InputStream inStream = new ByteArrayInputStream(testTermEntry.getBytes());

		Assert.assertNotNull(inStream);

		BufferedReader oboFile = new BufferedReader ( new InputStreamReader ( inStream ) );
		Ontology ontology = parser.parseOBO(oboFile, "so-xp/subsets/biosapiens",
                    "snippet from biosapiens protein feature ontology");
		Set<Term> keys = ontology.getTerms();

		Assert.assertTrue(keys.size() > 1);
		Assert.assertTrue(ontology.getTerm("SO:0000691").getAnnotation().containsProperty(NAMESPACE));
        Assert.assertEquals("sequence", ontology.getTerm("SO:0000691").getAnnotation().getProperty(NAMESPACE));
	}
}
