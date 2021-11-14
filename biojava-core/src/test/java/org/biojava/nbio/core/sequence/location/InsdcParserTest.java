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
 */
package org.biojava.nbio.core.sequence.location;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.util.Arrays;
import java.util.Collection;
import org.biojava.nbio.core.sequence.DataSource;
import org.biojava.nbio.core.sequence.location.InsdcParser.complexFeaturesAppendEnum;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author Jacek Grzebyta
 */

public class InsdcParserTest {

	private Logger log = LoggerFactory.getLogger(getClass());

	static Collection<String[]> data() {
		return Arrays.asList(new String[][] {

		});
	}

	/**
	 * Test for issue #254
	 *
	 * @throws Exception
	 */
	@ParameterizedTest
	@CsvSource({ "complement(CP001663.1:6463934..6465826),CP001663.1",
			"complement(NC_000932.1:69611..69724),NC_000932.1" })
	public void extractAccessionTest(String data, String expected) throws Exception {

		InsdcParser parser = new InsdcParser(DataSource.GENBANK);
		Location loc = parser.parse(data);

		if (!loc.isComplex()) {
			log.info("simple location: {}", data);
			log.debug("\taccession: '{}'  expected: '{}'", loc.getAccession().getID(), expected);
			assertEquals(expected, loc.getAccession().getID());
		}
	}

	@Test
	public void testParser() {
		String[] testStrings = { "J00194.1:100..202", "A00001.5:34..45", "43..129", "bond(55,110)",
				"bond(34,35),join(56..80),complement(45,73)",
				"order(complement(30,40),70..80),bond(34,35),join(56,80),complement(45..56)",
				"join(join(complement(30,40),complement(70..80)),bond(34,35),join(56,80),complement(45..56))",
				"complement(join(complement(2000..4000),complement(70..80)),bond(34,35),join(56,80),complement(45..56))",

		};
		InsdcParser p = new InsdcParser();
		p.setComplexFeaturesAppendMode(complexFeaturesAppendEnum.HIERARCHICAL);

		for (String s : testStrings) {
			Location l = p.parse(s);
			assertNotNull(l);
		}
	}
}
