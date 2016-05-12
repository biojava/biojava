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

import java.util.Arrays;
import java.util.Collection;
import org.biojava.nbio.core.sequence.DataSource;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author Jacek Grzebyta
 */
@RunWith(Parameterized.class)
public class InsdcParserTest {

	private Logger log = LoggerFactory.getLogger(getClass());
	private String data;
	private String expected;

	public InsdcParserTest(String data, String expected) {
		this.data = data;
		this.expected = expected;
	}

	@Parameterized.Parameters
	public static Collection<String[]> data() {
		return Arrays.asList(new String[][]{
			{"complement(CP001663.1:6463934..6465826)", "CP001663.1"},
			{"complement(NC_000932.1:69611..69724)", "NC_000932.1"}
		});
	}

	/**
	 * Test for issue #254
	 *
	 * @throws Exception
	 */
	@Test
	public void extractAccessionTest() throws Exception {
		log.info("test accession");
		log.debug("data: '{}'   expected: '{}'", data, expected);

		InsdcParser parser = new InsdcParser(DataSource.GENBANK);
		Location loc = parser.parse(data);

		if (!loc.isComplex()) {
			log.info("simple location: {}", data);
			log.debug("\taccession: '{}'  expected: '{}'", loc.getAccession().getID(), expected);
			Assert.assertEquals(expected, loc.getAccession().getID());
		}
	}
}
