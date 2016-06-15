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

import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * TODO: Temporary test is switched off. Currently results are messy:
 *
 * <code>
 * to test: complement(order(1,2..34,complement(34..45),A00001.5:34..45))	expected: 1..45(.)	received: 1..45(.)
 * to test: 1	expected: 1..1(+)	received: 1..1(+)
 * to test: 1..10	expected: 1..10(+)	received: 1..10(+)
 * to test: 1^2	expected: 1^2(+)	received: 1^2(+)
 * to test: complement(1..10)	expected: 1..10(-)	received: 1..10(-)
 * to test: join(1..2,7..8)	expected: 1..8(+)	received: 1..8(+)
 * to test: complement(join(1..2,7..8))	expected: 1..8(-)	received: 1..8(-)
 * to test: join(complement(1..2),complement(7..8))	expected: 1..8(-)	received: 1..8(-)
 * to test: join(1..2,join(4..5,complement(6..8))	expected: 1..8(.)	received: 1..8(.)
 * to test: join(5..10,1..3)	expected: 5..13(+ - circular)	received: 1..10(+)
 * </code>
 *
 * Serialisation to string should be fixed as well.
 *
 *
 * @author Jacek Grzebyta
 */
public class LocationParserTest {

	public static final InsdcParser PARSER = new InsdcParser();

	private Logger log = LoggerFactory.getLogger(getClass());

	@Test
	@Ignore
	public void basicLocationTests() {
		assertInsdcLoc("1", new SimpleLocation(1, 1, Strand.POSITIVE));

		assertInsdcLoc("1..10", new SimpleLocation(1, 10, Strand.POSITIVE));
		assertInsdcLoc("1^2", new SimpleLocation(
				new SimplePoint(1),
				new SimplePoint(2),
				Strand.POSITIVE, false, true));

		assertInsdcLoc("complement(1..10)", new SimpleLocation(1, 10, Strand.NEGATIVE));

		assertInsdcLoc("join(1..2,7..8)", new InsdcLocations.GroupLocation(
				new SimplePoint(1), new SimplePoint(8), Strand.POSITIVE,
				new SimpleLocation(1, 2, Strand.POSITIVE),
				new SimpleLocation(7, 8, Strand.POSITIVE)));

		assertInsdcLoc("complement(join(1..2,7..8))", new InsdcLocations.GroupLocation(
				new SimplePoint(1), new SimplePoint(8), Strand.NEGATIVE,
				new SimpleLocation(1, 2, Strand.NEGATIVE),
				new SimpleLocation(7, 8, Strand.NEGATIVE)));

		//Reverse relationship
		assertInsdcLoc("join(complement(1..2),complement(7..8))", new InsdcLocations.GroupLocation(
				new SimplePoint(1), new SimplePoint(8), Strand.NEGATIVE,
				new SimpleLocation(1, 2, Strand.NEGATIVE),
				new SimpleLocation(7, 8, Strand.NEGATIVE)));

		//Complex sub relations
		//should tests be designed for both modes?
		//PARSER.setComplexFeaturesAppendMode(InsdcParser.complexFeaturesAppendEnum.HIERARCHICAL);
		assertInsdcLoc("join(1..2,join(4..5,complement(6..8))", new InsdcLocations.GroupLocation(
				new SimplePoint(1), new SimplePoint(8), Strand.UNDEFINED,
				new SimpleLocation(1, 2, Strand.POSITIVE),
				new SimpleLocation(4, 8, Strand.UNDEFINED,
						new SimpleLocation(4, 5, Strand.POSITIVE),
						new SimpleLocation(6, 8, Strand.NEGATIVE))));

		assertInsdcLoc("join(5..10,1..3)", new InsdcLocations.GroupLocation(
				new SimplePoint(5), new SimplePoint(13), Strand.POSITIVE,
				true, //Circular genome
				new SimpleLocation(5, 10, Strand.POSITIVE),
				new SimpleLocation(1, 3, Strand.POSITIVE)));

		assertInsdcLoc("order(1..2,7..8)", new InsdcLocations.OrderLocation(
				new SimplePoint(1), new SimplePoint(8), Strand.POSITIVE,
				new SimpleLocation(1, 2, Strand.POSITIVE),
				new SimpleLocation(7, 8, Strand.POSITIVE)));
	}

	@Test
	@Ignore
	public void moreComplex() {
		assertInsdcLoc("complement(order(1,2..34,complement(34..45),A00001.5:34..45))",
				new InsdcLocations.OrderLocation(
						new SimplePoint(1), new SimplePoint(45), Strand.UNDEFINED,
						new SimpleLocation(1, 1, Strand.NEGATIVE),
						new SimpleLocation(2, 34, Strand.NEGATIVE),
						new SimpleLocation(34, 45, Strand.POSITIVE),
						new SimpleLocation(
								new SimplePoint(34), new SimplePoint(45),
								Strand.NEGATIVE,
								new AccessionID("A00001.5", PARSER.getDataSource()))));
	}

	public void assertInsdcLoc(String stringLoc, Location expected) {
		Location actual = PARSER.parse(stringLoc);
		log.info("to test: {}\texpected: {}\treceived: {}", stringLoc, expected.toString(), actual.toString());
		Assert.assertEquals("Asserting locations are the same", expected.toString(), actual.toString());
	}
}
