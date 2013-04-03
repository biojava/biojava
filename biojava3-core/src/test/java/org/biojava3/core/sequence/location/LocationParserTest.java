package org.biojava3.core.sequence.location;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.template.Location;
import org.junit.Assert;
import org.junit.Test;

public class LocationParserTest {

    public static final InsdcParser PARSER = new InsdcParser();

    @Test
    public void basicLocationTests() {
        assertInsdcLoc("1", new SimpleLocation(1, 1, Strand.POSITIVE));

        assertInsdcLoc("1..10", new SimpleLocation(1, 10, Strand.POSITIVE));
        assertInsdcLoc("1^2", new SimpleLocation(
                new SimplePoint(1),
                new SimplePoint(2),
                Strand.POSITIVE, false, true));

        assertInsdcLoc("complement(1..10)", new SimpleLocation(1, 10, Strand.NEGATIVE));

        assertInsdcLoc("join(1..2,7..8)", new SimpleLocation(
                new SimplePoint(1), new SimplePoint(8), Strand.POSITIVE,
                new SimpleLocation(1, 2, Strand.POSITIVE),
                new SimpleLocation(7, 8, Strand.POSITIVE)));

        assertInsdcLoc("complement(join(1..2,7..8))", new SimpleLocation(
                new SimplePoint(1), new SimplePoint(8), Strand.NEGATIVE,
                new SimpleLocation(1, 2, Strand.NEGATIVE),
                new SimpleLocation(7, 8, Strand.NEGATIVE)));

        //Reverse relationship
        assertInsdcLoc("join(complement(1..2),complement(7..8))", new SimpleLocation(
                new SimplePoint(1), new SimplePoint(8), Strand.NEGATIVE,
                new SimpleLocation(1, 2, Strand.NEGATIVE),
                new SimpleLocation(7, 8, Strand.NEGATIVE)));

        //Complex sub relations
        assertInsdcLoc("join(1..2,join(4..5,complement(6..8))", new SimpleLocation(
                new SimplePoint(1), new SimplePoint(8), Strand.UNDEFINED,
                new SimpleLocation(1, 2, Strand.POSITIVE),
                new SimpleLocation(4, 8, Strand.UNDEFINED,
                new SimpleLocation(4, 5, Strand.POSITIVE),
                new SimpleLocation(6, 8, Strand.NEGATIVE))));

        assertInsdcLoc("join(5..10,1..3)", new SimpleLocation(
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
        Assert.assertEquals("Asserting locations are the same", expected, actual);
    }
}
