package org.biojava3.core.sequence.location;

import static org.biojava3.core.sequence.Strand.UNDEFINED;
import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;
import org.biojava3.core.sequence.DNASequence;

import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.template.Sequence;
import org.junit.Test;

public class LocationTest {

    @Test
    public void testSubLocations() {
        List<SimpleLocation> expected = Arrays.asList(
                new SimpleLocation(1, 2, Strand.UNDEFINED),
                new SimpleLocation(3, 6, Strand.UNDEFINED),
                new SimpleLocation(7, 10, Strand.UNDEFINED),
                new SimpleLocation(11, 20, UNDEFINED));

        Location location = new SimpleLocation(1, 20, UNDEFINED,
                new SimpleLocation(1, 10, UNDEFINED,
                new SimpleLocation(1, 6, UNDEFINED,
                new SimpleLocation(1, 2, UNDEFINED), new SimpleLocation(3, 6, UNDEFINED)),
                new SimpleLocation(7, 10, UNDEFINED)),
                new SimpleLocation(11, 20, UNDEFINED));

        List<Location> actual = location.getRelevantSubLocations();

        assertEquals("Checking sublocations iterate as expected", toStr(expected), toStr(actual));

        String expectedDna = "AACCCCTTTTGGGGGGGGGG";
        Sequence<NucleotideCompound> s = new DNASequence(expectedDna + "CC");
        Sequence<NucleotideCompound> subSeq = location.getSubSequence(s);
        assertEquals("Checking subseq as expected", expectedDna, subSeq.getSequenceAsString());
    }

    @Test
    public void testBasicCircularLocation() {
        Location circularLocation = new SimpleLocation(
                new SimplePoint(3), new SimplePoint(52), Strand.POSITIVE, true,
                new SimpleLocation(3, 20, Strand.POSITIVE),
                new SimpleLocation(1, 20, Strand.POSITIVE),
                new SimpleLocation(1, 12, Strand.POSITIVE));
        assertEquals("Checking length as expected", circularLocation.getLength(), 50);

        String expectedDna = "CCCCTTTTGGGGGGGGGGAACCCCTTTTGGGGGGGGGGAACCCCTTTTGG";
        DNASequence s = new DNASequence("AACCCCTTTTGGGGGGGGGG");
        Sequence<NucleotideCompound> subSeq = circularLocation.getSubSequence(s);
        assertEquals("Checking subseq as expected", expectedDna, subSeq.getSequenceAsString());

    
        Location newCircularLocation = Location.Tools.circularLocation(
                3, 52, Strand.POSITIVE, 20);
        //Check this is the right set of coords even to use!
        Location negativeCoordsCircularLocation = Location.Tools.circularLocation(
                58,9, Strand.POSITIVE, 20);

        assertEquals("location objects should be equivalent", circularLocation, newCircularLocation);
        assertEquals("location objects should be equivalent even if they are on the wrong coord system", circularLocation.getSubLocations(), negativeCoordsCircularLocation.getSubLocations());
        assertEquals("Checking subseq as expected", expectedDna,
                newCircularLocation.getSubSequence(s).getSequenceAsString());
    }

    @Test
    public void testWithStrandSwitch() {
        DNASequence s = new DNASequence("AAAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGGAAATTTCCCG");
        Location location = new SimpleLocation(1, 50, Strand.UNDEFINED,
                new SimpleLocation(1, 10, Strand.POSITIVE), //AAAAAAAAAA
                new SimpleLocation(19, 23, Strand.NEGATIVE), //GGGAA
                new SimpleLocation(45, 50, Strand.POSITIVE)); //TTCCCG
        String expectedDna = "AAAAAAAAAAGGGAATTCCCG";
        Sequence<NucleotideCompound> subSeq = location.getRelevantSubSequence(s);
        assertEquals("Checking subseq as expected", expectedDna, subSeq.getSequenceAsString());
    }

    @Test
    public void testStrandFlip() {
        Location l = new SimpleLocation(3,17,Strand.POSITIVE);
        Location r = Location.Tools.location(18, 4, Strand.POSITIVE, 20);
        assertEquals("Locations should be the same even though they were expressed differently", l, r);
    }

    @Test(expected = IllegalStateException.class)
    public void badLocations() {
        new SimpleLocation(10, 1, Strand.UNDEFINED);
    }

    @Test
    public void modulateCircular() {
        int length = 20;
        assertEquals("Checking modulation", 12, Location.Tools.modulateCircularIndex(52, length));
        assertEquals("Checking modulation", 1, Location.Tools.modulateCircularIndex(21, length));
        assertEquals("Checking modulation", 3, Location.Tools.modulateCircularIndex(3, length));
    }

    @Test
    public void completePasses() {
        assertEquals("Checking passes", 4, Location.Tools.completeCircularPasses(52, 10));
        assertEquals("Checking passes", 2, Location.Tools.completeCircularPasses(36, 10));
    }

    private <L extends Location> String toStr(List<L> locations) {
        StringBuilder sb = new StringBuilder();
        for (L l : locations) {
            sb.append(l).append("|");
        }
        return sb.toString();
    }
}
