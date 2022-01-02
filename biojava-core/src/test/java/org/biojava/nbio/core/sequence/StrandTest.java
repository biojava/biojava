package org.biojava.nbio.core.sequence;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class StrandTest {

    @Test
    void reverse(){
        assertEquals(Strand.POSITIVE, Strand.NEGATIVE.getReverse());
        assertEquals(Strand.NEGATIVE, Strand.POSITIVE.getReverse());
        assertEquals(Strand.UNDEFINED, Strand.UNDEFINED.getReverse());
    }

    @Test
    void stringRepresentation() {
        assertEquals("+", Strand.POSITIVE.getStringRepresentation());
        assertEquals("-", Strand.NEGATIVE.getStringRepresentation());
        assertEquals(".", Strand.UNDEFINED.getStringRepresentation());
    }

    @Test
    void numberRepresentation() {
        assertEquals(1, Strand.POSITIVE.getNumericRepresentation());
        assertEquals(-1, Strand.NEGATIVE.getNumericRepresentation());
        assertEquals(0, Strand.UNDEFINED.getNumericRepresentation());
    }

}