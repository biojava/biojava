package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashSet;
import java.util.Set;

import org.junit.jupiter.api.RepeatedTest;
import org.junit.jupiter.api.Test;

class HashcoderTest {

    int seed = Hashcoder.SEED;
    @RepeatedTest(10)
    void hashcodeBool() {
        final int EXPECTED_TRUE = 712;
        final int EXPECTED_FALSE = 711;

        assertEquals(EXPECTED_TRUE, Hashcoder.hash(seed, true));
        assertFalse(EXPECTED_TRUE == Hashcoder.hash(seed, Boolean.TRUE));

        assertEquals(EXPECTED_FALSE, Hashcoder.hash(seed, false));
        assertFalse(EXPECTED_FALSE == Hashcoder.hash(seed, Boolean.FALSE));
    }

    @Test
    void hashcodeIntProducesUniqueValues() {
        Set<Integer> hashCodes = new HashSet<>();
        for (int i =0; i< 1000; i++) {
            hashCodes.add(Hashcoder.hash(seed, i));
        }
        assertEquals(1000, hashCodes.size());
    }

    @Test
    void extremeIntValuesHashed(){
        assertTrue(Hashcoder.hash(seed, Integer.MAX_VALUE) != 0);
        assertTrue(Hashcoder.hash(seed, Integer.MIN_VALUE) != 0);
    }

    @Test
    void hashcodeLongProducesUniqueValues() {
        Set<Integer> hashCodes = new HashSet<>();
        for (long i =0; i< 1000; i++) {
            hashCodes.add(Hashcoder.hash(seed, i));
        }
        assertEquals(1000, hashCodes.size());
      
    }

    @Test
    void extremeLongValuesHashed(){
        assertTrue(Hashcoder.hash(seed, Long.MAX_VALUE) != 0);
        assertTrue(Hashcoder.hash(seed, Long.MIN_VALUE) != 0);
    }

    @Test
    void hashcodeCharProducesUniqueValues() {
        final int NUM_CHAR=65535;
        Set<Integer> hashCodes = new HashSet<>();
        for (int i =0; i< Character.MAX_VALUE; i++) {
            char [] codes = Character.toChars(i);
            hashCodes.add(Hashcoder.hash(seed, codes[0]));
        }
        assertEquals(NUM_CHAR, hashCodes.size());
    }

    @Test
    void hashcodeFloatDifferentPrecisionSameHash() {
        Set<Integer> hashCodes = new HashSet<>();
        float [] floats = new float [] {1, 1.0f, 1.00f, 1.000f};
        for (float f: floats) {
            hashCodes.add(Hashcoder.hash(seed, f));
        }
        assertEquals(1, hashCodes.size());
    }

    @Test
    void hashcodeDoubleDifferentPrecisionSameHash() {
        Set<Integer> hashCodes = new HashSet<>();
        double [] doubles = new double [] {1, 1.0f, 1.00f, 1.000f};
        for (double d: doubles) {
            hashCodes.add(Hashcoder.hash(seed, d));
        }
        assertEquals(1, hashCodes.size());
    }

    @Test
    void hashcodeLong() {
        Set<Integer> hashCodes = new HashSet<>();
        double [] doubles = new double [] {1, 1.0f, 1.00f, 1.000f};
        for (double d: doubles) {
            hashCodes.add(Hashcoder.hash(seed, d));
        }
        assertEquals(1, hashCodes.size());
    }

    static class TestObject {
        // test spies
        boolean hashcodeInvoked = false;
        static int totalHashcodeInvocations = 0;

        public int hashCode() {
            totalHashcodeInvocations++;
            hashcodeInvoked = true;
            return 1;
        }   
    }

    @Test
    void hashCodeObjectInvokesObjectsHashcode(){
        TestObject test = new TestObject();
        Hashcoder.hash(seed, test);
        assertTrue(test.hashcodeInvoked);
    }

    @Test
    void hashCodeArrayInvokesHashcodeOnEachObject(){
        // 3 element array
        TestObject [] testObjects = new TestObject[]{new TestObject(),
             new TestObject(), new TestObject()};
        Hashcoder.hash(seed, testObjects);
        assertEquals(3, TestObject.totalHashcodeInvocations); 
    }

    @Test
    void hashcodeOfEmptyArrayReturnsSeed(){
        assertEquals(seed, Hashcoder.hash(seed, new TestObject[]{}));
    }
    
}
