package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertAll;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.RepeatedTest;
import org.junit.jupiter.api.Test;
/** 
*  
* @author Richard Adams
*/
class CRC64ChecksumTest {
    CRC64Checksum crc64 = null;
    private final String helloInCrc64Hex = "53C1111D27800000";
    private final Long helloInCrc64Decimal = 6035123792567599104L;

    @BeforeEach
    void before (){
         crc64 = new CRC64Checksum();
    }

    @Test
    @DisplayName("Default value is 0")
    void initialBehaviour() {
        assertEquals(0, crc64.getValue());
        assertEquals("0000000000000000", crc64.toString());
    }

    @RepeatedTest(10)
    void sameInputRepeatedlyGeneratesSameOutput(){
        crc64.update("hello");
        assertEquals(helloInCrc64Decimal, crc64.getValue());
        assertEquals(helloInCrc64Hex, crc64.toString());
    }

    @Test
    void afterResettingCrcIsZero(){
        crc64.update("hello");
        crc64.reset();
        assertEquals(0, crc64.getValue());
    }

    @Test
    void addingIncrementallyIsSameAsAllAtOnce(){
        crc64.update("h");
        crc64.update("e");
        crc64.update("l");
        crc64.update("l");
        crc64.update("o");
        assertEquals(helloInCrc64Hex, crc64.toString());
    }

    @Test
    void allbyteRange (){
        byte [] testBytes = new byte [] {1,2,3,4,5};
        crc64.update(testBytes, 0, testBytes.length);
        String allBytesHex = crc64.toString();
        crc64.reset();
        for (byte b: testBytes) {
            crc64.update(b);
        }
        assertEquals(allBytesHex, crc64.toString());
    }

    @Test 
    void partialByteRange (){
        byte [] testBytes = new byte [] {1,2,3,4,5};
        crc64.update(testBytes, 2, 1);
        String partialBytesHex = crc64.toString();
        crc64.reset();
        crc64.update(testBytes[2]);
        assertEquals(partialBytesHex, crc64.toString());
    }

    @Test 
    void partialByteRangeRejectsInvalidInput (){
        byte [] testBytes = new byte [] {1,2,3,4,5};
        assertAll(
            ()->assertThrows(IllegalArgumentException.class,
                ()->crc64.update(testBytes, -1, 0)),
            ()->assertThrows(IllegalArgumentException.class,
                ()->crc64.update(testBytes, 0, -1)),
            ()->assertThrows(IllegalArgumentException.class,
                ()->crc64.update(testBytes, 0, testBytes.length+1)),
            ()->assertThrows(IllegalArgumentException.class,
                ()->crc64.update(testBytes,  testBytes.length, 1))
        );
    }


    @Test
    void hexStringIsEqualToValue(){
        Long value = Long.parseLong(helloInCrc64Hex, 16);
        assertEquals(helloInCrc64Decimal, value);
    }
}