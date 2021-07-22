package org.biojava.nbio.core.util;
import static org.junit.Assert.assertThrows;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
class StringManipulationHelperTest {

    @Nested
    class PaddingTest {
	    @Test
        void padLeft() {
            assertEquals("     ",
            StringManipulationHelper.padLeft("",5));
            assertEquals("   xx",
                StringManipulationHelper.padLeft("xx",5));
            assertEquals("xxxxxx", StringManipulationHelper.padLeft("xxxxxx",5));
        }

        @Test
        void padRight() {
            assertEquals("     ",
            StringManipulationHelper.padRight("",5));
            assertEquals("xx   ",
                StringManipulationHelper.padRight("xx",5));
            assertEquals("xxxxxx", StringManipulationHelper.padRight("xxxxxx",5));
        }

        @ParameterizedTest
        @ValueSource(ints = {0,-1,-2})
        @DisplayName("invalid padding arguments throw IAE")
        void padInvalidValues(int invalidPadding) {
            assertThrows(IllegalArgumentException.class, 
              ()->StringManipulationHelper.padLeft(
               "anystring",invalidPadding));
            assertThrows(IllegalArgumentException.class, 
               ()->StringManipulationHelper.padRight(
                "anystring",invalidPadding));
        }
    
   }
   @Nested
   class InputStreamToString {

        @Test
        void basicString(){
           String singleLine = "hello";
           ByteArrayInputStream bais = new ByteArrayInputStream(singleLine.getBytes());
           assertEquals("hello\n", StringManipulationHelper.convertStreamToString(bais));
        }

        @ParameterizedTest
        @DisplayName("Newlines are converted to Unix newlines")
        @ValueSource(strings={"line1\r\nline2", "line1\nline2", "line1\rline2"})
        void multiLineConvertedToUnixNewLine(String multiline){
           ByteArrayInputStream bais = new ByteArrayInputStream(multiline.getBytes());
           assertEquals("line1\nline2\n", StringManipulationHelper.convertStreamToString(bais));
        }

        @Test
        void streamIsClosedAfterCompletion() throws IOException{
            // this is a stream that will throw IOException
            // if called after closing
            InputStream is = InputStream.nullInputStream();
            StringManipulationHelper.convertStreamToString(is);
            // attempt to read again after closing
            assertThrows(IOException.class, ()->is.read());
        }

        @Test
        void emptyStreamGeneratesEmptyString() {
            assertEquals("", StringManipulationHelper.convertStreamToString(
                new ByteArrayInputStream(new byte [0])));
        }
    }
}