package org.biojava.nbio.core.util;
import static org.biojava.nbio.core.util.StringManipulationHelper.equalsToIgnoreEndline;
import static org.junit.Assert.assertThrows;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

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
        // in java11 there is a NullInputStream for this
        class InputStreamTss extends InputStream {
        	boolean closed = false;
			@Override
			public int read() throws IOException {
				if (closed) {
					throw  new IOException();
				}
				return -1;
			}
			public void close() throws IOException {
				closed = true;
			}
        	
        }
        

        @Test
        void streamIsClosedAfterCompletion() throws IOException{
            // this is a stream that will throw IOException
            // if called after closing
            InputStream is =  new InputStreamTss();

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

    @Nested
    class equalsToIgnoreEndline{
        @Test
        void emptyOrNullStringsAreEqual() {
            assertTrue(equalsToIgnoreEndline("",""));
            assertTrue(equalsToIgnoreEndline(null, null));
        }

        @Test
        void emptyVsNullStringsAreNotEqual() {
            assertFalse(equalsToIgnoreEndline(null,""));
        }

        @Test
        @DisplayName("multiline strings with different line terminators are equal")
        void differentLineTerminatorsAreEqual() {
            assertTrue(equalsToIgnoreEndline("ab\ncd\nef","ab\r\ncd\r\nef"));
            assertTrue(equalsToIgnoreEndline("ab\r\ncd\nef","ab\rcd\ref"));
        }

        @Test
        @DisplayName("comparison is case-sensitive")
        void caseSensitive() {
            assertFalse(equalsToIgnoreEndline("ab\ncd\nef","ab\nCD\nef"));
        }

        @Test
        @DisplayName("multiline strings with different lengths are  unequal")
        void s2LongerThanS1() {
            assertFalse(equalsToIgnoreEndline("ab\ncd\nef","ab\ncd\nef\nextra-line"));
        }

        @Test
        @DisplayName("multiline strings with different lengths are  unequal")
        void s1LongerThanS2() {
            assertFalse(equalsToIgnoreEndline("ab\ncd\nef\nextra","ab\ncd\nef"));
        }
    }
    @Nested
    class JoinString{   
        List<String> empty = new ArrayList<>();
        List<String> items =  new ArrayList<>();
        void populateItems() {
            items.add("a");
            items.add("b");
            items.add("c");
        }
        @Test
        void join() {

            assertEquals("", StringManipulationHelper.join(empty,","));
            assertEquals("", StringManipulationHelper.join(null,","));
            items.add("a");
            assertEquals("a", StringManipulationHelper.join(items,","));
            items.add("b");
            items.add("c");
            assertEquals("a,b,c", StringManipulationHelper.join(items,","));
            assertEquals("abc", StringManipulationHelper.join(items,""));
        }
        @Test
        void delimiterCanBeAnyLength(){
            populateItems();
            assertEquals("a---b---c", StringManipulationHelper.join(items,"---"));      
        }
    }

    @Nested
    class EqualsToXml {

        String docType ="<!DOCTYPE " +
        "ex  [ <!ENTITY foo \"foo\"> <!ENTITY bar"+
        " \"bar\">]> <ex/>";
        @Test
        void isNotImplemented() {
            assertThrows(UnsupportedOperationException.class, 
            ()->StringManipulationHelper.equalsToXml(docType, docType));
        }
    }
}