package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.time.LocalDate;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class PrettyXMLWriterTest {

    StringWriter sw = null;
    XMLWriter xw = null;

    @BeforeEach
    void before() {
        sw = new StringWriter();
        xw = new PrettyXMLWriter(new PrintWriter(sw));
    }

    private static final String HTTP_TEST_NAMESPACE = "http://test-namespace";

    class TestObject {
        int a = 22;
        double b = 43.22;
        String text = "some text";
        LocalDate dt = LocalDate.of(2021, 8, 15);
        String timezone = "UTC";
    }

    @Test
    void simpleOutput() throws IOException {
        final String EXPECTED = "<to>\n" + "  <a>22</a>\n" + "  <b>43.22</b>\n" + "  <text>some text</text>\n"
                + "  <dt tz=\"UTC\">2021-08-15</dt>\n" + "</to>\n";
        TestObject to = new TestObject();

        xw.openTag("to");
        xw.openTag("a");
        xw.print(to.a + "");
        xw.closeTag("a");

        xw.openTag("b");
        xw.print(to.b + "");
        xw.closeTag("b");

        xw.openTag("text");
        xw.print(to.text);
        xw.closeTag("text");

        xw.openTag("dt");
        xw.attribute("tz", to.timezone);
        xw.print(to.dt.toString());
        xw.closeTag("dt");

        xw.closeTag("to");
        System.out.println(sw.toString());
        assertTrue(StringManipulationHelper.equalsToIgnoreEndline(EXPECTED, sw.toString()),
        	String.format("Strings are not equal (ignoring endline differences. expected: [%s], but it was:[%s]", EXPECTED, sw.toString()));
    }

    @Test
    void specialCharsAreEscaped() throws IOException {
        final String EXPECTED = "<dt>&#60;code&#62;some literal xml &#60;/code&#62;</dt>\n";
        xw.openTag("dt");
        xw.print("<code>some literal xml </code>");
        xw.closeTag("dt");

        assertTrue(StringManipulationHelper.equalsToIgnoreEndline(EXPECTED, sw.toString()),
            	String.format("Strings are not equal (ignoring endline differences. expected: [%s], but it was:[%s]", EXPECTED, sw.toString()));
    }

    @Test
    void namespacesAreAddedToElements() throws IOException {
        final String EXPECTED = "<ns1:dt xmlns:ns1=\"http://test-namespace\" ns1:myattr=\"1\">&#60;code&#62;some literal xml &#60;/code&#62;</ns1:dt>\n";

        xw.declareNamespace(HTTP_TEST_NAMESPACE, "test");
        // prefix
        xw.openTag(HTTP_TEST_NAMESPACE, "dt");
        xw.attribute(HTTP_TEST_NAMESPACE, "myattr", "1");
        xw.print("<code>some literal xml </code>");
        xw.closeTag(HTTP_TEST_NAMESPACE, "dt");

        assertTrue(StringManipulationHelper.equalsToIgnoreEndline(EXPECTED, sw.toString()),
            	String.format("Strings are not equal (ignoring endline differences. expected: [%s], but it was:[%s]", EXPECTED, sw.toString()));
    }
}
