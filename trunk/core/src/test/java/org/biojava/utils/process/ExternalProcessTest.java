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

/*
 *    ExternalProcessTest.java
 */
package org.biojava.utils.process;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.Properties;

import junit.extensions.RepeatedTest;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Test class for {@link org.biojava.utils.process.ExternalProcess}.
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * @version $Revision$
 */
public class ExternalProcessTest extends TestCase {

    /**
     * Input for the more command.
     */
    private static final String INPUT_MORE = "test";

    /**
     * Object to execute external processes.
     */
    private static ExternalProcess ep = null;

    /**
     * OS-specific string for the more command.
     */
    private static String commandMore = null;
    
    static {
        if (System.getProperty("os.name").toLowerCase().startsWith("windows")) {
            commandMore = "cmd /C more";
        } else {
            commandMore = "cat";
        }
        ep = new ExternalProcess();
    }

    /**
     * Main entry point. Runs a repeated test to stress the thread allocation 
     * and freeing of the <code>ExternalProcess</code> class.
     * @param args command line arguments.
     */
    public static void main(String[] args) {
        TestSuite suite = new TestSuite();
        suite.addTest(new TestSuite(ExternalProcessTest.class));
        junit.textui.TestRunner.run(new RepeatedTest(suite, 5000));
    }

    /**
     * Initializes the test.
     * @param name test name
     */
    public ExternalProcessTest(String name) {
        super(name);
    }

    /**
     * Test the input/output handling, using input/output streams and not 
     * readers/writers.
     * @throws IOException if execution fails.
     * @throws InterruptedException if execution fails.
     */
    public void testExecuteStreams() throws IOException,
            InterruptedException {
        ByteArrayInputStream input = new ByteArrayInputStream(
                INPUT_MORE.getBytes());
        ByteArrayOutputStream output = new ByteArrayOutputStream();
        ByteArrayOutputStream error = new ByteArrayOutputStream();
        ep.setCommands(commandMore);
        ep.setInputHandler(new SimpleInputHandler(input, "in"));
        ep.setOutputHandler(new SimpleOutputHandler(output, "out"));
        ep.setErrorHandler(new SimpleOutputHandler(error, "err"));
        int exitCode = ep.execute();
        assertTrue(exitCode == 0);
        assertEquals(INPUT_MORE, output.toString().trim());
        assertEquals("", error.toString());
    }

    /**
     * Test the input/output handling.
     * @throws IOException if execution fails.
     * @throws InterruptedException if execution fails.
     */
    public void testExecuteMore() throws IOException,
            InterruptedException {
        StringWriter outputString = new StringWriter();
        StringWriter errorString = new StringWriter();
        StringReader inputString = new StringReader(INPUT_MORE);
        ep.setCommands(commandMore);
        ep.setOutputHandler(new WriterOutputHandler(outputString, "out"));
        ep.setErrorHandler(new WriterOutputHandler(errorString, "err"));
        ep.setInputHandler(new ReaderInputHandler(inputString, "in"));
        int exitCode = ep.execute();
        assertTrue(exitCode == 0);
        assertEquals(INPUT_MORE, outputString.toString().trim());
        assertEquals("", errorString.toString());
    }

    /**
     * Tests the static execute method.
     * @throws Exception if execution fails.
     */
    public void testStaticExecute() throws Exception {
        StringWriter outputString = new StringWriter();
        StringWriter errorString = new StringWriter();
        int exitCode = ExternalProcess.execute(commandMore, INPUT_MORE, 
                outputString, errorString);
        assertTrue(exitCode == 0);
        assertEquals(INPUT_MORE, outputString.toString().trim());
        assertEquals("", errorString.toString());
    }

    /**
     * Tests the {@link ExternalProcess#resolveCommands(String, Properties)}
     * method.
     */
    public void testResolveCommands() {
        Properties props = new Properties();
        props.put("TEST", "test");
        String commands = "do %TEST% it";
        commands = ExternalProcess.resolveCommands(commands, props);
        assertEquals("do test it", commands);
    }
    
    /**
     * Tests the {@link ExternalProcess#joinCommands(Object[])} method.
     */
    public void testJoinCommands() {
        String[] commandList = new String[] {"do", "test", "it"};
        String commands = ExternalProcess.joinCommands(commandList);
        assertEquals("do test it", commands);
    }
}
