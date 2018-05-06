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

package org.biojava.nbio.core.sequence.io;

import java.io.File;
import java.net.URL;
import org.junit.*;

public class ABITracerTest {

    private String filePath = System.getProperty("user.dir") + "/src/test/resources/3730.ab1";

    public ABITracerTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of URL method, of class ABITracer.
     */
    @Test
    public void testURL() throws Exception {
        File file = new File(filePath);
        Assert.assertNotNull(file);
        URL url = file.toURI().toURL();
        Assert.assertNotNull(url);
        ABITrace tracer = new ABITrace(url);
        Assert.assertNotNull(tracer);
    }

    /**
     * Test of Local file method, of class ABITracer.
     */
    @Test
    public void testLocal() throws Exception {
        File file = new File(filePath);
        Assert.assertNotNull(file);
        ABITrace tracer = new ABITrace(file);
        Assert.assertNotNull(tracer);
    }
}
