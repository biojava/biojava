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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ABITracerTest {

    private final static Logger logger = LoggerFactory.getLogger(ABITracerTest.class);

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
     * Test of local method, of class ABITracer.
     */
    @Test
    public void testURL() throws Exception {
        URL url = new URL("https://github.com/biopython/biopython/blob/master/Tests/Abi/3730.ab1");
        Assert.assertNotNull(url);
        ABITrace tracer = new ABITrace(url);
        Assert.assertNotNull(tracer);
    }

    /**
     * Test of local method, of class ABITracer.
     */
    @Test
    public void testLocal() throws Exception {
        File file = new File("/3730.ab1");
        Assert.assertNotNull(file);
        ABITrace tracer = new ABITrace(file);
        Assert.assertNotNull(tracer);
    }
}
