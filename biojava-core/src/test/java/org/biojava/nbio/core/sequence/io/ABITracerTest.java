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
     * Test of process method, of class ABITracer.
     */
    @Test
    public void testProcess() throws Exception {
        logger.info("process");
        File file = new File("/3730.ab1");
        ABITrace tracer = new ABITrace(file);
        Assert.assertNotNull(tracer);
    }
}
