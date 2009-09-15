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
package org.biojava.utils.process;

import junit.framework.Test;
import junit.framework.TestSuite;

/** 
 * @author <a href="mailto:Martin.Szugat@GMX.net">Martin Szugat</a>
 * 
 */
public class AllTests {

    public static void main(String[] args) {
        junit.textui.TestRunner.run(AllTests.suite());
    }

    public static Test suite() {
        TestSuite suite = new TestSuite("Test for org.biojava.utils.process");
        //$JUnit-BEGIN$
        suite.addTestSuite(ExternalProcessTest.class);
        //$JUnit-END$
        return suite;
    }

}
