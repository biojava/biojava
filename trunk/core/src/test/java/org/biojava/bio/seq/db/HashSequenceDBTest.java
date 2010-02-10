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

package org.biojava.bio.seq.db;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Runs SequenceDB tests using HashSequenceDB.
 * 
 * @author Len Trigg
 */
public class HashSequenceDBTest extends AbstractSequenceDBTest {
    

    public HashSequenceDBTest(String name) {
        super(name);
    }

    protected SequenceDB getSequenceDB() throws Exception {
        return new HashSequenceDB();
    }
    
    public static Test suite() {
        return new TestSuite(HashSequenceDBTest.class);
    }
    
    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }    
}
