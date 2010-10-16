/*
 * CrossReferenceResolutionExceptionTest.java
 * JUnit based test
 *
 * Created on November 11, 2005, 3:48 PM
 */

package org.biojavax;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 *
 * @author Mark Schreiber
 */
public class CrossReferenceResolutionExceptionTest extends TestCase {
    
    CrossReferenceResolutionException ex;
    String mess;
    Exception cause;
    
    public CrossReferenceResolutionExceptionTest(String testName) {
        super(testName);
        mess = "message";
        cause = new Exception();
    }

    protected void setUp() throws Exception {
        ex = new CrossReferenceResolutionException(mess, cause);
    }

    protected void tearDown() throws Exception {
    }

    public static Test suite() {
        TestSuite suite = new TestSuite(CrossReferenceResolutionExceptionTest.class);
        
        return suite;
    }
    
    public void testGetMessage(){
        assertEquals(mess, ex.getMessage());
    }
    
    public void testGetCause(){
        assertEquals(cause, ex.getCause());
    }
}
