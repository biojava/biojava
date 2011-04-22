/**
 * 
 */
package org.biojava.bio.structure;

import org.biojava.bio.structure.PDBStatus;
import org.biojava.bio.structure.PDBStatus.Status;

import junit.framework.TestCase;

/**
 * @author Spencer Bliven <sbliven@ucsd.edu>
 *
 */
public class PDBStatusTest extends TestCase {

	/**
	 * Test {@link PDBStatus#getStatus(String)}.
	 * 
	 * <p>Uses the following PDBs:<br/>
	 * <pre>1HHB    OBSOLETE
	 *4HHB    CURRENT
	 *3HHB    CURRENT
	 *</pre>
	 */
	public void testGetStatus() {
		assertEquals(Status.OBSOLETE, PDBStatus.getStatus("1HHB"));
		assertEquals(Status.CURRENT, PDBStatus.getStatus("3HHB"));
		assertEquals(Status.CURRENT, PDBStatus.getStatus("4HHB"));
	}
	
	public void testGetReplacement() {
		//TODO should test a case with a longer/more complicated version tree
		assertEquals("4HHB",PDBStatus.getReplacement("1HHB",false));
		assertEquals("3HHB",PDBStatus.getReplacement("3HHB",false));
	}
	

	public void testGetCurrent() {
		//TODO should test a case with a longer/more complicated version tree
		assertEquals("4HHB",PDBStatus.getCurrent("1HHB"));
		assertEquals("3HHB",PDBStatus.getCurrent("3HHB"));
		assertEquals(null, PDBStatus.getCurrent("1CMW"));
	}
	
	public void testGetReplaces() {
		//TODO should test a case with a longer/more complicated version tree
		assertEquals("1HHB",PDBStatus.getReplaces("4HHB"));
		assertEquals("1HHB",PDBStatus.getReplaces("3HHB"));
		assertEquals(null, PDBStatus.getReplaces("1HHB"));
	}
	
}
