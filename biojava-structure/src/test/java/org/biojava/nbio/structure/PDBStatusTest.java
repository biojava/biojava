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
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.PDBStatus.Status;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;

/**
 * @author Spencer Bliven <sbliven@ucsd.edu>
 *
 */
public class PDBStatusTest {

	/**
	 * Test {@link PDBStatus#getStatus(String)}.
	 *
	 * <p>Uses the following PDBs:<br/>
	 * <pre>1HHB    OBSOLETE	replacedBy=4HHB
	 *4HHB    CURRENT	replaces=1HHB
	 *3HHB    CURRENT	replaces=1HHB
	 *</pre>
	 */
	@Test
	public void testGetStatus() throws IOException {
		Assert.assertEquals(Status.REMOVED, PDBStatus.getStatus("1HHB"));
		Assert.assertEquals(Status.CURRENT, PDBStatus.getStatus("3HHB"));
		Assert.assertEquals(Status.CURRENT, PDBStatus.getStatus("4HHB"));
	}

	@Test
	public void testGetStatusMultipleIds() throws IOException {
		String[] ids = {"1HHB", "3HHB", "4HHB"};
		Status[] statuses = PDBStatus.getStatus(ids);
		Assert.assertEquals(Status.REMOVED, statuses[0]);
		Assert.assertEquals(Status.CURRENT, statuses[1]);
		Assert.assertEquals(Status.CURRENT, statuses[2]);
	}

	@Test
	public void testGetCurrent() throws IOException {
		Assert.assertEquals("4HHB", PDBStatus.getCurrent("1HHB"));
		Assert.assertEquals("3HHB", PDBStatus.getCurrent("3HHB"));
		Assert.assertNull(PDBStatus.getCurrent("1CMW"));
		Assert.assertEquals("3ENI", PDBStatus.getCurrent("1KSA"));
		Assert.assertNull(PDBStatus.getCurrent("1CAT")); //rcsb.org returning null "id_code_replaced_by_latest"
		Assert.assertEquals("8CAT", PDBStatus.getCurrent("3CAT"));
		Assert.assertEquals("7CAT", PDBStatus.getCurrent("7CAT"));
	}
}
