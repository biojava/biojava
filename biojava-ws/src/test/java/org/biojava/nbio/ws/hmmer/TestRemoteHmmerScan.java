package org.biojava.nbio.ws.hmmer;

import java.util.SortedSet;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.junit.Test;
import static org.junit.Assert.*;

public class TestRemoteHmmerScan {
	
	/**
	 * Sequence for UniProt id P30340 (PDB 1SMT)
	 */
	private static final String TEST_SEQ = "MTKPVLQDGETVVCQGTHAAIASELQAIAPEVAQSLAEFFAVLADPNRLRLLSLLARSEL" + 
			"CVGDLAQAIGVSESAVSHQLRSLRNLRLVSYRKQGRHVYYQLQDHHIVALYQNALDHLQE" + 
			"CR";
	
	@Test
	public void testHmmerWs() throws Exception {

		ProteinSequence seq = new ProteinSequence(TEST_SEQ);

		// now we submit this sequence to the Hmmer web site
		RemoteHmmerScan hmmer = new RemoteHmmerScan();

		SortedSet<HmmerResult> results = hmmer.scan(seq);
		
		assertNotNull(results);
		// 2 results (domains) for P30340 (PDB 1smt) as of Jan 2018
		assertEquals(2, results.size());

		boolean gotSh2Domain = false;
		
		for (HmmerResult hmmerResult : results) {
			if (hmmerResult.getName().equals("HTH_5")) {
				gotSh2Domain = true;
			}
		}
		
		assertTrue("A HTH_5 domain should be present as one of the hmmer scan matches",gotSh2Domain);
	}
}
