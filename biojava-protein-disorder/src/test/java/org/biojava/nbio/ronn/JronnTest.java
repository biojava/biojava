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
package org.biojava.nbio.ronn;

import org.biojava.nbio.data.sequence.FastaSequence;
import org.biojava.nbio.ronn.Jronn.Range;
import org.junit.Test;

import static org.junit.Assert.assertEquals;


public class JronnTest {

	@Test
	public void verifyRanges() { 
  
	Range[]	ranges = Jronn.getDisorder(new FastaSequence("name", "LLRGRHLMNGTMIMRPWNFLNDHHFPKFFPHLIEQQAIWLADWWRKKHC" +
				"RPLPTRAPTMDQWDHFALIQKHWTANLWFLTFPFNDKWGWIWFLKDWTPGSADQAQRACTWFFCHGHDTN" +
				"CQIIFEGRNAPERADPMWTGGLNKHIIARGHFFQSNKFHFLERKFCEMAEIERPNFTCRTLDCQKFPWDDP" +
				"CSSTHSDCPKLEDLISFTETHGCSAADNADRPSQACHIGWAAMCEPTAMFMLMGSRCRCSFWPQNNAARHR" +
				"NFLIQIEMHSHLEHWIQTLHPQRPFLCNTWDDNWPICQFASQARGNSPDHHP"));
	assertEquals(4, ranges.length);
	assertEquals(53, ranges[0].from);
	assertEquals(59, ranges[0].to); 
	
	assertEquals(190, ranges[1].from);
	assertEquals(196, ranges[1].to);
	
	assertEquals(210, ranges[2].from);
	assertEquals(226, ranges[2].to);
	
	assertEquals(305, ranges[3].from);
	assertEquals(313, ranges[3].to);
	//System.out.println(Arrays.toString(ranges));
	}
}
