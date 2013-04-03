package org.biojava3.ronn;

import static org.junit.Assert.assertEquals;

import org.biojava3.data.sequence.FastaSequence;
import org.biojava3.ronn.Jronn.Range;
import org.junit.Test;


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
