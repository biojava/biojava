package org.biojava.bio.structure.align.symm.order;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.symm.CeSymm;
import org.biojava.bio.structure.align.symm.CeSymmTest;
import org.biojava.bio.structure.align.util.AtomCache;
import org.junit.Test;

/**
 * Originally part of {@link CeSymmTest}.
 * @author Spencer Bliven
 */
public class SequenceFunctionOrderDetectorTest {

	@Test
	public void testGetSymmetryOrder() throws IOException, StructureException, OrderDetectionFailedException {
		
		// List of alignments to try, along with proper symmetry
		Map<String,Integer> orderMap = new LinkedHashMap<String,Integer>();
		orderMap.put("1itb.A", 3); // b-trefoil, C3
//		orderMap.put("1tim.A", 8); // tim-barrel, C8
		orderMap.put("d1p9ha_", 1); // not rotational symmetry
		orderMap.put("d1jlya1", 3); // a very nice trefoil
//		orderMap.put("d1ijqa1", 6);
		orderMap.put("d1h70a_", 5);
		
		AtomCache cache = new AtomCache();
		
		for (String name : orderMap.keySet()) {
			CeSymm ce = new CeSymm();
			
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);
			
			AFPChain afpChain = ce.align(ca1, ca2);
			
			int order = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
			
			assertEquals("Wrong order for "+name,orderMap.get(name).intValue(), order);
		}
	}
	
}
