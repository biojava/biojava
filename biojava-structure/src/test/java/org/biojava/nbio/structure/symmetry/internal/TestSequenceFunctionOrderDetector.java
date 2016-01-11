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
package org.biojava.nbio.structure.symmetry.internal;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.RefinerFailedException;
import org.biojava.nbio.structure.symmetry.internal.SequenceFunctionOrderDetector;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.junit.Test;

/**
 * Originally part of {@link CeSymmTest}.
 * @author Spencer Bliven
 */
public class TestSequenceFunctionOrderDetector {

	@Test
	public void testGetSymmetryOrder() 
			throws IOException, StructureException, RefinerFailedException {
		// List of alignments to try, along with proper symmetry
		Map<String,Integer> orderMap = new HashMap<String,Integer>();
		orderMap.put("1itb.A",3); // b-trefoil, C3
		orderMap.put("1tim.A",2); // tim-barrel, C8
		//orderMap.put("d1p9ha_",-1); // not rotational symmetry
		orderMap.put("3HKE.A",2); // very questionable alignment
		orderMap.put("d1jlya1",3); // a very nice trefoil
		
		AtomCache cache = new AtomCache();
		
		for(String name : orderMap.keySet()) {
			CeSymm ce = new CeSymm();
			ce.getParameters().setRefineMethod(RefineMethod.NOT_REFINED);
			Atom[] ca1 = cache.getAtoms(name);
			
			ce.analyzeLevel(ca1);
			AFPChain afpChain = ce.getSelfAlignments().get(0);
			
			int order = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
			
			assertEquals("Wrong order for "+name,orderMap.get(name).intValue(), order);
		}
	}
	
}
