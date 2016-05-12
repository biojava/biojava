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
package org.biojava.nbio.structure.test.scop;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupIterator;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.scop.ScopInstallation;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

//import org.biojava.nbio.structure.scop.RemoteScopInstallation;



public class ScopTest {

	boolean debug = false;

	@Test
	public void testLocalScop() throws IOException, StructureException{ 

		if ( debug ){
			System.out.println("local");
		}
		long timeS = System.currentTimeMillis();
		ScopDatabase scop = new ScopInstallation();
		ScopDatabase defaultScop = ScopFactory.getSCOP();
		ScopFactory.setScopDatabase(scop);

		runSCOPTests();

		ScopFactory.setScopDatabase(defaultScop);

		long timeE = System.currentTimeMillis();

		if ( debug ){
			System.out.println(timeE-timeS);
		}
	}


	public void testRemoteScop(){

		if (debug) {
			System.out.println("remote");
		}
		long timeS = System.currentTimeMillis();

//		if (false){
//			ScopDatabase scop = new RemoteScopInstallation();
//			ScopDatabase defaultScop = ScopFactory.getSCOP();
//			ScopFactory.setScopDatabase(scop);
//
//			runSCOPTests();
//
//			ScopFactory.setScopDatabase(defaultScop);
//		}
		long timeE = System.currentTimeMillis();

		if ( debug ){
			System.out.println(timeE-timeS);
		}
	}


	private void runSCOPTests() throws IOException, StructureException {

		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75);

		List<ScopDomain> domains = scop.getDomainsForPDB("4HHB");

		assertTrue(domains.size() == 4);

		// test case sensitivity;
		List<ScopDomain> domains2 = scop.getDomainsForPDB("4hhb");
		assertTrue(domains2.size() == domains.size());

		//System.out.println(domains);


		String scop1m02 = "d1m02a_	1m02	A:	k.36.1.1	74353	cl=58788,cf=75796,sf=75797,fa=75798,dm=75799,sp=75800,px=74353";

		List<ScopDomain> domains1m02 = scop.getDomainsForPDB("1m02");
		assertTrue(domains1m02.size() == 1);
		ScopDomain d1 = domains1m02.get(0);

		assertNotNull(d1);

		assertEquals("The toString() methods for ScopDomains don't match the scop display",d1.toString(),scop1m02);


		List<ScopDomain> domains1cdg = scop.getDomainsForPDB("1CDG");
		assertTrue(domains1cdg.size() == 4);
		ScopDomain d2 = domains1cdg.get(0);
		assertEquals("Wrong SCOP Id", "d1cdga1", d2.getScopId());
		AtomCache cache = new AtomCache();
		
		Structure s = cache.getStructureForDomain(d2);
		/*
			The actual SCOP description is A:496-581.
			HET    MAL  A 688      23
			HET    MAL  A 689      23
			HET    MAL  A 690      23
			HET     CA  A 691       1
			HET     CA  A 692       1
			Thus, the two hetero-atoms are out of range.
		 */
		// t
		//checkRange(s,"A:496-581");
		// now with ligands!
		checkRange(s,"A:496-692");






		// check a domain with multiple ranges
		List<ScopDomain> domains1xzp = scop.getDomainsForPDB("1xzp");
		assertTrue(domains1xzp.size() ==4 );


		

		s = cache.getStructureForDomain(domains1xzp.get(0)); 

		Chain a = s.getPolyChainByPDB("A");

		// since biojava 5.0, there's no ligands in the polymer chains: thus we substract 3 SO4 molecules present in chain A
		assertEquals(176 -3,a.getAtomGroups().size());

		

		// check insertion codes

		List<ScopDomain> domains2bq6 = scop.getDomainsForPDB("2bq6");
		assertTrue(domains2bq6.size() == 2);
		ScopDomain target = scop.getDomainByScopID("d2bq6a1");

		assertNotNull(target);
		
		s = cache.getStructureForDomain(target);

		a = s.getPolyChainByPDB("A");
		assertEquals(a.getAtomGroups().size(),52);
		checkRange(s,"A:1A-49");

	}

	private void checkRange(Structure s, String range) {
		GroupIterator iter = new GroupIterator(s);
		Group g1 = iter.next();
		Group g2 =null;
		while (iter.hasNext()){
			g2 = iter.next();
		}
		assertNotNull(g1);
		assertNotNull(g2);
		String chainId = g1.getChain().getName();
		String rangeTest = chainId + ":"+ g1.getResidueNumber().toString()+"-"+ g2.getResidueNumber().toString();

		assertEquals("The expected range and the detected range don't match!", range, rangeTest);

	}


}
