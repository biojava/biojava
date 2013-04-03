package org.biojava.bio.structure.scop;

import java.util.List;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupIterator;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;

import junit.framework.TestCase;


public class ScopTest extends TestCase {

	boolean debug = false;

	public void testLocalScop(){

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

		if (false){
			ScopDatabase scop = new RemoteScopInstallation();
			ScopDatabase defaultScop = ScopFactory.getSCOP();
			ScopFactory.setScopDatabase(scop);

			runSCOPTests();

			ScopFactory.setScopDatabase(defaultScop);
		}
		long timeE = System.currentTimeMillis();

		if ( debug ){
			System.out.println(timeE-timeS);
		}
	}


	private void runSCOPTests(){

		ScopDatabase scop = ScopFactory.getSCOP();

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
		AtomCache cache = new AtomCache();
		try {
			Structure s = cache.getStructureForDomain(d2);

			//checkRange(s,"A:496-581");
			// now with ligands!
			checkRange(s,"A:496-692");



		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}


		// check a domain with multiple ranges
		List<ScopDomain> domains1xzp = scop.getDomainsForPDB("1xzp");
		assertTrue(domains1xzp.size() ==4 );


		try {
			Structure s = cache.getStructureForDomain(domains1xzp.get(0));
			Chain a = s.getChainByPDB("A");

			// now with ligands...
			assertEquals(a.getAtomGroups().size(),176);


		}catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}

		// check insertion codes

		List<ScopDomain> domains2bq6 = scop.getDomainsForPDB("2bq6");
		assertTrue(domains2bq6.size() == 2);
		ScopDomain target = scop.getDomainByScopID("d2bq6a1");

		assertNotNull(target);
		try {
			Structure s = cache.getStructureForDomain(target);

			Chain a = s.getChainByPDB("A");
			assertEquals(a.getAtomGroups().size(),52);
			checkRange(s,"A:1A-49");

		}catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}		
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
		String chainId = g1.getChain().getChainID();
		String rangeTest = chainId + ":"+ g1.getResidueNumber().toString()+"-"+ g2.getResidueNumber().toString();

		assertEquals("The expected range and the detected range don;t match!", rangeTest, range);

	}


}
