/**
 * 
 */
package org.biojava.nbio.structure.domain;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;

/**
 * @author Spencer Bliven
 *
 */
public class EcodDomainTest {
	private static EcodInstallation ecod;
	private static final String VERSION = "develop77";

	// Set up static ecod singleton
	static {
		ecod = new EcodInstallation();
		ecod.setVersion(VERSION);
	}

	@Test
	public void testDownloads() throws IOException {
		// Delete old VERSION
		File domainsFile = new File(ecod.getCacheLocation(),"ecod."+VERSION+".domains.txt");
		if( domainsFile.exists() ) {
			domainsFile.delete();
		}
		// Force download
		ecod.ensureDomainsFileInstalled();
		// Check for download
		assertTrue("No downloaded file at "+domainsFile.toString(),domainsFile.exists());
	}

	@Test
	public void testAllDomains() throws IOException {
		List<EcodDomain> domains = ecod.getAllDomains();
		assertEquals("Wrong number of domains",423779,domains.size());
	}

	@Test
	public void testByPDB() throws IOException {
		String pdbId;
		String[] expectedDomains;
		List<EcodDomain> domains;

		pdbId = "1lyw";
		expectedDomains = new String[] {"e1lyw.1","e1lyw.2","e1lyw.3","e1lyw.4"};
		domains = ecod.getDomainsForPDB(pdbId);

		matchNames(pdbId,expectedDomains,domains);

	}

	private void matchNames(String pdbId,String[] expected,List<EcodDomain> actual) {
		assertEquals("Wrong number of domains for "+pdbId, expected.length, actual.size());
		Set<String> exp = new HashSet<String>(Arrays.asList(expected));
		for(EcodDomain d : actual) {
			assertTrue("Unexpected domain "+d.getDomainId()+" in "+pdbId,exp.contains(d.getDomainId()));
		}
	}

	@Test
	public void testParsing() throws IOException {
		String ecodId;
		EcodDomain domain,expected;

		ecodId = "e1lyw.1";
		domain = ecod.getDomainsById(ecodId);
		expected = new EcodDomain(
				//				Long uid, String domainId, Boolean manual,
				20669l, "e1lyw.1", null,
				//				Integer xGroup, Integer hGroup, Integer tGroup, String pdbId,
				1,1,1,"1lyw",
				//				String chainId, String range, String architectureName,
				".", "A:3-97,B:106-346", "beta barrels",
				//				String xGroupName, String hGroupName, String tGroupName,
				//				String fGroupName, Boolean isAssembly, List<String> ligands
				"cradle loop barrel", "RIFT-related", "acid protease",
				"UNK_F_TYPE", false, Collections.singleton("EPE")
				);
		assertEquals(ecodId,expected,domain);

	}
}
