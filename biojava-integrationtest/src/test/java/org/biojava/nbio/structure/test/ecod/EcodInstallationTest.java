/*
 * BioJava development code
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
 */

package org.biojava.nbio.structure.test.ecod;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.core.util.ConcurrencyTools;
import org.biojava.nbio.structure.ecod.EcodDatabase;
import org.biojava.nbio.structure.ecod.EcodDomain;
import org.biojava.nbio.structure.ecod.EcodFactory;
import org.biojava.nbio.structure.ecod.EcodInstallation;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Spencer Bliven
 *
 */
public class EcodInstallationTest {

	private static final Logger logger = LoggerFactory.getLogger(EcodInstallationTest.class);
	private static final String VERSION = "develop78";

	static {
		//System.setProperty("Log4jContextSelector", "org.apache.logging.log4j.core.async.AsyncLoggerContextSelector");
	} 
	@Rule
	public TemporaryFolder tmpFolder = new TemporaryFolder();
	@Test
	public void testDownloads() throws IOException {
		// Use second installation with tmp location to avoid overwriting main cache
		EcodInstallation ecod2 = new EcodInstallation(tmpFolder.getRoot().getAbsolutePath(),VERSION);
		// Delete old VERSION
		File domainsFile = new File(ecod2.getCacheLocation(),"ecod."+VERSION+".domains.txt");
		if( domainsFile.exists() ) {
			domainsFile.delete();
		}
		// Force download
		ecod2.ensureDomainsFileInstalled();
		// Check for download
		assertTrue("No downloaded file at "+domainsFile.toString(),domainsFile.exists());
	}


	@Test
	public void testAllDomains() throws IOException {
		int expected;
		EcodDatabase ecod = EcodFactory.getEcodDatabase(VERSION);

		List<EcodDomain> domains = ecod.getAllDomains();
		expected = 423779; //version77
		expected = 423869; //version78
		assertEquals("Wrong number of domains",expected,domains.size());
	}

	@Test
	public void testByPDB() throws IOException {
		EcodDatabase ecod = EcodFactory.getEcodDatabase(VERSION);

		String pdbId;
		String[] expectedDomains;
		List<EcodDomain> domains;

		pdbId = "1lyw";
		expectedDomains = new String[] {"e1lyw.1","e1lyw.2","e1lyw.3","e1lyw.4"};
		domains = ecod.getDomainsForPdb(pdbId);

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
		EcodDatabase ecod = EcodFactory.getEcodDatabase(VERSION);

		String ecodId;
		EcodDomain domain,expected;

		ecodId = "e1lyw.1";
		domain = ecod.getDomainsById(ecodId);
		expected = new EcodDomain(
				//				Long uid, String domainId, Boolean manual,
				20669l, "e1lyw.1", false,
				//				Integer xGroup, Integer hGroup, Integer tGroup, Integer fGroup, String pdbId,
				1,1,1,2,"1lyw",
				//				String chainId, String range, String architectureName,
				".", "A:3-97,B:106-346", "beta barrels",
				//				String xGroupName, String hGroupName, String tGroupName,
				//				String fGroupName, Boolean isAssembly, List<String> ligands
				"cradle loop barrel", "RIFT-related", "acid protease",
				"EF00710",//"UNK_F_TYPE",
				20669l, Collections.singleton("EPE")
				);
		assertEquals(ecodId,expected,domain);

		ecodId = "e4v4fAA1";
		domain = ecod.getDomainsById(ecodId);
		assertNotNull(ecodId,domain);
		assertEquals(ecodId,ecodId,domain.getDomainId());
	}

	@Test
	public void testMultithreaded() throws IOException {
		final EcodInstallation ecod = (EcodInstallation) EcodFactory.getEcodDatabase(VERSION);
		ecod.clear();
		String[] ecodIds = new String[] {
				"e4s1gA1", "e4umoB1", "e4v0cA1", "e4v1af1", "e3j7yj1", "e4wfcA1","e4b0jP1",
		};
		List<Future<EcodDomain>> futureDomains = new ArrayList<Future<EcodDomain>>();
		for(final String ecodId : ecodIds) {
			Callable<EcodDomain> job = new Callable<EcodDomain>() {
				@Override
				public EcodDomain call() throws Exception {
					logger.info("Running "+ecodId);
					EcodDomain d = ecod.getDomainsById(ecodId);
					logger.info("Finished "+ecodId);
					return d;
				}
				@Override
				public String toString() {
					return "Job fetching ECOD "+ecodId;
				}
			};
			Future<EcodDomain> future = ConcurrencyTools.submit(job,ecodId);
			futureDomains.add(future);
		}
		int successful = 0;
		for(Future<EcodDomain> future : futureDomains) {
			try {
				EcodDomain domain = future.get(60, TimeUnit.SECONDS);
				if(domain != null) {
					successful++;
				}
			} catch (InterruptedException e) {
				logger.error("Job "+future+" interrupted",e);
			} catch (ExecutionException e) {
				logger.error("Job "+future+" error",e);
			} catch (TimeoutException e) {
				logger.error("Job "+future+" timed out",e);
			}

		}
		assertEquals(ecodIds.length, successful);
	}

	@Test
	public void testFilterByHierarchy() throws IOException {
		EcodDatabase ecod = EcodFactory.getEcodDatabase(VERSION);

		List<EcodDomain> filtered;
		Set<String> expected,actual;

		expected = new HashSet<String>(Arrays.asList(
				"e4il6R1 e4pj0R1 e4pj0r1 e4ub6R1 e4ub8R1".split(" ") ));
		filtered = ecod.filterByHierarchy("6106.1.1");
		actual = new HashSet<String>();
		for(EcodDomain d : filtered) {
			actual.add(d.getDomainId());
		}
		assertEquals(expected,actual);

		filtered = ecod.filterByHierarchy("6106.1");
		actual = new HashSet<String>();
		for(EcodDomain d : filtered) {
			actual.add(d.getDomainId());
		}
		assertEquals(expected,actual);

		filtered = ecod.filterByHierarchy("6106");
		actual = new HashSet<String>();
		for(EcodDomain d : filtered) {
			actual.add(d.getDomainId());
		}
		assertEquals(expected,actual);
	}

	@Test
	public void testVersion() throws IOException {
		EcodDatabase ecod3 = EcodFactory.getEcodDatabase("latest");
		String version = ecod3.getVersion();
		assertNotNull(version);
		assertNotEquals("latest", version);
	}

	/**
	 * Parses all known versions. Only fails due to exceptions, so manually check for warnings.
	 * Hierarchical field warnings are expected for versions prior to develop68.
	 * @throws IOException
	 */
	//@Ignore // Very slow parsing test
	@Test
	public void testAllVersions() throws IOException {
		// Fetch latest version
		EcodDatabase latest = EcodFactory.getEcodDatabase("latest");
		String latestVersionStr = latest.getVersion();
		int latestVersion = 0;
		Matcher match = Pattern.compile("develop([0-9]+)",Pattern.CASE_INSENSITIVE).matcher(latestVersionStr);
		if(match.matches())
			latestVersion = Integer.parseInt(match.group(1));
		latest = null;

		// List all versions
		int firstVersion = 45;
		int lastVersion = Math.max(78,latestVersion);
		List<String> versions = new ArrayList<String>();
		versions.add("latest");
		for(int version = firstVersion; version<= lastVersion;version++) {
			versions.add("develop"+version);
		}
		
		// Parse all versions
		for(String version : versions) {
			EcodInstallation ecod = (EcodInstallation)EcodFactory.getEcodDatabase(version);
			ecod.getAllDomains();
			System.out.println(version +" -> "+ ecod.getVersion());
			
			// Force garbage collection of all soft references
			// This shouldn't be required, but without it we get
			// 'OutOfMemoryError: GC overhead limit exceeded'.
			// Probably this is due to synchronization in EcodFactory blocking
			// the GC during parsing. -Spencer
			ecod = null;
			System.gc();
			try {
				@SuppressWarnings("unused")
				Object[] ignored = new Object[(int) Runtime.getRuntime().maxMemory()];
			} catch (Throwable e) {
				// Ignore OME
			}
		}
	}
}
