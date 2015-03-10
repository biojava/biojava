/**
 * 
 */
package org.biojava.nbio.structure.test.domain;

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

import org.biojava.nbio.core.util.ConcurrencyTools;
import org.biojava.nbio.structure.domain.EcodDomain;
import org.biojava.nbio.structure.domain.EcodInstallation;
import org.biojava.nbio.structure.io.util.FileDownloadUtils;
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
	private static EcodInstallation ecod;
	private static final String VERSION = "develop77";

	// Set up static ecod singleton
	static {
		ecod = new EcodInstallation();
		ecod.setVersion(VERSION);
	}

	static {
		//System.setProperty("Log4jContextSelector", "org.apache.logging.log4j.core.async.AsyncLoggerContextSelector");
		} 
	@Rule
	public TemporaryFolder tmpFolder = new TemporaryFolder();
	@Test
	public void testDownloads() throws IOException {
		// Use second installation with tmp location to avoid overwriting main cache
		EcodInstallation ecod2 = new EcodInstallation(tmpFolder.getRoot().getAbsolutePath());
		ecod2.setVersion(VERSION);
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

	@Test
	public void testMultithreaded() throws IOException {
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
}
