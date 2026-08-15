/**
 * BioJava development code
 *
 * This code may be freely distributed and modified under the terms of the GNU
 * Lesser General Public Licence. This should be distributed with the code. If
 * you do not have a copy, see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual authors. These
 * should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims, or to join the
 * biojava-l mailing list, visit the home page at:
 *
 * http://www.biojava.org/
 */
package org.biojava.nbio.structure.test.io.density;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.core.util.FileDownloadUtils;
import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.io.density.Ccp4Header;
import org.biojava.nbio.structure.io.density.DensityFileFormat;
import org.biojava.nbio.structure.io.density.DensityMapCache;
import org.biojava.nbio.structure.io.density.DensityMapKind;
import org.biojava.nbio.structure.io.density.DensityMapRequest;
import org.biojava.nbio.structure.io.density.DensityMapResult;
import org.biojava.nbio.structure.io.density.DensityMapSource;
import org.biojava.nbio.structure.io.density.NoDensityMapException;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 * Density fetching against the real services.
 * <p>
 * Deliberately frugal: the entries chosen keep a full run to a couple of
 * megabytes plus a few small metadata calls. In particular the cryo-EM path is
 * exercised with the size limit set so low that the 116 MB map is declined
 * before any of its body is transferred, which tests the whole resolution and
 * guard sequence without the download.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class DensityMapIntegrationTest {

	private File cacheRoot;
	private DensityMapCache cache;

	@Before
	public void setUp() throws IOException {
		cacheRoot = Files.createTempDirectory("bj-density-it").toFile();
		cache = new DensityMapCache(cacheRoot.getAbsolutePath());
	}

	@After
	public void tearDown() throws IOException {
		FileDownloadUtils.deleteDirectory(cacheRoot.toPath());
	}

	/** The default path for an X-ray entry: the smallest source answers first. */
	@Test
	public void fetchesAnXrayMapFromTheFirstSourceTried() throws IOException {
		DensityMapResult result = cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);

		assertEquals(DensityMapSource.RCSB_VOLUME_SERVER, result.getSource());
		assertEquals(DensityMapKind.TWO_FO_FC, result.getKind());
		assertTrue(result.isRenderable());
		assertFalse(result.isFromCache());
		assertTrue(result.getFileSizeBytes() > 1024);
		assertTrue("a .meta sidecar makes the result reconstructible offline",
				DensityMapResult.metaFileFor(result.getFile()).isFile());

		// second call must come from the cache without another download
		DensityMapResult again = cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);
		assertTrue(again.isFromCache());
		assertEquals(result.getFile(), again.getFile());
	}

	/**
	 * Both map kinds come out of one download, and the difference map is presented
	 * under the companion name that makes Jmol read the other data block.
	 */
	@Test
	public void bothKindsShareASingleDownload() throws IOException {
		DensityMapResult twoFoFc = cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);
		DensityMapResult foFc = cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.FO_FC);

		assertEquals(DensityMapKind.FO_FC, foFc.getKind());
		assertFalse("the difference map needs its own file name", twoFoFc.getFile().equals(foFc.getFile()));
		assertTrue("the marker has to be in the name for Jmol to select the FO-FC block",
				foFc.getFile().getName().contains("&diff=1"));
		assertEquals("both names must address the same bytes",
				twoFoFc.getFileSizeBytes(), foFc.getFileSizeBytes());
	}

	/** PDBe serves real CCP4 files, which the header check should recognise. */
	@Test
	public void pdbeServesAGenuineCcp4Map() throws IOException {
		cache.setSourceChain(DensityMapKind.TWO_FO_FC, Arrays.asList(DensityMapSource.PDBE_CCP4));
		DensityMapResult result = cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);

		assertEquals(DensityMapSource.PDBE_CCP4, result.getSource());
		assertEquals(DensityFileFormat.CCP4, result.getFormat());
		assertTrue("the CCP4 stamp should be present at byte 208", Ccp4Header.isCcp4(result.getFile()));
		assertTrue(FileDownloadUtils.validateFile(result.getFile()));
	}

	/**
	 * The whole cryo-EM route: resolve the EMDB entry, pick up the author contour
	 * level, and decline the full map on size without transferring it.
	 */
	@Test
	public void resolvesCryoEmEntriesAndHonoursTheSizeLimit() throws IOException {
		List<String> emdbIds = cache.getEmdbResolver().getEmdbIds(new PdbId("6hu9"));
		assertEquals(Arrays.asList("EMD-0262"), emdbIds);

		DensityMapResult result = cache.getDensityMap(new PdbId("6hu9"), DensityMapKind.AUTO);
		assertEquals(DensityMapKind.EM, result.getKind());
		assertEquals("EMD-0262", result.getEmdbId());
		assertNotNull("EM maps need the author contour level to be displayed properly",
				result.getRecommendedContourLevel());
		assertEquals(0.0263, result.getRecommendedContourLevel(), 1e-6);
		assertNotNull(result.getContourInSigma());

		// With only the full archive enabled and a tiny ceiling, the guard must fire
		// rather than pulling down 116 MB.
		DensityMapCache strict = new DensityMapCache(cacheRoot.getAbsolutePath());
		strict.setSourceChain(DensityMapKind.EM, Arrays.asList(DensityMapSource.EMDB_MAP));
		strict.setMaxDownloadBytes(1024);
		try {
			strict.getDensityMap(DensityMapRequest.builder(new PdbId("6hu9")).kind(DensityMapKind.EM).build());
			fail("the size guard should have declined the full EMDB map");
		} catch (NoDensityMapException e) {
			assertTrue(e.getAttempts().get(DensityMapSource.EMDB_MAP).contains("too large"));
		}
	}

	/** 4HHB was deposited in 1984 without structure factors, so nothing has a map for it. */
	@Test
	public void reportsWhyAnEntryHasNoDensity() throws IOException {
		cache.setSourceEnabled(DensityMapSource.WWPDB_MAP_COEFFICIENTS, true);
		try {
			cache.getDensityMap(new PdbId("4hhb"), DensityMapKind.AUTO);
			fail("4hhb has no deposited structure factors");
		} catch (NoDensityMapException e) {
			assertFalse(e.getAttempts().isEmpty());
			assertTrue(e.getAttempts().values().stream().anyMatch(r -> r.contains("404")));
		}
	}

	/**
	 * The wwPDB servers return the content MD5 as the ETag, so a coefficient
	 * download is checksum-verified without a separate hash file.
	 */
	@Test
	public void mapCoefficientsArriveWithAVerifiableChecksum() throws IOException {
		cache.setSourceEnabled(DensityMapSource.WWPDB_MAP_COEFFICIENTS, true);
		cache.setSourceChain(DensityMapKind.TWO_FO_FC, Arrays.asList(DensityMapSource.WWPDB_MAP_COEFFICIENTS));

		DensityMapResult result = cache.getDensityMap(DensityMapRequest.builder(new PdbId("1cbs"))
				.kind(DensityMapKind.TWO_FO_FC)
				.allowNonRenderableFormats(true)
				.build());

		assertEquals(DensityMapSource.WWPDB_MAP_COEFFICIENTS, result.getSource());
		assertFalse("structure factors are not a map and must not claim to be renderable",
				result.isRenderable());

		File hashFile = new File(result.getFile().getParentFile(), result.getFile().getName() + ".hash_MD5");
		assertTrue("an MD5 should have been recorded from the ETag", hashFile.isFile());
		assertTrue(FileDownloadUtils.validateFile(result.getFile()));

		// corrupt it and confirm the checksum actually catches it
		Files.write(result.getFile().toPath(), new byte[] {0, 1, 2, 3});
		assertFalse(FileDownloadUtils.validateFile(result.getFile()));
	}
}
