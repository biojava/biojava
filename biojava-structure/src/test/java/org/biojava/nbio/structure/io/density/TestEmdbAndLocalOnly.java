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
package org.biojava.nbio.structure.io.density;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

import org.biojava.nbio.core.util.FileDownloadUtils;
import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 * EMDB metadata parsing against captured responses, and the guarantee that
 * LOCAL_ONLY never opens a connection.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class TestEmdbAndLocalOnly {

	/**
	 * Nothing listens here, so any attempt to reach a server fails immediately
	 * rather than hanging or, worse, quietly succeeding.
	 */
	private static final String DEAD_SERVER = "http://localhost:1/";

	private File cacheRoot;

	@Before
	public void setUp() throws IOException {
		cacheRoot = Files.createTempDirectory("bj-density").toFile();
	}

	@After
	public void tearDown() throws IOException {
		FileDownloadUtils.deleteDirectory(cacheRoot.toPath());
		PdbeCcp4MapProvider.resetToDefaults();
		WwpdbMapCoefficientsProvider.resetToDefaults();
		EmdbMapProvider.resetToDefaults();
		EmdbEntryResolver.resetToDefaults();
	}

	private void copyResource(String resource, File target) throws IOException {
		target.getParentFile().mkdirs();
		try (InputStream in = getClass().getResourceAsStream(resource)) {
			assertNotNull("missing test resource " + resource, in);
			Files.copy(in, target.toPath());
		}
	}

	/**
	 * The map metadata carries the two numbers a viewer needs: the level the
	 * depositors recommend and the RMS deviation that converts it to sigma.
	 */
	@Test
	public void parsesEmdbMapMetadata() throws IOException {
		copyResource("emdb-map-EMD-0262.json",
				DensityCacheLayout.emdbMapInfoFile(cacheRoot, "EMD-0262"));

		EmdbEntryResolver resolver = new EmdbEntryResolver(cacheRoot);
		resolver.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		EmdbEntryInfo info = resolver.getEntryInfo("EMD-0262");

		assertNotNull(info);
		assertEquals("EMD-0262", info.getEmdbId());
		assertEquals(0.0263, info.getRecommendedContourLevel(), 1e-9);
		assertNotNull("sigma is needed to express the contour in sigma units", info.getSigma());
		// 119165 kB, which is well past the default download ceiling
		assertEquals(119165L * 1024L, info.getMapSizeBytes().longValue());
	}

	/**
	 * The size of a real EM map, and why the source order rather than the size
	 * guard is what keeps it from being fetched.
	 * <p>
	 * EMD-0262 is about 116 MB, which sits comfortably under the 256 MiB default
	 * ceiling: left to the ceiling alone it would be downloaded in full. What
	 * avoids that is putting the density servers ahead of the archive, which serve
	 * a few megabytes for the same entry. The ceiling is the backstop for the maps
	 * that run to gigabytes.
	 */
	@Test
	public void theFullEmMapIsLargeButWithinTheDefaultCeiling() throws IOException {
		copyResource("emdb-map-EMD-0262.json",
				DensityCacheLayout.emdbMapInfoFile(cacheRoot, "EMD-0262"));
		EmdbEntryResolver resolver = new EmdbEntryResolver(cacheRoot);
		resolver.setFetchBehavior(FetchBehavior.LOCAL_ONLY);

		long bytes = resolver.getEntryInfo("EMD-0262").getMapSizeBytes();
		assertTrue("expected a map of order 100 MB, got " + bytes, bytes > 100L * 1024 * 1024);
		assertTrue("the ceiling alone would not stop this download",
				bytes < DensityMapCache.DEFAULT_MAX_DOWNLOAD_BYTES);

		DensityMapCache cache = new DensityMapCache(cacheRoot.getAbsolutePath());
		List<DensityMapSource> em = cache.getSourceChain(DensityMapKind.EM);
		assertTrue("a density server must be tried before the full archive",
				em.indexOf(DensityMapSource.RCSB_VOLUME_SERVER) < em.indexOf(DensityMapSource.EMDB_MAP));
	}

	@Test
	public void contourConvertsToSigma() throws IOException {
		copyResource("emdb-map-EMD-0262.json",
				DensityCacheLayout.emdbMapInfoFile(cacheRoot, "EMD-0262"));
		EmdbEntryResolver resolver = new EmdbEntryResolver(cacheRoot);
		resolver.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		EmdbEntryInfo info = resolver.getEntryInfo("EMD-0262");

		DensityMapResult result = new DensityMapResult(new File("x.map"), DensityMapSource.EMDB_MAP,
				DensityFileFormat.CCP4_GZ, DensityMapKind.EM, new PdbId("6hu9"), "EMD-0262", "u", false,
				info.getRecommendedContourLevel(), info.getSigma());

		assertEquals(info.getRecommendedContourLevel() / info.getSigma(), result.getContourInSigma(), 1e-9);
	}

	/** A cached mapping is used whatever its age when downloads are off. */
	@Test
	public void localOnlyUsesACachedMappingWithoutAskingAnyServer() throws IOException {
		File mappingFile = DensityCacheLayout.emdbMappingFile(cacheRoot, new PdbId("6hu9"));
		mappingFile.getParentFile().mkdirs();
		Files.write(mappingFile.toPath(),
				("pdbId=6hu9\nemdbIds=EMD-0262\ncontourLevel=0.0263\nretrieved=1999-01-01T00:00:00Z\n")
						.getBytes(StandardCharsets.UTF_8));

		EmdbEntryResolver resolver = new EmdbEntryResolver(cacheRoot);
		resolver.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		// point every template at a dead port, so a lookup would fail loudly
		EmdbEntryResolver.setSearchUrlTemplate(DEAD_SERVER + "{pdbid_lc}");
		EmdbEntryResolver.setRcsbEntryUrlTemplate(DEAD_SERVER + "{pdbid_lc}");

		List<String> ids = resolver.getEmdbIds(new PdbId("6hu9"));
		assertEquals(1, ids.size());
		assertEquals("EMD-0262", ids.get(0));
	}

	@Test
	public void localOnlyReturnsNothingWhenNothingIsCached() {
		EmdbEntryResolver resolver = new EmdbEntryResolver(cacheRoot);
		resolver.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		EmdbEntryResolver.setSearchUrlTemplate(DEAD_SERVER + "{pdbid_lc}");
		EmdbEntryResolver.setRcsbEntryUrlTemplate(DEAD_SERVER + "{pdbid_lc}");

		assertTrue(resolver.getEmdbIds(new PdbId("6hu9")).isEmpty());
	}

	/**
	 * A cached map is served without any network access at all, and its metadata is
	 * rebuilt from the sidecar rather than from a server.
	 */
	@Test
	public void localOnlyServesACachedMapOffline() throws IOException {
		DensityMapCache cache = new DensityMapCache(cacheRoot.getAbsolutePath());
		cache.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		PdbeCcp4MapProvider.setServerBaseUrl(DEAD_SERVER);
		cache.setSourceEnabled(DensityMapSource.RCSB_VOLUME_SERVER, false);
		cache.setSourceEnabled(DensityMapSource.PDBE_VOLUME_SERVER, false);

		File cached = DensityCacheLayout.pdbMapFile(cacheRoot, new PdbId("1cbs"),
				DensityMapKind.TWO_FO_FC, DensityMapSource.PDBE_CCP4, DensityFileFormat.CCP4, null);
		cached.getParentFile().mkdirs();
		Files.write(cached.toPath(), fakeCcp4());

		DensityMapResult result = cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);

		assertEquals(DensityMapSource.PDBE_CCP4, result.getSource());
		assertTrue(result.isFromCache());
		assertEquals(cached, result.getFile());
		assertTrue("a sidecar should have been written for the recovered file",
				DensityMapResult.metaFileFor(cached).isFile());
	}

	@Test
	public void localOnlyRefusesWhenNothingIsCached() {
		DensityMapCache cache = new DensityMapCache(cacheRoot.getAbsolutePath());
		cache.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		PdbeCcp4MapProvider.setServerBaseUrl(DEAD_SERVER);

		try {
			cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);
			fail("expected NoDensityMapException");
		} catch (NoDensityMapException e) {
			assertFalse(e.getAttempts().isEmpty());
		} catch (IOException e) {
			fail("LOCAL_ONLY must not attempt any connection, but got " + e);
		}
	}

	private static byte[] fakeCcp4() {
		byte[] bytes = new byte[4096];
		byte[] stamp = Ccp4Header.MAP_STAMP.getBytes(StandardCharsets.US_ASCII);
		System.arraycopy(stamp, 0, bytes, Ccp4Header.MAP_STAMP_OFFSET, stamp.length);
		return bytes;
	}
}
