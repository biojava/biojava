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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
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
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

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

	@BeforeEach
	public void setUp() throws IOException {
		cacheRoot = Files.createTempDirectory("bj-density-it").toFile();
		cache = new DensityMapCache(cacheRoot.getAbsolutePath());
	}

	@AfterEach
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
		assertTrue(DensityMapResult.metaFileFor(result.getFile()).isFile(),
				"a .meta sidecar makes the result reconstructible offline");

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
		assertFalse(twoFoFc.getFile().equals(foFc.getFile()), "the difference map needs its own file name");
		assertTrue(foFc.getFile().getName().contains("&diff=1"),
				"the marker has to be in the name for Jmol to select the FO-FC block");
		assertEquals(twoFoFc.getFileSizeBytes(), foFc.getFileSizeBytes(),
				"both names must address the same bytes");
	}

	/** PDBe serves real CCP4 files, which the header check should recognise. */
	@Test
	public void pdbeServesAGenuineCcp4Map() throws IOException {
		cache.setSourceChain(DensityMapKind.TWO_FO_FC, Arrays.asList(DensityMapSource.PDBE_CCP4));
		DensityMapResult result = cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);

		assertEquals(DensityMapSource.PDBE_CCP4, result.getSource());
		assertEquals(DensityFileFormat.CCP4, result.getFormat());
		assertTrue(Ccp4Header.isCcp4(result.getFile()), "the CCP4 stamp should be present at byte 208");
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
		assertNotNull(result.getRecommendedContourLevel(),
				"EM maps need the author contour level to be displayed properly");
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
	 * Coefficients must arrive intact and be verifiable afterwards. Whether that
	 * verification includes a cryptographic digest depends on the server.
	 * <p>
	 * The divided archive paths on files.wwpdb.org and files.rcsb.org return the content
	 * MD5 as the ETag. The flat /validation/download/ endpoint this provider now uses
	 * returns neither an ETag nor a Content-Length, and neither does the beta archive on
	 * those two hosts, so no digest can be recorded there. The size sidecar is written
	 * from the bytes actually read, so it exists either way.
	 * <p>
	 * The digest is therefore asserted when the server offered one and skipped when it
	 * did not, rather than being required: requiring it would fail against the endpoint
	 * we use, and hard-coding the divided path would only work until the archive
	 * transition in July 2027.
	 */
	@Test
	public void mapCoefficientsArriveIntactAndVerifiable() throws IOException {
		cache.setSourceEnabled(DensityMapSource.WWPDB_MAP_COEFFICIENTS, true);
		cache.setSourceChain(DensityMapKind.TWO_FO_FC, Arrays.asList(DensityMapSource.WWPDB_MAP_COEFFICIENTS));

		DensityMapResult result = cache.getDensityMap(DensityMapRequest.builder(new PdbId("1cbs"))
				.kind(DensityMapKind.TWO_FO_FC)
				.allowNonRenderableFormats(true)
				.build());

		assertEquals(DensityMapSource.WWPDB_MAP_COEFFICIENTS, result.getSource());
		assertFalse(result.isRenderable(),
				"structure factors are not a map and must not claim to be renderable");

		// written from the observed byte count, so it is present whether or not the
		// server declared a length
		assertTrue(FileDownloadUtils.validateFile(result.getFile()),
				"a freshly downloaded file must validate against its own sidecars");

		File hashFile = new File(result.getFile().getParentFile(), result.getFile().getName() + ".hash_MD5");
		if (hashFile.isFile()) {
			String recorded = new String(Files.readAllBytes(hashFile.toPath()), StandardCharsets.UTF_8).trim();
			assertTrue(FileDownloadUtils.verifyHash(result.getFile(), FileDownloadUtils.Hash.MD5, recorded),
					"the recorded MD5 must match the file it describes");
		} else {
			System.out.println("No MD5 recorded for " + result.getSourceUrl()
					+ " - the server offered no usable ETag. Size validation still applies.");
		}

		// corrupt it and confirm validation actually catches it
		Files.write(result.getFile().toPath(), new byte[] {0, 1, 2, 3});
		assertFalse(FileDownloadUtils.validateFile(result.getFile()),
				"a truncated file must not validate");
	}
}
