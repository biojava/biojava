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
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.net.SocketTimeoutException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.core.util.HttpStatusException;
import org.biojava.nbio.structure.PdbId;
import org.junit.Before;
import org.junit.Test;

/**
 * The fallback chain, exercised with stub providers so that no server is
 * contacted.
 * <p>
 * The behaviour that matters most here is the difference between "this source
 * has nothing for the entry" and "this source could not be reached". The first
 * must move on to the next source; the second must abort, because reporting a
 * network outage as "no density exists" would be actively misleading.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class TestDensityFallbackChain {

	private DensityMapCache cache;
	private List<DensityMapSource> called;

	@Before
	public void setUp() {
		cache = new DensityMapCache(System.getProperty("java.io.tmpdir"));
		called = new ArrayList<>();
	}

	/** A provider that behaves however the test needs it to. */
	private class StubProvider implements DensityMapProvider {

		private final DensityMapSource source;
		private final DensityFileFormat format;
		private final IOException failure;
		private final boolean returnsNull;

		StubProvider(DensityMapSource source, DensityFileFormat format, IOException failure, boolean returnsNull) {
			this.source = source;
			this.format = format;
			this.failure = failure;
			this.returnsNull = returnsNull;
		}

		@Override
		public DensityMapSource getSource() {
			return source;
		}

		@Override
		public DensityFileFormat getFormat() {
			return format;
		}

		@Override
		public boolean supports(DensityMapKind kind) {
			return kind != DensityMapKind.AUTO;
		}

		@Override
		public DensityMapResult fetch(DensityMapRequest request) throws IOException {
			called.add(source);
			if (failure != null) {
				throw failure;
			}
			if (returnsNull) {
				return null;
			}
			return new DensityMapResult(new File("stub.map"), source, format, request.getKind(),
					request.getPdbId(), request.getEmdbId(), "stub://" + source, false, null, null);
		}
	}

	private StubProvider missing(DensityMapSource source) {
		return new StubProvider(source, DensityFileFormat.CCP4,
				new HttpStatusException(404, "stub://" + source, "Not Found"), false);
	}

	private StubProvider succeeds(DensityMapSource source) {
		return new StubProvider(source, DensityFileFormat.CCP4, null, false);
	}

	private void useOnly(DensityMapSource... sources) {
		cache.setSourceChain(DensityMapKind.TWO_FO_FC, Arrays.asList(sources));
		for (DensityMapSource s : DensityMapSource.values()) {
			cache.setSourceEnabled(s, Arrays.asList(sources).contains(s));
		}
	}

	@Test
	public void aMissingSourceFallsThroughToTheNext() throws IOException {
		cache.registerProvider(missing(DensityMapSource.RCSB_VOLUME_SERVER));
		cache.registerProvider(succeeds(DensityMapSource.PDBE_CCP4));
		useOnly(DensityMapSource.RCSB_VOLUME_SERVER, DensityMapSource.PDBE_CCP4);

		DensityMapResult result = cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);

		assertEquals(DensityMapSource.PDBE_CCP4, result.getSource());
		assertEquals(Arrays.asList(DensityMapSource.RCSB_VOLUME_SERVER, DensityMapSource.PDBE_CCP4), called);
	}

	@Test
	public void returningNullAlsoFallsThrough() throws IOException {
		cache.registerProvider(new StubProvider(DensityMapSource.RCSB_VOLUME_SERVER,
				DensityFileFormat.CCP4, null, true));
		cache.registerProvider(succeeds(DensityMapSource.PDBE_CCP4));
		useOnly(DensityMapSource.RCSB_VOLUME_SERVER, DensityMapSource.PDBE_CCP4);

		assertEquals(DensityMapSource.PDBE_CCP4,
				cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC).getSource());
	}

	@Test
	public void aTooLargeMapFallsThroughToASmallerSource() throws IOException {
		cache.registerProvider(new StubProvider(DensityMapSource.EMDB_MAP, DensityFileFormat.CCP4_GZ,
				new DensityMapTooLargeException("stub://big", 111426503L, 1024L), false));
		cache.registerProvider(succeeds(DensityMapSource.RCSB_VOLUME_SERVER));
		useOnly(DensityMapSource.EMDB_MAP, DensityMapSource.RCSB_VOLUME_SERVER);

		assertEquals(DensityMapSource.RCSB_VOLUME_SERVER,
				cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC).getSource());
	}

	/**
	 * A dropped connection is not evidence that the entry has no density, so the
	 * chain must stop rather than quietly try the rest and report "none available".
	 */
	@Test
	public void aTransportFailureAbortsTheChain() {
		cache.registerProvider(new StubProvider(DensityMapSource.RCSB_VOLUME_SERVER, DensityFileFormat.CCP4,
				new SocketTimeoutException("connection timed out"), false));
		cache.registerProvider(succeeds(DensityMapSource.PDBE_CCP4));
		useOnly(DensityMapSource.RCSB_VOLUME_SERVER, DensityMapSource.PDBE_CCP4);

		try {
			cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);
			fail("a transport failure should not be reported as a missing map");
		} catch (NoDensityMapException e) {
			fail("a transport failure must not be reported as NoDensityMapException");
		} catch (IOException expected) {
			assertEquals(Arrays.asList(DensityMapSource.RCSB_VOLUME_SERVER), called);
		}
	}

	/** A 5xx is a server problem, not an absent entry, so it must abort too. */
	@Test
	public void aServerErrorAbortsTheChain() {
		cache.registerProvider(new StubProvider(DensityMapSource.RCSB_VOLUME_SERVER, DensityFileFormat.CCP4,
				new HttpStatusException(503, "stub://x", "Service Unavailable"), false));
		cache.registerProvider(succeeds(DensityMapSource.PDBE_CCP4));
		useOnly(DensityMapSource.RCSB_VOLUME_SERVER, DensityMapSource.PDBE_CCP4);

		try {
			cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);
			fail("HTTP 503 should abort the chain");
		} catch (IOException expected) {
			assertTrue(expected instanceof HttpStatusException);
			assertEquals(503, ((HttpStatusException) expected).getStatusCode());
		}
	}

	@Test
	public void exhaustingEverySourceReportsWhatWasTried() {
		cache.registerProvider(missing(DensityMapSource.RCSB_VOLUME_SERVER));
		cache.registerProvider(missing(DensityMapSource.PDBE_CCP4));
		useOnly(DensityMapSource.RCSB_VOLUME_SERVER, DensityMapSource.PDBE_CCP4);

		try {
			cache.getDensityMap(new PdbId("4hhb"), DensityMapKind.TWO_FO_FC);
			fail("expected NoDensityMapException");
		} catch (NoDensityMapException e) {
			assertEquals(2, e.getAttempts().size());
			assertTrue(e.getAttempts().get(DensityMapSource.RCSB_VOLUME_SERVER).contains("404"));
			assertTrue(e.getMessage().contains("4HHB") || e.getMessage().contains("4hhb"));
		} catch (IOException e) {
			fail("expected NoDensityMapException but got " + e);
		}
	}

	@Test
	public void aDisabledSourceIsNotCalled() {
		cache.registerProvider(succeeds(DensityMapSource.PDBE_CCP4));
		useOnly(DensityMapSource.PDBE_CCP4);
		cache.setSourceEnabled(DensityMapSource.PDBE_CCP4, false);

		try {
			cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);
			fail("expected NoDensityMapException");
		} catch (NoDensityMapException e) {
			assertTrue(called.isEmpty());
			assertEquals("disabled", e.getAttempts().get(DensityMapSource.PDBE_CCP4));
		} catch (IOException e) {
			fail("expected NoDensityMapException but got " + e);
		}
	}

	/** Map coefficients are archival only, so they are off unless asked for. */
	@Test
	public void mapCoefficientsAreDisabledByDefault() {
		assertFalse(new DensityMapCache(System.getProperty("java.io.tmpdir"))
				.isSourceEnabled(DensityMapSource.WWPDB_MAP_COEFFICIENTS));
	}

	/**
	 * A viewer asks for renderable formats only; that must exclude the coefficients
	 * even when the source is enabled.
	 */
	@Test
	public void unrenderableFormatsAreSkippedWhenTheCallerCannotUseThem() {
		cache.registerProvider(new StubProvider(DensityMapSource.WWPDB_MAP_COEFFICIENTS,
				DensityFileFormat.MAP_COEFFICIENTS_CIF_GZ, null, false));
		useOnly(DensityMapSource.WWPDB_MAP_COEFFICIENTS);

		try {
			cache.getDensityMap(DensityMapRequest.builder(new PdbId("1cbs"))
					.kind(DensityMapKind.TWO_FO_FC)
					.allowNonRenderableFormats(false)
					.build());
			fail("expected NoDensityMapException");
		} catch (NoDensityMapException e) {
			assertTrue(called.isEmpty());
			assertTrue(e.getAttempts().get(DensityMapSource.WWPDB_MAP_COEFFICIENTS).contains("Fourier"));
		} catch (IOException e) {
			fail("expected NoDensityMapException but got " + e);
		}
	}

	@Test
	public void coefficientsAreUsedWhenTheCallerAllowsThem() throws IOException {
		cache.registerProvider(new StubProvider(DensityMapSource.WWPDB_MAP_COEFFICIENTS,
				DensityFileFormat.MAP_COEFFICIENTS_CIF_GZ, null, false));
		useOnly(DensityMapSource.WWPDB_MAP_COEFFICIENTS);

		DensityMapResult result = cache.getDensityMap(DensityMapRequest.builder(new PdbId("1cbs"))
				.kind(DensityMapKind.TWO_FO_FC)
				.allowNonRenderableFormats(true)
				.build());
		assertFalse(result.isRenderable());
	}

	/** AUTO tries the X-ray map first, then falls through to the EM one. */
	@Test
	public void autoFallsFromXrayToEm() {
		assertEquals(Arrays.asList(DensityMapKind.TWO_FO_FC, DensityMapKind.EM),
				DensityMapKind.AUTO.resolve());
		assertEquals(Arrays.asList(DensityMapKind.FO_FC), DensityMapKind.FO_FC.resolve());
	}

	/** X-ray and EM entries are served by deliberately different orders. */
	@Test
	public void defaultChainsPreferTheSmallestSource() {
		DensityMapCache fresh = new DensityMapCache(System.getProperty("java.io.tmpdir"));
		assertEquals(DensityMapSource.RCSB_VOLUME_SERVER,
				fresh.getSourceChain(DensityMapKind.TWO_FO_FC).get(0));
		assertEquals(DensityMapSource.RCSB_VOLUME_SERVER,
				fresh.getSourceChain(DensityMapKind.EM).get(0));
		// the full-resolution archive is the last resort for EM
		List<DensityMapSource> em = fresh.getSourceChain(DensityMapKind.EM);
		assertEquals(DensityMapSource.EMDB_MAP, em.get(em.size() - 1));
	}
}
