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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import org.biojava.nbio.core.util.FileDownloadUtils;
import org.biojava.nbio.core.util.HttpStatusException;
import org.biojava.nbio.structure.ExperimentalTechnique;
import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Downloads and caches electron density and cryo-EM maps, trying several sources
 * in turn until one produces a map.
 * <p>
 * Cached files live under the BioJava cache directory (<code>PDB_CACHE_DIR</code>),
 * laid out as described in {@link DensityCacheLayout}. Typical use is simply:
 * <pre>
 * DensityMapCache cache = new DensityMapCache();
 * DensityMapResult map = cache.getDensityMap(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);
 * File file = map.getFile();
 * </pre>
 * <p>
 * <b>Source order.</b> Sources are tried smallest-adequate-first, because they
 * differ enormously in size for the same entry and the smallest is usually
 * perfectly adequate to look at. For 1cbs a density-server slice is roughly a
 * tenth the size of the equivalent pair of CCP4 files; for the cryo-EM entry
 * behind EMD-0262 it is a few hundred kilobytes against 106&nbsp;MB. Anything
 * that cannot be displayed at all is tried last, and
 * {@link DensityMapSource#WWPDB_MAP_COEFFICIENTS} is disabled altogether by
 * default for that reason. Override with {@link #setSourceChain(DensityMapKind,
 * java.util.List)} or {@link #setSourceEnabled(DensityMapSource, boolean)}.
 * <p>
 * <b>Failures.</b> A source that has nothing for an entry is skipped and the next
 * is tried; if every source is exhausted, {@link NoDensityMapException} is thrown
 * carrying the reason from each one, so a caller can explain what happened.
 * Genuine transport failures abort the chain instead, so that a network outage is
 * never reported as "this entry has no density".
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class DensityMapCache {

	private static final Logger logger = LoggerFactory.getLogger(DensityMapCache.class);

	/**
	 * Default ceiling on a single download, 256&nbsp;MiB. Exceeding it is not an
	 * error: the chain moves on to a source that offers a smaller representation.
	 * Set to 0 to remove the limit.
	 */
	public static final long DEFAULT_MAX_DOWNLOAD_BYTES = 256L * 1024 * 1024;

	/** Order in which sources are tried for X-ray and neutron entries. */
	public static final List<DensityMapSource> DEFAULT_XRAY_SOURCE_CHAIN = Collections.unmodifiableList(Arrays.asList(
			DensityMapSource.RCSB_VOLUME_SERVER,
			DensityMapSource.PDBE_CCP4,
			DensityMapSource.PDBE_VOLUME_SERVER,
			DensityMapSource.WWPDB_MAP_COEFFICIENTS));

	/** Order in which sources are tried for cryo-EM entries. */
	public static final List<DensityMapSource> DEFAULT_EM_SOURCE_CHAIN = Collections.unmodifiableList(Arrays.asList(
			DensityMapSource.RCSB_VOLUME_SERVER,
			DensityMapSource.PDBE_VOLUME_SERVER,
			DensityMapSource.EMDB_MAP));

	private static DensityMapCache instance;

	private File cacheRoot;
	private FetchBehavior fetchBehavior = FetchBehavior.FETCH_FILES;
	private long maxDownloadBytes = DEFAULT_MAX_DOWNLOAD_BYTES;
	private final Map<DensityMapSource, DensityMapProvider> providers = new EnumMap<>(DensityMapSource.class);
	private final Set<DensityMapSource> disabled = EnumSet.of(DensityMapSource.WWPDB_MAP_COEFFICIENTS);
	private List<DensityMapSource> xraySourceChain = DEFAULT_XRAY_SOURCE_CHAIN;
	private List<DensityMapSource> emSourceChain = DEFAULT_EM_SOURCE_CHAIN;
	private EmdbEntryResolver emdbResolver;

	/**
	 * Creates a cache using the standard BioJava cache directory, i.e.
	 * <code>PDB_CACHE_DIR</code> falling back to <code>PDB_DIR</code>.
	 */
	public DensityMapCache() {
		this(new UserConfiguration().getCacheFilePath());
	}

	/**
	 * @param cachePath the directory to cache under
	 */
	public DensityMapCache(String cachePath) {
		this.cacheRoot = new File(FileDownloadUtils.expandUserHome(cachePath));
		this.emdbResolver = new EmdbEntryResolver(cacheRoot);
		buildDefaultProviders();
	}

	/**
	 * A lazily created shared instance, for callers that do not want to manage one.
	 *
	 * @return the shared cache
	 */
	public static synchronized DensityMapCache getInstance() {
		if (instance == null) {
			instance = new DensityMapCache();
		}
		return instance;
	}

	private void buildDefaultProviders() {
		providers.clear();
		register(new VolumeServerProvider(cacheRoot, VolumeServerProvider.Host.RCSB));
		register(new VolumeServerProvider(cacheRoot, VolumeServerProvider.Host.PDBE));
		register(new PdbeCcp4MapProvider(cacheRoot));
		register(new EmdbMapProvider(cacheRoot, emdbResolver));
		register(new WwpdbMapCoefficientsProvider(cacheRoot));
		applySettingsToProviders();
	}

	private void register(DensityMapProvider provider) {
		providers.put(provider.getSource(), provider);
	}

	private void applySettingsToProviders() {
		for (DensityMapProvider p : providers.values()) {
			if (p instanceof AbstractDensityMapProvider) {
				AbstractDensityMapProvider a = (AbstractDensityMapProvider) p;
				a.setCacheRoot(cacheRoot);
				a.setFetchBehavior(fetchBehavior);
				a.setMaxDownloadBytes(maxDownloadBytes);
			}
		}
		emdbResolver.setCacheRoot(cacheRoot);
		emdbResolver.setFetchBehavior(fetchBehavior);
	}

	/** @return the directory maps are cached under */
	public String getCachePath() {
		return cacheRoot.getAbsolutePath();
	}

	/** @param cachePath the directory to cache under */
	public void setCachePath(String cachePath) {
		this.cacheRoot = new File(FileDownloadUtils.expandUserHome(cachePath));
		applySettingsToProviders();
	}

	/** @return how aggressively cached maps are re-fetched */
	public FetchBehavior getFetchBehavior() {
		return fetchBehavior;
	}

	/**
	 * @param fetchBehavior how aggressively to re-fetch. {@link FetchBehavior#FETCH_REMEDIATED}
	 *        behaves as {@link FetchBehavior#FETCH_FILES}: the 2011 remediation date
	 *        is a coordinate-file concept with no meaning for maps.
	 */
	public void setFetchBehavior(FetchBehavior fetchBehavior) {
		this.fetchBehavior = fetchBehavior == null ? FetchBehavior.FETCH_FILES : fetchBehavior;
		applySettingsToProviders();
	}

	/** @return the ceiling on a single download in bytes, or 0 for no limit */
	public long getMaxDownloadBytes() {
		return maxDownloadBytes;
	}

	/** @param maxDownloadBytes the ceiling on a single download in bytes, or 0 for no limit */
	public void setMaxDownloadBytes(long maxDownloadBytes) {
		this.maxDownloadBytes = maxDownloadBytes;
		applySettingsToProviders();
	}

	/**
	 * The order sources are tried in for a given kind of map.
	 *
	 * @param kind the kind of map
	 * @return the source order
	 */
	public List<DensityMapSource> getSourceChain(DensityMapKind kind) {
		return kind == DensityMapKind.EM ? emSourceChain : xraySourceChain;
	}

	/**
	 * Overrides the order sources are tried in.
	 *
	 * @param kind {@link DensityMapKind#EM} to set the cryo-EM order, anything else
	 *        to set the X-ray order
	 * @param chain the new order
	 */
	public void setSourceChain(DensityMapKind kind, List<DensityMapSource> chain) {
		List<DensityMapSource> copy = Collections.unmodifiableList(new ArrayList<>(chain));
		if (kind == DensityMapKind.EM) {
			emSourceChain = copy;
		} else {
			xraySourceChain = copy;
		}
	}

	/**
	 * @param source the source
	 * @return whether it will be tried
	 */
	public boolean isSourceEnabled(DensityMapSource source) {
		return !disabled.contains(source);
	}

	/**
	 * Enables or disables a source.
	 * <p>
	 * {@link DensityMapSource#WWPDB_MAP_COEFFICIENTS} is disabled by default
	 * because it delivers structure factors rather than a map; enable it explicitly
	 * if you want the archival form and are prepared to run a Fourier transform.
	 *
	 * @param source the source
	 * @param enabled whether to try it
	 */
	public void setSourceEnabled(DensityMapSource source, boolean enabled) {
		if (enabled) {
			disabled.remove(source);
		} else {
			disabled.add(source);
		}
	}

	/**
	 * Replaces the provider used for a source. Mainly a testing seam, but also the
	 * way to plug in a site-local mirror.
	 *
	 * @param provider the provider to use
	 */
	public void registerProvider(DensityMapProvider provider) {
		register(provider);
		applySettingsToProviders();
	}

	/** @return the resolver used to map PDB entries to EMDB entries */
	public EmdbEntryResolver getEmdbResolver() {
		return emdbResolver;
	}

	/** @param resolver the resolver used to map PDB entries to EMDB entries */
	public void setEmdbResolver(EmdbEntryResolver resolver) {
		this.emdbResolver = resolver;
		applySettingsToProviders();
	}

	/**
	 * Fetches a density map, using the cache where possible.
	 *
	 * @param pdbId the entry
	 * @param kind the kind of map, or {@link DensityMapKind#AUTO} to take whatever
	 *        the entry has
	 * @return the map
	 * @throws NoDensityMapException if no enabled source has a map for this entry
	 * @throws IOException on transport failure
	 */
	public DensityMapResult getDensityMap(PdbId pdbId, DensityMapKind kind) throws IOException {
		return getDensityMap(DensityMapRequest.builder(pdbId).kind(kind).build());
	}

	/**
	 * Fetches a density map for a PDB or EMDB identifier.
	 *
	 * @param id a PDB identifier, or an EMDB identifier such as <code>EMD-0262</code>
	 * @param kind the kind of map wanted
	 * @return the map
	 * @throws NoDensityMapException if no enabled source has a map for this entry
	 * @throws IOException on transport failure
	 */
	public DensityMapResult getDensityMap(String id, DensityMapKind kind) throws IOException {
		return getDensityMap(DensityMapRequest.builder(id).kind(kind).build());
	}

	/**
	 * Fetches a density map for an already-loaded structure.
	 * <p>
	 * Knowing the structure lets the experimental method be read directly rather
	 * than guessed, so a cryo-EM entry goes straight to the EM sources instead of
	 * trying the X-ray ones first. No extra network request is involved.
	 *
	 * @param structure the structure
	 * @param kind the kind of map wanted
	 * @return the map
	 * @throws NoDensityMapException if no enabled source has a map for this entry
	 * @throws IOException on transport failure
	 */
	public DensityMapResult getDensityMap(Structure structure, DensityMapKind kind) throws IOException {
		PdbId pdbId = structure.getPdbId();
		if (pdbId == null) {
			throw new IOException("The structure has no PDB ID, so no density map can be looked up for it.");
		}
		DensityMapRequest request = DensityMapRequest.builder(pdbId).kind(kind).build();
		return getDensityMap(request, kindOrderFor(structure, kind));
	}

	/**
	 * Fetches a density map.
	 *
	 * @param request what is wanted
	 * @return the map
	 * @throws NoDensityMapException if no enabled source has a map for this entry
	 * @throws IOException on transport failure
	 */
	public DensityMapResult getDensityMap(DensityMapRequest request) throws IOException {
		return getDensityMap(request, request.getKind().resolve());
	}

	/**
	 * As {@link #getDensityMap(PdbId, DensityMapKind)} but returning an empty
	 * {@link Optional} rather than throwing when nothing is available. Genuine
	 * transport failures are still logged and swallowed, so use this only where
	 * "no map" and "could not reach the server" need not be told apart.
	 *
	 * @param pdbId the entry
	 * @param kind the kind of map wanted
	 * @return the map, if one could be obtained
	 */
	public Optional<DensityMapResult> findDensityMap(PdbId pdbId, DensityMapKind kind) {
		try {
			return Optional.of(getDensityMap(pdbId, kind));
		} catch (NoDensityMapException e) {
			logger.debug("{}", e.getMessage());
			return Optional.empty();
		} catch (IOException e) {
			logger.warn("Could not fetch a density map for {}: {}", pdbId.getId(), e.getMessage());
			return Optional.empty();
		}
	}

	/**
	 * Fetches both halves of the conventional X-ray pair: the 2mFo-DFc map and the
	 * mFo-DFc difference map.
	 *
	 * @param pdbId the entry
	 * @return whichever of the two could be obtained, in that order; possibly empty
	 */
	public List<DensityMapResult> getDifferenceMapPair(PdbId pdbId) {
		List<DensityMapResult> results = new ArrayList<>(2);
		findDensityMap(pdbId, DensityMapKind.TWO_FO_FC).ifPresent(results::add);
		findDensityMap(pdbId, DensityMapKind.FO_FC).ifPresent(results::add);
		return results;
	}

	/**
	 * Looks a map up in the cache without contacting any server.
	 *
	 * @param pdbId the entry
	 * @param kind the kind of map
	 * @param source the source it would have come from
	 * @return the cached map, or <code>null</code> if it is not cached
	 */
	public DensityMapResult getCached(PdbId pdbId, DensityMapKind kind, DensityMapSource source) {
		DensityMapProvider provider = providers.get(source);
		if (provider == null) {
			return null;
		}
		File file = DensityCacheLayout.pdbMapFile(cacheRoot, pdbId, kind, source, provider.getFormat(), null);
		return file.isFile() ? DensityMapResult.readMeta(file) : null;
	}

	/**
	 * Removes every cached density file for an entry, including the sidecars.
	 *
	 * @param pdbId the entry
	 * @return the number of files deleted
	 */
	public int deleteDensityMaps(PdbId pdbId) {
		String id = DensityCacheLayout.shortIdOrFull(pdbId).toLowerCase();
		File dir = DensityCacheLayout.pdbMapFile(cacheRoot, pdbId, DensityMapKind.TWO_FO_FC,
				DensityMapSource.PDBE_CCP4, DensityFileFormat.CCP4, null).getParentFile();
		File[] files = dir == null ? null : dir.listFiles((d, name) -> name.startsWith(id + "_"));
		if (files == null) {
			return 0;
		}
		int deleted = 0;
		for (File f : files) {
			if (f.delete()) {
				deleted++;
			}
		}
		return deleted;
	}

	private DensityMapResult getDensityMap(DensityMapRequest request, List<DensityMapKind> kinds) throws IOException {
		Map<DensityMapSource, String> attempts = new LinkedHashMap<>();

		for (DensityMapKind kind : kinds) {
			DensityMapRequest kindRequest = request.withKind(kind);

			if (kind == DensityMapKind.EM && kindRequest.getEmdbId() == null) {
				String emdbId = resolveEmdbId(kindRequest);
				if (emdbId == null) {
					attempts.put(DensityMapSource.EMDB_MAP, "no associated EMDB entry");
					continue;
				}
				kindRequest = kindRequest.withEmdbId(emdbId);
			}

			List<DensityMapSource> chain = kindRequest.getSourceChain() != null
					? kindRequest.getSourceChain() : getSourceChain(kind);

			for (DensityMapSource source : chain) {
				if (!isSourceEnabled(source)) {
					attempts.put(source, "disabled");
					continue;
				}
				DensityMapProvider provider = providers.get(source);
				if (provider == null) {
					attempts.put(source, "no provider registered");
					continue;
				}
				if (!provider.supports(kind)) {
					continue; // structurally impossible; not worth reporting
				}
				if (!kindRequest.isAllowNonRenderableFormats() && !provider.getFormat().isJmolLoadable()) {
					attempts.put(source, "cannot be displayed without a Fourier transform");
					continue;
				}

				try {
					DensityMapResult result = provider.fetch(kindRequest);
					if (result != null) {
						return withContourLevel(result);
					}
					attempts.put(source, "no map for this entry");
				} catch (HttpStatusException e) {
					if (!e.isNotFound()) {
						// A server error or an authentication problem is a real failure,
						// not evidence that the entry has no density.
						throw e;
					}
					attempts.put(source, "HTTP " + e.getStatusCode());
				} catch (DensityMapTooLargeException e) {
					attempts.put(source, "too large (" + e.getSizeBytes() + " bytes)");
				} catch (IOException e) {
					throw new IOException("Failed to fetch a density map for "
							+ (request.getPdbId() == null ? request.getEmdbId() : request.getPdbId().getId())
							+ " from " + source + ": " + e.getMessage(), e);
				}
			}
		}

		throw new NoDensityMapException(request.getPdbId(), request.getKind(), attempts);
	}

	/**
	 * Fills in the author-recommended contour level for an EM map when the source
	 * that supplied the file did not know it.
	 * <p>
	 * A density server returns voxels and nothing else, but an EM map is
	 * conventionally displayed at the level its depositors chose rather than at a
	 * multiple of sigma, so a viewer needs that number whichever source the map came
	 * from. It costs one small metadata request, cached thereafter.
	 */
	private DensityMapResult withContourLevel(DensityMapResult result) {
		if (result.getKind() != DensityMapKind.EM || result.getRecommendedContourLevel() != null
				|| result.getEmdbId() == null || emdbResolver == null) {
			return result;
		}
		EmdbEntryInfo info = emdbResolver.getEntryInfo(result.getEmdbId());
		if (info == null || (info.getRecommendedContourLevel() == null && info.getSigma() == null)) {
			return result;
		}
		DensityMapResult enriched = new DensityMapResult(result.getFile(), result.getSource(), result.getFormat(),
				result.getKind(), result.getPdbId(), result.getEmdbId(), result.getSourceUrl(), result.isFromCache(),
				info.getRecommendedContourLevel(), info.getSigma());
		enriched.writeMeta();
		return enriched;
	}

	private String resolveEmdbId(DensityMapRequest request) {
		if (request.getPdbId() == null || emdbResolver == null) {
			return null;
		}
		List<String> ids = emdbResolver.getEmdbIds(request.getPdbId());
		return ids.isEmpty() ? null : ids.get(0);
	}

	/**
	 * Chooses which kinds to try, and in what order, using the structure's declared
	 * experimental method. Reading it from the structure costs nothing; note that
	 * the resolution field is deliberately not consulted, since BioJava parses it
	 * incorrectly for some cryo-EM entries (biojava/biojava#1000).
	 */
	private List<DensityMapKind> kindOrderFor(Structure structure, DensityMapKind kind) {
		if (kind != DensityMapKind.AUTO) {
			return kind.resolve();
		}
		Set<ExperimentalTechnique> techniques = structure.getPDBHeader() == null
				? null : structure.getPDBHeader().getExperimentalTechniques();
		if (techniques != null && techniques.contains(ExperimentalTechnique.ELECTRON_MICROSCOPY)
				&& !ExperimentalTechnique.isCrystallographic(techniques)) {
			return Arrays.asList(DensityMapKind.EM, DensityMapKind.TWO_FO_FC);
		}
		return kind.resolve();
	}
}
