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
import java.net.URL;
import java.net.URLConnection;
import java.util.Date;

import org.biojava.nbio.core.util.FileDownloadUtils;
import org.biojava.nbio.core.util.HttpStatusException;
import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.io.LocalPDBDirectory;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Shared machinery for the concrete density map providers: cache lookup,
 * download with validation, size limiting and format sanity checking.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public abstract class AbstractDensityMapProvider implements DensityMapProvider {

	private static final Logger logger = LoggerFactory.getLogger(AbstractDensityMapProvider.class);

	/**
	 * Anything smaller than this is not a usable map. Density servers answer with a
	 * short body rather than an error status when they have nothing, so a size
	 * floor is a necessary part of deciding whether a download succeeded.
	 */
	public static final long MIN_DENSITY_FILE_SIZE = 1024L;

	private File cacheRoot;
	private FetchBehavior fetchBehavior = FetchBehavior.FETCH_FILES;
	private long maxDownloadBytes = DensityMapCache.DEFAULT_MAX_DOWNLOAD_BYTES;

	/**
	 * @param cacheRoot the BioJava cache directory that the density directory sits in
	 */
	protected AbstractDensityMapProvider(File cacheRoot) {
		this.cacheRoot = cacheRoot;
	}

	/** @return the BioJava cache directory */
	public File getCacheRoot() {
		return cacheRoot;
	}

	/** @param cacheRoot the BioJava cache directory */
	public void setCacheRoot(File cacheRoot) {
		this.cacheRoot = cacheRoot;
	}

	/** @return how aggressively this provider re-fetches */
	public FetchBehavior getFetchBehavior() {
		return fetchBehavior;
	}

	/** @param fetchBehavior how aggressively to re-fetch */
	public void setFetchBehavior(FetchBehavior fetchBehavior) {
		this.fetchBehavior = fetchBehavior == null ? FetchBehavior.FETCH_FILES : fetchBehavior;
	}

	/** @return the download size limit in bytes, or 0 for no limit */
	public long getMaxDownloadBytes() {
		return maxDownloadBytes;
	}

	/** @param maxDownloadBytes the download size limit in bytes, or 0 for no limit */
	public void setMaxDownloadBytes(long maxDownloadBytes) {
		this.maxDownloadBytes = maxDownloadBytes;
	}

	@Override
	public boolean supports(DensityMapKind kind) {
		return kind != null && kind != DensityMapKind.AUTO;
	}

	/**
	 * The identifier spelling to put into a URL: four characters and lower case
	 * where the entry has a short form.
	 *
	 * @param pdbId the entry
	 * @return the identifier for use in a URL
	 */
	protected String urlId(PdbId pdbId) {
		return DensityCacheLayout.shortIdOrFull(pdbId).toLowerCase();
	}

	/**
	 * The effective fetch behaviour for a request, preferring the request's own
	 * setting over this provider's.
	 *
	 * @param request the request
	 * @return the behaviour to apply
	 */
	protected FetchBehavior effectiveFetchBehavior(DensityMapRequest request) {
		return request.getFetchBehavior() == null ? fetchBehavior : request.getFetchBehavior();
	}

	/**
	 * The effective size limit for a request, preferring the request's own setting.
	 *
	 * @param request the request
	 * @return the limit in bytes, or 0 for no limit
	 */
	protected long effectiveMaxBytes(DensityMapRequest request) {
		return request.getMaxDownloadBytes() < 0 ? maxDownloadBytes : request.getMaxDownloadBytes();
	}

	/**
	 * The cache directory to use for a request, preferring the request's override.
	 *
	 * @param request the request
	 * @return the cache directory
	 */
	protected File effectiveCacheRoot(DensityMapRequest request) {
		return request.getCacheDir() == null ? cacheRoot : request.getCacheDir();
	}

	/**
	 * Obtains a map, serving it from the cache when the fetch behaviour permits and
	 * downloading it otherwise.
	 *
	 * @param request what was asked for
	 * @param url where to fetch from
	 * @param target where to cache it
	 * @param kind the concrete kind of map being fetched
	 * @param emdbId the EMDB entry, if this is an EM map; otherwise <code>null</code>
	 * @param contourLevel the author-recommended contour level, or <code>null</code>
	 * @param sigma the map RMS deviation, or <code>null</code>
	 * @return the result, or <code>null</code> if the source has nothing for this entry
	 * @throws IOException on transport failure
	 */
	protected DensityMapResult obtain(DensityMapRequest request, URL url, File target, DensityMapKind kind,
			String emdbId, Double contourLevel, Double sigma) throws IOException {

		FetchBehavior behavior = effectiveFetchBehavior(request);

		if (isCacheUsable(target, url, behavior)) {
			// The kind always comes from the request, never from the sidecar: a single
			// cached file can hold more than one kind of map, so the sidecar's kind
			// records what was asked for first, not what the caller wants now.
			DensityMapResult cached = DensityMapResult.readMeta(target);
			Double cachedContour = contourLevel != null || cached == null ? contourLevel
					: cached.getRecommendedContourLevel();
			Double cachedSigma = sigma != null || cached == null ? sigma : cached.getSigma();
			DensityMapResult result = new DensityMapResult(target, getSource(), getFormat(), kind,
					request.getPdbId(), emdbId, url.toString(), true, cachedContour, cachedSigma);
			if (cached == null) {
				// The file is good but its description was lost; write one rather than
				// downloading several megabytes again.
				result.writeMeta();
			}
			return result;
		}

		if (behavior == FetchBehavior.LOCAL_ONLY) {
			throw new HttpStatusException(404, url.toString(),
					"not in the local cache and downloads are disabled (FetchBehavior.LOCAL_ONLY)");
		}

		File dir = target.getAbsoluteFile().getParentFile();
		if (!dir.isDirectory() && !dir.mkdirs()) {
			throw new IOException("Could not create density cache directory " + dir);
		}

		enforceSizeLimit(url, effectiveMaxBytes(request), request);

		logger.info("Fetching {} density map for {} from {}", kind,
				request.getPdbId() == null ? emdbId : request.getPdbId().getId(), url);
		FileDownloadUtils.downloadFileWithValidation(url, target, null, FileDownloadUtils.Hash.UNKNOWN,
				FileDownloadUtils.ETagPolicy.USE_IF_HEX_DIGEST);

		if (!isPlausibleMap(target)) {
			// Do not leave a bad file behind to be picked up as a cache hit later.
			deleteWithSidecars(target);
			throw new IOException("The content downloaded from " + url + " is not a usable "
					+ getFormat() + " density map.");
		}

		DensityMapResult result = new DensityMapResult(target, getSource(), getFormat(), kind,
				request.getPdbId(), emdbId, url.toString(), false, contourLevel, sigma);
		result.writeMeta();
		return result;
	}

	/**
	 * Whether a cached file may be used as-is.
	 *
	 * @param target the cached file
	 * @param url where it came from
	 * @param behavior the fetch behaviour in force
	 * @return <code>true</code> if the cached file should be served
	 */
	protected boolean isCacheUsable(File target, URL url, FetchBehavior behavior) {
		if (behavior == FetchBehavior.FORCE_DOWNLOAD) {
			return false;
		}
		if (!target.isFile() || target.length() < MIN_DENSITY_FILE_SIZE) {
			return false;
		}
		if (!FileDownloadUtils.validateFile(target)) {
			logger.info("Cached density map [{}] failed validation and will be re-downloaded.", target);
			return false;
		}
		if (!isPlausibleMap(target)) {
			logger.info("Cached file [{}] is not a usable {} map and will be re-downloaded.", target, getFormat());
			return false;
		}
		if (behavior == FetchBehavior.LOCAL_ONLY) {
			return true;
		}
		if (behavior == FetchBehavior.FETCH_IF_OUTDATED) {
			Date serverDate = LocalPDBDirectory.getLastModifiedTime(url);
			if (serverDate == null) {
				// Density servers generate their responses on the fly and send no
				// Last-Modified header. Treating that as "outdated" would re-download
				// on every single call, so an unknown timestamp keeps the cache.
				logger.debug("No server timestamp for {}; keeping the cached copy.", url);
				return true;
			}
			return target.lastModified() >= serverDate.getTime();
		}
		// FETCH_FILES, and FETCH_REMEDIATED which has no meaning for maps.
		return true;
	}

	/**
	 * Checks that a file looks like the format this provider delivers. Only the
	 * CCP4/MRC formats carry a recognisable stamp; the others are accepted on size
	 * alone.
	 *
	 * @param file the file to check
	 * @return <code>true</code> if the file is plausibly a map of this format
	 */
	protected boolean isPlausibleMap(File file) {
		if (file == null || file.length() < MIN_DENSITY_FILE_SIZE) {
			return false;
		}
		DensityFileFormat format = getFormat();
		if (format == DensityFileFormat.CCP4 || format == DensityFileFormat.CCP4_GZ) {
			return Ccp4Header.isCcp4Quietly(file);
		}
		return true;
	}

	/**
	 * Refuses a download whose declared size exceeds the limit, before any of the
	 * body is transferred.
	 *
	 * @param url the resource
	 * @param maxBytes the limit, or 0 for no limit
	 * @param request the request being served, used for reporting
	 * @throws DensityMapTooLargeException if the resource is too large
	 */
	protected void enforceSizeLimit(URL url, long maxBytes, DensityMapRequest request)
			throws DensityMapTooLargeException {
		if (maxBytes <= 0) {
			return;
		}
		long size = declaredSize(url);
		if (size > maxBytes) {
			reportTooLarge(url, size, maxBytes, request);
			throw new DensityMapTooLargeException(url.toString(), size, maxBytes);
		}
	}

	/**
	 * Announces that a map was skipped for being too large.
	 * <p>
	 * This is reported through the logger and on both standard output and standard
	 * error. A user who asked for a map and silently got a coarser one from another
	 * source deserves to be told why, and a library log configuration that discards
	 * warnings should not be able to hide it.
	 *
	 * @param url the resource that was skipped
	 * @param size its size in bytes
	 * @param maxBytes the limit in bytes
	 * @param request the request being served
	 */
	protected void reportTooLarge(URL url, long size, long maxBytes, DensityMapRequest request) {
		String entry = request.getPdbId() != null ? request.getPdbId().getId()
				: (request.getEmdbId() != null ? request.getEmdbId() : "?");
		String message = String.format(
				"BioJava density: skipping %s for %s (%,d bytes exceeds the %,d byte limit). "
				+ "Trying a smaller representation from another source; "
				+ "raise DensityMapCache.setMaxDownloadBytes() to allow it.",
				url, entry, size, maxBytes);
		logger.warn(message);
		System.out.println(message);
		System.err.println(message);
	}

	/**
	 * The <code>Content-Length</code> a server declares for a resource.
	 *
	 * @param url the resource
	 * @return the size in bytes, or a negative value if the server did not say
	 */
	protected long declaredSize(URL url) {
		try {
			URLConnection connection = FileDownloadUtils.prepareURLConnection(url.toString(), 30000);
			if (connection instanceof java.net.HttpURLConnection) {
				((java.net.HttpURLConnection) connection).setRequestMethod("HEAD");
			}
			connection.connect();
			try {
				return connection.getContentLengthLong();
			} finally {
				if (connection instanceof java.net.HttpURLConnection) {
					((java.net.HttpURLConnection) connection).disconnect();
				}
			}
		} catch (IOException e) {
			logger.debug("Could not determine the size of {}: {}", url, e.getMessage());
			return -1;
		}
	}

	/**
	 * Deletes a cached map along with its validation and metadata sidecars.
	 *
	 * @param target the cached map
	 */
	protected void deleteWithSidecars(File target) {
		File dir = target.getAbsoluteFile().getParentFile();
		if (dir == null) {
			return;
		}
		File[] siblings = dir.listFiles((d, name) -> name.startsWith(target.getName()));
		if (siblings == null) {
			return;
		}
		for (File f : siblings) {
			if (!f.delete()) {
				logger.debug("Could not delete [{}]", f);
			}
		}
	}
}
