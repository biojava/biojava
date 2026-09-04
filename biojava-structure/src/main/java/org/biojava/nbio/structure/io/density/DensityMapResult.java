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
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.time.Instant;
import java.util.Properties;

import org.biojava.nbio.structure.PdbId;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A density map that was successfully obtained, together with everything a
 * caller needs to know in order to use it.
 * <p>
 * Which source answered matters to the caller and is therefore part of the
 * result: the file might be a full CCP4 grid, a downsampled volume slice, or
 * &mdash; if non-renderable formats were allowed &mdash; a set of structure
 * factors that must be Fourier-transformed before anything can be drawn. See
 * {@link #isRenderable()}.
 * <p>
 * Every field is persisted to a <code>.meta</code> sidecar beside the cached
 * file, so a result can be reconstructed later without contacting any server.
 * That is what lets the cache serve
 * {@link org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior#LOCAL_ONLY}
 * requests fully offline.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class DensityMapResult {

	private static final Logger logger = LoggerFactory.getLogger(DensityMapResult.class);

	/** Extension of the metadata sidecar written beside every cached map. */
	public static final String META_EXT = ".meta";

	private final File file;
	private final DensityMapSource source;
	private final DensityFileFormat format;
	private final DensityMapKind kind;
	private final PdbId pdbId;
	private final String emdbId;
	private final String sourceUrl;
	private final boolean fromCache;
	private final Double recommendedContourLevel;
	private final Double sigma;

	/**
	 * @param file the cached map file
	 * @param source which service supplied it
	 * @param format the file format
	 * @param kind what the values mean; never {@link DensityMapKind#AUTO}
	 * @param pdbId the PDB entry, may be <code>null</code>
	 * @param emdbId the EMDB entry, may be <code>null</code>
	 * @param sourceUrl the URL it came from
	 * @param fromCache whether it was already on disk rather than freshly downloaded
	 * @param recommendedContourLevel the author-recommended contour level in absolute
	 *        map units, or <code>null</code> if unknown
	 * @param sigma the RMS deviation of the map, or <code>null</code> if unknown
	 */
	public DensityMapResult(File file, DensityMapSource source, DensityFileFormat format, DensityMapKind kind,
			PdbId pdbId, String emdbId, String sourceUrl, boolean fromCache,
			Double recommendedContourLevel, Double sigma) {
		this.file = file;
		this.source = source;
		this.format = format;
		this.kind = kind;
		this.pdbId = pdbId;
		this.emdbId = emdbId;
		this.sourceUrl = sourceUrl;
		this.fromCache = fromCache;
		this.recommendedContourLevel = recommendedContourLevel;
		this.sigma = sigma;
	}

	/** @return the cached map file */
	public File getFile() {
		return file;
	}

	/** @return which service supplied the map */
	public DensityMapSource getSource() {
		return source;
	}

	/** @return the file format */
	public DensityFileFormat getFormat() {
		return format;
	}

	/** @return what the map values mean; never {@link DensityMapKind#AUTO} */
	public DensityMapKind getKind() {
		return kind;
	}

	/** @return the PDB entry, or <code>null</code> */
	public PdbId getPdbId() {
		return pdbId;
	}

	/** @return the EMDB entry, or <code>null</code> */
	public String getEmdbId() {
		return emdbId;
	}

	/** @return the URL the map was fetched from */
	public String getSourceUrl() {
		return sourceUrl;
	}

	/** @return <code>true</code> if the file was already cached rather than downloaded now */
	public boolean isFromCache() {
		return fromCache;
	}

	/**
	 * Whether a viewer can contour this file as it stands.
	 *
	 * @return <code>false</code> for map coefficients, which need a Fourier
	 *         transform first
	 * @see DensityFileFormat#isJmolLoadable()
	 */
	public boolean isRenderable() {
		return format != null && format.isJmolLoadable();
	}

	/**
	 * The contour level recommended by the depositors, in absolute map units.
	 * Normally only available for EMDB maps, where it is the conventional way to
	 * contour rather than a multiple of sigma.
	 *
	 * @return the level, or <code>null</code> if unknown
	 */
	public Double getRecommendedContourLevel() {
		return recommendedContourLevel;
	}

	/**
	 * @return the RMS deviation of the map values, or <code>null</code> if unknown
	 */
	public Double getSigma() {
		return sigma;
	}

	/**
	 * The recommended contour expressed in multiples of sigma, for viewers that
	 * prefer to work that way.
	 *
	 * @return the level divided by sigma, or <code>null</code> if either is unknown
	 *         or sigma is zero
	 */
	public Double getContourInSigma() {
		if (recommendedContourLevel == null || sigma == null || sigma == 0.0) {
			return null;
		}
		return recommendedContourLevel / sigma;
	}

	/** @return the size of the cached file in bytes, or 0 if it is missing */
	public long getFileSizeBytes() {
		return file == null ? 0 : file.length();
	}

	/**
	 * The metadata sidecar file for a cached map.
	 *
	 * @param mapFile the cached map
	 * @return the sidecar path, which need not exist
	 */
	public static File metaFileFor(File mapFile) {
		return new File(mapFile.getAbsoluteFile().getParentFile(), mapFile.getName() + META_EXT);
	}

	/**
	 * Writes this result's metadata beside the cached file, so that it can be
	 * reconstructed later without any network access.
	 */
	public void writeMeta() {
		Properties p = new Properties();
		put(p, "source", source);
		put(p, "format", format);
		put(p, "kind", kind);
		put(p, "pdbId", pdbId == null ? null : pdbId.getId());
		put(p, "emdbId", emdbId);
		put(p, "url", sourceUrl);
		put(p, "downloaded", Instant.now().toString());
		put(p, "bytes", getFileSizeBytes());
		put(p, "contourLevel", recommendedContourLevel);
		put(p, "sigma", sigma);
		File meta = metaFileFor(file);
		try (OutputStream out = Files.newOutputStream(meta.toPath())) {
			p.store(out, "BioJava density map metadata");
		} catch (IOException e) {
			// Losing the sidecar costs us the offline description, not the map itself.
			logger.warn("Could not write density metadata [{}]: {}", meta, e.getMessage());
		}
	}

	/**
	 * Reconstructs a result from a cached file and its metadata sidecar.
	 *
	 * @param mapFile the cached map
	 * @return the reconstructed result, or <code>null</code> if the sidecar is
	 *         missing or unusable
	 */
	public static DensityMapResult readMeta(File mapFile) {
		File meta = metaFileFor(mapFile);
		if (!meta.isFile()) {
			return null;
		}
		Properties p = new Properties();
		try (InputStream in = Files.newInputStream(meta.toPath())) {
			p.load(in);
		} catch (IOException e) {
			logger.warn("Could not read density metadata [{}]: {}", meta, e.getMessage());
			return null;
		}
		try {
			String pdb = p.getProperty("pdbId");
			return new DensityMapResult(mapFile,
					DensityMapSource.valueOf(p.getProperty("source")),
					DensityFileFormat.valueOf(p.getProperty("format")),
					DensityMapKind.valueOf(p.getProperty("kind")),
					pdb == null || pdb.isEmpty() ? null : new PdbId(pdb),
					emptyToNull(p.getProperty("emdbId")),
					p.getProperty("url"),
					true,
					parseDouble(p.getProperty("contourLevel")),
					parseDouble(p.getProperty("sigma")));
		} catch (RuntimeException e) {
			logger.warn("Density metadata [{}] is unusable: {}", meta, e.getMessage());
			return null;
		}
	}

	private static void put(Properties p, String key, Object value) {
		p.setProperty(key, value == null ? "" : String.valueOf(value));
	}

	private static String emptyToNull(String s) {
		return s == null || s.isEmpty() ? null : s;
	}

	private static Double parseDouble(String s) {
		if (s == null || s.isEmpty()) {
			return null;
		}
		try {
			return Double.valueOf(s);
		} catch (NumberFormatException e) {
			return null;
		}
	}

	@Override
	public String toString() {
		return String.format("%s %s map from %s (%s, %d bytes)%s",
				pdbId == null ? emdbId : pdbId.getId(), kind, source, format, getFileSizeBytes(),
				isRenderable() ? "" : " [needs an FFT before it can be displayed]");
	}
}
