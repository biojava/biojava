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

import java.io.IOException;

/**
 * Thrown when a density map exceeds the configured download size limit.
 * <p>
 * Cryo-EM primary maps in particular can be very large &mdash; hundreds of
 * megabytes is common and gigabyte maps exist &mdash; while a downsampled slice
 * of the same map from a density server is usually a fraction of a percent of
 * that size and quite adequate for display. Rather than failing, the fallback
 * chain treats this as "try the next source", which is exactly why the density
 * servers are tried ahead of the full-resolution archives.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class DensityMapTooLargeException extends IOException {

	private static final long serialVersionUID = 1L;

	private final long sizeBytes;
	private final long limitBytes;

	/**
	 * @param url the resource that was too large
	 * @param sizeBytes its size, or a negative value if only a lower bound is known
	 * @param limitBytes the configured limit
	 */
	public DensityMapTooLargeException(String url, long sizeBytes, long limitBytes) {
		super(String.format("%s is %s, which exceeds the %s download limit.",
				url, describe(sizeBytes), describe(limitBytes)));
		this.sizeBytes = sizeBytes;
		this.limitBytes = limitBytes;
	}

	private static String describe(long bytes) {
		if (bytes < 0) {
			return "of unknown size";
		}
		if (bytes < 1024) {
			return bytes + " B";
		}
		double value = bytes;
		String[] units = {"kB", "MB", "GB", "TB"};
		int unit = -1;
		while (value >= 1024 && unit < units.length - 1) {
			value /= 1024;
			unit++;
		}
		return String.format("%.1f %s", value, units[unit]);
	}

	/**
	 * @return the size of the resource in bytes, or a negative value if unknown
	 */
	public long getSizeBytes() {
		return sizeBytes;
	}

	/**
	 * @return the configured limit in bytes
	 */
	public long getLimitBytes() {
		return limitBytes;
	}
}
