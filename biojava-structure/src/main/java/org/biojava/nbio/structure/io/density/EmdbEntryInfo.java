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

/**
 * The few facts about an EMDB entry's primary map that matter when fetching and
 * displaying it.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class EmdbEntryInfo {

	private final String emdbId;
	private final Double recommendedContourLevel;
	private final Double sigma;
	private final Long mapSizeBytes;

	/**
	 * @param emdbId the entry identifier, canonical form
	 * @param recommendedContourLevel the author-recommended contour level in absolute
	 *        map units, or <code>null</code>
	 * @param sigma the RMS deviation of the map values, or <code>null</code>
	 * @param mapSizeBytes the size of the primary map file, or <code>null</code>
	 */
	public EmdbEntryInfo(String emdbId, Double recommendedContourLevel, Double sigma, Long mapSizeBytes) {
		this.emdbId = emdbId;
		this.recommendedContourLevel = recommendedContourLevel;
		this.sigma = sigma;
		this.mapSizeBytes = mapSizeBytes;
	}

	/** @return the entry identifier in canonical <code>EMD-1234</code> form */
	public String getEmdbId() {
		return emdbId;
	}

	/**
	 * The contour level the depositors recommend, in absolute map units. EM maps
	 * are conventionally displayed at this level rather than at a multiple of
	 * sigma, so it is the right default for a viewer.
	 *
	 * @return the level, or <code>null</code> if the entry does not state one
	 */
	public Double getRecommendedContourLevel() {
		return recommendedContourLevel;
	}

	/**
	 * @return the RMS deviation of the map values, which converts between absolute
	 *         and sigma-relative contour levels, or <code>null</code> if unknown
	 */
	public Double getSigma() {
		return sigma;
	}

	/**
	 * @return the size of the full primary map in bytes, or <code>null</code> if
	 *         unknown. Used to decline a download before it starts.
	 */
	public Long getMapSizeBytes() {
		return mapSizeBytes;
	}

	@Override
	public String toString() {
		return String.format("%s[contour=%s, sigma=%s, %s bytes]",
				emdbId, recommendedContourLevel, sigma, mapSizeBytes);
	}
}
