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
 * The file format a density map was delivered in.
 * <p>
 * The important distinction here is {@link #isJmolLoadable()}. Most of these
 * formats are sampled grids that a viewer can contour directly; map
 * <i>coefficients</i> are not, they are structure factors that require a Fourier
 * transform first. Handing the latter to a viewer produces nothing at all, so
 * the difference has to be visible to callers.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public enum DensityFileFormat {

	/** A CCP4/MRC sampled map, as served pre-computed by PDBe. */
	CCP4(".ccp4", true),

	/** A gzipped CCP4/MRC map, the form EMDB distributes its primary maps in. */
	CCP4_GZ(".map.gz", true),

	/** A BinaryCIF volume slice from a Mol* density server. */
	BCIF_VOLUME(".bcif", true),

	/** A text CIF volume slice from a Mol* density server. */
	CIF_VOLUME(".cif", true),

	/**
	 * Structure-factor amplitudes and phases in mmCIF, as published with the wwPDB
	 * validation reports. <b>Not a density grid.</b> A Fourier transform (for
	 * example <code>gemmi sf2map</code>, or CCP4's <code>fft</code> after
	 * <code>cif2mtz</code>) is needed before these can be displayed.
	 */
	MAP_COEFFICIENTS_CIF_GZ(".cif.gz", false);

	private final String extension;
	private final boolean jmolLoadable;

	DensityFileFormat(String extension, boolean jmolLoadable) {
		this.extension = extension;
		this.jmolLoadable = jmolLoadable;
	}

	/**
	 * @return the file extension used for cached files of this format, including
	 *         the leading dot
	 */
	public String getExtension() {
		return extension;
	}

	/**
	 * Whether a viewer can contour this file as it stands.
	 *
	 * @return <code>false</code> for {@link #MAP_COEFFICIENTS_CIF_GZ}, which needs
	 *         an FFT first; <code>true</code> for the sampled grid formats
	 */
	public boolean isJmolLoadable() {
		return jmolLoadable;
	}

	/**
	 * Whether files of this format are gzip-compressed on disk. Jmol detects gzip
	 * from the magic bytes, so such files do not need decompressing before display.
	 *
	 * @return <code>true</code> for the compressed formats
	 */
	public boolean isCompressed() {
		return this == CCP4_GZ || this == MAP_COEFFICIENTS_CIF_GZ;
	}
}
