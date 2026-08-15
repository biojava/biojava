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
 * A remote service that density data can be fetched from.
 * <p>
 * These differ enormously in size for the same entry, which is why more than one
 * is supported. For PDB entry 1cbs a density-server slice at the coarsest detail
 * level is about 210&nbsp;kB and contains both the 2Fo-Fc and Fo-Fc maps, where
 * the two pre-computed CCP4 files come to about 2.1&nbsp;MB together. For cryo-EM
 * the gap is far wider: the primary map of EMD-0262 is about 106&nbsp;MB, against
 * roughly 480&nbsp;kB for the equivalent density-server slice.
 * <p>
 * Note that RCSB's own <code>edmaps.rcsb.org</code> service, which used to serve
 * DSN6 and MTZ files, was shut down in October 2024. The map coefficients
 * published with the wwPDB validation reports replaced it, but those are
 * structure factors rather than a sampled grid; see
 * {@link DensityFileFormat#MAP_COEFFICIENTS_CIF_GZ}.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public enum DensityMapSource {

	/**
	 * RCSB's Mol* density server at <code>maps.rcsb.org</code>, which serves
	 * downsampled BinaryCIF volume slices for both X-ray and EM entries.
	 */
	RCSB_VOLUME_SERVER("rcsbvs"),

	/**
	 * PDBe's Mol* density server, equivalent to {@link #RCSB_VOLUME_SERVER}.
	 */
	PDBE_VOLUME_SERVER("pdbevs"),

	/**
	 * PDBe's pre-computed full CCP4 maps. This is the source Jmol itself uses for
	 * its built-in map-loading shortcuts.
	 */
	PDBE_CCP4("pdbe"),

	/**
	 * The primary map of an EMDB entry, at full resolution. Can be very large.
	 */
	EMDB_MAP("emdb"),

	/**
	 * Map coefficients from the wwPDB validation reports. Archival only: these
	 * cannot be displayed without an FFT.
	 */
	WWPDB_MAP_COEFFICIENTS("wwpdb");

	private final String fileToken;

	DensityMapSource(String fileToken) {
		this.fileToken = fileToken;
	}

	/**
	 * The short token used to distinguish this source in a cached file name, so
	 * that maps of the same kind from different sources never collide.
	 *
	 * @return the token
	 */
	public String getFileToken() {
		return fileToken;
	}
}
