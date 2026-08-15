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

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * The kind of density map being requested, i.e. what the values in the grid mean.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public enum DensityMapKind {

	/**
	 * The 2mFo-DFc "best" map: the electron density itself, showing how the model
	 * fits the experimental data. This is what is normally meant by "the electron
	 * density" of an X-ray structure.
	 */
	TWO_FO_FC("2fofc"),

	/**
	 * The mFo-DFc difference map, showing density that the model does not account
	 * for (positive) and model that has no density to support it (negative). It is
	 * conventionally displayed as a signed pair of surfaces.
	 */
	FO_FC("fofc"),

	/**
	 * The primary map of a cryo-EM or cryo-ET reconstruction. Strictly this is a
	 * Coulomb potential map rather than an electron density map, and it is
	 * conventionally contoured at an absolute author-recommended level rather than
	 * in multiples of sigma.
	 */
	EM("em"),

	/**
	 * Not a map in itself: asks for whichever map the entry actually has. Resolves
	 * to {@link #TWO_FO_FC} first and then {@link #EM}, which covers X-ray and
	 * cryo-EM entries without the caller having to know which it is holding.
	 */
	AUTO(null);

	private final String fileToken;

	DensityMapKind(String fileToken) {
		this.fileToken = fileToken;
	}

	/**
	 * The short token used to distinguish this kind in a cached file name.
	 *
	 * @return the token, or <code>null</code> for {@link #AUTO}, which is never
	 *         itself cached
	 */
	public String getFileToken() {
		return fileToken;
	}

	/**
	 * Expands this kind into the concrete kinds to try, in order.
	 *
	 * @return a single-element list for a concrete kind, or the X-ray-then-EM order
	 *         for {@link #AUTO}
	 */
	public List<DensityMapKind> resolve() {
		if (this == AUTO) {
			return Arrays.asList(TWO_FO_FC, EM);
		}
		return Collections.singletonList(this);
	}

	/**
	 * Whether this kind is conventionally displayed as a signed pair of surfaces,
	 * one positive and one negative.
	 *
	 * @return <code>true</code> only for {@link #FO_FC}
	 */
	public boolean isDifferenceMap() {
		return this == FO_FC;
	}
}
