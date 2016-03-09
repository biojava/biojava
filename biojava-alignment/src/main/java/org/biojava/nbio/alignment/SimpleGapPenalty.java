/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on June 9, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment;

import java.io.Serializable;

import org.biojava.nbio.alignment.template.GapPenalty;

/**
 * Implements a data structure for the gap penalties used during a sequence alignment routine.
 *
 * @author Mark Chapman
 */
public class SimpleGapPenalty implements GapPenalty, Serializable {

	private static final long serialVersionUID = 3945671344135815456L;

	private static int dgop = 10, dgep = 1;

	/**
	 * Sets the default gap extension penalty.
	 *
	 * @param gep the default gap extension penalty
	 */
	public static void setDefaultExtensionPenalty(int gep) {
		dgep = gep;
	}

	/**
	 * Sets the default gap open penalty.
	 *
	 * @param gop the default gap open penalty
	 */
	public static void setDefaultOpenPenalty(int gop) {
		dgop = gop;
	}

	private GapPenalty.Type type;
	private int gop, gep;

	/**
	 * Creates a new set of gap penalties using the defaults.
	 */
	public SimpleGapPenalty() {
		this(dgop, dgep);
	}

	/**
	 * Creates a new set of gap penalties.
	 *
	 * @param gop the gap open penalty; should be nonnegative
	 * @param gep the gap extension penalty; should be nonnegative
	 */
	public SimpleGapPenalty(int gop, int gep) {
		this.gop = -Math.abs(gop);
		this.gep = -Math.abs(gep);
		setType();
	}

	/**
	 * <strong>Returns the negative of the extension penalty passed to the constructor.</strong>
	 */
	@Override
	public int getExtensionPenalty() {
		return gep;
	}

	/**
	 * <strong>Returns the negative of the opening penalty passed to the constructor.</strong>
	 */
	@Override
	public int getOpenPenalty() {
		return gop;
	}

	@Override
	public Type getType() {
		return type;
	}

	/**
	 * @param gep Should be nonnegative
	 */
	@Override
	public void setExtensionPenalty(int gep) {
		this.gep = -Math.abs(gep);
		setType();
	}

	/**
	 * @param gop Should be nonnegative
	 */
	@Override
	public void setOpenPenalty(int gop) {
		this.gop = -Math.abs(gop);
		setType();
	}

	// helper method to set the type given the open and extension penalties
	private void setType() {
		type = (gop == 0) ? GapPenalty.Type.LINEAR : ((gep == 0) ? GapPenalty.Type.CONSTANT : GapPenalty.Type.AFFINE);
	}

}
