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
 * Created on 2013-05-23
 *
 */
package org.biojava.nbio.structure.symmetry.core;

/**
 * Holds the results of quaternary symmetry perception.
 *
 * @author Peter Rose
 * @author Aleix Lafita
 *
 */
public class QuatSymmetryResults {

	private Subunits subunits;
	private RotationGroup rotationGroup;
	private HelixLayers helixLayers;
	private SymmetryPerceptionMethod method;
	private boolean local = false;

	public QuatSymmetryResults(Subunits subunits, RotationGroup rotationGroup,
			SymmetryPerceptionMethod method) {
		this.subunits = subunits;
		this.rotationGroup = rotationGroup;
		this.method = method;
	}

	public QuatSymmetryResults(Subunits subunits, HelixLayers helixLayers,
			SymmetryPerceptionMethod method) {
		this.subunits = subunits;
		this.helixLayers = helixLayers;
		this.method = method;
	}

	/**
	 * Returns protein subunit information that was used to determine symmetry
	 * information
	 *
	 * @return
	 */
	public Subunits getSubunits() {
		return subunits;
	}

	/**
	 * Returns rotation group (point group) information representing rotational
	 * quaternary symmetry, see
	 * http://en.wikipedia.org/wiki/Rotation_group_SO(3)
	 *
	 * @return rotation group
	 */
	public RotationGroup getRotationGroup() {
		return rotationGroup;
	}

	/**
	 * Returns helix layers (layer lines) as a list of helices that describe a
	 * helical structure
	 */
	public HelixLayers getHelixLayers() {
		return helixLayers;
	}

	/**
	 * Returns the method used for symmetry perception.
	 *
	 * @return method
	 */
	public SymmetryPerceptionMethod getMethod() {
		return method;
	}

	/**
	 * Returns the symmetry group. For point groups returns the point group
	 * symbol and for helical symmetry returns "H".
	 * 
	 * @return symmetry group symbol
	 */
	public String getSymmetry() {
		if (helixLayers != null && helixLayers.size() > 0) {
			return "H";
		} else if (rotationGroup != null && rotationGroup.getOrder() > 0) {
			return rotationGroup.getPointGroup();
		}
		return "";
	}

	public QuatSymmetryScores getScores() {
		if (helixLayers != null && helixLayers.size() > 0) {
			return helixLayers.getScores();
		} else if (rotationGroup != null && rotationGroup.getOrder() > 0) {
			return rotationGroup.getScores();
		}
		return new QuatSymmetryScores();
	}

	public int getNucleicAcidChainCount() {
		return subunits.getNucleicAcidChainCount();
	}

	@Override
	public String toString() {
		return "QuatSymmetryResults [symmetry=" + getSymmetry()
				+ ", stoichiometry=" + subunits.getStoichiometry()
				+ ", method=" + method + ", local=" + local + "]";
	}

	public boolean isLocal() {
		return local;
	}

	public void setLocal(boolean local) {
		this.local = local;
	}

}
