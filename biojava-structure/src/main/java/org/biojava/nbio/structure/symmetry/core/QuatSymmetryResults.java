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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.cluster.Subunit;
import org.biojava.nbio.structure.cluster.SubunitCluster;
import org.biojava.nbio.structure.cluster.SubunitClusterUtils;

/**
 * Holds the results of quaternary symmetry perception obtained with
 * {@link QuatSymmetryDetector}.
 *
 * @author Peter Rose
 * @author Aleix Lafita
 *
 */
public class QuatSymmetryResults {

	// Optional: the query Structure, if any
	private Structure structure;

	// Information about the clustering process
	private List<SubunitCluster> clusters;
	private boolean local = false;

	// Cached properties
	private String stoichiometry;
	private List<Subunit> subunits;

	// Information about the symmetry
	private SymmetryPerceptionMethod method;
	private RotationGroup rotationGroup;
	private HelixLayers helixLayers;
	private boolean pseudosymmetric = false;

	// TODO we should unify rotational and roto-translational results

	/**
	 * Constructor for rotational symmetries.
	 * 
	 * @param clusters
	 *            List of SubunitCluster used to calculate symmetry
	 * @param rotationGroup
	 * @param method
	 */
	public QuatSymmetryResults(List<SubunitCluster> clusters,
			RotationGroup rotationGroup, SymmetryPerceptionMethod method) {

		this.clusters = clusters;
		this.stoichiometry = SubunitClusterUtils
				.getStoichiometryString(clusters);
		
		subunits = new ArrayList<Subunit>();
		for (SubunitCluster c : clusters) {
			subunits.addAll(c.getSubunits());
		}
			
		this.rotationGroup = rotationGroup;
		this.method = method;
	}

	/**
	 * Constructor for roto-translational symmetries.
	 * 
	 * @param clusters
	 *            List of SubunitCluster used to calculate symmetry
	 * @param helixLayers
	 * @param method
	 */
	public QuatSymmetryResults(List<SubunitCluster> clusters,
			HelixLayers helixLayers, SymmetryPerceptionMethod method) {

		this.clusters = clusters;
		this.stoichiometry = SubunitClusterUtils
				.getStoichiometryString(clusters);

		this.helixLayers = helixLayers;
		this.method = method;
	}

	/**
	 * Returns the List of SubunitCluster used to calculate symmetry.
	 *
	 * @return an unmodifiable view of the original List
	 */
	public List<SubunitCluster> getSubunitClusters() {
		return Collections.unmodifiableList(clusters);
	}
	
	/**
	 * Returns the List of Subunits used to calculate symmetry.
	 *
	 * @return an unmodifiable view of the List
	 */
	public List<Subunit> getSubunits() {
		return Collections.unmodifiableList(subunits);
	}
	
	/**
	 * Return the number of Subunits involved in the symmetry.
	 * 
	 * @return the number of Subunits
	 */
	public int getSubunitCount() {
		return subunits.size();
	}
	
	/**
	 * @return rotation group (point group) information representing rotational
	 *         quaternary symmetry.
	 */
	public RotationGroup getRotationGroup() {
		return rotationGroup;
	}

	/**
	 * @return helix layers (layer lines) as a list of helices that describe a
	 *         helical structure.
	 */
	public HelixLayers getHelixLayers() {
		return helixLayers;
	}

	/**
	 * @return the method used for symmetry perception.
	 */
	public SymmetryPerceptionMethod getMethod() {
		return method;
	}

	/**
	 * @return the symmetry group symbol. For point groups returns the point
	 *         group symbol and for helical symmetry returns "H".
	 */
	public String getSymmetry() {
		if (helixLayers != null && helixLayers.size() > 0) {
			return "H";
		} else if (rotationGroup != null && rotationGroup.getOrder() > 0) {
			return rotationGroup.getPointGroup();
		}
		return "";
	}

	/**
	 * @return the quaternary scores as an object
	 */
	public QuatSymmetryScores getScores() {
		if (helixLayers != null && helixLayers.size() > 0) {
			return helixLayers.getScores();
		} else if (rotationGroup != null && rotationGroup.getOrder() > 0) {
			return rotationGroup.getScores();
		}
		return new QuatSymmetryScores();
	}

	public String getStoichiometry() {
		return stoichiometry;
	}

	public boolean isPseudoStoichiometric() {
		return SubunitClusterUtils.isPseudoStoichiometric(clusters);
	}

	/**
	 * A local result means that only a subset of the original Subunits was used
	 * for symmetry determination.
	 * 
	 * @return true if local result, false otherwise
	 */
	public boolean isLocal() {
		return local;
	}

	/**
	 * A local result means that only a subset of the original Subunits was used
	 * for symmetry determination.
	 * 
	 * @param local
	 *            true if local result, false otherwise
	 */
	void setLocal(boolean local) {
		this.local = local;
	}

	/**
	 * A symmetry result is pseudosymmetric when using pseudostoichiometry
	 * extended the symmetry over the
	 * 
	 * @return true is pseudosymmetric, false otherwise
	 */
	public boolean isPseudosymmetric() {
		return pseudosymmetric;
	}

	/**
	 * A symmetry result is pseudosymmetric when using pseudostoichiometry
	 * extended the symmetry over the
	 * 
	 * @param pseudosymmetric true if pseudosymmetric, false otherwise
	 */
	public void setPseudosymmetric(boolean pseudosymmetric) {
		this.pseudosymmetric = pseudosymmetric;
	}
	
	public Structure getStructure() {
		return structure;
	}

	public void setStructure(Structure structure) {
		this.structure = structure;
	}

	@Override
	public String toString() {
		return "QuatSymmetryResults [stoichiometry: " + getStoichiometry()
				+ ", symmetry: " + getSymmetry() + ", pseudo-stoichiometric: "
				+ isPseudoStoichiometric() + ", local: " + local + ", method: "
				+ method + "]";
	}

}
