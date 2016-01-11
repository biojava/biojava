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
 *
 */
public class QuatSymmetryResults {
	private Subunits subunits = null;
	private RotationGroup rotationGroup = null;
	private HelixLayers helixLayers = null;
	private String method = null;
	private double sequenceIdentityThreshold = 0;
	private boolean local = false;
	private boolean preferredResult = false;
	
	public QuatSymmetryResults(Subunits subunits, RotationGroup rotationGroup, String method) {
	//	SymmetryDeviation sd = new SymmetryDeviation(subunits, rotationGroup);
	//	rotationGroup.setSymmetryDeviation(sd.getSymmetryDeviation());
		this.subunits = subunits;
		this.rotationGroup = rotationGroup;
		this.method = method;
	}
	
	public QuatSymmetryResults(Subunits subunits, HelixLayers helixLayers, String method) {
	//	SymmetryDeviation sd = new SymmetryDeviation(subunits, helixLayers);
	//	helixLayers.setSymmetryDeviation(sd.getSymmetryDeviation());
		this.subunits = subunits;
		this.helixLayers = helixLayers;
		this.method = method;	
	}
	
	/**
	 * Returns protein subunit information that was used to determine symmetry information
	 * 
	 * @return
	 */
	public Subunits getSubunits() {
		return subunits;
	}
	
	/**
	 * Returns rotation group (point group) information representing rotational quaternary symmetry,
	 * see http://en.wikipedia.org/wiki/Rotation_group_SO(3)
	 * 
	 * @return rotation group
	 */
	public RotationGroup getRotationGroup() {
		return rotationGroup;
	}

	/*
	 * Returns helix layers (layer lines) as a list of helices that describe a helical structure
	 */
	public HelixLayers getHelixLayers() {
		return helixLayers;
	}

	/**
	 * Returns name of method used for symmetry perception.
	 * 
	 * @return method
	 */
	public String getMethod() {
		return method;
	}

	
    /**
     * Returns the symmetry group. For point groups returns the point group symbol
     * and for helical symmetry returns "H".
     * @return symmetry symbol
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
	
	/**
	 * Returns the average Calpha trace RMSD for all symmetry operations
	 * @return
	 * @deprecated use {@link getScores()} instead.  
	 */
	@Deprecated
	public double getAverageTraceRmsd() {
		if (helixLayers != null && helixLayers.size() > 0) {
			return helixLayers.getScores().getRmsd();
		} else if (rotationGroup != null && rotationGroup.getOrder() > 0) {
			return rotationGroup.getScores().getRmsd();
		}
		return 0;
	}
	
	/**
	 * Returns the average Calpha trace Tm for all symmetry operations
	 * @return
	 * @deprecated use {@link getScores()} instead.  
	 */
	@Deprecated
	public double getAverageTraceTmScoreMin() {
		if (helixLayers != null && helixLayers.size() > 0) {
			return helixLayers.getScores().getTm();
		} else if (rotationGroup != null && rotationGroup.getOrder() > 0) {
			return rotationGroup.getScores().getTm();
		}
		return 0;
	}
	
	public int getNucleicAcidChainCount() {
		return subunits.getNucleicAcidChainCount();
	}
	
	public double getSequenceIdentityThreshold() {
		return sequenceIdentityThreshold;
	}

	public void setSequenceIdentityThreshold(double sequenceIdentityThreshold) {
		this.sequenceIdentityThreshold = sequenceIdentityThreshold;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		
	    sb.append("Stoichiometry         : ");
	    sb.append(getSubunits().getStoichiometry());
	    sb.append("\n");
	    sb.append("Pseudostoichiometry   : ");
	    sb.append(getSubunits().isPseudoStoichiometric());
	    sb.append("\n");
	    sb.append("Pseudosymmetry        : ");
	    sb.append(getSubunits().isPseudoSymmetric());
	    sb.append("\n");
	    sb.append("Min sequence identity : ");
	    sb.append(Math.round(getSubunits().getMinSequenceIdentity()*100));
	    sb.append("\n");
	    sb.append("Max sequence identity : ");
	    sb.append(Math.round(getSubunits().getMaxSequenceIdentity()*100));
	    sb.append("\n");
	    sb.append("Symmetry              : ");
	    sb.append(getSymmetry());				
	    sb.append("\n");
	    sb.append("Symmetry RMSD         : ");
	    sb.append((float) getAverageTraceRmsd());
	    sb.append("\n");
	    sb.append("Symmetry TmScoreMin   : ");
	    sb.append((float) getAverageTraceTmScoreMin());
	    sb.append("\n");
	    sb.append("Prefered result       : ");
	    sb.append(isPreferredResult());
	    sb.append("\n");
	    
	    return sb.toString();
	}
	
	/**
	 * Return true 
	 * @return
	 */
	public boolean isLocal() {
		return local;
	}

	public void setLocal(boolean local) {
		this.local = local;
	}

	public boolean isPreferredResult() {
		return preferredResult;
	}

	public void setPreferredResult(boolean preferredResult) {
		this.preferredResult = preferredResult;
	}
}
