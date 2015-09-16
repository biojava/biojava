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
 */

package org.biojava.nbio.structure.symmetry.core;

import org.biojava.nbio.structure.symmetry.geometry.MomentsOfInertia;
import org.biojava.nbio.structure.symmetry.geometry.SuperPosition;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import java.util.*;
import java.util.Map.Entry;


/**
 * A bean to represent info about the set of subunits being considered for a
 * QuatSymmetryDetector alignment.
 * @author Peter Rose
 */
public class Subunits {
    private List<Point3d[]> caCoords = Collections.emptyList();
    private List<Integer> sequenceClusterIds = Collections.emptyList();
    private List<Boolean> pseudoStoichiometry = Collections.emptyList();
    private List<Double> minSequenceIdentity = Collections.emptyList();
    private List<Double> maxSequenceIdentity = Collections.emptyList();
    private List<String> chainIds = Collections.emptyList();
    private List<Integer> modelNumbers =  Collections.emptyList();
    private List<Integer> folds = Collections.emptyList();
    private List<Point3d> originalCenters = new ArrayList<Point3d>();
    private List<Point3d> centers = new ArrayList<Point3d>();
    private List<Vector3d> unitVectors = new ArrayList<Vector3d>();
    private int nucleicAcidChainCount = 0;
    private boolean pseudoSymmetric = false;

    private Point3d centroid;
    private MomentsOfInertia momentsOfInertia = new MomentsOfInertia();

    /**
     * All inputs should contain one element per subunit.
     * @param caCoords CA coordinates of all subunits
     * @param sequenceClusterIds ID of the cluster that each subunit belongs to
     * @param pseudoStoichiometry Whether pseudosymmetry was used when clustering the subunit
     * @param minSequenceIdentity Minimum sequence identity to other cluster members
     * @param maxSequenceIdentity Maximum sequence identity to other cluster members
     * @param folds Valid symmetry orders for this stoichiometry
     * @param chainIds Chain ID for the subunit
     * @param modelNumbers Model number for the subunit
     */
    public Subunits(List<Point3d[]> caCoords, List<Integer> sequenceClusterIds, List<Boolean> pseudoStoichiometry, List<Double> minSequenceIdentity, List<Double> maxSequenceIdentity, List<Integer> folds, List<String> chainIds, List<Integer> modelNumbers) {
    	this.caCoords = caCoords;
    	this.sequenceClusterIds = sequenceClusterIds;
    	this.pseudoStoichiometry = pseudoStoichiometry;
    	this.minSequenceIdentity = minSequenceIdentity;
    	this.maxSequenceIdentity = maxSequenceIdentity;
    	this.folds = folds;
    	this.chainIds = chainIds;
    	this.modelNumbers = modelNumbers;
    }

    public List<Point3d[]> getTraces() {
        return caCoords;
    }

    public int getSubunitCount() {
        run();
        if (centers == null) {
        	return 0;
        }
        return centers.size();
    }
    
    public List<Integer> getSequenceClusterIds() {
    	return sequenceClusterIds;
    }
    
    public boolean isPseudoStoichiometric() {
    	for (Boolean b: pseudoStoichiometry) {
    		if (b) {
    			return true;
    		}
    	}
    	return false;
    }
    
    public boolean isPseudoSymmetric() {
		return pseudoSymmetric;
	}

	public void setPseudoSymmetric(boolean pseudoSymmetric) {
		this.pseudoSymmetric = pseudoSymmetric;
	}

	public double getMinSequenceIdentity() {
    	double minId = 1.0;
    	for (double seqId: minSequenceIdentity) {
    		minId = Math.min(seqId,  minId);
    	}
    	return minId;
    }
    
    public double getMaxSequenceIdentity() {
    	double maxId = 1.0;
    	for (double seqId: maxSequenceIdentity) {
    		maxId = Math.min(seqId,  maxId);
    	}
    	return maxId;
    }
      
    public List<String> getChainIds() {
    	return chainIds;
    }
    
    public List<Integer> getModelNumbers() {
    	return modelNumbers;
    }
    
    public List<Integer>getFolds() {
    	return folds;
    }
    
    public String getStoichiometry() {	
		// count number of members in each cluster
    	Map<Integer, Integer> map = new TreeMap<Integer, Integer>();
		for (Integer id: sequenceClusterIds) {
			Integer value = map.get(id);
			if (value == null) {
				value = new Integer(1);
			} else {
				value++;
			}
			map.put(id, value);
		}
		
		// build formula string
    	String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		StringBuilder formula = new StringBuilder();
		for (Entry<Integer, Integer> entry: map.entrySet()) {
			String key = "?";
			int id = entry.getKey();
			if (id < alpha.length()) {
				key = alpha.substring(id, id+1);
			}
			formula.append(key);
			if (entry.getValue() > 1) {
				formula.append(entry.getValue());
			}
		}

        return formula.toString();
    }
    
    public int getCalphaCount() {
    	int count = 0;
    	for (Point3d[] trace: caCoords) {
    		count += trace.length;
    	}
    	return count;
    }
    
    public int getLargestSubunit() {
    	int index = -1;
    	int maxLength = 0;
    	for (int i = 0; i < caCoords.size(); i++) {
    		int length = caCoords.get(i).length;
    		if (length > maxLength) {
    			index = i;
    		}
    	}
    	return index;
    }

    public List<Point3d> getCenters() {
        run();
        return centers;
    }
    
    public List<Vector3d> getUnitVectors() {
        run();
        return unitVectors;
    }

    public List<Point3d> getOriginalCenters() {
        run();
        return originalCenters;
    }

    public Point3d getCentroid() {
        run();
        return centroid;
    }

    public MomentsOfInertia getMomentsOfInertia() {
    	run();
    	return momentsOfInertia;
    }
    
    /**
	 * @return the nucleicAcidChainCount
	 */
	public int getNucleicAcidChainCount() {
		run();
		return nucleicAcidChainCount;
	}

	/**
	 * @param nucleicAcidChainCount the nucleicAcidChainCount to set
	 */
	public void setNucleicAcidChainCount(int nucleicAcidChainCount) {
		this.nucleicAcidChainCount = nucleicAcidChainCount;
	}

	public boolean overlaps(Subunits subunits) {
    	Set<String> set1 = getSignatures(this);
    	Set<String> set2 = getSignatures(subunits);
    	set1.retainAll(set2);
    	return set1.size() > 0;
    }
    
    public boolean contains(Subunits subunits) {
    	Set<String> set1 = getSignatures(this);
    	Set<String> set2 = getSignatures(subunits);
    	return set1.containsAll(set2);
    }
    
    private static Set<String> getSignatures(Subunits subunits) {
    	Set<String> set = new HashSet<String>(subunits.getSubunitCount());
    	for (int i = 0; i < subunits.getSubunitCount(); i++) {
    		set.add(subunits.getChainIds().get(i) + "_" + subunits.getModelNumbers().get(i));
    	}
    	return set;
    }
    
    private void run() {
        if (centers.size() > 0) {
            return;
        }
        calcOriginalCenters();
        calcCentroid();
        calcCenters();
        calcMomentsOfIntertia();
    }

    private void calcOriginalCenters() {
        for (Point3d[] trace: caCoords) {
            Point3d com = SuperPosition.centroid(trace);
            originalCenters.add(com);
        }
    }

    private void calcCentroid() {
        Point3d[] orig = originalCenters.toArray(new Point3d[originalCenters.size()]);
        centroid = SuperPosition.centroid(orig);
    }

    private void calcCenters() {
        for (Point3d p: originalCenters) {
            Point3d c = new Point3d(p);
            c.sub(centroid);
            centers.add(c);
            Vector3d v = new Vector3d(c);
            v.normalize();
            unitVectors.add(v);
        }
    }
    
    public Point3d getLowerBound() {
    	Point3d lower = new Point3d();
    	for (Point3d p: centers) {
    		if (p.x < lower.x) {
    			lower.x = p.x;
    		}
    		if (p.y < lower.y) {
    			lower.y = p.y;
    		}
    		if (p.z < lower.z) {
    			lower.z = p.z;
    		}
    	}
    	return lower;
    }
    
    public Point3d getUpperBound() {
    	Point3d upper = new Point3d();
    	for (Point3d p: centers) {
    		if (p.x > upper.x) {
    			upper.x = p.x;
    		}
    		if (p.y > upper.y) {
    			upper.y = p.y;
    		}
    		if (p.z > upper.z) {
    			upper.z = p.z;
    		}
    	}
    	return upper;
    }
    
    private void calcMomentsOfIntertia() {
    	for (Point3d[] trace: caCoords) {
    		for (Point3d p: trace) {
    			momentsOfInertia.addPoint(p, 1.0f);
    		}
    	}
    }
}
