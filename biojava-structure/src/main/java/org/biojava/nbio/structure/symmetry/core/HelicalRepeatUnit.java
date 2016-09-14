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

import org.biojava.nbio.structure.geometry.CalcPoint;
import javax.vecmath.Point3d;

import java.util.*;
import java.util.Map.Entry;

public class HelicalRepeatUnit {
	private QuatSymmetrySubunits subunits = null;
	private List<Point3d> repeatUnitCenters = new ArrayList<Point3d>();
	private List<Point3d[]> repeatUnits = new ArrayList<Point3d[]>();
	private List<List<Integer>> repeatUnitIndices = new ArrayList<List<Integer>>();
	private Map<Integer[], Integer> interactingNeighbors = Collections.emptyMap();

public HelicalRepeatUnit(QuatSymmetrySubunits subunits) {
	this.subunits = subunits;
}

public List<Point3d> getRepeatUnitCenters() {
	if (repeatUnitCenters.isEmpty()) {
		run();
	}
	return repeatUnitCenters;
}

public List<Point3d[]> getRepeatUnits() {
	if (repeatUnits.isEmpty()) {
		run();
	}
	return repeatUnits;
}

public List<List<Integer>> getRepeatUnitIndices() {
	return repeatUnitIndices;
}

public Map<Integer[], Integer> getInteractingRepeatUnits() {
	if (interactingNeighbors.isEmpty()) {
		run();
	}
	return interactingNeighbors;
}

private void run() {
	this.repeatUnitCenters = calcRepeatUnitCenters();
	if (this.repeatUnitCenters.size() == 0) {
		return;
	}
	this.repeatUnits = calcRepeatUnits();
	this.interactingNeighbors = findInteractingNeigbors();
}

private List<Point3d> calcRepeatUnitCenters() {
	
	// TODO why do we use models here? it should not matter. Setting to 0 all
	List<Integer> models = new ArrayList<Integer>(subunits.getSubunitCount());
	for (int s = 0; s <subunits.getSubunitCount(); s++)
		models.add(0);
	Set<Integer> uniqueModels = new HashSet<Integer>(Arrays.asList(1));
	
	int modelCount = uniqueModels.size();
	List<Integer> folds = this.subunits.getFolds();
	int maxFold = folds.get(folds.size()-1);

	List<Point3d> repeatCenters = new ArrayList<Point3d>();
	List<Point3d> centers = subunits.getCenters();

//	if (modelCount == maxFold && subunits.getSubunitCount() > 3) {
	if (maxFold%modelCount == 0 && modelCount > 1 && subunits.getSubunitCount() > 3) {
//		System.out.println("calcRepeatUnitCenters case 1");
		for (int i = 0; i < modelCount; i++) {
			List<Integer> subunitIndices = new ArrayList<Integer>();
			Point3d p = new Point3d();
			int count = 0;
//			System.out.println("Models: " + models.size());
			for (int j = 0; j < models.size(); j++) {
				if (models.get(j) == i) {
					p.add(centers.get(j));
					subunitIndices.add(j);
					count++;
				}
			}
//			System.out.println("count: " + count);
			p.scale(1.0/count);
//			System.out.println("Orig Repeat unit: " + p);
			repeatCenters.add(p);
			repeatUnitIndices.add(subunitIndices);
		}
	} else {
//		System.out.println("calcRepeatUnitCenters case21");
		// TODO need to group subunits into repeating groups
		// Case of 3B5U: A14, but seems to form (A2)*7 and symmetry related subunits don't have direct contact
		List<Integer> sequenceClusterIds = subunits.getClusterIds();
		for (int i = 0; i < subunits.getSubunitCount(); i++) {
			List<Integer> subunitIndices = new ArrayList<Integer>(1);
			if (sequenceClusterIds.get(i) == 0) {
				repeatCenters.add(new Point3d(centers.get(i)));
//				System.out.println("Orig Repeat unit: " + centers.get(i));
				subunitIndices.add(i);
				repeatUnitIndices.add(subunitIndices);
			}
		}
	}

	// helixes should have at least 3 repeat centers
//	System.out.println("Number of repeat centers: " + repeatCenters.size());
	if (repeatCenters.size() < 3) {
		repeatCenters.clear();
	}

	return repeatCenters;
}

private List<Point3d[]> calcRepeatUnits() {
	
	// TODO why do we use models here? it should not matter. Setting to 0 all
	List<Integer> models = new ArrayList<Integer>(
			subunits.getSubunitCount());
	for (int s = 0; s < subunits.getSubunitCount(); s++)
		models.add(0);
	Set<Integer> uniqueModels = new HashSet<Integer>(Arrays.asList(1));
		
	int modelCount = uniqueModels.size();
	List<Integer> folds = this.subunits.getFolds();
	int maxFold = folds.get(folds.size()-1);

	List<Point3d[]> repeatTraces = new ArrayList<Point3d[]>();
	List<Point3d[]> traces = subunits.getTraces();

//	if (modelCount == maxFold && subunitCount > 3) {
	if (maxFold%modelCount == 0 && modelCount > 1 && subunits.getSubunitCount() > 3) {
		for (int i = 0; i < modelCount; i++) {
			List<Point3d> coords = new ArrayList<Point3d>();
			for (int j = 0; j < models.size(); j++) {
				if (models.get(j) == i) {
					coords.addAll(Arrays.asList(traces.get(j)));
				}
			}
			Point3d[] x = new Point3d[coords.size()];
			coords.toArray(x);
//			repeatTraces.add(x); // make sure we don't accidently change the original coordinates
			repeatTraces.add(CalcPoint.clonePoint3dArray(x));
		}
	} else {
		List<Integer> sequenceClusterIds = subunits.getClusterIds();
		for (int i = 0; i < subunits.getSubunitCount(); i++) {
			if (sequenceClusterIds.get(i) == 0) {
				Point3d[] x = subunits.getTraces().get(i);
				repeatTraces.add(CalcPoint.clonePoint3dArray(x));
			}
		}
	}

//	for (int i = 0; i < repeatTraces.size(); i++) {
//		System.out.println("Repeat " + i);
//		System.out.println(Arrays.toString(repeatTraces.get(i)));
//	}
	return repeatTraces;
}

private Map<Integer[], Integer> findInteractingNeigbors() {
	Map<Integer[], Integer>  contactMap = new HashMap<Integer[], Integer>();

	Map<Integer, List<Integer[]>> distanceMap = findClosestPairs(8);
	for (List<Integer[]> pairs: distanceMap.values())
	for (Integer[] pair: pairs) {
		Integer contacts = calcContactNumber(repeatUnits.get(pair[0]), repeatUnits.get(pair[1]));
//		System.out.println("contacts: " + pair[0] + "-" + pair[1] + ": " + contacts);
		if (contacts > 0) {
			contactMap.put(pair, contacts);
		}
	}

	return contactMap;
}

private Map<Integer, List<Integer[]>> findClosestPairs(int maxNeighbors) {
	Map<Integer, List<Integer[]>>  reducedMap = new TreeMap<Integer, List<Integer[]>>();

	Map<Integer, List<Integer[]>>  distanceMap = new TreeMap<Integer, List<Integer[]>>();
	int nCenters = repeatUnitCenters.size();
//	System.out.println("repeatUnitCenters: " + repeatUnitCenters);

	for (int i = 0; i < nCenters-1; i++) {
		for (int j = i+1; j < nCenters; j++) {
			float dist = (float)repeatUnitCenters.get(i).distance(repeatUnitCenters.get(j));
			// round to 2 digits precision
//			System.out.println("dist pair: " + i + "-" + j + ": " + dist);
			dist *= 100;
			int intDist = Math.round(dist);
			List<Integer[]> pairs = distanceMap.get(intDist);
			// save only one representative pair for each distance
			if (pairs == null) {
				pairs = new ArrayList<Integer[]>();
			}
			Integer[] pair = new Integer[2];
			pair[0] = i;
			pair[1] = j;
			pairs.add(pair);
			distanceMap.put(intDist, pairs);
		}
	}

	int count = 0;
	for (Entry<Integer, List<Integer[]>> entry: distanceMap.entrySet()) {
		if (! (reducedMap.containsKey(entry.getKey())) ) {
			reducedMap.put(entry.getKey(), entry.getValue());
//			System.out.println("dist pair: " + entry.getValue() + ": " + entry.getKey());
			count++;
			if (count == maxNeighbors) {
				break;
			}
		}
	}
	distanceMap.clear();

	return reducedMap;
}

private static int calcContactNumber(Point3d[] a, Point3d[] b) {
	int contacts = 0;
	for (Point3d pa : a) {
		for (Point3d pb : b) {
			if (pa.distance(pb) < 10) {
				contacts++;
			}
		}
	}
	return contacts;
}
}
