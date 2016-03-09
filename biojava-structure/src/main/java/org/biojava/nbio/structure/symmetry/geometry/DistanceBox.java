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
package org.biojava.nbio.structure.symmetry.geometry;

import javax.vecmath.Point3d;
import java.util.*;

/**
 *
 * @author Peter
 */
public class DistanceBox<T> {
	private Map<Long, List<T>> map;
	private Map<Long, List<T>> layerMap;
	private boolean modified;
	private double inverseBinWidth;
	// relative locations of 26 boxes surrounding a box at (0, 0, 0)
	private static final long[] offset = new long[] {
		0 + ( 0 * 10000) + ( 0 * 1000000000L),
	   -1 + (-1 * 10000) + (-1 * 1000000000L),
	   -1 + (-1 * 10000) + ( 0 * 1000000000L),
	   -1 + (-1 * 10000) + ( 1 * 1000000000L),
	   -1 + ( 0 * 10000) + (-1 * 1000000000L),
	   -1 + ( 0 * 10000) + ( 0 * 1000000000L),
	   -1 + ( 0 * 10000) + ( 1 * 1000000000L),
	   -1 + ( 1 * 10000) + (-1 * 1000000000L),
	   -1 + ( 1 * 10000) + ( 0 * 1000000000L),
	   -1 + ( 1 * 10000) + ( 1 * 1000000000L),
		0 + (-1 * 10000) + (-1 * 1000000000L),
		0 + (-1 * 10000) + ( 0 * 1000000000L),
		0 + (-1 * 10000) + ( 1 * 1000000000L),
		0 + ( 0 * 10000) + (-1 * 1000000000L),
		0 + ( 0 * 10000) + ( 1 * 1000000000L),
		0 + ( 1 * 10000) + (-1 * 1000000000L),
		0 + ( 1 * 10000) + ( 0 * 1000000000L),
		0 + ( 1 * 10000) + ( 1 * 1000000000L),
		1 + (-1 * 10000) + (-1 * 1000000000L),
		1 + (-1 * 10000) + ( 0 * 1000000000L),
		1 + (-1 * 10000) + ( 1 * 1000000000L),
		1 + ( 0 * 10000) + (-1 * 1000000000L),
		1 + ( 0 * 10000) + ( 0 * 1000000000L),
		1 + ( 0 * 10000) + ( 1 * 1000000000L),
		1 + ( 1 * 10000) + (-1 * 1000000000L),
		1 + ( 1 * 10000) + ( 0 * 1000000000L),
		1 + ( 1 * 10000) + ( 1 * 1000000000L)
	};

	private List<T> tempBox = new ArrayList<T>(offset.length);

	/** Creates a new instance of DistanceBox */
	public DistanceBox(double binWidth) {
		map = new HashMap<Long, List<T>>();
		layerMap = new HashMap<Long, List<T>>();
		this.inverseBinWidth = 1.0f/binWidth;
		this.modified = true;
	}

	public void addPoint(Point3d point, T object) {
		long i = (long) StrictMath.rint(point.x * inverseBinWidth);
		long j = (long) StrictMath.rint(point.y * inverseBinWidth);
		long k = (long) StrictMath.rint(point.z * inverseBinWidth);
		long location = i + (j * 10000L) + (k * 1000000000L);

		List<T> box = map.get(location);

		if (box == null) {
			box = new ArrayList<T>();
			map.put(location, box);
		}

		box.add(object);
		modified = true;
	}

	public List<T> getNeighborsWithCache(Point3d point) {
		if (modified) {
			layerMap.clear();
			modified = false;
		}

		long i = (long) StrictMath.rint(point.x * inverseBinWidth);
		long j = (long) StrictMath.rint(point.y * inverseBinWidth);
		long k = (long) StrictMath.rint(point.z * inverseBinWidth);
		long location = i + (j * 10000L) + (k * 1000000000L);

		List<T> box = layerMap.get(location);

		if (box == null) {
			box = getBoxTwo(location);
			layerMap.put(location, box);
		}

		return box;
	}

	public List<T> getNeighbors(Point3d point) {
		if (modified) {
			layerMap.clear();
			modified = false;
		}

		long i = (long) StrictMath.rint(point.x * inverseBinWidth);
		long j = (long) StrictMath.rint(point.y * inverseBinWidth);
		long k = (long) StrictMath.rint(point.z * inverseBinWidth);
		long location = i + (j * 10000L) + (k * 1000000000L);

		List<T> box = getBoxTwo(location);
		return box;
	}

	public List<T> getIntersection(DistanceBox<T> distanceBox) {
		List<T> intersection = new ArrayList<T>();
		HashSet<Long> checkedLocations = new HashSet<Long>();

		for (Iterator<Long> iter = map.keySet().iterator(); iter.hasNext();) {
			long location = iter.next();
			boolean overlap = false;
			for (int i = 0, n = offset.length; i < n; i++) {
				long loc = location + offset[i];
				if (distanceBox.contains(loc)) {
					overlap = true;
					break;
				}
			}
			if (overlap) {
				for (int i = 0, n = offset.length; i < n; i++) {
					long loc = location + offset[i];
					if (checkedLocations.contains(loc))
						continue;

					checkedLocations.add(loc);
					if (contains(loc)) {
						intersection.addAll(map.get(loc));
					}
				}
			}
		}
		return intersection;
	}

	private List<T> getBoxTwo(long location) {
		tempBox.clear();
		for (int i = 0, n = offset.length; i < n; i++) {
			List<T> box = map.get(location + offset[i]);
			if (box != null) {
				tempBox.addAll(box);
			}
		}
		// ensure that boxTwo has no empty element by copying from tempBox of defined size
		List<T> boxTwo = null;
		if (tempBox.size() == 0) {
			boxTwo = Collections.emptyList();
		} else if (tempBox.size() == 1) {
			boxTwo = Collections.singletonList(tempBox.get(0));
		} else {
			boxTwo = new ArrayList<T>(tempBox);
		}
		return boxTwo;
	}

	private boolean contains(long location) {
		return map.containsKey(location);
	}
}
