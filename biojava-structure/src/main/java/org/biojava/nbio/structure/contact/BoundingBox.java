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
package org.biojava.nbio.structure.contact;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import java.io.Serializable;
import java.util.Arrays;


/**
 * A bounding box for short cutting some geometrical calculations.
 *
 * See http://en.wikipedia.org/wiki/Bounding_volume
 *
 * @author Jose Duarte
 *
 */
public class BoundingBox implements Serializable {

	private static final long serialVersionUID = 1L;
	private static final Logger logger = LoggerFactory.getLogger(StructureInterfaceList.class);


	public double xmin;
	public double xmax;
	public double ymin;
	public double ymax;
	public double zmin;
	public double zmax;

	public BoundingBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
		this.xmin = xmin;
		this.xmax = xmax;
		this.ymin = ymin;
		this.ymax = ymax;
		this.zmin = zmin;
		this.zmax = zmax;
	}

	public BoundingBox(BoundingBox bb) {
		this.xmin = bb.xmin;
		this.xmax = bb.xmax;
		this.ymin = bb.ymin;
		this.ymax = bb.ymax;
		this.zmin = bb.zmin;
		this.zmax = bb.zmax;
	}

	/**
	 * Constructs a BoundingBox by calculating maxs and mins of given array of atoms.
	 * @param atoms
	 */
	public BoundingBox (Point3d[] atoms) {

		if (atoms.length==0) logger.error("Error! Empty list of atoms");

		xmax = atoms[0].x;
		xmin = xmax;
		ymax = atoms[0].y;
		ymin = ymax;
		zmax = atoms[0].z;
		zmin = zmax;

		for(int i=1;i<atoms.length;i++) {
			if(atoms[i].x > xmax) xmax = atoms[i].x;
			else if(atoms[i].x < xmin) xmin = atoms[i].x;

			if(atoms[i].y > ymax) ymax = atoms[i].y;
			else if(atoms[i].y < ymin) ymin = atoms[i].y;

			if(atoms[i].z > zmax) zmax = atoms[i].z;
			else if(atoms[i].z < zmin) zmin = atoms[i].z;
		}

	}


	/** Returns the dimensions of this bounding box.
	 *
	 * @return a double array (x,y,z) with the dimensions of the box.
	 */
	public double[] getDimensions(){

		double[] dim = new double[3];

		dim[0] = xmax-xmin;
		dim[1] = ymax-ymin;
		dim[2] = zmax-zmin;

		return dim;

	}

	/**
	 * Given a set of bounding boxes returns a bounding box that bounds all of them.
	 * @param boxes
	 */
	public BoundingBox(BoundingBox[] boxes) {

		if (boxes.length==0) logger.error("Error! Empty list of bounding boxes");

		xmax = boxes[0].xmax;
		xmin = boxes[0].xmin;
		ymax = boxes[0].ymax;
		ymin = boxes[0].ymin;
		zmax = boxes[0].zmax;
		zmin = boxes[0].zmin;

		for (int i=1;i<boxes.length;i++) {
			if(boxes[i].xmax > xmax) xmax = boxes[i].xmax;
			else if(boxes[i].xmin < xmin) xmin = boxes[i].xmin;
			if(boxes[i].ymax > ymax) ymax = boxes[i].ymax;
			else if(boxes[i].ymin < ymin) ymin = boxes[i].ymin;
			if(boxes[i].zmax > zmax) zmax = boxes[i].zmax;
			else if(boxes[i].zmin < zmin) zmin = boxes[i].zmin;
		}

	}

	private class Bound implements Comparable<Bound> {
		int cardinal;
		double value;
		public Bound(int cardinal,double value) {
			this.cardinal = cardinal;
			this.value = value;
		}
		@Override
		public int compareTo(Bound o) {
			return Double.compare(this.value,o.value);
		}
		@Override
		public String toString() {
			return "["+cardinal+","+value+"]";
		}
	}

	/**
	 * Returns true if this bounding box overlaps given one, i.e. they are within
	 * one cutoff distance in one of their 3 dimensions.
	 * @param cutoff
	 * @return
	 */
	public boolean overlaps(BoundingBox o, double cutoff) {
		if (this==o) return true;
		// x dimension
		if (!areOverlapping(xmin,xmax,o.xmin,o.xmax,cutoff)) {
			return false;
		}
		// y dimension
		if (!areOverlapping(ymin,ymax,o.ymin,o.ymax,cutoff)) {
			return false;
		}
		// z dimension
		if (!areOverlapping(zmin,zmax,o.zmin,o.zmax,cutoff)) {
			return false;
		}
		return true;
	}

	private boolean areOverlapping(double imin, double imax, double jmin, double jmax, double cutoff) {

		Bound[] bounds = {new Bound(0,imin), new Bound(1,imax),
				new Bound(2,jmin), new Bound(3,jmax)};

		Arrays.sort(bounds);

		if ((bounds[0].cardinal==0 && bounds[1].cardinal==1)) {
			if ((bounds[2].value-bounds[1].value)>cutoff) {
				return false;
			}
		} else if (bounds[0].cardinal==2 && bounds[1].cardinal==3) {
			if ((bounds[2].value-bounds[1].value)>cutoff) {
				return false;
			}
		}

		return true;

	}
	
	/**
	 * Check if a given point falls within this box
	 * @param atom
	 * @return
	 */
	public boolean contains(Point3d atom) {
		double x = atom.x;
		double y = atom.y;
		double z = atom.z;
		return xmin <= x && x <= xmax
				&& ymin <= y && y <= ymax
				&& zmin <= z && z <= zmax;
	}

	public void translate(Vector3d translation) {
		xmin+=translation.x;
		xmax+=translation.x;
		ymin+=translation.y;
		ymax+=translation.y;
		zmin+=translation.z;
		zmax+=translation.z;
	}

	/**
	 * Returns an array of size 2 with min and max values of given double array
	 * @param array
	 * @return
	 */
	public double[] getMinMax(double[] array) {
		double[] minmax = new double[2];

		double max = Double.MIN_VALUE;
		double min = Double.MAX_VALUE;

		for(double value : array) {
			if(value > max) max = value;
			if(value < min) min = value;
		}

		minmax[0] = min;
		minmax[1] = max;
		return minmax;
	}

	@Override
	public String toString() {
		return String.format("[(%7.2f,%7.2f),(%7.2f,%7.2f),(%7.2f,%7.2f)]", xmin,xmax,ymin,ymax,zmin,zmax);
	}
}
