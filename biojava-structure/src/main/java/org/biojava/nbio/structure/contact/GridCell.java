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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;


/**
 * A grid cell to be used in contact calculation via spatial hashing algorithm.
 *
 * @author Jose Duarte
 *
 */
public class GridCell {


	/**
	 * Shared empty array so that cells that never receive indices (e.g. the j indices when only one
	 * set of atoms was added to the grid) don't allocate anything at all.
	 */
	private static final int[] EMPTY = new int[0];

	/**
	 * Capacity of the index arrays on first insertion. Cell occupancy depends on the cutoff (the cell
	 * side is the cutoff), ranging from a handful of atoms for small cutoffs to a few tens for large
	 * ones, so we start small and grow geometrically.
	 */
	private static final int INITIAL_CAPACITY = 8;

	private Grid grid;

	/**
	 * The indices of the i atoms in this cell, held as a primitive array to avoid the boxing (and the
	 * pointer chasing it entails) of a Collection of Integers: these are read in the innermost loop of
	 * the contact calculation. Only the first {@link #numIindices} elements are meaningful.
	 */
	private int[] iIndices;
	private int numIindices;

	/**
	 * The indices of the j atoms in this cell. See {@link #iIndices}.
	 */
	private int[] jIndices;
	private int numJindices;

	public GridCell(Grid parent){
		iIndices = EMPTY;
		jIndices = EMPTY;
		this.grid = parent;
	}

	public void addIindex(int serial){
		if (numIindices == iIndices.length) {
			iIndices = Arrays.copyOf(iIndices, numIindices == 0 ? INITIAL_CAPACITY : numIindices * 2);
		}
		iIndices[numIindices++] = serial;
	}

	public void addJindex(int serial){
		if (numJindices == jIndices.length) {
			jIndices = Arrays.copyOf(jIndices, numJindices == 0 ? INITIAL_CAPACITY : numJindices * 2);
		}
		jIndices[numJindices++] = serial;
	}

	public int getNumIindices() {
		return numIindices;
	}

	public int getNumJindices() {
		return numJindices;
	}

	/**
	 * Calculates all distances of atoms within this cell returning those that are within the given cutoff
	 * as a list of Contacts containing the indices of the pair and the calculated distance.
	 *
	 * If {@link Grid#getJAtoms()} is null, distances are within the iAtoms only
	 * @return
	 */
	public List<Contact> getContactsWithinCell(){

		List<Contact> contacts = new ArrayList<>();

		Point3d[] iAtoms = grid.getIAtoms();
		Point3d[] jAtoms = grid.getJAtoms();
		// we compare squared distances to the squared cutoff, so that the expensive square root is
		// only computed for the pairs that are actually in contact (the large majority are not)
		double cutoffSq = grid.getCutoffSq();

		if (jAtoms==null) {
			for (int a=0; a<numIindices; a++) {
				int i = iIndices[a];
				Point3d atomI = iAtoms[i];
				for (int b=0; b<numIindices; b++) {
					int j = iIndices[b];
					if (j>i) {
						double distanceSq = atomI.distanceSquared(iAtoms[j]);
						if (distanceSq<cutoffSq) contacts.add(new Contact(i, j, Math.sqrt(distanceSq)));
					}
				}
			}

		} else {
			for (int a=0; a<numIindices; a++) {
				int i = iIndices[a];
				Point3d atomI = iAtoms[i];
				for (int b=0; b<numJindices; b++) {
					int j = jIndices[b];
					double distanceSq = atomI.distanceSquared(jAtoms[j]);
					if (distanceSq<cutoffSq) contacts.add(new Contact(i, j, Math.sqrt(distanceSq)));
				}
			}
		}

		return contacts;
	}

	/**
	 * Calculates all distances of atoms between this cell and the given cell returning those that are
	 * within the given cutoff as a list of Contacts containing the indices of the pair and the calculated distance.
	 *
	 * @param otherCell
	 * @return
	 */
	public List<Contact> getContactsToOtherCell(GridCell otherCell){

		List<Contact> contacts = new ArrayList<>();

		Point3d[] iAtoms = grid.getIAtoms();
		Point3d[] jAtoms = grid.getJAtoms();
		// we compare squared distances to the squared cutoff, so that the expensive square root is
		// only computed for the pairs that are actually in contact (the large majority are not)
		double cutoffSq = grid.getCutoffSq();


		if (jAtoms==null) {

			int[] otherIndices = otherCell.iIndices;
			int otherNum = otherCell.numIindices;
			for (int a=0; a<numIindices; a++) {
				int i = iIndices[a];
				Point3d atomI = iAtoms[i];
				for (int b=0; b<otherNum; b++) {
					int j = otherIndices[b];
					if (j>i) {
						double distanceSq = atomI.distanceSquared(iAtoms[j]);
						if (distanceSq<cutoffSq) contacts.add(new Contact(i, j, Math.sqrt(distanceSq)));
					}
				}
			}

		} else {

			int[] otherIndices = otherCell.jIndices;
			int otherNum = otherCell.numJindices;
			for (int a=0; a<numIindices; a++) {
				int i = iIndices[a];
				Point3d atomI = iAtoms[i];
				for (int b=0; b<otherNum; b++) {
					int j = otherIndices[b];
					double distanceSq = atomI.distanceSquared(jAtoms[j]);
					if (distanceSq<cutoffSq) contacts.add(new Contact(i, j, Math.sqrt(distanceSq)));
				}
			}

		}

		return contacts;
	}

	/**
	 * Tests whether any atom in this cell has a contact with the specified query atom
	 * @param iAtoms the first set of atoms to which the iIndices correspond
	 * @param jAtoms the second set of atoms to which the jIndices correspond, or null
	 * @param query test point
	 * @param cutoff
	 * @return
	 */
	public boolean hasContactToAtom(Point3d[] iAtoms, Point3d[] jAtoms, Point3d query, double cutoff) {
		// only the comparison matters here, so we can stay in squared distance space and avoid square roots altogether
		double cutoffSq = cutoff * cutoff;
		for (int a=0; a<numIindices; a++) {
			double distanceSq = iAtoms[iIndices[a]].distanceSquared(query);
			if( distanceSq<cutoffSq)
				return true;
		}
		if (jAtoms!=null) {
			for (int a=0; a<numJindices; a++) {
				double distanceSq = jAtoms[jIndices[a]].distanceSquared(query);
				if( distanceSq<cutoffSq)
					return true;
			}
		}
		return false;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return String.format("GridCell [%d iAtoms,%d jAtoms]", numIindices, numJindices);
	}


}
