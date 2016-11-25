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
import java.util.List;

import javax.vecmath.Point3d;


/**
 * A grid cell to be used in contact calculation via spatial hashing algorithm.
 *
 * @author Jose Duarte
 *
 */
public class GridCell {


	private Grid grid;
	private ArrayList<Integer> iIndices;
	private ArrayList<Integer> jIndices;

	public GridCell(Grid parent){
		iIndices = new ArrayList<Integer>();
		jIndices = new ArrayList<Integer>();
		this.grid = parent;
	}

	public void addIindex(int serial){
		iIndices.add(serial);
	}

	public void addJindex(int serial){
		jIndices.add(serial);
	}

	public int getNumIindices() {
		return iIndices.size();
	}

	public int getNumJindices() {
		return jIndices.size();
	}

	/**
	 * Calculates all distances of atoms within this cell returning those that are within the given cutoff
	 * as a list of Contacts containing the indices of the pair and the calculated distance.
	 * 
	 * If {@link Grid#getJAtoms()} is null, distances are within the iAtoms only
	 * @return
	 */
	public List<Contact> getContactsWithinCell(){

		List<Contact> contacts = new ArrayList<Contact>();
		
		Point3d[] iAtoms = grid.getIAtoms();
		Point3d[] jAtoms = grid.getJAtoms();
		double cutoff = grid.getCutoff();

		if (jAtoms==null) {
			for (int i:iIndices) {
				for (int j:iIndices) {
					if (j>i) {
						double distance = iAtoms[i].distance(iAtoms[j]);
						if (distance<cutoff) contacts.add(new Contact(i, j, distance));
					}
				}
			}

		} else {
			for (int i:iIndices) {
				for (int j:jIndices) {
					double distance = iAtoms[i].distance(jAtoms[j]);
					if (distance<cutoff) contacts.add(new Contact(i, j, distance));
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
	 * @param iAtoms the first set of atom coordinates to which the iIndices correspond
	 * @param jAtoms the second set of atom coordinates to which the jIndices correspond, if null distances are within the iAtoms only
	 * @param cutoff
	 * @return
	 */
	public List<Contact> getContactsToOtherCell(GridCell otherCell){

		List<Contact> contacts = new ArrayList<Contact>();

		Point3d[] iAtoms = grid.getIAtoms();
		Point3d[] jAtoms = grid.getJAtoms();
		double cutoff = grid.getCutoff();

		
		if (jAtoms==null) {

			for (int i:iIndices) {
				for (int j:otherCell.iIndices) {
					if (j>i) {
						double distance = iAtoms[i].distance(iAtoms[j]);
						if (distance<cutoff) contacts.add(new Contact(i, j, distance));
					}
				}
			}

		} else {

			for (int i:iIndices) {
				for (int j:otherCell.jIndices) {
					double distance = iAtoms[i].distance(jAtoms[j]);
					if (distance<cutoff) contacts.add(new Contact(i, j, distance));
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
		for( int i : iIndices ) {
			double distance = iAtoms[i].distance(query);
			if( distance<cutoff)
				return true;
		}
		if (jAtoms!=null) {
			for( int i : jIndices ) {
				double distance = jAtoms[i].distance(query);
				if( distance<cutoff)
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
		return String.format("GridCell [%d iAtoms,%d jAtoms]",iIndices.size(),jIndices==null?"-":jIndices.size());
	}

	
}
