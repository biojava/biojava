package org.biojava.bio.structure.contact;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;


/**
 * A grid cell to be used in contact calculation via geometric hashing algorithm.
 * 
 * @author duarte_j
 *
 */
public class GridCell {
	
	
	private ArrayList<Integer> iIndices;
	private ArrayList<Integer> jIndices;
	
	public GridCell(){
		iIndices = new ArrayList<Integer>();
		jIndices = new ArrayList<Integer>();
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
	 * as a list of AtomContacts
	 * @param iAtoms the first set of atoms to which the iIndices correspond
	 * @param jAtoms the second set of atoms to which the jIndices correspond, if null distances are within the iAtoms only
	 * @param cutoff
	 * @return
	 */
	public List<AtomContact> getContactsWithinCell(Atom[] iAtoms, Atom[] jAtoms, double cutoff){
		
		List<AtomContact> contacts = new ArrayList<AtomContact>();
		
		if (jAtoms==null) {
			for (int i:iIndices) {
				for (int j:iIndices) {
					if (j>i) {
						double distance = Calc.getDistance(iAtoms[i], iAtoms[j]);
						if (distance<cutoff) contacts.add(new AtomContact(new Pair<Atom>(iAtoms[i],iAtoms[j]),distance));
					}
				}
			}
			
		} else {
			for (int i:iIndices) {
				for (int j:jIndices) {
					double distance = Calc.getDistance(iAtoms[i], jAtoms[j]);
					if (distance<cutoff) contacts.add(new AtomContact(new Pair<Atom>(iAtoms[i],jAtoms[j]),distance));
				}
			}
		}

		return contacts;
	}
	
	/**
	 * Calculates all distances of atoms between this cell and the given cell returning those that are
	 * within the given cutoff as a list of AtomContacts
	 * @param otherCell
	 * @param iAtoms the first set of atoms to which the iIndices correspond
	 * @param jAtoms the second set of atoms to which the jIndices correspond, if null distances are within the iAtoms only
	 * @param cutoff
	 * @return
	 */
	public List<AtomContact> getContactsToOtherCell(GridCell otherCell , Atom[] iAtoms, Atom[] jAtoms, double cutoff){
		
		List<AtomContact> contacts = new ArrayList<AtomContact>();
		
		if (jAtoms==null) {
			
			for (int i:iIndices) {
				for (int j:otherCell.iIndices) {			
					if (j>i) {
						double distance = Calc.getDistance(iAtoms[i], iAtoms[j]);
						if (distance<cutoff) contacts.add(new AtomContact(new Pair<Atom>(iAtoms[i],iAtoms[j]),distance));
					}
				}
			}
			
		} else {
			
			for (int i:iIndices) {
				for (int j:otherCell.jIndices) {
					double distance = Calc.getDistance(iAtoms[i], jAtoms[j]);
					if (distance<cutoff) contacts.add(new AtomContact(new Pair<Atom>(iAtoms[i],jAtoms[j]),distance));
				}
			}
			
		}
		
		return contacts;
	}
	
}
