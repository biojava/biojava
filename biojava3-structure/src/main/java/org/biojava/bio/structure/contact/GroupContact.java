package org.biojava.bio.structure.contact;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Group;

/**
 * A pair of residues that are in contact
 * @author duarte_j
 *
 */
public class GroupContact {

	private Pair<Group> pair;
	
	private List<AtomContact> atomContacts;
	
	public GroupContact() {
		atomContacts = new ArrayList<AtomContact>();
	}
	
	public void addAtomContact(AtomContact atomContact) {
		atomContacts.add(atomContact);
	}

	public Pair<Group> getPair() {
		return pair;
	}
	
	public void setPair(Pair<Group> pair) {
		this.pair = pair;
	}
	
	public double getMinDistance() {
		if (atomContacts.size()==0) return 0;
		
		double minDistance = Double.MAX_VALUE;
		for (AtomContact atomContact:atomContacts) {
			if (atomContact.getDistance()<minDistance) 
				minDistance = atomContact.getDistance();
		}
		return minDistance;
	}
	
	@Override
	public String toString() {
		return pair.getFirst().getResidueNumber()+","+pair.getSecond().getResidueNumber();
	}
}
