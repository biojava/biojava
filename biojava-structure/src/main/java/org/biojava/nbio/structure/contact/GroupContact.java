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

import org.biojava.nbio.structure.Group;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * A pair of residues that are in contact
 * @author duarte_j
 *
 */
public class GroupContact implements Serializable {


	private static final long serialVersionUID = 1L;

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

	public int getNumAtomContacts() {
		return atomContacts.size();
	}

	public List<AtomContact> getAtomContacts() {
		return atomContacts;
	}

	/**
	 * Returns the list of atom contacts in this GroupContact that are within the given distance.
	 * @param distance
	 * @return
	 */
	public List<AtomContact> getContactsWithinDistance(double distance) {

		List<AtomContact> list = new ArrayList<AtomContact>();
		for (AtomContact contact:this.atomContacts) {
			if (contact.getDistance()<distance) {
				list.add(contact);
			}
		}
		return list;
	}

	@Override
	public String toString() {
		return pair.getFirst().getResidueNumber()+","+pair.getSecond().getResidueNumber();
	}
}
