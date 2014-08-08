package org.biojava.bio.structure.contact;

import java.util.HashMap;
import java.util.Iterator;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;


public class GroupContactSet implements Iterable<GroupContact>{ 

	private HashMap<Pair<ResidueNumber>, GroupContact> contacts;
	
	public GroupContactSet() {
		contacts = new HashMap<Pair<ResidueNumber>, GroupContact>();
	}
	
	public GroupContactSet(AtomContactSet atomContacts) {
		contacts = new HashMap<Pair<ResidueNumber>, GroupContact>();
		atoms2groups(atomContacts);
	}
	
	private void atoms2groups(AtomContactSet atomContacts) {
		
		
		for (AtomContact atomContact:atomContacts) {

			Pair<Atom> atomPair = atomContact.getPair();

			Group iResidue = atomPair.getFirst().getGroup();
			Group jResidue = atomPair.getSecond().getGroup();
			
			// we skip the self-residue contacts
			if (iResidue.equals(jResidue)) continue;
			
			Pair<Group> residuePair = new Pair<Group> (iResidue, jResidue);
			Pair<ResidueNumber> pair = new Pair<ResidueNumber>(iResidue.getResidueNumber(), jResidue.getResidueNumber());
			
			if (!contacts.containsKey(pair)) {
				
				GroupContact groupContact = new GroupContact();
				groupContact.setPair(residuePair);
				groupContact.addAtomContact(atomContact);
				
				contacts.put(pair, groupContact);
				
			} else {
				
				GroupContact groupContact = contacts.get(pair);
				
				groupContact.addAtomContact(atomContact);
				
			}
			
		}
	}

	public void add(GroupContact groupContact) {
		contacts.put(getResNumberPairFromContact(groupContact),groupContact);
	}
	
	public boolean hasContact(Group group1, Group group2) {
		return contacts.containsKey(
				new Pair<ResidueNumber>(group1.getResidueNumber(),group2.getResidueNumber()));
	}
	
	/**
	 * Returns the corresponding GroupContact or null if no contact exists between the 2 given groups
	 * @param group1
	 * @param group2
	 * @return
	 */
	public GroupContact getContact(Group group1, Group group2) {
		return contacts.get(
				new Pair<ResidueNumber>(group1.getResidueNumber(),group2.getResidueNumber()));
	}
	
	public int size() {
		return contacts.size();
	}

	@Override
	public Iterator<GroupContact> iterator() {
		return contacts.values().iterator();
	}
	
	private Pair<ResidueNumber> getResNumberPairFromContact(GroupContact groupContact) {
		return new Pair<ResidueNumber>(groupContact.getPair().getFirst().getResidueNumber(),groupContact.getPair().getSecond().getResidueNumber());
		
	}
}
