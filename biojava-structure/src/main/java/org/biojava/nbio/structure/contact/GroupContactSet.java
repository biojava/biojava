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

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * A set of residue-residue contacts 
 * 
 * @author duarte_j
 *
 */
public class GroupContactSet implements Iterable<GroupContact>{ 

	private HashMap<Pair<ResidueNumber>, GroupContact> contacts;
	
	/**
	 * A cached HashSet to be used only if hasContact(Pair<ResidueIdentifier>) is called
	 */
	private HashSet<Pair<ResidueIdentifier>> residueIdContacts;
	
	public GroupContactSet() {
		contacts = new HashMap<Pair<ResidueNumber>, GroupContact>();
	}
	
	/**
	 * Constructs a <code>GroupContactSet</code> by collapsing the given <code>AtomContactSet</code> into
	 * residue-residue (group-group) contacts. 
	 * @param atomContacts
	 */
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
	
	/**
	 * Tell whether the given group pair is a contact in this GroupContactSet,
	 * the comparison is done by matching residue numbers and chain identifiers
	 * @param group1
	 * @param group2
	 * @return
	 */
	public boolean hasContact(Group group1, Group group2) {
		return hasContact(group1.getResidueNumber(),group2.getResidueNumber());
	}
	
	/**
	 * Tell whether the given pair is a contact in this GroupContactSet,
	 * the comparison is done by matching residue numbers and chain identifiers
	 * @param resNumber1
	 * @param resNumber2
	 * @return
	 */
	public boolean hasContact(ResidueNumber resNumber1, ResidueNumber resNumber2) {
		return contacts.containsKey(new Pair<ResidueNumber>(resNumber1, resNumber2));
	}
	
	/**
	 * Tell whether the given pair is a contact in this GroupContactSet,
	 * in a chain-identifier independent way: contacts happening between different copies of 
	 * the same Compound(Entity) will be considered equal as long as they have the same
	 * residue numbers.
	 * @param resId1
	 * @param resId2
	 * @return
	 */
	public boolean hasContact(ResidueIdentifier resId1, ResidueIdentifier resId2) {
		if (residueIdContacts == null) {
			initResidueIdContacts();
		}
		
		return residueIdContacts.contains(new Pair<ResidueIdentifier>(resId1, resId2));		
	}
	
	private void initResidueIdContacts() {
		residueIdContacts = new HashSet<Pair<ResidueIdentifier>>();
		for (Pair<ResidueNumber> pairResNum:contacts.keySet()) {
			ResidueNumber resNumFirst = pairResNum.getFirst();
			ResidueNumber resNumSecond = pairResNum.getSecond();
			residueIdContacts.add(new Pair<ResidueIdentifier>(
					new ResidueIdentifier(resNumFirst.getSeqNum(), resNumFirst.getInsCode()),
					new ResidueIdentifier(resNumSecond.getSeqNum(), resNumSecond.getInsCode())
					) ); 
		}
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
