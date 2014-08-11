package org.biojava.bio.structure.contact;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;

import org.biojava.bio.structure.Atom;

/**
 * A set of atom-atom contacts to hold the results of intra and inter-chain contact calculations
 * 
 * @author duarte_j
 *
 */
public class AtomContactSet implements Iterable<AtomContact> {

	private HashMap<Pair<AtomIdentifier>, AtomContact> contacts;
	
	public AtomContactSet() {
		this.contacts = new HashMap<Pair<AtomIdentifier>,AtomContact>();
	}
	
	public void add(AtomContact contact) {
		this.contacts.put(getAtomIdPairFromContact(contact), contact);
	}
	
	public void addAll(Collection<AtomContact> list) {
		for (AtomContact contact:list) {
			this.contacts.put(getAtomIdPairFromContact(contact), contact);
		}
	}
	
	public boolean hasContact(Atom atom1, Atom atom2) {
		return contacts.containsKey(new Pair<AtomIdentifier>(
					new AtomIdentifier(atom1.getPDBserial(),atom1.getGroup().getChainId()),
					new AtomIdentifier(atom2.getPDBserial(),atom2.getGroup().getChainId()) ));
	}

	/**
	 * Returns the corresponding AtomContact or null if no contact exists between the 2 given atoms
	 * @param atom1
	 * @param atom2
	 * @return
	 */
	public AtomContact getContact(Atom atom1, Atom atom2) {
		return contacts.get(new Pair<AtomIdentifier>(
				new AtomIdentifier(atom1.getPDBserial(),atom1.getGroup().getChainId()),
				new AtomIdentifier(atom2.getPDBserial(),atom2.getGroup().getChainId()) ));		
	}
	
	public int size() {
		return contacts.size();
	}
	
	@Override
	public Iterator<AtomContact> iterator() {
		return contacts.values().iterator();
	}
	
	private Pair<AtomIdentifier> getAtomIdPairFromContact(AtomContact contact) {
		Pair<AtomIdentifier> pair = new Pair<AtomIdentifier>(
				new AtomIdentifier(contact.getPair().getFirst().getPDBserial(),contact.getPair().getFirst().getGroup().getChainId()), 
				new AtomIdentifier(contact.getPair().getSecond().getPDBserial(),contact.getPair().getSecond().getGroup().getChainId()));
		
		return pair;
	}
}
