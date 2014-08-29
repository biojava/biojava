package org.biojava.bio.structure.contact;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.structure.Atom;


/**
 * A set of atom-atom contacts to hold the results of intra and inter-chain contact calculations
 * 
 * @author duarte_j
 *
 */
public class AtomContactSet implements Iterable<AtomContact> {

	private HashMap<Pair<AtomIdentifier>, AtomContact> contacts;
	private double cutoff;
	
	public AtomContactSet(double cutoff) {
		this.cutoff = cutoff;
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
	
	/**
	 * Returns true if at least 1 contact from this set is within the given distance.
	 * Note that if the distance given is larger than the distance cutoff used to 
	 * calculate the contacts then nothing will be found. 
	 * @param distance
	 * @return
	 * @throws IllegalArgumentException if given distance is larger than distance cutoff 
	 * used for calculation of contacts 
	 */
	public boolean hasContactsWithinDistance(double distance) {
		
		if (distance>=cutoff) 
			throw new IllegalArgumentException("Given distance "+
					String.format("%.2f", distance)+" is larger than contacts' distance cutoff "+
					String.format("%.2f", cutoff));
		
		for (AtomContact contact:this.contacts.values()) {
			if (contact.getDistance()<distance) {
				return true;
			}
		}
		return false;		
	}
	
	/**
	 * Returns the list of contacts from this set that are within the given distance.
	 * @param distance
	 * @return
	 * @throws IllegalArgumentException if given distance is larger than distance cutoff 
	 * used for calculation of contacts  
	 */
	public List<AtomContact> getContactsWithinDistance(double distance) {
		
		if (distance>=cutoff) 
			throw new IllegalArgumentException("Given distance "+
					String.format("%.2f", distance)+" is larger than contacts' distance cutoff "+
					String.format("%.2f", cutoff));

		List<AtomContact> list = new ArrayList<AtomContact>();
		for (AtomContact contact:this.contacts.values()) {
			if (contact.getDistance()<distance) {
				list.add(contact);
			}
		}
		return list;				
	}
}
