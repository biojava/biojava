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

import java.io.Serializable;
import java.util.*;


/**
 * A set of atom-atom contacts to hold the results of intra and inter-chain contact calculations
 * <p>
 * Contacts are keyed by the <i>ordered</i> pair of {@link AtomIdentifier}s of the 2 atoms, i.e. the
 * pair (a,b) and the pair (b,a) are 2 different keys. Thus look-ups ({@link #hasContact(Atom, Atom)},
 * {@link #getContact(Atom, Atom)}) must give the 2 atoms in the same order in which the contact was
 * calculated. The order produced by the calculation ({@link Grid}) is:
 * <ul>
 * <li>for contacts within a single set of atoms (e.g. the intra-chain contacts of
 * <code>StructureTools.getAtomsInContact(Chain, double)</code>), the atom that comes first in the
 * atom array is the first member of the pair</li>
 * <li>for contacts between 2 sets of atoms (e.g. the inter-chain contacts of
 * <code>StructureTools.getAtomsInContact(Chain, Chain, double, boolean)</code>, or a
 * {@link StructureInterface}), the atom belonging to the first set is the first member of the pair</li>
 * </ul>
 * <p>
 * Note that the order is not the order of PDB serials or of any other property of the atoms
 * themselves: it is only the order in which the atoms were given to the calculation.
 *
 * @author duarte_j
 *
 */
public class AtomContactSet implements Serializable, Iterable<AtomContact> {


	private static final long serialVersionUID = 1L;

	/**
	 * The default load factor of a {@link HashMap}, needed to size the map from an expected number of
	 * entries.
	 */
	private static final float DEFAULT_LOAD_FACTOR = 0.75f;

	private HashMap<Pair<AtomIdentifier>, AtomContact> contacts;
	private double cutoff;

	public AtomContactSet(double cutoff) {
		this.cutoff = cutoff;
		this.contacts = new HashMap<>();
	}

	/**
	 * Creates an AtomContactSet sized to hold the given number of contacts, so that the underlying map
	 * doesn't have to be repeatedly resized and rehashed as contacts are added. Contact calculations
	 * produce hundreds of thousands of contacts for a large structure, where the repeated rehashing
	 * that growing from the default capacity entails is a significant part of the cost.
	 * @param cutoff the distance cutoff
	 * @param expectedSize the number of contacts expected to be added. Only affects performance: an
	 *                     over-estimate merely leaves the map larger than it needs to be.
	 */
	public AtomContactSet(double cutoff, int expectedSize) {
		this.cutoff = cutoff;
		this.contacts = new HashMap<>((int) (expectedSize / DEFAULT_LOAD_FACTOR) + 1);
	}

	/**
	 * Adds the given contact to this set, keyed by the ordered pair of its 2 atoms. If a contact
	 * with the same ordered pair of atoms is already present it is replaced.
	 * @param contact the contact to add
	 */
	public void add(AtomContact contact) {
		this.contacts.put(getAtomIdPairFromContact(contact), contact);
	}

	/**
	 * Adds all given contacts to this set, see {@link #add(AtomContact)}.
	 * @param list the contacts to add
	 */
	public void addAll(Collection<AtomContact> list) {
		for (AtomContact contact:list) {
			this.contacts.put(getAtomIdPairFromContact(contact), contact);
		}
	}

	/**
	 * Tells whether a contact exists between the 2 given atoms, <i>in the given order</i>.
	 * <p>
	 * The 2 atoms have to be passed in the same order in which the contacts of this set were
	 * calculated, otherwise this returns false even if the 2 atoms are within the distance cutoff.
	 * See the class documentation for the ordering convention. If the order is not known, both orders
	 * have to be queried.
	 * @param atom1 the first atom of the pair
	 * @param atom2 the second atom of the pair
	 * @return true if the 2 atoms are in contact in the given order, false otherwise
	 * @see #getContact(Atom, Atom)
	 */
	public boolean hasContact(Atom atom1, Atom atom2) {
		return hasContact(
					new AtomIdentifier(atom1.getPDBserial(),atom1.getGroup().getChainId()),
					new AtomIdentifier(atom2.getPDBserial(),atom2.getGroup().getChainId()) );
	}

	/**
	 * Tells whether a contact exists between the 2 given atom identifiers, <i>in the given order</i>,
	 * see {@link #hasContact(Atom, Atom)}.
	 * @param atomId1 the identifier of the first atom of the pair
	 * @param atomId2 the identifier of the second atom of the pair
	 * @return true if the 2 atoms are in contact in the given order, false otherwise
	 */
	public boolean hasContact(AtomIdentifier atomId1, AtomIdentifier atomId2) {
		return contacts.containsKey(new Pair<AtomIdentifier>(atomId1,atomId2));
	}

	/**
	 * Returns the contact between the 2 given atoms <i>in the given order</i>, or null if there is
	 * no such contact in this set.
	 * <p>
	 * As in {@link #hasContact(Atom, Atom)} the order of the 2 atoms matters: they have to be passed
	 * in the same order in which the contacts of this set were calculated, otherwise null is returned
	 * even if the 2 atoms are within the distance cutoff. See the class documentation for the
	 * ordering convention.
	 * @param atom1 the first atom of the pair
	 * @param atom2 the second atom of the pair
	 * @return the contact between the 2 atoms in the given order, or null if there is none
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
		Pair<AtomIdentifier> pair = new Pair<>(
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

		List<AtomContact> list = new ArrayList<>();
		for (AtomContact contact:this.contacts.values()) {
			if (contact.getDistance()<distance) {
				list.add(contact);
			}
		}
		return list;
	}
}
