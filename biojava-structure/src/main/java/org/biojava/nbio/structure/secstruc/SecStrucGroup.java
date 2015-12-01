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
package org.biojava.nbio.structure.secstruc;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.HetatomImpl;

/**
 * Helper Group for the secondary structure prediction algorithm. It provides a
 * faster interface for retrieving the important Atom types and it is applicable
 * to any Group implementation that fulfills the {@link Group#hasAminoAtoms()}
 * condition.
 * 
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
class SecStrucGroup extends HetatomImpl {

	private static final long serialVersionUID = 313490286720467714L;

	private Atom N;
	private Atom CA;
	private Atom C;
	private Atom O;
	private Atom H;

	private Group original;

	public SecStrucGroup() {
		super();
	}

	@Override
	public String toString() {

		StringBuffer str = new StringBuffer("SecStrucGroup ");
		str.append(residueNumber);
		str.append(" ");
		str.append(pdb_name);
		str.append(" ");
		str.append(pdb_flag);
		if (pdb_flag) {
			str.append(" atoms: ");
			str.append(atoms.size());
		}

		return str.toString();
	}

	public Group getOriginal() {
		return original;
	}

	public void setOriginal(Group original) {
		this.original = original;
	}

	public Atom getC() {
		return C;
	}

	public void setC(Atom c) {
		addAtom(c);
		C = c;
	}

	public Atom getCA() {
		return CA;
	}

	public void setCA(Atom ca) {
		addAtom(ca);
		CA = ca;
	}

	public Atom getH() {
		return H;
	}

	public void setH(Atom h) {
		addAtom(h);
		H = h;
	}

	public Atom getN() {
		return N;
	}

	public void setN(Atom n) {
		addAtom(n);
		N = n;
	}

	public Atom getO() {
		return O;
	}

	public void setO(Atom o) {
		addAtom(o);
		O = o;
	}
}