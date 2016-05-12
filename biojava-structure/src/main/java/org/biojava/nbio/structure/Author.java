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

package org.biojava.nbio.structure;

import java.io.Serializable;

/**
 * Describes author attributes for author information in a PDB file.
 * @author Jules Jacobsen
 */
public class Author implements Serializable{

	private static final long serialVersionUID = 4840370515056666418L;
	private String surname = "";
	private String initials = "";


	public String getInitials() {
		return initials;
	}

	public void setInitials(String initials) {
		this.initials = initials;
	}

	public String getSurname() {
		return surname;
	}

	public void setSurname(String surname) {
		this.surname = surname;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final Author other = (Author) obj;
		if ((this.surname == null) ? (other.surname != null) : !this.surname.equals(other.surname)) {
			return false;
		}
		return !((this.initials == null) ? (other.initials != null) : !this.initials.equals(other.initials));
	}

	@Override
	public int hashCode() {
		int hash = 3;
		hash = 19 * hash + (this.surname != null ? this.surname.hashCode() : 0);
		hash = 19 * hash + (this.initials != null ? this.initials.hashCode() : 0);
		return hash;
	}

	@Override
	public String toString() {
		return initials + surname;
	}



}
