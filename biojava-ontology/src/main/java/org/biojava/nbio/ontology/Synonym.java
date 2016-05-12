/*
 *                  BioJava development code
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
 * Created on Jan 24, 2008
 *
 */

package org.biojava.nbio.ontology;

import java.util.Comparator;


public class Synonym implements Comparable<Synonym>{


	public final static int UNKNOWN_SCOPE = -1;
	public final static int RELATED_SYNONYM = 0;
	public final static int EXACT_SYNONYM = 1;
	public final static int NARROW_SYNONYM = 2;
	public final static int BROAD_SYNONYM = 3;

	int scope;
	String category;
	String name;

	@Override
	public String toString(){
		String txt = "Synonym: name:"+name+ " category:" + category + " scope: " +scope;
		return txt;
	}

	public final static Comparator<Synonym> COMPARATOR = new Comparator<Synonym>() {
		@Override
		public int compare(Synonym a, Synonym b) {
			if (a == null && b == null)
				return 0;
			else if (a == null)
				return -1;
			else if (b == null)
				return 1;
			else {
				if ((a.getCategory() == null) && (b.getCategory() == null))
					return 0;
				else if ( a.getCategory()==null)
					return -1;
				else if ( b.getCategory()==null)
					return 1;

				return a.getCategory().compareToIgnoreCase(
						b.getCategory());
			}
		}
	};

	public Synonym() {
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getCategory() {
		return category;
	}
	public void setCategory(String category) {
		this.category = category;
	}
	public int getScope() {
		return scope;
	}
	public void setScope(int scope) {
		this.scope = scope;
	}
	@Override
	public int compareTo(Synonym o) {
		return COMPARATOR.compare(this, o);
	}
}
