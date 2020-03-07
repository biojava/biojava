/*
 * BioJava development code
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
 * Author: Daniel Asarnow
 * Date:   2012-6-23
 */

package org.biojava.nbio.structure.cath;

/** The categories found within CATH.
 *
 * The CATH node types are:
 * 'C' (class), 'A' (architecture), 'T' (topology), 'H' (homologous superfamily),
 * 'S' (sequence family, S35), 'O' (orthologous sequence family, S60),
 * 'L' ("like" sequence family, S95), 'I' (identical, S100) and
 * 'D' (domain, S100 count).
 *
 * @author Daniel Asarnow
 */
public enum CathCategory {
	Class,
	Architecture,
	Topolgy,
	Homology,
	SequenceFamily,
	OrthologousSequenceFamily,
	LikeSequenceFamily,
	IdenticalSequenceFamily,
	DomainCounter;

	static final String lut = "CATHSOLID";

	public static CathCategory fromString(String type) {
        switch (type) {
            case "C":
                return Class;
            case "A":
                return Architecture;
            case "T":
                return Topolgy;
            case "H":
                return Homology;
            case "S":
                return SequenceFamily;
            case "O":
                return OrthologousSequenceFamily;
            case "L":
                return LikeSequenceFamily;
            case "I":
                return IdenticalSequenceFamily;
//        } else if ( type.equals("D") ) {
            default:
                return DomainCounter;
        }
	}

	@Override
	public String toString() {
		switch (this) {
			case Class:
				return "C";
			case Architecture:
				return "A";
			case Topolgy:
				return "T";
			case Homology:
				return "H";
			case SequenceFamily:
				return "S";
			case OrthologousSequenceFamily:
				return "O";
			case LikeSequenceFamily:
				return "L";
			case IdenticalSequenceFamily:
				return "I";
//            case DomainCounter:
			default:
				return "D";
		}
	}

	public static CathCategory fromCathCode(String code) {
		int count = 0;
		int idx = 0;
		while ((idx = code.indexOf(".",idx)) != -1) {
			count++;
			idx++;
		}
		return fromString(lut.substring(count,count+1));
	}

}
