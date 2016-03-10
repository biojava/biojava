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
package org.biojava.nbio.structure.io.util;

/**
 * <b>Internal use only. Do not use this class.</b>
 *
 * A class for temporary storing of data when parsing PDB Files.
 *
 * Why is this here? Because making private inner classes in PDBFileParser would
 * make that hulking giant of a class even more gigantic.
 *
 * @author Ulysse Carion
 */
public class PDBTemporaryStorageUtils {
	/**
	 * Temporary data storage for LINK records. This is necessary because LINK
	 * records precede the atoms they correspond to in a PDB file, so we must
	 * store the information encoded in a LINK record until we actually know
	 * about the atoms a LINK refers to.
	 *
	 * @author Ulysse Carion
	 */
	public static class LinkRecord {
		private String name1;
		private String altLoc1;
		private String resName1;
		private String chainID1;
		private String resSeq1;
		private String iCode1;

		private String name2;
		private String altLoc2;
		private String resName2;
		private String chainID2;
		private String resSeq2;
		private String iCode2;

		private String sym1;
		private String sym2;

		public LinkRecord(String name1, String altLoc1, String resName1,
				String chainID1, String resSeq1, String iCode1, String name2,
				String altLoc2, String resName2, String chainID2,
				String resSeq2, String iCode2, String sym1, String sym2) {
			this.name1 = name1;
			this.altLoc1 = altLoc1;
			this.resName1 = resName1;
			this.chainID1 = chainID1;
			this.resSeq1 = resSeq1;
			this.iCode1 = iCode1;
			this.name2 = name2;
			this.altLoc2 = altLoc2;
			this.resName2 = resName2;
			this.chainID2 = chainID2;
			this.resSeq2 = resSeq2;
			this.iCode2 = iCode2;
			this.sym1 = sym1;
			this.sym2 = sym2;
		}

		public String getName1() {
			return name1;
		}

		public String getAltLoc1() {
			return altLoc1;
		}

		public String getResName1() {
			return resName1;
		}

		public String getChainID1() {
			return chainID1;
		}

		public String getResSeq1() {
			return resSeq1;
		}

		public String getiCode1() {
			return iCode1;
		}

		public String getName2() {
			return name2;
		}

		public String getAltLoc2() {
			return altLoc2;
		}

		public String getResName2() {
			return resName2;
		}

		public String getChainID2() {
			return chainID2;
		}

		public String getResSeq2() {
			return resSeq2;
		}

		public String getiCode2() {
			return iCode2;
		}

		public String getSym1() {
			return sym1;
		}

		public String getSym2() {
			return sym2;
		}

		@Override
		public String toString() {
			String s = "[LINK Record:\n";

			s += "Atom 1:\n";
			s += "\tName: " + name1 + "\n";
			s += "\tAlt Loc: " + altLoc1 + "\n";
			s += "\tRes name: " + resName1 + "\n";
			s += "\tChain ID: " + chainID1 + "\n";
			s += "\tRes Seq: " + resSeq1 + "\n";
			s += "\tIns Code: " + iCode1 + "\n";

			s += "Atom 2:\n";
			s += "\tName: " + name2 + "\n";
			s += "\tAlt Loc: " + altLoc2 + "\n";
			s += "\tRes name: " + resName2 + "\n";
			s += "\tChain ID: " + chainID2 + "\n";
			s += "\tRes Seq: " + resSeq2 + "\n";
			s += "\tIns Code: " + iCode2 + "\n";

			s += "]";

			return s;
		}
	}
}
