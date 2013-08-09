package org.biojava.bio.structure.io.util;

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
	}
}
