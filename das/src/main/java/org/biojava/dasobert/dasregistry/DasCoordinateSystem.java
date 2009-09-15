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
 * Created on 15.04.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.dasobert.dasregistry;

/**
 * a Bean to be returned via SOAP. It takes care of the DAS - coordinate Systems
 * 
 * @author Andreas Prlic
 */
public class DasCoordinateSystem {

	String name;// name also seems to be Authority?
	String category;
	String organism_name;
	int ncbi_tax_id;
	String uniqueId;
	String version;
	String testCode;

	public DasCoordinateSystem() {
		uniqueId = "";
		name = "";
		category = "";
		organism_name = "";
		ncbi_tax_id = 0;
		version = "";
		testCode = "";
	}

	public boolean equals(DasCoordinateSystem other) {
		boolean match = true;
		//System.out.println("comparing in coordinate system " + this.toString() );
		//System.out.println("to other toString()="+other.toString());
		// URI has piority
//		if ( !uniqueId.equalsIgnoreCase(
//		other.getUniqueId())){
//			System.out.println("failed on first test unique id in database="+uniqueId+" is not equal to "+other.uniqueId);
//			return false;
//		}

		
		//added by jw for registry
		if (!organism_name.equals(other.getOrganismName())) {
			//System.out.println("mismatch in name |"+organism_name+"| other=|"+other.getOrganismName()+"|");
			match = false;
			return match;
		}
		if (ncbi_tax_id != other.getNCBITaxId()) {
			//System.out.println("mismatch in ncbi tax id " + ncbi_tax_id +
			//" != " + other.getNCBITaxId());
			match = false;
			return match;
		}
		if (!version.equals(other.getVersion())) {
			//System.out.println("mismatch in version");
			match = false;
			return match;
		}
		if (!category.equals(other.getCategory())) {
			//System.out.println("mismatch in category");
			match = false;
			return match;
		}
		if (!name.equals(other.getName())) {
			//System.out.println("mismatch in name");
			match = false;
			return match;
		}
		
		
		//if(match){
		//System.out.println(" match: " + match);
		//}
		//organism_name = "";
		
		
		

		return match;
	}

	public int hashCode() {
		int h = 7;

		h = 31 * h + (null == name ? 0 : name.hashCode());
		h = 31 * h + (null == category ? 0 : category.hashCode());

		return h;
	}

	public Object clone() {
		DasCoordinateSystem d = new DasCoordinateSystem();
		d.setTestCode(testCode);
		d.setCategory(category);
		d.setName(name);
		d.setNCBITaxId(ncbi_tax_id);
		d.setUniqueId(getUniqueId());
		d.setOrganismName(getOrganismName());
		d.setVersion(getVersion());
		return d;
	}

	public String getTestCode() {
		return testCode;
	}

	public void setTestCode(String testCode) {
		if (testCode == null)
			testCode = "";
		this.testCode = testCode.trim();
	}

	public void setUniqueId(String id) {
		uniqueId = id.trim();
	}

	public String getUniqueId() {
		return uniqueId;
	}

	/**
	 * set the name / authority for this coordiante system e.g. UniProt, PDB,
	 * Ensembl, etc.
	 * 
	 * @param n
	 *            the name
	 */
	public void setName(String n) {
		name = n.trim();
	}

	/**
	 * get the name / authority for this coordiante system e.g. UniProt, PDB,
	 * Ensembl, etc.
	 * 
	 * @return the name / authority of this coordinate system
	 * 
	 */
	public String getName() {
		return name;
	}

	public void setCategory(String c) {
		category = c.trim();
	}

	/**
	 * returns the type of the coordinate system. e.g if it is a Chromosomal,
	 * Protein Sequence, Protein Structure, etc.
	 * 
	 * @return the category
	 */
	public String getCategory() {
		return category;
	}

	public void setOrganismName(String t) {
		this.organism_name = t.trim();
	}

	public String getOrganismName() {
		return organism_name;
	}

	public void setNCBITaxId(int id) {
		ncbi_tax_id = id;
	}

	public int getNCBITaxId() {
		return ncbi_tax_id;
	}

	public String getVersion() {
		return version;
	}

	public void setVersion(String version) {
		if (version == null)
			version = "";
		this.version = version.trim();
	}

	public String toString() {
		String nam = name;
		if (!version.equals(""))
			nam += "_" + version;

		if (organism_name.equals(""))
			return nam + "," + category;
		else
			return nam + "," + category + "," + organism_name;
	}

	public static DasCoordinateSystem fromString(String rawString) {
		String[] spl = rawString.split(",");
		DasCoordinateSystem dcs = new DasCoordinateSystem();
		if (spl.length == 2) {
			dcs.setName(spl[0]);
			dcs.setCategory(spl[1]);
		}
		if (spl.length == 3) {
			dcs.setName(spl[0]);
			dcs.setCategory(spl[1]);
			dcs.setOrganismName(spl[2]);
		}
		return dcs;
	}

	
	public String getAuthority() {
		return getName();
	}

	public void setAuthority(String authority) {
		this.name = authority.trim();
	}

}
