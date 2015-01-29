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
 * Created on August 13, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment.io;

import java.util.HashSet;
import java.util.Set;

import org.biojava3.alignment.io.StockholmStructure.DatabaseReference;

/**
 * Stores all the content parsed from the #=GS lines
 * 
 * @since 3.0.5
 * @author Amr AL-Hossary
 * @author Marko Vaz
 * 
 */
class StockholmSequenceAnnotation {
	private String accessionNumber;
	private CharSequence description;

	private Set<DatabaseReference> dbReferences;
	private String organism;
	private String organismClassification;
	private String look;

	
	public String getDescription() {
		return description.toString();
	}
	public void setDescription(CharSequence description) {
		this.description = description;
	}
	public void addToDescription(CharSequence description) {
		if (this.description == null) {
			this.description = new StringBuffer(description);
		} else if (this.description instanceof StringBuffer){
			((StringBuffer) this.description).append(description);
		}else {
			this.description = new StringBuffer(this.description).append(description);
		}
	}

	public Set<DatabaseReference> getDbReferences() {
		return dbReferences;
	}

	public void setDbReferences(Set<DatabaseReference> dbReferences) {
		this.dbReferences = dbReferences;
	}
	/**
	 * @param dbReference the string without the initial annotation identifier ( #=GS DR )
	 */
	public void addDBReference(String dbReferenceRepresentingString) {
		if (this.dbReferences == null) {
			this.dbReferences = new HashSet<DatabaseReference>();
		} 
		dbReferences.add(new DatabaseReference(dbReferenceRepresentingString));
	}

	
	
	
	public String getAccessionNumber() {
		return accessionNumber;
	}
	public void setAccessionNumber(String accessionNumber) {
		this.accessionNumber = accessionNumber;
	}
	public String getOrganism() {
		return organism;
	}
	public void setOrganism(String organism) {
		this.organism = organism;
	}
	public String getOrganismClassification() {
		return organismClassification;
	}
	public void setOrganismClassification(String organismClassification) {
		this.organismClassification = organismClassification;
	}
	public String getLook() {
		return look;
	}
	public void setLook(String look) {
		this.look = look;
	}

}