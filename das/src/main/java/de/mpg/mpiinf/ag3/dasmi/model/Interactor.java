/* 
 * 
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA, or see the FSF site: http://www.fsf.org.
 * 
 */
package de.mpg.mpiinf.ag3.dasmi.model;

import java.util.List;
import java.util.ArrayList;

import de.mpg.mpiinf.ag3.dasmi.Constants;

/**
 * The interactor, i.e. a protein, domain etc. Contains fields to describe the basic entities of
 * an interactor, like name or database accesion, and optional details for further description
 * 
 * @author Hagen Blankenburg, Max Planck Institute for Informatics
 *
 */
public class Interactor {
	
	private int intId;
	private String id = null;
	private String name = null;
	private String dbSource = null;
	private String dbSourceCvId = null;
	private String dbVersion = null;
	private String dbAccessionId = null;
	private String dbCoordSys = null;
	private String description = null;
	
	private Sequence sequence = null;
	private List <Detail>details = new ArrayList<Detail>();
	
	/**
	 * Basic initialization
	 *
	 */
	public Interactor(){
		this.intId = Constants.INVALID_INT;
		this.name = Constants.INVALID_STRING;
		this.dbSource = Constants.INVALID_STRING;
		this.dbSourceCvId = Constants.INVALID_STRING;
		this.dbVersion = Constants.INVALID_STRING;
		this.dbAccessionId = Constants.INVALID_STRING;
		this.dbCoordSys = Constants.INVALID_STRING;
		this.description = Constants.INVALID_STRING;
	}
	

	/**
	 * @return the id
	 */
	public String getId() {
		return id;
	}

	/**
	 * @param id the id to set
	 */
	public void setId(String id) {
		this.id = id;
	}
	
	/**
	 * @return the details
	 */
	public List getDetails() {
		return details;
	}



	/**
	 * @param details the details to set
	 */
	public void setDetails(List<Detail> details) {
		this.details = details;
	}
	
	
	/**
	 * Adds a detail to the list
	 * @param detail the detail to add
	 */
	public void addDetail(Detail detail){
		this.details.add(detail);
	}



	/**
	 * @return the sequence
	 */
	public Sequence getSequence() {
		return sequence;
	}



	/**
	 * @param sequence the sequence to set
	 */
	public void setSequence(Sequence sequence) {
		this.sequence = sequence;
	}


	/**
	 * @return the dbAccessionId
	 */
	public String getDbAccessionId() {
		return dbAccessionId;
	}

	/**
	 * @param dbAccessionId the dbAccessionId to set
	 */
	public void setDbAccessionId(String dbAccessionId) {
		this.dbAccessionId = dbAccessionId;
	}

	/**
	 * @return the dbCoordSys
	 */
	public String getDbCoordSys() {
		return dbCoordSys;
	}

	/**
	 * @param dbCoordSys the dbCoordSys to set
	 */
	public void setDbCoordSys(String dbCoordSys) {
		this.dbCoordSys = dbCoordSys;
	}

	/**
	 * @return the dbSource
	 */
	public String getDbSource() {
		return dbSource;
	}

	/**
	 * @param dbSource the dbSource to set
	 */
	public void setDbSource(String dbSource) {
		this.dbSource = dbSource;
	}

	/**
	 * @return the dbSourceCvId
	 */
	public String getDbSourceCvId() {
		return dbSourceCvId;
	}

	/**
	 * @param dbSourceCvId the dbSourceCvId to set
	 */
	public void setDbSourceCvId(String dbSourceCvId) {
		this.dbSourceCvId = dbSourceCvId;
	}

	/**
	 * @return the dbVersion
	 */
	public String getDbVersion() {
		return dbVersion;
	}

	/**
	 * @param dbVersion the dbVersion to set
	 */
	public void setDbVersion(String dbVersion) {
		this.dbVersion = dbVersion;
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * @return the intId
	 */
	public int getIntId() {
		return intId;
	}
	
	/**
	 * 
	 * @return the description
	 */
	public String getDescription(){
		return this.description;
	}
	
	
	/**
	 * 
	 * @param description the description to set
	 */
	public void setDescription(String description){
		this.description = description;
	}
	
	/**
	 * Checks whether two interactors are equal, ie. have the same accession, source,
	 * version and coord sys
	 * @param check The interactor to compare with
	 * @return flag
	 * @returntrue if both interactors are the same, false otherwise
	 */
	public boolean equals(Interactor check){
		//TODO: if  the equals method is implemented
		// one also should implement the hashCode method
		// see also http://www.geocities.com/technofundo/tech/java/equalhash.html
		
		
		if (this.dbAccessionId == check.dbAccessionId &&
			this.dbCoordSys == check.dbCoordSys &&
			this.dbSource == check.dbSource &&
			this.dbSourceCvId == check.dbSourceCvId &&
			this.dbVersion == check.dbVersion){
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * Basic to string
	 * @return the output string
	 */
	public String toString(){
		String help = "Name : " + this.name + " - "
		+ "dbSource: " + this.dbSource + " - "
		+ "dbSourceCvId: " + this.dbSourceCvId + " - "
		+ "AccessionId: " + this.dbAccessionId + " - "
		+ "Version: " + this.dbVersion + "\n";
		if (details != null){
		help += "\t\tDetails: \n";
			for (int i = 0; i < details.size(); i++){
				help += "\t\t\t" + details.get(i).toString() + "\n";
			}
		}
		help += "\n";
		return help;
	}
	
	
	 /**
     * Hashes the properties of the interactor.
     * @return the hased string
     */
    public String hashify(){
    	//TODO: should this actually be the hashCode method? AP
    	// then one should return hash.hashCode() ...
    	
    	String hash = this.dbAccessionId + this.dbSource + this.dbVersion;
    	
    	return hash.toLowerCase();
    }
	

}
