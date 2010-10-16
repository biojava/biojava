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

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

import de.mpg.mpiinf.ag3.dasmi.Constants;

/**
 * An interaction element with several participants and futher interaction details
 * 
 * @author Hagen Blankenburg, Max Planck Institute for Informatics
 *
 */
public class Interaction {
	
	private List<String> sources = null; // contains the names of all sources that report this interaction

	private int id;
	private String name = null;
	private String dbSource = null;
	private String dbSourceCvId = null;
	private String dbVersion = null;
	private String dbAccessionId = null;
	private List<Participant> participants = new ArrayList<Participant>();
	private List<Detail> details = new ArrayList<Detail>();
	
	
	/**
	 * Basic initialzation
	 *
	 */
	public Interaction(){
		this.id = Constants.INVALID_INT;
		this.sources = new ArrayList<String>();
		this.name = Constants.INVALID_STRING;
		this.dbSource = Constants.INVALID_STRING;
		this.dbSourceCvId = Constants.INVALID_STRING;
		this.dbVersion = Constants.INVALID_STRING;
		this.dbAccessionId = Constants.INVALID_STRING;
	}
		
	
	/**
	 * 
	 * @return the id
	 */
	public int getId(){
		return this.id;
	}
	
	/**
	 * 
	 * @param id the id to set
	 */
	public void setId(int id){
		this.id = id;
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
	 * @return the details
	 */
	public List<Detail> getDetails() {
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
	 * @return the interactors
	 */
	public List<Participant> getParticipants() {
		return participants;
	}
	
	
	/**
	 * @param interactors the interactors to set
	 */
	public void setParticipants(List<Participant> participants) {
		this.participants = participants;
	}
	
	
	/**
	 * Adds an interactor ref to the list
	 * @param ref
	 */
	public void addParticipant(Participant participant){
		this.participants.add(participant);
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
	 * 
	 * @return the sources
	 */
	public List<String> getSources(){
		return this.sources;
	}
	
	/**
	 * @param sources the sources to set
	 */
	public void setSources(List<String> sources){
		this.sources = sources;
	}
	
	/**
	 * Adds a source to the sources list. If it is already present, nothing happens
	 * @param sourceId The id of the source to add
	 */
	public void addSource(String sourceId){
		sourceId = sourceId.toLowerCase();
		if (!this.sources.contains(sourceId)){
			this.sources.add(sourceId);
		}
	}
	
	/**
	 * Removes a source from the sources list
	 * @param sourceId the id of the source that is to be removed
	 */
	public void removeSource(String sourceId){
		sourceId = sourceId.toLowerCase();
		this.sources.remove(sourceId);	
	}
	
	/**
	 * Checks whether a certain source is alredy present
	 * @param sourceId The if of the source to check
	 * @return True if it is already contained, false otherwise
	 */
	public boolean containsSource(String sourceId){
		sourceId = sourceId.toLowerCase();
		Iterator<String> iter = sources.iterator();
		while (iter.hasNext()){
			String s = iter.next();
			if (sourceId.equals(s)){
				return true;
			}
		}
		return false;
	}
	
	
	/**
	 * Basic drawing
	 */
	public String toString(){
		String help = "Name: " + this.name + " - " 
		+ "Source: " + this.dbSource + " -  "
		+ "SourceCvId: " + this.dbSourceCvId + " - "
		+ "AccessionId: " + this.dbAccessionId + " - "
		+ "Version: " + this.dbVersion + " - "
		+ "\n"
		+ "Participants: \n ";
		if (participants != null){
			for (int i = 0; i < participants.size(); i++){
				help += "\t" + participants.get(i).toString();
			}
		}
		help += "Details: \n";
		if (details != null){
			for (int i = 0; i < details.size(); i++){
				help += "\t" + details.get(i).toString() + "\n";
			}
		}
		help+="\n";
		return help;
	}
	
}
