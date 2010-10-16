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

import java.util.ArrayList;
import java.util.List;

/**
 * An interactor taking part in a certain interaction is called a participant. This distinction is 
 * necessary to give the interactor interaction specific properties, like experimental role
 * 
 * @author Hagen Blankenburg, Max Planck Institute for Informatics
 *
 */
public class Participant {
	
	private String id = null;
	private Interactor interactor = null;
	private List <Detail>details = new ArrayList<Detail>();
	
	/**
	 * Empty construcotr
	 *
	 */
	public Participant(){}
	
	/**
	 * Basic initialization
	 * @param interactor 
	 */
	public Participant(Interactor interactor) {
		this.interactor = interactor;
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
	 * @param interactor 
	 */
	public void setInteractor(Interactor interactor) {
		this.interactor = interactor;
	}
	
	/**
	 * @return the Interactor
	 * 
	 */
	public Interactor getInteractor(){
		return this.interactor;
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
	 * Basic to string
	 * @return the string
	 */
	public String toString(){
		String help = null;
		if (interactor != null){
			help = "Interactor:\n\t\t" +  this.interactor.toString();
			if (details != null){
				help += "\tDetails:\n";
				for (int i = 0; i < details.size(); i++){
					help += "\t\t" + details.get(i).toString() +"\n";
				}
			}
		}
		return help;
	}
	
	
	

}
