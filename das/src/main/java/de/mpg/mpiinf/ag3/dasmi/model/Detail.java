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
 * A detail element, contains additional information on an interactor (e.g. domains), 
 * an interaction (e.g. confidence score) or a participant (e.g. experimental role)
 * 
 * @author Hagen Blankenburg, Max Planck Institute for Informatics
 *
 */
public class Detail {
	
	private List<String> sources = null; // contains the names of all sources that report this interaction
	private String property = null;
	private String value = null;
	private String propertyCvId = null;
	private String valueCvId = null;
	private Range range = null;
	
	
	/**
	 * Basic initialization
	 *
	 */
	public Detail(){
		this.sources = new ArrayList<String>();
		this.property = Constants.INVALID_STRING;
		this.propertyCvId = Constants.INVALID_STRING;
		this.value = Constants.INVALID_STRING;
		this.valueCvId = Constants.INVALID_STRING;
	}
	
	/**
	 * A new detail with the property value combination
	 * @param property The property of the detail
	 * @param value The value of the detail
	 */
	public Detail(String property, String value){
		this.sources = new ArrayList<String>();
		this.property = property;
		this.propertyCvId = Constants.INVALID_STRING;
		this.value = value;
		this.valueCvId = Constants.INVALID_STRING;
	}

	/**
	 * @return the property
	 */
	public String getProperty() {
		return property;
	}

	/**
	 * @param property the property to set
	 */
	public void setProperty(String property) {
		this.property = property;
	}

	/**
	 * @return the propertyCvId
	 */
	public String getPropertyCvId() {
		return propertyCvId;
	}

	/**
	 * @param propertyCvId the propertyCvId to set
	 */
	public void setPropertyCvId(String propertyCvId) {
		this.propertyCvId = propertyCvId;
	}

	/**
	 * @return the range
	 */
	public Range getRange() {
		return range;
	}

	/**
	 * @param range the range to set
	 */
	public void setRange(Range range) {
		this.range = range;
	}

	/**
	 * @return the value
	 */
	public String getValue() {
		return value;
	}

	/**
	 * @param value the value to set
	 */
	public void setValue(String value) {
		this.value = value;
	}

	/**
	 * @return the valueCvId
	 */
	public String getValueCvId() {
		return valueCvId;
	}

	/**
	 * @param valueCvId the valueCvId to set
	 */
	public void setValueCvId(String valueCvId) {
		this.valueCvId = valueCvId;
	}
		
	/**
	 * 
	 * @return
	 */
	public List<String> getSources(){
		return this.sources;
	}
	
	/**
	 * 
	 * @param sources
	 */
	public void setSources(List<String> sources){
		this.sources = sources;
	}
		
	/**
	 * Adds the id of a source to the sources list
	 * @param sourceId The id of the source to be added
	 */
	public void addSource(String sourceId){
		sourceId = sourceId.toLowerCase();
		if (!this.sources.contains(sourceId)){
			this.sources.add(sourceId);
		}
	}
		
	/**
	 * Removes the source with the attributed id from the
	 * source list. 
	 * @param sourceId id of the source that should be removed
	 */
	public void removeSource(String sourceId){
		sourceId = sourceId.toLowerCase();
		this.sources.remove(sourceId);	
	}
		
	/**
	 * Determinds if two details are the same, i.e., have the same 
	 * property and value. 
	 * @param comp The Detail to compare with
	 * @return True if they are equal, false otherwise
	 */
	public boolean isEqual(Detail comp){
		if (comp.getProperty().equalsIgnoreCase(this.property)
			&& comp.getValue().equals(this.value)){
			return true;
		}
		return false;
	}
	
	/**
	 * Checks if a certain source is already contained in the sources list
	 * @param source The source to check 
	 * @return True if the source is already contained, false otherwise
	 */
	public boolean containsSource(String source){
		Iterator<String> iter = sources.iterator();
		while (iter.hasNext()){
			String s = iter.next();
			if (source.equals(s)){
				return true;
			}
		}
		return false;
	}
		
	/**
	 * Simple printout function
	 */
	public String toString(){
		String help = "Property: " + this.property + 
					 " - PropCvId: " + this.propertyCvId +
					 " - Value: " + this.value +
					 " - ValCvId: " + this.valueCvId;
		if (range != null){
			help += "\n\t\t\t\t" + range.toString() + "\n";
		}
		return help;
		
	}

}
