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

import de.mpg.mpiinf.ag3.dasmi.Constants;

/**
 * Class representing a sequence. The type of the sequence is determined by the interactor type.
 * 
 * @author Hagen Blankenburg, Max Planck Institute for Informatics
 *
 */
public class Sequence {
	
	private int start;
	private int end;
	private String sequence; 
	
	/**
	 * Basic initialization
	 *
	 */
	public Sequence(){
		this.sequence = Constants.INVALID_STRING;
		this.start = Constants.INVALID_INT;
		this.end = Constants.INVALID_INT;
	}
	

	/**
	 * @return the end
	 */
	public int getEnd() {
		return end;
	}
	
	/**
	 * @param end the end to set
	 */
	public void setEnd(int end) {
		this.end = end;
	}
	
	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}
	
	/**
	 * @param start the start to set
	 */
	public void setStart(int start) {
		this.start = start;
	}
	
	/**
	 * 
	 * @return the sequence
	 */
	public String getSequence(){
		return this.sequence;
	}
	
	/**
	 * 
	 * @param sequence the sequence to set
	 */
	public void setSequence(String sequence){
		this.sequence = sequence;
	}
	
	
	/**
	 * Returns a specific part of the sequence 
	 * @param from the start of the sequence chunk
	 * @param till the end of the sequence chunk
	 * @return the Sequence
	 */
	public String getSequence(int from, int till){
		if (from > 0 && from < end){
			return sequence.substring(from); 
		}else{
			return null;
		} 
	}

}
