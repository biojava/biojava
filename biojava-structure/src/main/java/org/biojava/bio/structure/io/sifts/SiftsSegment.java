/**
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
 * Created on Feb 22, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.io.sifts;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class SiftsSegment implements Serializable{

	
	/**
	 * 
	 */
	private static final long serialVersionUID = -8005129863256307153L;
	String segId;
	String start;
	String end;
	
	List<SiftsResidue> residues;
	
	public SiftsSegment(){
		this(null,null,null);
	}
	
	public SiftsSegment(String segId, String start, String end) {
		this.segId = segId;
		this.start = start;
		this.end = end;
		residues = new ArrayList<SiftsResidue>();
	}

	public String getSegId() {
		return segId;
	}

	public void setSegId(String segId) {
		this.segId = segId;
	}

	public String getStart() {
		return start;
	}

	public void setStart(String start) {
		this.start = start;
	}

	public String getEnd() {
		return end;
	}

	public void setEnd(String end) {
		this.end = end;
	}

	public void addResidue(SiftsResidue pos) {
		residues.add(pos);
		
	}

	public List<SiftsResidue> getResidues(){
		return residues;
	}
	
	public void setResidues(List<SiftsResidue> residues){
		this.residues = residues;
	}

	@Override
	public String toString() {
		return "SiftsSegment [segId=" + segId + ", start=" + start + ", end="
				+ end + ", residues=" + residues + "]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((end == null) ? 0 : end.hashCode());
		result = prime * result
				+ ((residues == null) ? 0 : residues.hashCode());
		result = prime * result + ((segId == null) ? 0 : segId.hashCode());
		result = prime * result + ((start == null) ? 0 : start.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SiftsSegment other = (SiftsSegment) obj;
		if (end == null) {
			if (other.end != null)
				return false;
		} else if (!end.equals(other.end))
			return false;
		if (residues == null) {
			if (other.residues != null)
				return false;
		} else if (!residues.equals(other.residues))
			return false;
		if (segId == null) {
			if (other.segId != null)
				return false;
		} else if (!segId.equals(other.segId))
			return false;
		if (start == null) {
			if (other.start != null)
				return false;
		} else if (!start.equals(other.start))
			return false;
		return true;
	}
	
	
}
