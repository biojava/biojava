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

public class SiftsEntity implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 750353252427491487L;
	String type;
	String entityId;
	
	List<SiftsSegment> segments;
	
	public SiftsEntity(){
		this(null,null);
	}
	
	public SiftsEntity(String type, String entityId) {
		this.type = type;
		this.entityId = entityId;
		segments = new ArrayList<SiftsSegment>();
	}

	public void addSegment(SiftsSegment s) {
		segments.add(s);
		
	}
	
	public List<SiftsSegment> getSegments(){
		return segments;
	}
	
	public void setSegments(List<SiftsSegment> segments){
		this.segments = segments;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public String getEntityId() {
		return entityId;
	}

	public void setEntityId(String entityId) {
		this.entityId = entityId;
	}

	@Override
	public String toString() {
		return "SiftsEntity [type=" + type + ", entityId=" + entityId
				+ ", segments=" + segments + "]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((entityId == null) ? 0 : entityId.hashCode());
		result = prime * result
				+ ((segments == null) ? 0 : segments.hashCode());
		result = prime * result + ((type == null) ? 0 : type.hashCode());
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
		SiftsEntity other = (SiftsEntity) obj;
		if (entityId == null) {
			if (other.entityId != null)
				return false;
		} else if (!entityId.equals(other.entityId))
			return false;
		if (segments == null) {
			if (other.segments != null)
				return false;
		} else if (!segments.equals(other.segments))
			return false;
		if (type == null) {
			if (other.type != null)
				return false;
		} else if (!type.equals(other.type))
			return false;
		return true;
	}
	
	

}
