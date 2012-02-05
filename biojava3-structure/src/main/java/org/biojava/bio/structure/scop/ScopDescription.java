package org.biojava.bio.structure.scop;

import java.io.Serializable;
import java.io.StringWriter;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

/** Contains data from
 * dir.des.scop.txt_1.75
 * 
 * <p>e.g 
 * <pre>
 * SunID	Cat	Class   	Name	Description
 * -----	---	-----   	----	-----------
 * 26154	px	b.47.1.2	d1nrs.1	1nrs L:,H:
 * 125030	px	b.47.1.2	d1zgia1	1zgi A:1A-245
 * </pre>
 * 
 * @author Andreas Prlic
 *
 */
@XmlRootElement(name = "ScopDescription", namespace ="http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class ScopDescription implements Serializable,Cloneable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 8579808155176839161L;
	int sunID;
	ScopCategory category;
	String classificationId;
	String name;
	String description;


	public String toString(){
		StringWriter buf = new StringWriter();

        buf.append(String.valueOf(sunID));
		buf.append("\t");
		buf.append(category.toString());
		buf.append("\t");
		buf.append(classificationId);
		buf.append("\t");
		buf.append(name);
		buf.append("\t");
		buf.append(description);

		return buf.toString();
	}


	public int getSunID()
	{
		return sunID;
	}
	public void setSunID(int sunID)
	{
		this.sunID = sunID;
	}
	public ScopCategory getCategory()
	{
		return category;
	}
	public void setCategory(ScopCategory category)
	{
		this.category = category;
	}
	public String getClassificationId()
	{
		return classificationId;
	}
	public void setClassificationId(String classificationId)
	{
		this.classificationId = classificationId;
	}
	public String getName()
	{
		return name;
	}
	public void setName(String name)
	{
		this.name = name;
	}
	public String getDescription()
	{
		return description;
	}
	public void setDescription(String description)
	{
		this.description = description;
	}

	// Methods to return parts of the classificationID

	/**
	 * Return a portion of the classificationID corresponding to the specified
	 * category (class, fold, superfamily, family).
	 * 
	 * <p>Example: for SCOP family "b.5.1.1",
	 * getClassificationId(ScopCategory.Superfamily) => "b.5.1"
	 */
	public String getClassificationId(ScopCategory category) {
		if(classificationId == null || classificationId.isEmpty()) {
			return null;
		}

		int numParts = 0;
		switch(category) {
		case Family:      numParts++;
		case Superfamily: numParts++;
		case Fold:        numParts++;
		case Class:       numParts++; break;
		default:
			throw new IllegalArgumentException("Only Class, Fold, Superfamily, and Family are supported.");
		}

		int endChar = -1;
		for(int i = 0;i<numParts-1;i++) {
			endChar = classificationId.indexOf('.', endChar+1);
			if(endChar<0) {
				// Not enough items in the classification for this category
				return null;
			}
		}
		endChar = classificationId.indexOf('.', endChar+1);
		if(endChar<0) {
			// category goes to the end
			return classificationId;
		}
		else {
			return classificationId.substring(0, endChar);
		}

	}

	/**
	 * @return
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
		+ ((category == null) ? 0 : category.hashCode());
		result = prime
		* result
		+ ((classificationId == null) ? 0 : classificationId.hashCode());
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		result = prime * result + sunID;
		return result;
	}


	/**
	 * Compares the fields sunID, category, classificationId, and name for equality
	 * 
	 * @param obj
	 * @return
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		ScopDescription other = (ScopDescription) obj;
		if (category == null) {
			if (other.category != null) {
				return false;
			}
		} else if (!category.equals(other.category)) {
			return false;
		}
		if (classificationId == null) {
			if (other.classificationId != null) {
				return false;
			}
		} else if (!classificationId.equals(other.classificationId)) {
			return false;
		}
		if (name == null) {
			if (other.name != null) {
				return false;
			}
		} else if (!name.equals(other.name)) {
			return false;
		}
		if (sunID != other.sunID) {
			return false;
		}
		return true;
	}


	@Override
	protected Object clone() throws CloneNotSupportedException {
		ScopDescription n = new ScopDescription();
		
		n.setCategory(getCategory());
		n.setClassificationId(getClassificationId());
		n.setDescription(getDescription());
		n.setName(getName());
		n.setSunID(getSunID());
		return n;
	}

	
	
}
