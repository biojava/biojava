package org.biojava.bio.structure.scop;

import java.io.Serializable;
import java.io.StringWriter;
import java.util.List;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;


/** Container for the information for a domain. Contains a line in the file
 * dir.cla.scop.txt_1.75
 * 
 * e.g d1dlwa_	1dlw	A:	a.1.1.1	14982	cl=46456,cf=46457,sf=46458,fa=46459,dm=46460,sp=46461,px=14982
 * @author Andreas Prlic
 *
 */
@XmlRootElement(name = "ScopDomain", namespace ="http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class ScopDomain implements Serializable, Cloneable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 5890476209571654301L;
	
	String scopId;
	String pdbId;
	List<String> ranges;
	String classificationId;
	Integer sunid;

	int classId;
	int foldId;
	int superfamilyId;
	int familyId;
	int domainId;
	int speciesId;
	int px;




	@Override
	public String toString() {
		StringWriter buf = new StringWriter();
		buf.append(scopId);
		buf.append("\t") ;
		buf.append(pdbId);
		buf.append( "\t");

		int rangePos = 0;
		for (String range: ranges){
			rangePos++;

			buf.append(range);

			if ( ( ranges.size()> 1 ) && (rangePos < ranges.size()))
				buf.append(",");
		}
		buf.append("\t") ;
		buf.append(classificationId);
		buf.append("\t") ;
        buf.append(String.valueOf(sunid));
		buf.append("\t") ;
		
		buf.append("cl=");
        buf.append(String.valueOf(classId));
		buf.append(",cf=");
        buf.append(String.valueOf(foldId));
		buf.append(",sf=");
		buf.append(String.valueOf(superfamilyId));
		buf.append(",fa=");
		buf.append(String.valueOf(familyId));        
		buf.append(",dm=");
        buf.append(String.valueOf(domainId));
		buf.append(",sp=");
        buf.append(String.valueOf(speciesId));
		buf.append(",px=");
        buf.append(String.valueOf(px));
		
		
		return buf.toString();
	}

	public String getScopId() {
		return scopId;
	}
	public void setScopId(String scopId) {
		this.scopId = scopId;
	}
	public String getPdbId() {
		return pdbId;
	}
	public void setPdbId(String pdbId) {
		this.pdbId = pdbId;
	}
	public List<String> getRanges() {
		return ranges;
	}
	public void setRanges(List<String> ranges) {
		this.ranges = ranges;
	}
	public String getClassificationId() {
		return classificationId;
	}
	public void setClassificationId(String classificationId) {
		this.classificationId = classificationId;
	}
	public Integer getSunid() {
		return sunid;
	}
	public void setSunid(Integer sunid) {
		this.sunid = sunid;
	}
	public int getClassId() {
		return classId;
	}
	public void setClassId(int classId) {
		this.classId = classId;
	}
	public int getFoldId() {
		return foldId;
	}
	public void setFoldId(int foldId) {
		this.foldId = foldId;
	}
	public int getSuperfamilyId() {
		return superfamilyId;
	}
	public void setSuperfamilyId(int superfamilyId) {
		this.superfamilyId = superfamilyId;
	}
	public int getFamilyId() {
		return familyId;
	}
	public void setFamilyId(int familyId) {
		this.familyId = familyId;
	}
	public int getDomainId() {
		return domainId;
	}
	public void setDomainId(int domainId) {
		this.domainId = domainId;
	}
	public int getSpeciesId() {
		return speciesId;
	}
	public void setSpeciesId(int speciesId) {
		this.speciesId = speciesId;
	}
	public int getPx() {
		return px;
	}
	public void setPx(int px) {
		this.px = px;
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
	
		
		ScopDomain n = new ScopDomain();
		n.setClassId(getClassId());
		n.setClassificationId(getClassificationId());
		n.setDomainId(getDomainId());
		n.setFamilyId(getFamilyId());
		n.setFoldId(getFoldId());
		n.setPdbId(getPdbId());
		n.setPx(getPx());
		n.setRanges(getRanges());
		n.setScopId(getScopId());
		n.setSpeciesId(getSpeciesId());
		n.setSunid(getSunid());
		n.setSuperfamilyId(getSuperfamilyId());
		
		
		return n;
		
		
	}

	


}
