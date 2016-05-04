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
 */
package org.biojava.nbio.structure.scop;

import java.io.IOException;
import java.io.Serializable;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

import org.biojava.nbio.structure.ResidueRange;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.SubstructureIdentifier;
import org.biojava.nbio.structure.align.util.AtomCache;


/** Container for the information for a domain. Contains a line in the file
 * dir.cla.scop.txt_1.75
 *
 * e.g d1dlwa_	1dlw	A:	a.1.1.1	14982	cl=46456,cf=46457,sf=46458,fa=46459,dm=46460,sp=46461,px=14982
 *
 * Instantiated using {@link ScopDatabase#getDomainByScopID(String)}
 * @author Andreas Prlic
 *
 */
@XmlRootElement(name = "ScopDomain", namespace ="http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class ScopDomain implements Serializable, Cloneable, StructureIdentifier {

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
		StringBuilder buf = new StringBuilder();
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

		super.clone();

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

	/**
	 * Returns the chains this domain is defined over; contains more than 1 element only if this domains is a multi-chain domain.
	 */
	public Set<String> getChains() {
		Set<String> chains = new HashSet<String>();
		List<ResidueRange> rrs = ResidueRange.parseMultiple(getRanges());
		for (ResidueRange rr : rrs) chains.add(rr.getChainName());
		return chains;
	}

	@Override
	public String getIdentifier() {
		return getScopId();
	}

	public List<ResidueRange> getResidueRanges() {
		return ResidueRange.parseMultiple(ranges);
	}

	@Override
	public SubstructureIdentifier toCanonical() {
		return new SubstructureIdentifier(getPdbId(), ResidueRange.parseMultiple(getRanges()));
	}

	@Override
	public Structure reduce(Structure input) throws StructureException {
		return toCanonical().reduce(input);
	}

	@Override
	public Structure loadStructure(AtomCache cache) throws StructureException,
	IOException {
		return cache.getStructureForPdbId(pdbId);
	}

}
