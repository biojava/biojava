/*
 * BioJava development code
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
 * Author: Daniel Asarnow
 * Date:   2012-7-23
 */

package org.biojava.bio.structure.cath;

/**
 * @author Daniel Asarnow
 */

import java.io.Serializable;
import java.util.Date;
import java.util.List;

/**
 * A class which represents a single CATH domain.
 */
public class CathDomain implements Serializable {

    public static final long serialVersionUID = 1L;

    /**
     * The CATH domain code. Always 7 characters in length, combining the PDB and chain letter with the number of the domain within CATH.
     * Example: 1aoiA00
     * If the chain letter '0', domain refers to an entire PDB entry.
     */
    String domainName; // 7 characters 1oaiA00

    /**
     * The class number of this domain.
     */
    Integer classId; // C

    /**
     * The architecture number of this domain.
     */
    Integer architectureId; // A

    /**
     * The topology number of this domain.
     */
    Integer topologyId; // T

    /**
     * The homologous superfamily number of this domain.
     */
    Integer homologyId; // H

    /**
     * The sequence family (35% identity) number of this domain.
     */
    Integer sequenceFamilyId; // S

    /**
     * The "orthologous" sequence family (60% identity) number of this domain.
     */
    Integer orthologousSequenceFamilyId; // O

    /**
     * The "Like" sequence family (95% identity) number of this domain.
     */
    Integer likeSequenceFamilyId; // L

    /**
     * The identical sequence family (100% identity) number of this domain.
     */

    Integer identicalSequenceFamilyId; // I

    /**
     * The count of this domain among the identical sequence family members.
     */
    Integer domainCounter; // D

    /**
     * The domain length..
     */
    Integer length;

    /**
     * The resolution of the domain structure. Nominally in Angstroms,
     * the values 999.000 and 1000.000 signify NMR structures and obsolete structures, respectively.
     */
    Double resolution;

    /**
     * The format and version of the CathDomainDescriptionFile.
     */
    String format;

    /**
     * The CATH version.
     */
    String version;

    Date date;

    /**
     * The so-called name field holds a potentially long description of the domain.
     */
    String name;

    /**
     * Complete source organism listing.
     */
    String source;

    /**
     * FASTA header.
     */
    String sequenceHeader;

    /**
     * FASTA sequence.
     */
    String sequence;

    /**
     * List of all sub-domain segments.
     */
    List<CathSegment> segments;

    /**
     * A (potentially long) comment. Usually empty.
     */
    String comment;

    public String getDomainName() {
        return domainName;
    }

    public void setDomainName(String domainName) {
        this.domainName = domainName;
    }

    public String getPdbId() {
        return domainName.substring(0, 4) +
                (!domainName.substring(4, 5).equals("0") ? "." + domainName.substring(4, 5) : ""); //TODO ask about
    }

    public Integer getDomainId() {
        return Integer.parseInt(domainName.substring(5));
    }

    public Integer getClassId() {
        return classId;
    }

    public void setClassId(Integer classId) {
        this.classId = classId;
    }

    public Integer getArchitectureId() {
        return architectureId;
    }

    public void setArchitectureId(Integer architectureId) {
        this.architectureId = architectureId;
    }

    public Integer getTopologyId() {
        return topologyId;
    }

    public void setTopologyId(Integer topologyId) {
        this.topologyId = topologyId;
    }

    public Integer getHomologyId() {
        return homologyId;
    }

    public void setHomologyId(Integer homologyId) {
        this.homologyId = homologyId;
    }

    public Integer getSequenceFamilyId() {
        return sequenceFamilyId;
    }

    public void setSequenceFamilyId(Integer sequenceFamilyId) {
        this.sequenceFamilyId = sequenceFamilyId;
    }

    public Integer getOrthologousSequenceFamilyId() {
        return orthologousSequenceFamilyId;
    }

    public void setOrthologousSequenceFamilyId(Integer orthologousSequenceFamilyId) {
        this.orthologousSequenceFamilyId = orthologousSequenceFamilyId;
    }

    public Integer getLikeSequenceFamilyId() {
        return likeSequenceFamilyId;
    }

    public void setLikeSequenceFamilyId(Integer likeSequenceFamilyId) {
        this.likeSequenceFamilyId = likeSequenceFamilyId;
    }

    public Integer getIdenticalSequenceFamilyId() {
        return identicalSequenceFamilyId;
    }

    public void setIdenticalSequenceFamilyId(Integer identicalSequenceFamilyId) {
        this.identicalSequenceFamilyId = identicalSequenceFamilyId;
    }

    public Integer getDomainCounter() {
        return domainCounter;
    }

    public void setDomainCounter(Integer domainCounter) {
        this.domainCounter = domainCounter;
    }

    public Integer getLength() {
        return length;
    }

    public void setLength(Integer length) {
        this.length = length;
    }

    public Double getResolution() {
        return resolution;
    }

    public void setResolution(Double resolution) {
        this.resolution = resolution;
    }

    public void setCATH(String cathCode) {
        String[] token = cathCode.split("[.]");
        setClassId(Integer.parseInt(token[0]));
        setArchitectureId(Integer.parseInt(token[1]));
        setTopologyId(Integer.parseInt(token[2]));
        setHomologyId(Integer.parseInt(token[3]));
    }

    public String getCATH() {
        return Integer.toString(getClassId()) + "." +
                Integer.toString(getArchitectureId()) + "." +
                Integer.toString(getTopologyId()) + "." +
                Integer.toString(getHomologyId());
    }

    public void setSOLID(String cathCode) {
        String[] token = cathCode.split("[.]");
        setSequenceFamilyId(Integer.parseInt(token[0]));
        setOrthologousSequenceFamilyId(Integer.parseInt(token[1]));
        setLikeSequenceFamilyId(Integer.parseInt(token[2]));
        setIdenticalSequenceFamilyId(Integer.parseInt(token[3]));
        setDomainCounter(Integer.parseInt(token[4]));
    }

    public String getSOILD() {
        return Integer.toString(getSequenceFamilyId()) + "." +
                Integer.toString(getOrthologousSequenceFamilyId()) + "." +
                Integer.toString(getLikeSequenceFamilyId()) + "." +
                Integer.toString(getIdenticalSequenceFamilyId()) + "." +
                Integer.toString(getDomainCounter());
    }

    public Integer getClassificationId(CathCategory cathCategory) {
        switch (cathCategory) {
            case Class:
                return getClassId();
            case Architecture:
                return getArchitectureId();
            case Topolgy:
                return getTopologyId();
            case Homology:
                return getHomologyId();
            case SequenceFamily:
                return getSequenceFamilyId();
            case OrthologousSequenceFamily:
                return getOrthologousSequenceFamilyId();
            case LikeSequenceFamily:
                return getLikeSequenceFamilyId();
            case IdenticalSequenceFamily:
                return getIdenticalSequenceFamilyId();
            case DomainCounter:
                return getDomainCounter();
            default:
                return null;
        }
    }

    public String getFormat() {
        return format;
    }

    public void setFormat(String format) {
        this.format = format;
    }

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public Date getDate() {
        return date;
    }

    public void setDate(Date date) {
        this.date = date;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getSource() {
        return source;
    }

    public void setSource(String source) {
        this.source = source;
    }

    public String getSequenceHeader() {
        return sequenceHeader;
    }

    public void setSequenceHeader(String sequenceHeader) {
        this.sequenceHeader = sequenceHeader;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public List<CathSegment> getSegments() {
        return segments;
    }

    public void setSegments(List<CathSegment> segments) {
        this.segments = segments;
    }

    public String getComment() {
        return comment;
    }

    public void setComment(String comment) {
        this.comment = comment;
    }

	@Override
	public String toString() {
		return "CathDomain [domainName=" + domainName + ", classId=" + classId
				+ ", architectureId=" + architectureId + ", topologyId="
				+ topologyId + ", homologyId=" + homologyId
				+ ", sequenceFamilyId=" + sequenceFamilyId
				+ ", orthologousSequenceFamilyId="
				+ orthologousSequenceFamilyId + ", likeSequenceFamilyId="
				+ likeSequenceFamilyId + ", identicalSequenceFamilyId="
				+ identicalSequenceFamilyId + ", domainCounter="
				+ domainCounter + ", length=" + length + ", resolution="
				+ resolution + ", format=" + format + ", version=" + version
				+ ", date=" + date + ", name=" + name + ", source=" + source
				+ ", sequenceHeader=" + sequenceHeader + ", sequence="
				+ sequence + ", segments=" + segments + ", comment=" + comment
				+ "]";
	}

//    @Override
//    public String toString() {
//    	StringBuffer buf = new StringBuffer();
//    	buf.append("CathDomain " + domainName + " ");
//    }

}
