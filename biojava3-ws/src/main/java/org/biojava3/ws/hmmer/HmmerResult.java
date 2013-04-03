package org.biojava3.ws.hmmer;

import java.io.Serializable;
import java.util.SortedSet;

/** The results of a Hmmer search for a single sequence
 * 
 * @author Andreas Prlic
 * @since 3.0.3
 */
public class HmmerResult implements Comparable<HmmerResult>, Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -6016026193090737943L;
	
	String desc ;
	Float score;
	Float evalue;
	Double pvalue;
	String acc;
	Integer dcl;
	String name;
	Integer ndom;
	Integer nreported;
	
	SortedSet<HmmerDomain>domains;
	
	public SortedSet<HmmerDomain> getDomains() {
		return domains;
	}
	public void setDomains(SortedSet<HmmerDomain> domains) {
		this.domains = domains;
	}
	public String getDesc() {
		return desc;
	}
	public void setDesc(String desc) {
		this.desc = desc;
	}
	public Float getScore() {
		return score;
	}
	public void setScore(Float score) {
		this.score = score;
	}
	public Float getEvalue() {
		return evalue;
	}
	public void setEvalue(Float evalue) {
		this.evalue = evalue;
	}
	public Double getPvalue() {
		return pvalue;
	}
	public void setPvalue(Double pvalue) {
		this.pvalue = pvalue;
	}
	public String getAcc() {
		return acc;
	}
	public void setAcc(String acc) {
		this.acc = acc;
	}
	public Integer getDcl() {
		return dcl;
	}
	public void setDcl(Integer dcl) {
		this.dcl = dcl;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public Integer getNdom() {
		return ndom;
	}
	public void setNdom(Integer ndom) {
		this.ndom = ndom;
	}
	public Integer getNreported() {
		return nreported;
	}
	public void setNreported(Integer nreported) {
		this.nreported = nreported;
	}
	@Override
	public String toString() {
		return "HmmerResult [acc=" + acc + ", desc=" + desc + ", score=" + score + ", evalue="
				+ evalue + ", pvalue=" + pvalue + ", dcl="
				+ dcl + ", name=" + name + ", ndom=" + ndom + ", nreported="
				+ nreported + ", domains=" + domains + "]";
	}
	
	
	@Override
	public int compareTo(HmmerResult o) {
		// 	sort  by the start position of the first domain
		
		if ( emptyDomains(this) && emptyDomains(o)){
			return 0;
		}
		
		if ( ! emptyDomains(this) && emptyDomains(o))
			return -1;
		
		if ( emptyDomains(this) && (! emptyDomains(o)))
			return 1;
		
		// ok when we are here, both domains are not empty
		
		HmmerDomain me = this.getDomains().first();
		HmmerDomain other = o.getDomains().first();
		
		//System.out.println(" domains: " + me.getHmmAcc() + " " + other.getHmmAcc()+ " " + me.getSqFrom().compareTo(other.getSqFrom()));
		
		return(me.getSqFrom().compareTo(other.getSqFrom()));
	}
	private boolean emptyDomains(HmmerResult o) {
		if ( o.getDomains() == null || o.getDomains().size() == 0)
			return true;
		return false;
	}
	
	
	/** Get the overlap between two HmmerResult objects
	 * 
	 * @param other
	 * @return 0 if no overlap, otherwise the length of the overlap
	 */
	public int getOverlapLength(HmmerResult other){
		
		int overlap = 0;
		for ( HmmerDomain d1 : getDomains()){
			for (HmmerDomain d2 : other.getDomains()){
				overlap += getOverlap(d1, d2);
			}
		}
		return overlap;
		
	}
	
	private int getOverlap(HmmerDomain one, HmmerDomain other){
		int xs = one.getSqFrom();
		int ys = one.getSqTo();
		int as = other.getSqFrom();
		int bs = other.getSqTo();

		int overlap = 0;
		//1:
		
		if ((( xs< as)  && ( as<ys)) || ((xs < bs) && ( bs <= ys)) || (as<xs && ys<bs)) {
			
			//2:

			if ( xs < as) {
				if ( ys < bs) 
					overlap = ys-as;
				else
					overlap = bs-as;
			} else {
				if  ( ys < bs) 
					overlap = ys -xs;
				else 
					overlap = bs - xs;

			} 

		}

		return overlap;
	}
	
}
