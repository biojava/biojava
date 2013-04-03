package org.biojava3.ws.hmmer;

import java.io.Serializable;

/** Provides the details of a domain hit 
 * 
 * @author Andreas Prlic
 * @since 3.0.3
 */
public class HmmerDomain implements Comparable<HmmerDomain>, Serializable{

	
	/**
	 * 
	 */
	private static final long serialVersionUID = 8004302800150892757L;
	
	Integer sqFrom;
	Integer sqTo;
	Integer aliLenth;
	Integer simCount;
	Integer hmmFrom;
	Integer hmmTo;
	String hmmName;
	String hmmDesc;
	String hmmAcc;
	public Integer getSqFrom() {
		return sqFrom;
	}
	public void setSqFrom(Integer sqFrom) {
		this.sqFrom = sqFrom;
	}
	public Integer getSqTo() {
		return sqTo;
	}
	public void setSqTo(Integer sqTo) {
		this.sqTo = sqTo;
	}
	public Integer getAliLenth() {
		return aliLenth;
	}
	public void setAliLenth(Integer aliLenth) {
		this.aliLenth = aliLenth;
	}
	public Integer getSimCount() {
		return simCount;
	}
	public void setSimCount(Integer simCount) {
		this.simCount = simCount;
	}
	public Integer getHmmFrom() {
		return hmmFrom;
	}
	public void setHmmFrom(Integer hmmFrom) {
		this.hmmFrom = hmmFrom;
	}
	public Integer getHmmTo() {
		return hmmTo;
	}
	public void setHmmTo(Integer hmmTo) {
		this.hmmTo = hmmTo;
	}
	public String getHmmName() {
		return hmmName;
	}
	public void setHmmName(String hmmName) {
		this.hmmName = hmmName;
	}
	public String getHmmDesc() {
		return hmmDesc;
	}
	public void setHmmDesc(String hmmDesc) {
		this.hmmDesc = hmmDesc;
	}
	public String getHmmAcc() {
		return hmmAcc;
	}
	public void setHmmAcc(String hmmAcc) {
		this.hmmAcc = hmmAcc;
	}
	@Override
	public String toString() {
		return "HmmerDomain [hmmAcc=" + hmmAcc + ", hmmDesc=" + hmmDesc 
				+ "sqFrom=" + sqFrom + ", sqTo=" + sqTo
				+ ", aliLenth=" + aliLenth + ", simCount=" + simCount
				+ ", hmmFrom=" + hmmFrom + ", hmmTo=" + hmmTo + ", hmmName="
				+ hmmName  + "]" ;
				
	}
	@Override
	public int compareTo(HmmerDomain o) {
		if (emptyDomain(this) &&emptyDomain(o))
			return 0;
		
		if ( ! emptyDomain(this) &&emptyDomain(o))
			return -1;
		if (emptyDomain(this) && (! emptyDomain(o)))
			return 1;
		
		return (this.getSqFrom().compareTo(o.getSqFrom()));
	}
	private boolean emptyDomain(HmmerDomain o) {
		
		if  ( o.getSqFrom() == null)
			return true;
		return false;
	}
	
	/*
	 * 
	 * that's the full data that is currently responded from the JSON api:
	[{"ievalue":"0.021","cevalue":"1.0e-05",
		"alimline":"+++Wc+ ++++ +GW+ ++ +","jali":62,
		"alicsline":0,"aliIdCount":4,"aliSimCount":17,
		"aliM":55,"alisqto":"62","aliL":164,"alimemsize":
			"kenWcrvradGatGWiyqslL\u0000+++Wc+ ++++ +GW+ ++ +\u0000NGEWCEAQTKNGQGWVPSNYI\u000089**************98765\u0000000010555\u0000PF06347.7\u0000Bacterial SH3 domain\u0000seq\u0000\u0000",
			"alihmmacc":"PF06347.7","oasc":"0.88",
			"aliaseq":"NGEWCEAQTKNGQGWVPSNYI","aliN":21,
			"iali":42,"alihindex":"10555","is_reported":"1",
			"alimodel":"kenWcrvradGatGWiyqslL","alippline":"89**************98765",
			"aliSim":0.80952380952381,"alisqacc":"","jenv":"64",
			"alihmmname":"SH3_4","alihmmdesc":"Bacterial SH3 domain",
			"alihmmto":"53","alirfline":0,
			"bitscore":14.0345134735107,"bias":"0.20",
			"alisqfrom":"42",
			"alisqname":"1","aliappline":0,
			"alihmmfrom":"33","alimem":46912854337136,
			"aliId":0.19047619047619,"is_included":"0",
			"alisqdesc":"","ienv":"36"}] */
					
					
					
}
