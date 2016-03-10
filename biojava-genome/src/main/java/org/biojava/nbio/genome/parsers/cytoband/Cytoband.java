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
 * created at 20 Feb 2014
 * Author: ap3
 */

package org.biojava.nbio.genome.parsers.cytoband;

import java.io.Serializable;

public class Cytoband implements Serializable, Comparable<Cytoband> {


	 /**
	 *
	 */
	private static final long serialVersionUID = 2805976387404499650L;
	String chromosome;
	Integer start;
	Integer end;
	String locus;
	StainType type;



	public String getChromosome() {
		return chromosome;
	}



	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}



	public Integer getStart() {
		return start;
	}



	public void setStart(Integer start) {
		this.start = start;
	}



	public Integer getEnd() {
		return end;
	}



	public void setEnd(Integer end) {
		this.end = end;
	}



	public StainType getType() {
		return type;
	}



	public void setType(StainType type) {
		this.type = type;
	}



	@Override
	public int compareTo(Cytoband o) {

		if ( this.chromosome.equals( o.chromosome)) {
			return this.start.compareTo(o.start);
		} else {



			Short s1 = null;
			try {
				s1 = Short.parseShort(chromosome.substring(3));
			} catch (NumberFormatException ex){}
			Short s2 = null;
			try {
				s2 = Short.parseShort(o.chromosome.substring(3));
			}catch (NumberFormatException ex){}

			if (s1 == null || s2 == null){
				return this.chromosome.compareTo(o.chromosome);
			} else {
				return s1.compareTo(s2);
			}
		}

	}



	public String getLocus() {
		return locus;
	}



	public void setLocus(String locus) {
		this.locus = locus;
	}



	@Override
	public String toString() {
		return "Cytoband [chromosome=" + chromosome + ", start=" + start
				+ ", end=" + end + ", locus=" + locus + ", type=" + type + "]";
	}







}
