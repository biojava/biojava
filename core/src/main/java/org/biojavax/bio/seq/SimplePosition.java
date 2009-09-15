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

package org.biojavax.bio.seq;


/**
 * A simple implementation of the Position interface.
 * @author Richard Holland
 * @since 1.5
 */
public class SimplePosition implements Position {
    
    private boolean fs;
    private boolean fe;
    private int s;
    private int e;
    private String t;
    
    /**
     * Constructs a point position, with no fuzzy start or
     * end. (eg. 1, 2, or 3).
     * @param p the point position
     */
    public SimplePosition(int p) {
        this(false,false,p);
    }
    
    /**
     * Constructs a range position, with no fuzzy start or
     * end. (eg. 1..2, 2..5, or 3..8).
     * @param s the start position
     * @param e the end position
     */
    public SimplePosition(int s, int e) {
        this(false,false,s,e,null);
    }
    
    /**
     * Constructs a point position, with optionally fuzzy start and
     * end. (eg. <1 or 3> or 2 or even <5>).
     * @param fs fuzzy start?
     * @param fe fuzzy end?
     * @param p the point position
     */
    public SimplePosition(boolean fs, boolean fe, int p) {
        this(fs,fe,p,p,null);
    }
    
    /**
     * Constructs a range position, with optionally fuzzy start and
     * end. (eg. <1.2 or 1^3> or 2.2 or even <5^6>). The type of the
     * range is given, it should normally be one of the two defined
     * in the Position interface, but its up to you.
     * @param fs fuzzy start?
     * @param fe fuzzy end?
     * @param s the start of the position
     * @param e the end of the position
     * @param t the type of the position
     */
    public SimplePosition(boolean fs, boolean fe, int s, int e, String t) {
        this.fs = fs;
        this.fe = fe;
        this.s = s;
        this.e = e;
        this.t = t;
    }
    
    //Hibernate only - futureproofing
    protected SimplePosition() {}
    
    /**
     * {@inheritDoc}
     */
    public boolean getFuzzyStart() { return this.fs; }
    
    /**
     * {@inheritDoc}
     */
    public boolean getFuzzyEnd() { return this.fe; }
    
    /**
     * {@inheritDoc}
     */
    public int getStart() { return this.s; }
    
    /**
     * {@inheritDoc}
     */
    public int getEnd()  { return this.e; }
    
    /**
     * {@inheritDoc}
     */
    public String getType() { return this.t; }
    
    // Hibernate requirement - not for public use - futureproofing
    void setFuzzyStart(boolean fs) { this.fs = fs; }
    
    // Hibernate requirement - not for public use - futureproofing
    void setFuzzyEnd(boolean fe) { this.fe = fe; }
    
    // Hibernate requirement - not for public use - futureproofing
    void setStart(int s) { this.s = s; }
    
    // Hibernate requirement - not for public use - futureproofing
    void setEnd(int e) { this.e = e; }
    
    // Hibernate requirement - not for public use - futureproofing
    void setType(String t) { this.t = t; }
    
    /** 
     * {@inheritDoc}
     */
    public Position translate(int distance) {
        return new SimplePosition(this.fs,this.fe,this.s+distance,this.e+distance,this.t);
    }
        
    /** 
     * {@inheritDoc}
     * Two positions are equal if they share all parameters in common, 
     * eg. fuzzy start+end, start, end, type.
     */
    public boolean equals(Object o) {
        if (!(o instanceof Position)) return false;
        if (o==this) return true;
        Position them = (Position)o;
        if (this.getFuzzyStart() != them.getFuzzyStart()) return false;
        if (this.getFuzzyEnd() != them.getFuzzyEnd()) return false;
        if (this.getStart()!=them.getStart()) return false;
        if (this.getEnd()!=them.getEnd()) return false;
        if (this.getType()!=null || them.getType()!=null) {
            if (this.getType()!=null && them.getType()!=null) {
                if (!this.getType().equals(them.getType())) return false;
            } else return false;
        }
        return true;
    }  
    
    /** 
     * {@inheritDoc}
     */
    public String toString() {
        StringBuffer sb = new StringBuffer();
        if (this.getFuzzyStart()) sb.append("<");
        sb.append(this.s);
        if (s!=e) {
            sb.append(this.t);
            sb.append(this.e);
        }
        if (this.getFuzzyEnd()) sb.append(">");
        return sb.toString();
    }
    
    // Hibernate requirement - not for public use - futureproofing
    private Integer id;
    
    /**
     * Gets the Hibernate ID. Should be used with caution.
     * @return the Hibernate ID, if using Hibernate.
     */
    public Integer getId() { return this.id; }
    
    /**
     * Sets the Hibernate ID. Should be used with caution.
     * @param id the Hibernate ID, if using Hibernate.
     */
    public void setId(Integer id) { this.id = id;}
    
}
