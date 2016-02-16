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
 * Created on 28.04.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.FileConvert;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


/**
 * Implementation of an Atom of a PDB file.
 * currently the coordinates of an atom are represented by a double[3] array.
 * @author Andreas Prlic
 * @since 1.4
 * @version %I% %G%
 */
public class AtomImpl implements Atom, Serializable, PDBRecord {

	private static final long serialVersionUID = -2258364127420562883L;
	
	private String name     ;
	private Element element;
	private double[] coords ;
	private int pdbserial   ;
	private short charge		;

	private float occupancy ;
	private float tempfactor;

	private char altLoc ;
	private Group parent;

	private List<Bond> bonds;

	public AtomImpl () {
		name       = null        ;
		element    = Element.R;
		coords     = new double[3];
		occupancy  = 0.0f       ;
		tempfactor = 0.0f       ;
		altLoc 	   = 0;
		parent     = null;
		bonds      = Collections.emptyList();
		charge     = 0				;
	}
	
	/** 
	 * {@inheritDoc}
	 */
	@Override
	public void   setName(String s) { name = s ;}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getName()         { return name ;}

	/** 
	 * {@inheritDoc}
	 */
	@Override
	public void setPDBserial(int i) { pdbserial = i    ; }

	/** 
	 * {@inheritDoc}
	 */
	@Override
	public int  getPDBserial()      { return pdbserial ; }

	/** 
	 * {@inheritDoc}
	 */
	@Override
	public void     setCoords( double[] c ) { coords = c   ; }

	/** 
	 * {@inheritDoc}
	 */
	@Override
	public double[] getCoords()            { return coords ; }

	@Override
	public void setX(double x) {
		coords[0] = x ;
	}
	@Override
	public void setY(double y) {
		coords[1] = y ;
	}
	@Override
	public void setZ(double z) {
		coords[2] = z ;
	}

	/** 
	 * {@inheritDoc} 
	 */
	@Override
	public double getX() { return coords[0]; }

	/** 
	 * {@inheritDoc}
	 */
	@Override
	public double getY() { return coords[1]; }

	/** 
	 * {@inheritDoc}
	 */
	@Override
	public double getZ() { return coords[2]; }

	/** 
	 * Set alternate Location.
	 * @see #getAltLoc
	 */
	@Override
	public void setAltLoc(Character c) {
		// after changing altLoc from Character to char, we do this to keep the interface the same as it used to be - JD 2016-01-27
		if (c==null) 
			altLoc = 0;
		else 
			altLoc = c ;
	}
	
	/** 
	 * Get alternate Location.
	 * @return a Character object representing the alt loc value
	 * @see #setAltLoc
	 */
	@Override
	public Character getAltLoc() {
		// after changing altLoc from Character to char, we do this to keep the interface the same as it used to be - JD 2016-01-27
		if (altLoc==0 ) return null;
		return altLoc ;
	}

	@Override
	public String toString() {
		return name + " " + element + " " + pdbserial + " " + coords[0] + " " + coords[1] + " " + coords[2];
	}

	@Override
	public void   setOccupancy(float occu){ 
		occupancy = occu ;
	}
	
	@Override
	public float getOccupancy(){ 
		return occupancy; 
	}

	@Override
	public void   setTempFactor(float temp) { 
		tempfactor = temp ;
	}
	
	@Override
	public float getTempFactor() { 
		return tempfactor; 
	}

	/** returns and identical copy of this  object .
	 * @return  and identical copy of this  object
	 */
	@Override
	public Object clone() {
		AtomImpl n = new AtomImpl();
		n.setOccupancy(getOccupancy());
		n.setTempFactor(getTempFactor());
		n.altLoc = altLoc; // since char is a primitive we can do this (to avoid going through getter/setter that check for nulls)
		n.setCharge(getCharge());
		double[] coords = getCoords();
		n.setX(coords[0]);
		n.setY(coords[1]);
		n.setZ(coords[2]);
		n.setPDBserial(getPDBserial());
		n.setName(getName());
		n.setElement(getElement());
		
		// TODO bonds are not cloned here, do we need to clone them? - JD 2016-01-27
		
		return n ;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void setGroup(Group parent){
		this.parent = parent;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Group getGroup(){
		return parent;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Element getElement() {
		return element;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void setElement(Element e) {
		this.element = e;

	}

	@Override
	public String toPDB() {

		return FileConvert.toPDB(this);
	}

	@Override
	public void toPDB(StringBuffer buf) {
		FileConvert.toPDB(this,buf);

	}

	@Override
	public List<Bond> getBonds() {
		return bonds;
	}

	@Override
	public void addBond(Bond bond) {
		if (bonds.isEmpty()) {
			// most atoms have a maximum of 3 heavy atom neighbors, so use this as the default size
			bonds = new ArrayList<Bond>(3);
		}
		bonds.add(bond);
	}
	
	@Override
	public short getCharge() {
		// Get the charge
		return charge;
	}
	
	@Override
	public void setCharge(short inputCharge) {
		// Set the charge
		charge = inputCharge;
		
	}
}
