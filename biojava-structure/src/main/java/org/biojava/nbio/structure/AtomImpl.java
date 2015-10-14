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
	String name     ;
	Element element;
	double[] coords ;
	String pdbline  ;
	int pdbserial   ;

	double occupancy ;
	double tempfactor;

	Character altLoc ;
	Group parent;
	long id;

	private List<Bond> bonds;

	public AtomImpl () {
		name     = null        ;
		element = Element.R;
		coords   = new double[3];
		pdbline  = ""          ;
		occupancy  = 0.0       ;
		tempfactor = 0.0       ;
		altLoc = ' ';
		altLoc = null;
		parent = null;
		bonds = Collections.emptyList();
	}
	/** Get the Hibernate database ID.
	 *
	 * @return the id
	 * @see #setId(long)
	 */
	public long getId() {
		return id;
	}

	/** Set the Hibernate database ID.
	 *
	 * @param id the hibernate id
	 * @see #getId()
	 */
	public void setId(long id) {
		this.id = id;
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

	/** set alternate Location.
	 * @see #getAltLoc
	 */
	@Override
	public void setAltLoc(Character c) {
		altLoc = c ;
	}
	/** get alternate Location.
	 * @return a Character object representing the alt loc value
	 * @see #setAltLoc
	 */
	@Override
	public Character getAltLoc() {
		return altLoc ;
	}

	/** string representation. */
	@Override
	public String toString() {
		return name + " " + element + " " + pdbserial + " " + coords[0] + " " + coords[1] + " " + coords[2];
	}

	@Override
	public void   setOccupancy(double occu){ occupancy = occu ;} ;
	@Override
	public double getOccupancy(){ return occupancy; } ;

	@Override
	public void   setTempFactor(double temp){ tempfactor = temp ;} ;
	@Override
	public double getTempFactor(){ return tempfactor; } ;

	/** returns and identical copy of this  object .
	 * @return  and identical copy of this  object
	 */
	@Override
	public Object clone() {
		AtomImpl n = new AtomImpl();
		n.setOccupancy(getOccupancy());
		n.setTempFactor(getTempFactor());
		n.setAltLoc(getAltLoc());
		double[] coords = getCoords();
		n.setX(coords[0]);
		n.setY(coords[1]);
		n.setZ(coords[2]);
		n.setPDBserial(getPDBserial());
		n.setName(getName());
		n.setElement(getElement());
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
}
