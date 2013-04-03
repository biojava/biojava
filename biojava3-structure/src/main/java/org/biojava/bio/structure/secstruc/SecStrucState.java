/*
 *                    PDB web development code
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
 *
 * Created on Aug 5, 2009
 * Created by ap3
 *
 */

package org.biojava.bio.structure.secstruc;

public class SecStrucState{
	public boolean isBend() {
		return bend;
	}

	public void setBend(boolean bend) {
		this.bend = bend;
	}

	double phi;
	double psi;
	double omega;

	float kappa;
	
	HBond accept1;
	HBond accept2;
	HBond donor1;
	HBond donor2;

	boolean[] turn;
	boolean bend;
	
	SecStrucType secStruc;
	SecStrucType threeState;

	public SecStrucState(){
		phi = 360;
		psi = 360 ;
		omega = 360;
		accept1 = new HBond();
		accept2 = new HBond();
		donor1  = new HBond();
		donor2  = new HBond();
		turn = new boolean[3];
		turn[0] = false;
		turn[1] = false;
		turn[2] = false;
		bend = false;
		secStruc = SecStrucType.coil;
		threeState = SecStrucType.coil;
		kappa = 360;
	}

	public float getKappa() {
		return kappa;
	}

	public void setKappa(float kappa) {
		this.kappa = kappa;
	}

	public String toString(){
		StringBuffer buf = new StringBuffer();
		buf.append(secStruc.toString()+threeState.toString() + " a1:"+accept1 + " a2:" + accept2 + " d1:" + donor1 + " d2:" + donor2);
		return buf.toString();
	}

	public boolean[] getTurn() {
		return turn;
	}

	public void setTurn(boolean[] turn) {
		this.turn = turn;
	}

	public HBond getAccept1() {
		return accept1;
	}
	public void setAccept1(HBond accept1) {
		this.accept1 = accept1;
	}
	public HBond getAccept2() {
		return accept2;
	}
	public void setAccept2(HBond accept2) {
		this.accept2 = accept2;
	}
	public HBond getDonor1() {
		return donor1;
	}
	public void setDonor1(HBond donor1) {
		this.donor1 = donor1;
	}
	public HBond getDonor2() {
		return donor2;
	}
	public void setDonor2(HBond donor2) {
		this.donor2 = donor2;
	}
	public double getPhi() {
		return phi;
	}
	public void setPhi(double phi) {
		this.phi = phi;
	}
	public double getPsi() {
		return psi;
	}
	public void setPsi(double psi) {
		this.psi = psi;
	}


	public SecStrucType getSecStruc() {
		return secStruc;
	}

	public void setSecStruc(SecStrucType secStruc) {
		this.secStruc = secStruc;
	}

	public SecStrucType getThreeState() {
		return threeState;
	}

	public void setThreeState(SecStrucType threeState) {
		this.threeState = threeState;
	}


	public double getOmega() {
		return omega;
	}

	public void setOmega(double omega) {
		this.omega = omega;
	}




}
