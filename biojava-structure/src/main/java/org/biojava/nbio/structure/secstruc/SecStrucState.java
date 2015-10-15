package org.biojava.nbio.structure.secstruc;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.StructureTools;

/**
 * This class extends the basic container for secondary structure
 * annotation, including all the information used in the DSSP algorithm.
 * 
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public class SecStrucState extends SecStrucInfo {
	
	private double phi;
	private double psi;
	private double omega;
	private float kappa;
	
	private HBond accept1; //from CO of partner to NH of this
	private HBond accept2; //this is the donor of accept partner
	private HBond donor1; //from CO of this to NH of partner
	private HBond donor2; //this is the acceptor of donor partner

	//Symbols: starting '>', ending '<', or both 'X'. 
	//Number means bracketed n-turn residue without h-bond
	private char[] turn;
	private boolean bend;
	
	private BetaBridge bridge1;
	private BetaBridge bridge2;

	public SecStrucState(Group g, String ass, SecStrucType t){
		super(g, ass, t);
		
		phi = 360;
		psi = 360;
		omega = 360;
		
		accept1 = new HBond();
		accept2 = new HBond();
		donor1  = new HBond();
		donor2  = new HBond();
		
		bridge1 = null;
		bridge2 = null;
		
		turn = new char[3];
		turn[0] = ' ';
		turn[1] = ' ';
		turn[2] = ' ';
		bend = false;
		
		kappa = 360;
	}
	
	public boolean isBend() {
		return bend;
	}

	public void setBend(boolean bend) {
		this.bend = bend;
	}

	public float getKappa() {
		return kappa;
	}

	public void setKappa(float kappa) {
		this.kappa = kappa;
	}

	public char[] getTurn() {
		return turn;
	}

	/**
	 * Set the turn column corresponding to 3,4 or 5 helix patterns.
	 * If starting > or ending < was set and the opposite is being set,
	 * the value will be converted to X. If a number was set, it will be
	 * overwritten by the new character.
	 * @param c character in the column
	 * @param t turn of the helix {3,4,5}
	 */
	public void setTurn(char c, int t) {
		if (turn[t-3] == 'X') return;
		else if (turn[t-3] == '<' && c == '>' 
				|| turn[t-3] == '>' && c == '<'){
			turn[t-3] = 'X';
		} else if (turn[t-3] == '<' || turn[t-3] == '>') return;
		else turn[t-3] = c;
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

	public double getOmega() {
		return omega;
	}

	public void setOmega(double omega) {
		this.omega = omega;
	}
	
	public BetaBridge getBridge1() {
		return bridge1;
	}

	public void setBridge(BetaBridge bridge) {
		if (bridge1 == null) bridge1 = bridge;
		else if (bridge1.equals(bridge)) return;
		else if (bridge2 == null) bridge2 = bridge;
		else if (bridge2.equals(bridge)) return; 
		else throw new IllegalStateException(
				"Cannot store more than 2 beta Bridge in one residue");
	}

	public BetaBridge getBridge2() {
		return bridge2;
	}

	public String printDSSPline(int index) {
		
		StringBuffer buf = new StringBuffer();
		
		//#
		if (index < 9) buf.append("    ");
		else if (index < 99) buf.append("   ");
		else if (index < 999) buf.append("  ");
		else buf.append(" ");
		buf.append(index + 1);
		
		//RESIDUE
		int resnum = parent.getResidueNumber().getSeqNum();
		if (resnum < 10) buf.append("    ");
		else if (resnum < 100) buf.append("   ");
		else buf.append("  ");
		buf.append(resnum);
		Character insCode = parent.getResidueNumber().getInsCode();
		if (insCode != null) buf.append(insCode);
		else buf.append(" ");
		buf.append(parent.getChainId());
		if (parent.getChainId().length() == 1) buf.append(" ");
		
		//AA
		char aaLetter = StructureTools.get1LetterCode(parent.getPDBName());
		buf.append(aaLetter+"  ");

		//STRUCTURE
		buf.append(type).append(" ");
		
		for (int t=0; t<3; t++){
			buf.append(turn[t]);
		}
		
		buf.append("  ");
		
		if (isBend()) buf.append('S');
		else buf.append(" ");
		
		buf.append(" ");

		int bp1 = 0;
		if (bridge1 != null){
			if (bridge1.partner1 != index) bp1 = bridge1.partner1+1;
			else bp1 = bridge1.partner2+1;
		}
		//TODO a clever way to do this?
		if (bp1 < 10) buf.append("   "+bp1);
		else if (bp1 < 100) buf.append("  "+bp1);
		else if (bp1 < 1000) buf.append(" "+bp1);
		else buf.append(bp1);
		
		int bp2 = 0;
		if (bridge2 != null){
			if (bridge2.partner1 != index) bp2 = bridge2.partner1+1;
			else bp2 = bridge2.partner2+1;
		}
		if (bp2 < 10) buf.append("   "+bp2);
		else if (bp2 < 100) buf.append("  "+bp2);
		else if (bp2 < 1000) buf.append(" "+bp2);
		else buf.append(bp2);
		
		//beta-sheet label TODO
		buf.append(" ");
		
		//ACC TODO
		buf.append("     ");
		
		//N-H-->O
		int p1 = accept1.getPartner();
		double e1 =  (accept1.getEnergy() / 1000.0);
		if (e1 < 0.0) p1 -= index;
		buf.append(String.format( "%6d,%4.1f",p1,e1));

		//O-->H-N
		int p2 = donor1.getPartner();
		double e2 = (donor1.getEnergy() / 1000.0);
		if (e2 < 0.0) p2 -= index;
		buf.append(String.format( "%6d,%4.1f",p2,e2));

		//N-H-->O
		int p3 = accept2.getPartner();
		double e3 =  (accept2.getEnergy() / 1000.0);
		if (e3 < 0.0) p3 -= index;
		buf.append(String.format( "%6d,%4.1f",p3,e3));

		//O-->H-N
		int p4 = donor2.getPartner();
		double e4 = (donor2.getEnergy() / 1000.0);
		if (e4 < 0.0) p4 -= index;
		buf.append(String.format( "%6d,%4.1f",p4,e4));
		
		//TCO TODO
		buf.append("        ");
		
		//KAPPA
		buf.append(String.format("%6.1f",kappa));
		
		//ALPHA TODO
		buf.append("      ");
		
		//PHI PSI
		buf.append(String.format("%6.1f %6.1f ", phi, psi));
		
		//X-CA Y-CA Z-CA
		Atom ca = parent.getAtom("CA");
		buf.append(String.format("%6.1f %6.1f %6.1f", 
				ca.getX(), ca.getY(), ca.getZ()));
		
		return buf.toString();
	}

}
