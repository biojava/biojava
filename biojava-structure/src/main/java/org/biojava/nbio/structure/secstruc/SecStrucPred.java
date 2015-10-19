package org.biojava.nbio.structure.secstruc;

import org.biojava.nbio.structure.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

/** 
 * Calculate and assign the secondary structure (SS) to the 
 * Groups of a Structure object. This object also stores the result
 * of the prediction.
 * <p>
 * The rules for SS calculation are the ones defined by DSSP:
 * Kabsch,W. and Sander,C. (1983) Biopolymers 22, 2577-2637.
 * original DSSP article see at:
 * <a href="http://www.cmbi.kun.nl/gv/dssp/dssp.pdf">dssp.pdf</a>. 
 * Some parts are also taken from: T.E.Creighton, Proteins - 
 * Structure and Molecular Properties, 2nd Edition, Freeman 1994.
 * 
 * @author Andreas Prlic
 * @author Aleix Lafita
 * 
 */
public class SecStrucPred {
	
	/** 
	 * DSSP assigns helices one residue shorter at each end, because the
	 * residues at (i-1) and (i+n+1) are not assigned helix type although
	 * they contain a consistent turn (H-bond). If this parameter
	 * is true, the helices will be the length of the original DSSP
	 * convention. If it is false, they will be two residue longer.
	 */
	private static final boolean DSSP_HELICES = true;

	private static final Logger logger = 
			LoggerFactory.getLogger(SecStrucPred.class);

	/** min distance between two residues */
	public static final double MINDIST = 0.5;

	/** min distance of two CA atoms if H-bonds are allowed to form */
	public static final double CA_MIN_DIST = 9.0;
	
	/** max distance CA atoms in peptide bond (backbone discontinuity) */
	public static final double MAX_PEPTIDE_BOND_LENGTH = 2.5;

	/** Minimal H-bond energy in cal/mol */
	public static final int HBONDLOWENERGY  = -9900;

	/** higher limit for H-bond energy */
	public static final double HBONDHIGHENERGY = -500.0;

	/** constant for electrostatic energy
	 * <pre>
	 *      f  *  q1 *   q2  *  scale
	 * Q = -332 * 0.42 * 0.20 * 1000.0
	 *</pre>
	 *
	 * q1 and q2 are partial charges which are placed on the C,O
	 * (+q1,-q1) and N,H (-q2,+q2)
	 */
	public static final double Q = -27888.0;

	private SecStrucGroup[] groups;
	private List<Ladder> ladders;
	private List<BetaBridge> bridges;

	public SecStrucPred(){
		ladders = new ArrayList<Ladder>();
		bridges = new ArrayList<BetaBridge>();
	}

	/** 
	 * Predicts the secondary structure of this Structure object,
	 * using a DSSP implementation.
	 *
	 * @param s Structure to predict the SS
	 * @param assign sets the SS information to the Groups of s
	 * @return a List of SS annotation objects
	 */
	public List<SecStrucState> predict(Structure s, boolean assign) 
			throws StructureException {

		groups = initGroupArray(s);

		if (groups.length < 5) {
			// not enough groups to do anything
			throw new StructureException("Not enough backbone groups in the"
					+ " Structure to calculate the secondary structure ("
					+ groups.length+" given, minimum 5)" );
		}

		calculateHAtoms();
		calculateHBonds();
		calculateDihedralAngles();

		calculateTurns();
		buildHelices();
		
		detectBends();
		detectStrands();
		
		List<SecStrucState> secstruc = new ArrayList<SecStrucState>();
		for (SecStrucGroup sg : groups){
			SecStrucState ss = (SecStrucState) 
					sg.getProperty(Group.SEC_STRUC);
			//Add to return list and assign to original if flag is true
			secstruc.add(ss);
			if (assign) sg.getOriginal().setProperty(Group.SEC_STRUC, ss);
		}
		return secstruc;	
	}

	private void detectStrands() {

		//Find all the beta bridges of the structure
		for (int i = 1; i < groups.length-1; i++) findBridges(i);

		//Create Ladders
		createLadders();
		
		//Detect beta bulges between ladders
		connectLadders();

		//AND store SS assignments for Sheets, Strands and Bridges
		updateSheets();
	}
	
	private void createLadders(){
		
		for (BetaBridge b : bridges){
			boolean found = false;
			for (Ladder ladder : ladders){
				if (shouldExtendLadder(ladder, b)) {
					found = true;
					ladder.to++; //we go forward in this direction
					switch(b.type){
					case parallel:
						ladder.lto++; //increment second strand
						break;
					case antiparallel:
						ladder.lfrom--; //decrement second strand
						break;
					}
					break;
				}
			}
			if (!found){
				//Create new ladder with a single Bridge
				Ladder l = new Ladder();
				l.from = b.partner1;
				l.to = b.partner1;
				l.lfrom = b.partner2;
				l.lto = b.partner2;
				l.btype = b.type;
				ladders.add(l);
			}
		}
	}


	private void updateSheets() {
		
		logger.debug(" got " +ladders.size() + "  ladders!");
		
		for (Ladder ladder : ladders){
			logger.debug(ladder.toString());

			for (int lcount = ladder.from; lcount <= ladder.to; lcount++) {

				SecStrucState state = getSecStrucState(lcount);
				SecStrucType stype = state.getType();

				int diff = ladder.from - lcount;
				int l2count = ladder.lfrom - diff ;

				SecStrucState state2 = getSecStrucState(l2count);
				SecStrucType stype2 = state2.getType();

				if ( ladder.from != ladder.to ) {
					setSecStrucType(lcount, SecStrucType.extended);
					setSecStrucType(l2count, SecStrucType.extended);
				}
				else {
					if ( !stype.isHelixType() && 
							( !stype.equals(SecStrucType.extended)))
						setSecStrucType(lcount,SecStrucType.bridge);

					if ( ! stype2.isHelixType() &&
							(! stype2.equals(SecStrucType.extended)))
						setSecStrucType(l2count,SecStrucType.bridge);
				}
			}

			// Check if two ladders are connected. both sides are 'E'

			if (ladder.connectedTo == 0) continue;
			Ladder conladder = ladders.get(ladder.connectedTo);

			if (ladder.btype.equals(BridgeType.antiparallel)) {
				/* set one side */
				for (int lcount = ladder.from; lcount <= conladder.to;
						lcount++) {
					setSecStrucType(lcount, SecStrucType.extended);

				}
				/* set other side */
				for (int lcount = conladder.lto;
						lcount <= ladder.lfrom;
						lcount++) {
					setSecStrucType(lcount, SecStrucType.extended);
				}

			} else {
				/* set one side */
				for ( int lcount = ladder.from;
						lcount <= conladder.to;
						lcount++) {

					setSecStrucType(lcount, SecStrucType.extended);
				}
				/* set other side */
				for ( int lcount =  ladder.lfrom;
						lcount <= conladder.lto;
						lcount++) {

					setSecStrucType(lcount, SecStrucType.extended);
				}
			}
		}
	}

	private void connectLadders() {

		for (int i = 0 ; i < ladders.size(); i++) {
			for ( int j = i ; j < ladders.size(); j++){
				Ladder l1 = ladders.get(i);
				Ladder l2 = ladders.get(j);
				if (hasBulge(l1,l2)) {
					l1.connectedTo = j;
					l2.connectedFrom = i;
					logger.debug("Bulge from " + i + " to " + j);
				}
			}
		}


	}

	/**
	 * For beta structures, we define explicitly: a bulge-linked 
	 * ladder consists of two (perfect) ladder or bridges of the 
	 * same type connected by at most one extra residue on one 
	 * strand and at most four extra residues on the other strand,
	 * all residues in bulge-linked ladders are marked "E,"
	 * including the extra residues.
	 */
	private boolean hasBulge(Ladder l1, Ladder l2) {
		
		boolean bulge = ((l1.btype.equals(l2.btype)) &&
				(l2.from - l1.to < 6) &&
				(l1.to < l2.from) &&
				(l2.connectedTo == 0));

		if (!bulge) return bulge;

		switch(l1.btype){
		case parallel:
			bulge = ( (l2.lfrom - l1.lto > 0) &&
					((( l2.lfrom -l1.lto < 6) &&
							(l2.from - l1.to < 3)) ||
							( l2.lfrom - l1.lto <3)));
			
			break;
			
		case antiparallel:
			bulge = ( (l1.lfrom - l2.lto > 0) &&
					(((l1.lfrom -l2.lto < 6) &&
							( l2.from - l1.to < 3)) ||
							(l1.lfrom - l2.lto < 3)));
			
			break;
		}
		
		return bulge;
	}

	private void registerBridge(int i, int j, BridgeType btype) {
				
		BetaBridge bridge = new BetaBridge(i,j,btype);
		
		getSecStrucState(i).setBridge(bridge);
		getSecStrucState(j).setBridge(bridge);

		bridges.add(bridge);
	}

	/**
	 * Conditions to extend a ladder with a given beta Bridge:
	 * <li>The bridge and ladder are of the same type.
	 * <li>The smallest bridge residue is sequential to the first
	 * 		strand ladder.
	 * <li>The second bridge residue is either sequential (parallel)
	 * 		or previous (antiparallel) to the second strand of the ladder
	 * </li>
	 * @param ladder the ladder candidate to extend
	 * @param b the beta bridge that would extend the ladder
	 * @return true if the bridge b extends the ladder
	 */
	private boolean shouldExtendLadder(Ladder ladder, BetaBridge b) {

		//Only extend if they are of the same type
		boolean sameType = b.type.equals(ladder.btype);
		if (!sameType) return false;
		
		//Only extend if residue 1 is sequential to ladder strand
		boolean sequential = (b.partner1 == ladder.to+1);
		if (!sequential) return false;
		
		switch(b.type){
		case parallel:
			//Residue 2 should be sequential to second strand
			if (b.partner2 == ladder.lto+1) return true;
			break;
		case antiparallel:
			//Residue 2 should be previous to second strand
			if (b.partner2 == ladder.lfrom-1) return true;
			break;
		}
		return false;
	}

	/**
	 * Two nonoverlapping stretches of three residues each, i-1,i,i+1 and
	 * j-1,j,j+1, form either a parallel or antiparallel bridge, depending on
	 * which of two basic patterns is matched. We assign a bridge between
	 * residues i and j if there are two H bonds characteristic of beta-
	 * structure; in particular:
	 * <p>
	 * Parallel Bridge(i,j) =: [Hbond(i-1,j) and Hbond(j,i+1)] 
	 * 							or [Hbond(j-1,i) and Hbond(i,j+1)]
	 * <p>
	 * Antiparallel Bridge(i,j) =: [Hbond(i,j) and Hbond(j,i)] 
	 * 								or [Hbond(i-1,j+1) and Hbond(j-1,i+1)]
	 */
	private void findBridges(int i) {
		
		for (int j = i+3; j < groups.length-1; j++){

			BridgeType btype = null;

			if ((isBonded(i-1,j) && isBonded(j,i+1)) ||
					(isBonded(j-1,i) && isBonded(i,j+1))) {
				btype = BridgeType.parallel;
			}
			else if ((isBonded(i,j) && isBonded(j,i)) ||
					(isBonded(i-1,j+1) && (isBonded(j-1,i+1)))) {
				btype = BridgeType.antiparallel;
			}
			
			if (btype != null){
				registerBridge(i, j, btype);
			}
		}

	}

	private void detectBends() {

		for (int i = 2 ; i < groups.length-2 ;i++){
						
			//Check if all atoms form peptide bonds (backbone discontinuity)
			boolean bonded = true;
			for (int k=0; k<4; k++){
				int index = i+k-2;
				Atom C = groups[index].getC();
				Atom N = groups[index+1].getN();
				//Peptide bond C-N
				if (Calc.getDistance(C, N) > MAX_PEPTIDE_BOND_LENGTH){
					bonded = false;
					break;
				}
			}
			if (!bonded) continue;
			
			SecStrucGroup im2 = groups[i-2];
			SecStrucGroup g = groups[i];
			SecStrucGroup ip2 = groups[i+2];

			Atom caim2 = im2.getCA();
			Atom cag   = g.getCA();
			Atom caip2 = ip2.getCA();
			
			//Create vectors ( Ca i to Ca i-2 ) ; ( Ca i to CA i + 2 )
			Atom caminus2 = Calc.subtract(caim2,cag);
			Atom caplus2  = Calc.subtract(cag,caip2);

			double angle = Calc.angle(caminus2, caplus2);

			SecStrucState state = getSecStrucState(i); 
			state.setKappa((float) angle);

			//Angles = 360 should be discarded
			if (angle > 70.0 && angle < 359.99) {
				setSecStrucType(i, SecStrucType.bend);
				state.setBend(true);
			}
		}
	}

	private void calculateDihedralAngles() throws StructureException {

		// dihedral angles
		// phi: C-N-CA-C
		// psi: N-CA-C-N
		// Chi1: N-CA-CB-CG, N-CA-CB-OG(SER),N-CA-CB-OG1(Thr),
		// N-CA-CB-CG1(ILE/VAL), N-CA-CB-SG(CYS)
		// Omega: CA-C-N-CA

		for (int i=0 ; i < groups.length-1 ;  i++){

			SecStrucGroup a = groups[i];
			SecStrucGroup b = groups[i+1];

			Atom a_N   = a.getN();
			Atom a_CA  = a.getCA();
			Atom a_C  = a.getC();

			Atom b_N  = b.getN();
			Atom b_CA = b.getCA();
			Atom b_C  = b.getC();

			double phi = Calc.torsionAngle(a_C,b_N,b_CA,b_C);
			double psi = Calc.torsionAngle(a_N,a_CA,a_C,b_N);
			double omega = Calc.torsionAngle(a_CA,a_C,b_N,b_CA);

			SecStrucState state1 = (SecStrucState) 
					a.getProperty(Group.SEC_STRUC);
			SecStrucState state2 = (SecStrucState) 
					b.getProperty(Group.SEC_STRUC);
			
			state2.setPhi(phi);
			state1.setPsi(psi);
			state1.setOmega(omega);
		}
	}

	@Override
	public String toString() {
		return printDSSP();
	}
	
	/**
	 * Generate a DSSP file format ouput String of this SS prediction.
	 * @return String in DSSP output file format
	 */
	public String printDSSP() {
		
		StringBuffer buf = new StringBuffer();
		String nl = System.getProperty("line.separator");
		
		//Header Line
		buf.append("==== Secondary Structure Definition by BioJava"
				+ " DSSP implementation, Version October 2015 ===="+nl);
		
		//First line with column definition
		buf.append("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     "
				+ "N-H-->O    O-->H-N    N-H-->O    O-->H-N    "
				+ "TCO  KAPPA ALPHA  PHI    PSI    "
				+ "X-CA   Y-CA   Z-CA ");

		for (int i =0 ; i < groups.length ;i++){
			buf.append(nl);
			SecStrucState ss = getSecStrucState(i);
			buf.append(ss.printDSSPline(i));
		}

		return buf.toString();
	}
	
	/**
	 * Generate a summary of this SS prediction with information about 
	 * the three types of helix turns in different row sequences.
	 * <p>
	 * This is similar to the summary output of Jmol, and useful to visualize
	 * the helix patterns.
	 * 
	 * @return String helix summary
	 */
	public String printHelixSummary() {
		
		StringBuffer g = new StringBuffer(); //3-10 helix
		StringBuffer h = new StringBuffer(); //alpha helix
		StringBuffer i = new StringBuffer(); //pi-helix
		StringBuffer ss = new StringBuffer(); //SS summary
		StringBuffer aa = new StringBuffer(); //AA one-letter
		String nl = System.getProperty("line.separator");
		
		g.append(	"3 turn: ");
		h.append(	"4 turn: ");
		i.append(	"5 turn: ");
		ss.append(	"SS:     ");
		aa.append(	"AA:     ");
		
		for (int k = 0; k < groups.length; k++){
			
			SecStrucState state = getSecStrucState(k);
			g.append(state.getTurn()[0]);
			h.append(state.getTurn()[1]);
			i.append(state.getTurn()[2]);
			ss.append(state.getType());
			aa.append(StructureTools.get1LetterCode(groups[k].getPDBName()));
		}
		
		return g.toString()+nl+h.toString()+nl+
				i.toString()+nl+ss.toString()+nl+aa.toString();
	}
	
	/**
	 * Generate a FASTA sequence with the SS annotation letters in the
	 * aminoacid sequence order. 
	 * @return String in FASTA sequence format
	 */
	public String printFASTA() {
		
		StringBuffer buf = new StringBuffer();
		String nl = System.getProperty("line.separator");
		buf.append(">"+groups[0].getChain().getStructure().getIdentifier()+nl);
		
		for (int g = 0; g < groups.length; g++){
			buf.append(getSecStrucState(g).getType());
		}
		return buf.toString();
	}
	
	@Override
	public boolean equals(Object o){
		
		if (!(o instanceof SecStrucPred)) return false;
		else {
			SecStrucPred ss = (SecStrucPred) o;
			if (groups.length != ss.groups.length) return false;
			
			for (int g=0; g<groups.length; g++){
				SecStrucInfo g1 = getSecStrucState(g);
				SecStrucInfo g2 = ss.getSecStrucState(g);
				if (!g1.equals(g2)) return false;
			}
			return true;
		}
	}

	private static SecStrucGroup[] initGroupArray(Structure s) {
		List<SecStrucGroup> groupList = new ArrayList<SecStrucGroup>();
		
		for ( Chain c : s.getChains()){

			for (Group g : c.getAtomGroups()){
				
				//We can also calc secstruc if it is a modified amino acid
				if ( g.hasAminoAtoms()) {

					SecStrucGroup sg = new SecStrucGroup();
					sg.setResidueNumber(g.getResidueNumber());
					sg.setPDBFlag(true);
					sg.setPDBName(g.getPDBName());
					sg.setChain(g.getChain());

					Atom N = g.getAtom(StructureTools.N_ATOM_NAME);
					Atom CA =  g.getAtom(StructureTools.CA_ATOM_NAME);
					Atom C = g.getAtom(StructureTools.C_ATOM_NAME);
					Atom O =  g.getAtom(StructureTools.O_ATOM_NAME);
					if ( N == null || CA == null || C == null || O == null)
						continue;

					sg.setN((Atom)   N.clone());
					sg.setCA((Atom) CA.clone());
					sg.setC((Atom)   C.clone());
					sg.setO((Atom)  O.clone());
					sg.setOriginal(g);
					
					SecStrucState state = new SecStrucState(sg, 
							SecStrucInfo.BIOJAVA_ASSIGNMENT, 
							SecStrucType.coil);
					
					sg.setProperty(Group.SEC_STRUC, state);
					groupList.add(sg);
				}
			}
		}
		return groupList.toArray(new SecStrucGroup[groupList.size()]);
	}

	/** 
	 * Calculate the coordinates of the H atoms. They are usually
	 * missing in the PDB files as only few experimental methods allow
	 * to resolve their location.
	 */
	private void calculateHAtoms() throws StructureException {

		for ( int i = 0 ; i < groups.length-1  ; i++) {

			SecStrucGroup a  = groups[i];
			SecStrucGroup b  = groups[i+1];

			if ( !b.hasAtom("H") ) {
				//Atom H = calc_H(a.getC(), b.getN(), b.getCA());
				Atom H = calcSimple_H(a.getC(), a.getO(), b.getN());
				b.setH(H);
			}
		}
	}

	/** 
	 * Calculate the HBonds between different groups.
	 * see Creighton page 147 f
	 */
	private void calculateHBonds() throws StructureException {

		if (groups.length < 5) return;

		for (int i=0 ; i < groups.length ; i++){

			SecStrucGroup one = groups[i];

			for (int j=i+1 ; j<groups.length ; j++){

				SecStrucGroup two = groups[j];
				//TODO use contacts package to speed up n^2 distance
				//if distance too big - for sure no HBonds - sppedup
				double dist = Calc.getDistance(one.getCA(),two.getCA());
				if (dist >= CA_MIN_DIST) continue;

				checkAddHBond(i,j);

				//"backwards" hbonds are not allowed
				if (j!=(i+1)) checkAddHBond(j,i);
			}
		}
	}

	private void checkAddHBond(int i, int j){
		
		SecStrucGroup one = groups[i];

		if (one.getPDBName().equals("PRO")){
			logger.debug("Ignore: PRO " + one.getResidueNumber());
			return;
		}
		if (!one.hasAtom("H")) {
			logger.debug("Residue "+one.getResidueNumber()+" has no H");
			return;
		}

		SecStrucGroup two = groups[j];
		
		double energy = 0;
		
		try {
			energy = calculateHBondEnergy(one,two);
		} catch (Exception e){
			logger.warn("Energy calculation failed", e);
			return;
		}
		logger.debug("Energy between positions ("+i+","+j+"): "+energy);

		trackHBondEnergy(i,j,energy);
	}

	/**
	 * Calculate HBond energy of two groups in cal/mol
	 * see Creighton page 147 f
	 * <p>
	 * Jeffrey, George A., An introduction to hydrogen bonding, 
	 * Oxford University Press, 1997.
	 * categorizes hbonds with donor-acceptor distances of
	 * 2.2-2.5 &aring; as "strong, mostly covalent",
	 * 2.5-3.2 &aring; as "moderate, mostly electrostatic",
	 * 3.2-4.0 &aring; as "weak, electrostatic".
	 * Energies are given as 40-14, 15-4, and <4 kcal/mol respectively.
	 */
	private static double calculateHBondEnergy(SecStrucGroup one, 
			SecStrucGroup two) throws StructureException {

		Atom N = one.getN();
		Atom H = one.getH();

		Atom O = two.getO();
		Atom C = two.getC();

		double dno = Calc.getDistance(O,N);
		double dhc = Calc.getDistance(C,H);
		double dho = Calc.getDistance(O,H);
		double dnc = Calc.getDistance(C,N);

		logger.debug("     cccc: " + one.getResidueNumber() + 
				" " + one.getPDBName() + " " +two.getResidueNumber()+ 
				" " + two.getPDBName() + String.format(" O ("+
				O.getPDBserial()+")..N ("+ N.getPDBserial()+
				"):%4.1f  |  ho:%4.1f - hc:%4.1f + nc:%4.1f - no:%4.1f ", 
				dno,dho,dhc,dnc,dno));

		//there seems to be a contact!
		if ( (dno < MINDIST) || (dhc < MINDIST) || 
				(dnc < MINDIST) || (dno < MINDIST)) {
			return HBONDLOWENERGY;
		}

		double e1 = Q / dho - Q / dhc;
		double e2 = Q / dnc - Q / dno;

		double energy = e1 + e2;

		logger.debug(String.format("      N (%d) O(%d): %4.1f : %4.2f ",
				N.getPDBserial(),O.getPDBserial(), (float) dno, energy));
		
		//Avoid too strong energy
		if (energy > HBONDLOWENERGY) return energy;

		return HBONDLOWENERGY ;
	}

	/**
	 * Store Hbonds in the Groups.
	 * DSSP allows two HBonds per aminoacids to allow bifurcated bonds.
	 */
	private  void trackHBondEnergy(int i, int j, double energy) {
		
		if (groups[i].getPDBName().equals("PRO")) {
			logger.debug("Ignore: PRO " + groups[i].getResidueNumber());
			return;
		}

		SecStrucState stateOne = getSecStrucState(i);
		SecStrucState stateTwo = getSecStrucState(j);

		double acc1e = stateOne.getAccept1().getEnergy();
		double acc2e = stateOne.getAccept2().getEnergy();

		double don1e = stateTwo.getDonor1().getEnergy();
		double don2e = stateTwo.getDonor2().getEnergy();

		//Acceptor: N-H-->O
		if (energy < acc1e) {
			logger.debug(energy +"<"+acc1e);
			stateOne.setAccept2(stateOne.getAccept1());

			HBond bond = new HBond();
			bond.setEnergy(energy);
			bond.setPartner(j);

			stateOne.setAccept1(bond);

		} else if ( energy < acc2e ) {
			logger.debug(energy +"<"+acc2e);
			
			HBond bond = new HBond();
			bond.setEnergy(energy);
			bond.setPartner(j);

			stateOne.setAccept2(bond);
		}

		//The other side of the bond: donor O-->N-H
		if (energy <  don1e) {
			logger.debug(energy +"<"+don1e);
			stateTwo.setDonor2(stateTwo.getDonor1());

			HBond bond = new HBond();
			bond.setEnergy(energy);
			bond.setPartner(i);

			stateTwo.setDonor1(bond);

		} else if ( energy < don2e ) {
			logger.debug(energy +"<"+don2e);

			HBond bond = new HBond();
			bond.setEnergy(energy);
			bond.setPartner(i);

			stateTwo.setDonor2(bond);
		}
	}

	/** 
	 * Detect helical turn patterns.
	 */
	private void calculateTurns(){

		for (int i = 0 ; i< groups.length; i++){
			for (int turn = 3; turn <= 5; turn++) {
				
				if (i+turn >= groups.length) continue;

				//Check for H bond from NH(i+n) to CO(i)
				if (isBonded(i, i+turn)) {
					logger.debug("Turn at ("+i+","+(i+turn)+") turn "+turn);
					getSecStrucState(i).setTurn('>', turn);
					getSecStrucState(i+turn).setTurn('<', turn);
					//Bracketed residues get the helix number
					for (int j=i+1; j<i+turn; j++){
						Integer t = turn;
						char helix = t.toString().charAt(0);
						getSecStrucState(j).setTurn(helix, turn);
					}
				}
			}
		}
	}

	/** 
	 * Test if two groups are forming an H-Bond. The bond tested is
	 * from the CO of group i to the NH of group j. Acceptor (i) and
	 * donor (j). The donor of i has to be j, and the acceptor of j 
	 * has to be i.
	 * DSSP defines H-Bonds if the energy < -500 cal/mol.
	 * 
	 * @param one group one
	 * @param two group two
	 * @return flag if the two are forming an Hbond
	 */
	private boolean isBonded(int i, int j) {

		SecStrucState one = getSecStrucState(i);
		SecStrucState two = getSecStrucState(j);

		double don1e = one.getDonor1().getEnergy();
		double don2e = one.getDonor2().getEnergy();
		double acc1e = two.getAccept1().getEnergy();
		double acc2e = two.getAccept2().getEnergy();
		
		int don1p = one.getDonor1().getPartner();
		int don2p = one.getDonor2().getPartner();
		int acc1p = two.getAccept1().getPartner();
		int acc2p = two.getAccept2().getPartner();

		//Either donor from i is j, or accept from j is i
		boolean hbond = (don1p == j && don1e < HBONDHIGHENERGY) ||
				(don2p == j && don2e < HBONDHIGHENERGY) ||
				(acc1p == i && acc1e < HBONDHIGHENERGY) ||
				(acc2p == i && acc2e < HBONDHIGHENERGY);
				
		if (hbond){				
			logger.debug("*** H-bond from CO of " + i + " to NH of " + j);
			return true;
		}
		return false ;
	}

	/**
	 * Use unit vectors NC and NCalpha Add them. Calc unit vector and
	 * substract it from N.
	 * C coordinates are from amino acid i-1
	 * N, CA atoms from amino acid i
	 *
	 * @link http://openbioinformatics.blogspot.com/
	 * 		2009/08/how-to-calculate-h-atoms-for-nitrogens.html
	 */
	@SuppressWarnings("unused")
	private static Atom calc_H(Atom C, Atom N, Atom CA)
			throws StructureException {

		Atom nc  = Calc.subtract(N,C);
		Atom nca = Calc.subtract(N,CA);

		Atom u_nc  = Calc.unitVector(nc)   ;
		Atom u_nca = Calc.unitVector(nca);

		Atom added = Calc.add(u_nc,u_nca);

		Atom U = Calc.unitVector(added);

		// according to Creighton distance N-H is 1.03 +/- 0.02A
		Atom H = Calc.add(N,U);

		H.setName("H");
		// this atom does not have a pdbserial number ...
		return H;

	}

	private static Atom calcSimple_H(Atom c, Atom o, Atom n) 
			throws StructureException {

		Atom h = Calc.subtract(c,o);
		double dist = Calc.getDistance(o,c);
		//System.out.println(dist);
		double x = n.getX() + h.getX() / dist;
		double y = n.getY() + h.getY() / dist;
		double z = n.getZ() + h.getZ() / dist;

		h.setX(x);
		h.setY(y);
		h.setZ(z);

		h.setName("H");
		return h;
	}

	private void buildHelices(){

		//Alpha-helix (i+4), 3-10-helix (i+3), Pi-helix (i+5)
		checkSetHelix(4, SecStrucType.helix4);
		checkSetHelix(3, SecStrucType.helix3);
		checkSetHelix(5, SecStrucType.helix5);

		checkSetTurns();
	}

	private void checkSetTurns() {
		
		SecStrucType type = SecStrucType.turn;
		
		for (int idx = 0; idx < 3; idx++) {
			for (int i = 0; i < groups.length-1; i++) {
				
				SecStrucState state = getSecStrucState(i);
				char[] turn = state.getTurn();
			
				//Any turn opening matters
				if (turn[idx] == '>' || turn[idx] == 'X') {
					//Mark following n residues as turn
					for (int k=1; k<idx+3; k++){
						setSecStrucType(i+k, type);
					}
					if (!DSSP_HELICES) {
						setSecStrucType(i, type);
						setSecStrucType(i+idx+3, type);
					}
				}
			}
		}
	}

	/**
	 * A minimal helix is defined by two consecutive n-turns.
	 * For example, a 4-helix, of minimal length 4 from residues 
	 * i to (i+3), requires turns (of type 4) at residues (i-1) and i.
	 * <p>
	 * Note that the orignal DSSP implementation does not assign
	 * helix type to residue (i-1) and residue (i+n+1), although 
	 * they contain a helix turn. As they state in the original paper,
	 * "the helices are one residue shorter than they would be according
	 * to rule 6.3 of IUPAC-IUB".
	 * 
	 * @param n
	 * @param type
	 */
	private void checkSetHelix(int n, SecStrucType type){

		int idx = n - 3;
		logger.debug("Set helix " + type + " " + n + " " + idx);
		
		for (int i = 1; i < groups.length-n; i++) {

			SecStrucState state = getSecStrucState(i);
			SecStrucState previousState = getSecStrucState(i-1);
			
			//Check that no other helix was assgined to this range
			if (state.getType().compareTo(type) < 0) continue;
			if (getSecStrucState(i+1).getType().compareTo(type) < 0) continue;

			char turn = state.getTurn()[idx];
			char pturn = previousState.getTurn()[idx];

			//Two consecutive n-turns present to define a n-helix
			if ((turn=='>' || turn=='X') && (pturn=='>' || pturn=='X')) {
				//Mark following n residues as turn
				for (int k=0; k<n; k++){
					setSecStrucType(i+k, type);
				}
				if (!DSSP_HELICES) {
					setSecStrucType(i-1, type);
					setSecStrucType(i+n, type);
				}
			}
		}
	}

	/**
	 * Set the new type only if it has more preference than the
	 * current residue SS type.
	 * @param pos
	 * @param type
	 */
	private void setSecStrucType(int pos, SecStrucType type){
		SecStrucState ss = getSecStrucState(pos);
		if (type.compareTo(ss.getType()) < 0) ss.setType(type);
	}

	private SecStrucState getSecStrucState(int pos){
		Group g = groups[pos];
		SecStrucState state = (SecStrucState) g.getProperty(Group.SEC_STRUC);
		return state;
	}

}