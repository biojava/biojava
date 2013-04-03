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
 * Created on 26.04.2004
 * @author Andreas Prlic
 * @since 1.5
 *
 */

package org.biojava.bio.structure.secstruc;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.PDBParseException;




/** Assign secondary structure to a Structure object.
 * tries to use the rules as defined by DSSP
 * Kabsch,W. and Sander,C. (1983) Biopolymers 22, 2577-2637.
 * original DSSP article see at : <a href="http://www.cmbi.kun.nl/gv/dssp/dssp.pdf">dssp.pdf</a>
 *
 * some bits also taken from
 * T.E.Creighton
 * Proteins - Structure and Molecular Properties
 * 2nd Edition, Freeman 1994
 */
public class SecStruc {
	private static final boolean debug = false;

	/** the minimal distance between two residues */
	public static double MINDIST       = 0.5      ;

	/** the minimal distance of two CA atoms if H-bonds are allowed to form.
	 *
	 */
	public static int CA_MIN_DIST = 9;

	/** Minimal H-bond energy in cal / mol */
	public static int HBONDLOWENERGY  = -9900   ;

	/** higher limit for H-bond energy */
	public static double HBONDHIGHENERGY = -500.0     ;

	/** constant for electrostatic energy
	 * <pre>
	 * 	     f  *   q1 *   q2 * scale
	 * Q = -332 * 0.42 * 0.20 * 1000.0
	 *</pre>
	 *
	 * q1 and q2 are partial charges which are placed on the C,O
	 * (+q1,-q1) and N,H (-q2,+q2)
	 */

	public static double Q            = -27888.0 ;
	//public static double Q            = ( -332 * 0.42 * 0.2 * 1000 ); // -27888.0

	private SecStrucGroup[] groups;

	List<DistEn> distVsEnergy;
	List<Ladder> ladders;

	public SecStruc(){
		distVsEnergy = new ArrayList<DistEn>();
		ladders = new ArrayList<Ladder>();

	}

	public static void main(String[] args){
		try {



			AtomCache cache = new AtomCache();
			cache.setAutoFetch(false);
			cache.setPath("/Users/andreas/WORK/PDB/");

			// a small one
			Structure s = cache.getStructure("5pti");
			//Structure s = pdbr.getStructureById("1bsp");
			//Structure s = pdbr.getStructureById("1co7");

			//Structure s = pdbr.getStructureById("1buz");

			//BiojavaJmol jmol = new BiojavaJmol();


			SecStruc sec = new SecStruc();
			sec.assign(s);
			//jmol.setStructure(s);
			//jmol.evalString("select *.HC2 ; spacefill 0.4 ; color green ; select 2,5 ; set display selected;");
			//System.out.println(s.toPDB());

			System.out.println(sec);

		} catch (Exception e){
			e.printStackTrace();
		}
	}

	/** assigns the secondary structure to the groups in this Structure object
	 * and set the results in the group properties.
	 *
	 * @param s
	 */
	public void assign(Structure s)
			throws StructureException {

		groups = initGroupArray(s);

		if ( groups.length < 5) {
			// not enough groups to do anything
			throw new StructureException("not enough groups in structure to calculate secondary structure ("+ groups.length+")" );
		}


		calculateHAtoms();

		/*for (int j=0 ; j<3;j++ ) {
         Group g = groups[j];
         System.out.println(g);
         for (int i=0 ; i< g.size();i++){
            Atom a = g.getAtom(i);
            System.out.println(a);
         }
      }*/

		calculateHBonds();
		calculateTurns();
		calculateDihedralAngles();
		buildHelices();

		detectBends();
		detectStrands();


	}

	private void detectStrands() {

		for (int i =1 ; i < groups.length -1 ;i++){
			testBridge(i);
		}

		// detect beta bulges
		connectLadders();

		// and store results for Sheets and Strands
		updateSheets();
	}


	private void updateSheets() {
		if ( debug)
			System.out.println(" got "  +ladders.size() + "  ladders!");
		for (Ladder ladder : ladders){

			if ( debug)
				System.out.println(ladder);

			for ( int lcount = ladder.from;
					lcount <= ladder.to;
					lcount++  ) {


				SecStrucState state = getSecStrucState(lcount);
				SecStrucType stype = state.getSecStruc();

				int diff = ladder.from - lcount;
				int l2count = ladder.lfrom - diff ;

				SecStrucState state2 = getSecStrucState(l2count);
				SecStrucType stype2 = state2.getSecStruc();

				if ( ladder.from != ladder.to ) {
					if ( ! stype.isHelixType())
						setSecStrucType(lcount, SecStrucType.extended);

					if ( ! stype2.isHelixType()){
						setSecStrucType(l2count,SecStrucType.extended);
					}			    	  				
				}
				else {
					if ( ! stype.isHelixType() && 
							( ! stype.equals(SecStrucType.extended))) {
						setSecStrucType(lcount,SecStrucType.bridge);
					}

					if ( ! stype2.isHelixType() &&
							(! stype2.equals(SecStrucType.extended))) {
						setSecStrucType(l2count,SecStrucType.bridge);
					}	
				}
			}

			// Check if two ladders are connected. both sides are 'E'  

			if ( ladder.connectedTo == 0) 
				continue;

			Ladder conladder = ladders.get(ladder.connectedTo);

			if ( ladder.getBtype().equals(BridgeType.antiparallel)) {
				/* set one side */
				for ( int lcount = ladder.from;
						lcount <= conladder.to;
						lcount++) {
					testSetExtendedSecStrucState(lcount);

				}
				/* set other side */
				for (int lcount = conladder.lto;
						lcount <= ladder.lfrom;
						lcount++) {
					testSetExtendedSecStrucState(lcount);
				}

			} else {
				/* set one side */
				for ( int lcount = ladder.from;
						lcount <= conladder.to;
						lcount++) {

					testSetExtendedSecStrucState(lcount);
				}
				/* set other side */
				for ( int lcount =  ladder.lfrom;
						lcount <= conladder.lto;
						lcount++) {

					testSetExtendedSecStrucState(lcount);
				}

			}

		}

	}

	private void testSetExtendedSecStrucState(int lcount) {
		SecStrucState state = getSecStrucState(lcount);
		SecStrucType stype = state.getSecStruc();
		if ( ! stype.isHelixType()){
			setSecStrucType(lcount,SecStrucType.extended);
		}
	}

	private void connectLadders() {

		for ( int i = 0 ; i < ladders.size(); i++) {
			for ( int j = i ; j < ladders.size() ; j++){
				Ladder l1 = ladders.get(i);
				Ladder l2 = ladders.get(j);
				if ( hasBulge(l1,l2) ) {

					l1.connectedTo = j;
					l2.connectedFrom = i;
					System.out.println("BUlge from " + i + " to " + j);
				}
			}
		}


	}

	private boolean hasBulge(Ladder l1, Ladder l2) {
		boolean bulge =  ( (l1.getBtype().equals(l2.getBtype()) ) &&
				( l2.from - l1.to < 6) &&
				( l1.to < l2.from) &&
				(l2.connectedTo == 0) );

		//	    		 ( ( ladder1->type == ladder2->type) &&
		//	  	       ( ladder2->from - ladder1->to < 6 ) &&
		//	  	       ( ladder1->to < ladder2->from) &&
		//	  	       ( ladder2->conto == 0 ) ) ;

		if ( ! bulge )
			return bulge;

		if ( l1.getBtype().equals(BridgeType.parallel)) {
			bulge = ( (l2.lfrom - l1.lto > 0) &&
					((( l2.lfrom -l1.lto < 6) &&
							(l2.from - l1.to < 3)) ||
							( l2.lfrom - l1.lto <3)));
			//case parallel:
			//			 bulge = ( ( ladder2->lfrom - ladder1->lto > 0 ) &&
			//				   ((( ladder2->lfrom - ladder1->lto < 6 ) &&
			//				    ( ladder2->from  - ladder1->to  < 3 )) ||
			//				   ( ladder2->lfrom - ladder1->lto < 3 )) );
			//			 break;
		} else {
			bulge = ( (l1.lfrom - l2.lto > 0) &&
					(((l1.lfrom -l2.lto < 6) &&
							( l2.from - l1.to <3  )) ||
							(l1.lfrom - l2.lto < 3))
					);
		}
		return bulge;


		//	       case antiparallel:
		//		 bulge = ( ( ladder1->lfrom - ladder2->lto > 0 ) &&
		//			   ((( ladder1->lfrom - ladder2->lto < 6 ) &&
		//			     ( ladder2->from - ladder1->to  < 3 )) ||
		//			    ( ladder1->lfrom - ladder2->lto < 3 ))   );
		//		 break;
		//	       }
		//	     }

	}

	private void registerBridge(int start, int end, BridgeType btype){
		if ( start > end)
			return;

		boolean found = false;
		for (Ladder ladder : ladders){
			if (shouldExtendLadder(ladder, start,end,btype)) {
				// extend laddder
				found = true;
				ladder.to++;
				if ( btype.equals(BridgeType.parallel)) {
					ladder.lto++;
				} else {
					ladder.lfrom--;
				}
				break;

			}
		}
		if ( ! found){
			// create new ladder!
			Ladder l = new Ladder();
			l.setFrom(start);
			l.setTo(end);
			l.setBtype(btype);
			l.setLfrom(start);
			l.setLto(end);
			ladders.add(l);
		}

	}

	private boolean shouldExtendLadder(Ladder ladder, int start, int end,
			BridgeType btype) {

		return ( (btype.equals(ladder.getBtype()) ) &&
				( start  == ladder.getTo() +1) &&
				(
						(( end == ladder.getLto() + 1) &&
								( btype.equals(BridgeType.parallel) )) ||
								(( end == ladder.getLfrom() - 1 ) &&
										(btype.equals(BridgeType.parallel) ))										
						) );



	}

	private void testBridge(int i) {

		int currpos = i + 3;
		for ( int foundNrBridges = 0 ; foundNrBridges < 2 && ( currpos < groups.length-1); currpos++){



			BridgeType btype = null;

			if ( (isBonded(i+1, currpos) && isBonded(currpos, i-1) ) ||
					( isBonded(currpos +1 , i) && isBonded(i, currpos-1)) ) {
				// parallel
				btype = BridgeType.parallel;
			}
			else if ( (isBonded(i, currpos) && isBonded(currpos, i)) ||
					(isBonded(i+1, currpos -1 ) && (isBonded(currpos + 1 , i - 1)))) {
				// antiparallel
				btype = BridgeType.antiparallel;
			}
			if (btype != null){
				foundNrBridges++;
				registerBridge(i, currpos, btype);
			}
		}

	}

	private void detectBends() {
		if ( groups.length < 5)
			return;

		for (int i =2 ; i < groups.length -2 ;i++){

			/* Create vectors ( Ca i to Ca i-2 ) ; ( Ca i to CA i +2 ) */

			SecStrucGroup im2 = groups[i-2];
			SecStrucGroup g = groups[i];
			SecStrucGroup ip2 = groups[i+2];


			Atom caim2 = im2.getCA();
			Atom cag   = g.getCA();
			Atom caip2 = ip2.getCA();

			try {
				Atom caminus2 = Calc.subtract(caim2,cag);

				Atom caplus2  = Calc.subtract(cag,caip2);

				double angle    = Calc.angle(caminus2,caplus2);

				SecStrucState state = getSecStrucState(i); 

				state.setKappa((float)angle);

				if (angle > 70.0) {
					if ( state.getSecStruc().equals(SecStrucType.coil)) 
						state.setSecStruc(SecStrucType.bend);

					state.setBend(true);
					//d[i].bend = 'S';
				}
			} catch (Exception e){
				e.printStackTrace();
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

			SecStrucState state1 = (SecStrucState) a.getProperty("secstruc");
			SecStrucState state2 = (SecStrucState) b.getProperty("secstruc");
			state2.setPhi(phi);
			state1.setPsi(psi);
			state1.setOmega(omega);

		}

	}

	/*
	 *
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		StringBuffer buf = new StringBuffer();
		String nl = System.getProperty("line.separator");
		buf.append("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI     PSI     Omega    X-CA   Y-CA   Z-CA ");
		buf.append(nl);

		for (int i =0 ; i < groups.length ;i++){
			Group g = groups[i];
			SecStrucState state = (SecStrucState) g.getProperty("secstruc");
			//	System.out.println("XX"+i+" "+g.getResidueNumber().toString() + " " + g.getPDBName() + iter.getCurrentChain().getName() + " " + state);
			buf.append(i + 1).append(" ");
			buf.append(g.getChainId()).append(" ");
			buf.append(g.getPDBName()).append(" ");
			buf.append(g.getResidueNumber().toString()).append("\t");

			boolean[] turns = state.getTurn();
			for (int t=0;t<3;t++){
				if (turns[t]) {
					buf.append('>');
				} else {
					buf.append(' ');
				}
			}

			// tmp filler
			buf.append(state.getSecStruc().type);
			buf.append(" ");
			if ( state.isBend())
				buf.append('S');
			else 
				buf.append(" ");
			buf.append("                 ");

			int p1 = state.getAccept1().getPartner();
			if ( p1 != 0)
				p1 -= i;
			double e1 =  (state.getAccept1().getEnergy() / 1000.0);
			buf.append(String.format( "%6d,%4.1f",p1,e1));

			int p2 = state.getDonor1().getPartner();
			if ( p2 != 0)
				p2 -= i;
			double e2 = (state.getDonor1().getEnergy() / 1000.0);
			buf.append(String.format( "%6d,%4.1f",p2,e2 ));

			int p3 = state.getAccept1().getPartner() ;
			if ( p3 != 0)
				p3 -= i;
			double e3 =  (state.getAccept2().getEnergy() / 1000.0);
			buf.append(String.format( "%6d,%4.1f",p3,e3));

			int p4 = state.getDonor2().getPartner();
			if ( p4 != 0)
				p4 -= i;
			double e4 = (state.getDonor2().getEnergy() / 1000.0);
			buf.append(String.format( "%6d,%4.1f",p4,e4 ));

			double kappa = state.getKappa();

			double phi = state.getPhi();
			double psi = state.getPsi();
			double omega = state.getOmega();
			buf.append("        ");
			buf.append(String.format("%6.1f",kappa));
			buf.append("       ");

			buf.append(String.format("%6.1f %6.1f %6.1f", phi,psi, omega));

			buf.append(nl);
		}

		return buf.toString();
	}

	private static SecStrucGroup[] initGroupArray(Structure s) {
		List<SecStrucGroup> groupList = new ArrayList<SecStrucGroup>();
		//GroupIterator iter = new GroupIterator(s);
		for ( Chain c : s.getChains()){

			for (Group g : c.getAtomGroups()){
				//System.out.println(g);
				//			 we can also calc secstruc if hetatom is a modified amino acid.
				if ( g.hasAminoAtoms()) {

					SecStrucGroup sg = new SecStrucGroup();
					sg.setResidueNumber(g.getResidueNumber());
					sg.setPDBFlag(true);
					try {
						sg.setPDBName(g.getPDBName());
					} catch (PDBParseException e){
						e.printStackTrace();
					}
					sg.setParent(g.getChain());

					try {

						sg.setN((Atom)   g.getAtomByPDBname(" N  ").clone());
						sg.setCA((Atom)  g.getAtomByPDBname(" CA ").clone());
						sg.setC((Atom)   g.getAtomByPDBname(" C  ").clone());
						sg.setO((Atom)   g.getAtomByPDBname(" O  ").clone());
						sg.setOriginal(g);
						// create H in calc_H


					} catch (StructureException e){
						e.printStackTrace();
						// one of these atoms is not found!
						continue;
					}

					SecStrucState state = new SecStrucState();
					Map<String,Object> m = sg.getProperties();
					if ( m == null) {
						m = new HashMap<String, Object>();
						sg.setProperties(m);
					}

					m.put("secstruc",state);


					groupList.add(sg);
				} else {
					//System.out.println("not an amino group");
				}
			}
		}

		return (SecStrucGroup[]) groupList.toArray(new SecStrucGroup[groupList.size()]);
	}



	/** calculate the coordinates for the H atoms. They are usually
	 * missing in the PDB files as only few experimental methods allow
	 * to resolve their location
	 */
	private  void calculateHAtoms()
			throws StructureException
			{

		for ( int i = 0 ; i < groups.length-1  ; i++) {

			SecStrucGroup a  = groups[i];
			SecStrucGroup b  = groups[i+1];

			if ( ! b.hasAtom("H") ) {

				//System.out.println(cur);
				// calculate the coordinate for the H atom
				//Atom H = calc_H(a.getC(), b.getN(), b.getCA());

				// alternative:
				Atom H = calcSimple_H(a.getC(), a.getO(),b.getN());

				b.setH(H);

				/*System.out.println("added H for " + i + " " + H);
            for ( int aa = 0 ; aa < b.size() ; aa++){
               Atom at = b.getAtom(aa);
               System.out.println(aa + " " + at.getFullName() + " "+ Calc.getDistance(at,H));
            }*/


			}

		}
			}

	/** calculate the HBonds between different groups ...
	 * see Creighton page 147 f
	 *
	 */
	private void calculateHBonds()
			throws StructureException
			{
		//System.out.println("groups length: " + groups.length);

		if ( groups.length < 5)
			return;

		// skip the first residue , unable to calc H for it ...
		for (int i=1 ; i < groups.length ;  i++){

			SecStrucGroup one = groups[i];

			if ( ! one.hasAtom("H")) {
				System.out.println(" no H at " + i);
				continue;
			}

			for ( int j = i+1 ; j < groups.length ; j++){

				SecStrucGroup two = groups[j];

				// check if distance is  too large.
				// if too big - for sure no HBonds ...
				double dist = Calc.getDistance(one.getCA(),two.getCA());

				// speed up...
				if ( dist >= CA_MIN_DIST  )
					continue;
				//System.out.println("calc " + i + " " + j + " "+  dist);

				checkAddHBond(i,j);

				// "backwards" hbonds are not allowed
				if ( j != (i+1) ) {

					checkAddHBond(j,i);

				}
				//System.out.println(" ");
			}
		}


			}


	private void checkAddHBond(int i, int j){
		SecStrucGroup one = groups[i];

		if (one.getPDBName().equals("PRO")){
			if (debug)
				System.out.println("     ignore: PRO " +     one.getResidueNumber().toString());
			return ;
		}

		SecStrucGroup two = groups[j];
		if ( ! two.hasAtom("H")) {
			System.err.println("two has no H " + j);
			return;
		}



		double energy = 0;
		try {
			energy = calculateHBondEnergy(one,two);
		} catch (Exception e){
			e.printStackTrace();
			return;
		}
		//System.out.println(" " + energy);

		trackHBondEnergy(i,j,energy);

	}

	/** calculate HBond energy of two groups in cal/mol ...
	 * see Creighton page 147 f
	 *
	 * Jeffrey, George A., An introduction to hydrogen bonding, Oxford University Press, 1997.
	 * categorizes hbonds with donor-acceptor distances of
	 * 2.2-2.5 &aring; as "strong, mostly covalent",
	 * 2.5-3.2 &aring; as "moderate, mostly electrostatic",
	 * 3.2-4.0 &aring; as "weak, electrostatic".
	 *  Energies are given as 40-14, 15-4, and <4 kcal/mol respectively.
	 *
	 */


	public  double calculateHBondEnergy(SecStrucGroup one, SecStrucGroup two )
			throws StructureException{


		//System.out.println("calcHBondEnergy" + one + "|" + two);

		Atom N = one.getN();
		Atom H = one.getH();

		Atom O = two.getO();
		Atom C = two.getC();

		double dno = Calc.getDistance(O,N);
		double dhc = Calc.getDistance(C,H);
		double dho = Calc.getDistance(O,H);
		double dnc = Calc.getDistance(C,N);

		if  ( debug ){

			System.out.println("     cccc: " + one.getResidueNumber().toString() + " " + one.getPDBName() + " " +two.getResidueNumber().toString() + " " + two.getPDBName() +
					String.format(" O ("+O.getPDBserial()+")..N ("+ N.getPDBserial()+"):%4.1f  |  ho:%4.1f - hc:%4.1f + nc:%4.1f - no:%4.1f " , dno,dho,dhc,dnc,dno));

		}
		//System.out.println( cn > ch && oh < 3.0f);

		double contact = MINDIST ;

		//		 there seems to be a contact!
		if ( (dno < contact) || (dhc < contact) || (dnc < contact) || (dno < contact)) {
			//System.out.println("!!! contact " + one + " " + two);
			return HBONDLOWENERGY ;
		}

		double e1 = Q / dho  - Q / dhc ;
		double e2 = Q / dnc  - Q / dno ;

		double energy = e1 + e2;

		if ( debug )
			System.out.println(String.format("      N (%d) O(%d): %4.1f : %4.2f ",N.getPDBserial(),O.getPDBserial(), (float)dno , energy));

		// bond too weak
		//if ( energy > HBONDHIGHENERGY)
		//	return energy;

		// test to avoid bond too strong
		if ( energy > HBONDLOWENERGY)
			return energy;


		return HBONDLOWENERGY ;

	}

	/**
	 * calculate distance between two atoms.
	 *
	 * @param a  an Atom object
	 * @param b  an Atom object
	 * @return a double
	 * @throws StructureException ...
	 */
	public static BigDecimal getPreciseDistance(Atom a, Atom b)
			throws StructureException
			{
		double x = a.getX() - b.getX();
		double y = a.getY() - b.getY();
		double z = a.getZ() - b.getZ();

		double s  = x * x  + y * y + z * z;
		BigSqrt sqrt = new BigSqrt();
		BigDecimal d = new BigDecimal(s);
		BigDecimal dist = sqrt.sqrt(d);

		return dist ;
			}



	/** store Hbonds inamino acids
	 * DSSP allows two HBonds / aminoacids to allow bifurcated bonds ...
	 */
	private  void trackHBondEnergy(int i, int j, double energy) {

		Group one = groups[i];
		Group two = groups[j];

		if ( one.getPDBName().equals("PRO"))
			return;

		SecStrucState stateOne = (SecStrucState) one.getProperty("secstruc");
		SecStrucState stateTwo = (SecStrucState) two.getProperty("secstruc");



		double acc1e = stateOne.getAccept1().getEnergy();
		double acc2e = stateOne.getAccept2().getEnergy();

		double don1e = stateTwo.getDonor1().getEnergy();
		double don2e = stateTwo.getDonor2().getEnergy();

		//if ( energy < 0)
		//   System.out.println("--- tracking Hbond " + i + " " + j + " " + energy + " (accept: " + acc1e + " " + acc2e +") (donor: "+don1e + " " + don2e+")");

		if (energy <  acc1e) {
			//System.out.println(energy +"<"+acc1e) ;
			stateOne.setAccept2(stateOne.getAccept1());

			HBond bond = new HBond();
			bond.setEnergy(energy);
			bond.setPartner(j);

			stateOne.setAccept1(bond);

		} else if ( energy < acc2e ) {
			//System.out.println(energy +"<"+acc2e) ;
			HBond bond = new HBond();
			bond.setEnergy(energy);
			bond.setPartner(j);

			stateOne.setAccept2(bond);
		}

		// and now the other side of the bond ..


		if (energy <  don1e) {
			stateTwo.setDonor2(stateTwo.getDonor1());

			HBond bond = new HBond();
			bond.setEnergy(energy);
			bond.setPartner(i);

			stateTwo.setDonor1(bond);

		} else if ( energy < don2e ) {

			//System.out.println(energy +"<"+don2e) ;

			HBond bond = new HBond();
			bond.setEnergy(energy);
			bond.setPartner(i);

			stateTwo.setDonor2(bond);

		}

		//System.out.println(stateOne);
		//one.setProperty("secstruc", stateOne);
		//two.setProperty("secstruc", stateTwo);
		/*groups[i] = one;
      groups[j] = two;
		 */
	}

	/** detect helical turn patterns
	 *
	 *
	 */
	private void calculateTurns(){

		int l = groups.length;
		for (int i = 0 ; i< l ; i++){

			for ( int turn = 3; turn <= 5 ; turn++ ) {
				if (i+turn >= l)
					continue;

				if ( isBonded(i+turn, i)){
					//System.out.println("is bondend " + (i+turn) + i );
					for ( int j=i;j<i+turn+1;j++){
						//System.out.println("turn at i:" + i + " j:" + j + " turn" + turn);
						SecStrucGroup group = groups[j];
						SecStrucState state = (SecStrucState) group.getProperty("secstruc");
						boolean[] turns = state.getTurn();
						turns[turn-3] = true;

					}
				}
			}
		}
	}


	/** test if two groups are forming an H-Bond
	 * DSSP defines H-Bonds if the energy < -500 cal /mol
	 * @param one group one
	 * @param two group two
	 * @return flag if the two are forming an Hbond
	 */

	private boolean isBonded(int i, int j) {

		Group one = groups[i];
		//Group two = groups[j];

		SecStrucState stateOne = (SecStrucState)one.getProperty("secstruc");
		//SecStrucState stateTwo = (SecStrucState)two.getProperty("secstruc");

		//System.out.println("*** bonded? " + i + " " + j + " " + stateOne);
		double acc1e    = stateOne.getAccept1().getEnergy();
		double acc2e    = stateOne.getAccept2().getEnergy();

		int partnerAcc1 = stateOne.getAccept1().getPartner();
		int partnerAcc2 = stateOne.getAccept2().getPartner();

		if (
				( ( partnerAcc1 == j ) && (acc1e < HBONDHIGHENERGY) )
				||
				( ( partnerAcc2 == j ) && (acc2e < HBONDHIGHENERGY) )
				) {
			//System.out.println("*** yes is bonded " + i + " " + j);
			return true ;
		}
		return false ;
	}



	/**
	 * Use unit vectors NC and NCalpha Add them. Calc unit vector and
	 * substract it from N.
	 * C coordinates are from amino acid i-1
	 * N, CA atoms from amino acid i
	 *
	 * see also:
	 * @link{http://openbioinformatics.blogspot.com/2009/08/how-to-calculate-h-atoms-for-nitrogens.html}
	 *
	 *
	 */
	@SuppressWarnings("unused")
	private static Atom calc_H(Atom C, Atom N, Atom CA)
			throws StructureException
			{



		Atom nc  = Calc.subtract(N,C);
		Atom nca = Calc.subtract(N,CA);

		Atom u_nc  = Calc.unitVector(nc)   ;
		Atom u_nca = Calc.unitVector(nca);

		Atom added = Calc.add(u_nc,u_nca);

		Atom U     = Calc.unitVector(added);

		// according to Creighton distance N-H is 1.03 +/- 0.02A
		Atom H = Calc.add(N,U);

		H.setName("H");
		H.setFullName(" H  ");
		// this atom does not have a pdbserial number ...
		return H;


			}


	private static Atom calcSimple_H(Atom c,Atom o, Atom n) throws StructureException{

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
		h.setFullName(" H  ");
		return h;


	}

	public SecStrucGroup[] getGroups(){
		return groups;
	}

	private void buildHelices(){



		// check for minimum size

		if ( groups.length < 5 )

			return;




		/** test if two groups are forming an H-Bond

		 * DSSP defines H-Bonds if the energy < -500 cal /mol

		 * @param one group one

		 * @param two group two

		 * @return flag if the two are forming an Hbond

		 */


		// Alpha Helix (i+4)



		checkSetHelix(4, SecStrucType.helix4);

		checkSetHelix(3, SecStrucType.helix3);

		checkSetHelix(5, SecStrucType.helix5);

		checkSetTurns();

	}

	private void checkSetTurns() {
		for (int i =0 ; i < groups.length -3 ;i++){

			Group g = groups[i];

			SecStrucState state = (SecStrucState) g.getProperty("secstruc");

			SecStrucType type = state.getSecStruc();

			if ( type.isHelixType() )
				continue;

			boolean[] turns = state.getTurn();
			for ( int t = 0 ; t < 3 ; t ++){
				if ( turns[t])
				{

					for ( int l = i+1 ; l < i+t+3; l++) {
						//System.out.println("turn: " + i + " " + type);
						if ( l >= groups.length)
							break;
						SecStrucType typel =getSecStrucType(l);
						if ( typel.equals(SecStrucType.coil))
							setSecStrucType(l, SecStrucType.turn);

					}
				}
			}
		}
	}

	private void checkSetHelix( int prange,SecStrucType type){

		int range = prange - 3;
		//System.out.println("set helix " + type + " " + prange + " " + range);
		for (int i =1 ; i < groups.length -range -1 ;i++){

			Group g = groups[i];

			SecStrucState state = (SecStrucState) g.getProperty("secstruc");


			Group prevG = groups[i-1];

			SecStrucState prevState = (SecStrucState) prevG.getProperty("secstruc");

			Group nextG = groups[i+1];

			SecStrucState nextState = (SecStrucState) nextG.getProperty("secstruc");

			boolean[] turns = state.getTurn();

			boolean[] pturns = prevState.getTurn();

			boolean[] nturns = nextState.getTurn();

			// DSSP sets helices one amino acid too short....

			if ( turns[range] && pturns[range] && nturns[range]) {

				//Check if no secstruc assigned

				boolean empty = true;

				for ( int curr =i; curr <= i+range ;curr++){
					//System.out.println("   " + i + " range: " + prange + " " + range);
					SecStrucType cstate = getSecStrucType(curr);

					if (  cstate.isHelixType()){

						empty = false;

						break;

					}

				}


				// none is assigned yet, set to helix type

				if ( empty ) {

					for ( int curr =i; curr <= i+range ;curr++){

						setSecStrucType(curr,type);

					}

				}

			}

		}

	}

	private void setSecStrucType(int pos, SecStrucType state){

		Group g = groups[pos];

		SecStrucState s= (SecStrucState) g.getProperty("secstruc");

		s.setSecStruc(state);

	}

	private SecStrucType getSecStrucType(int pos){

		SecStrucState s = getSecStrucState(pos);

		return s.getSecStruc();

	}


	private SecStrucState getSecStrucState(int pos){
		Group g = groups[pos];
		SecStrucState state = (SecStrucState) g.getProperty("secstruc");
		return state;
	}





}









