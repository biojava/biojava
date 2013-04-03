/*
 *                  BioJava development code
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
 * Created on May 21, 2006
 *
 */
package org.biojava.bio.structure.align;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.ce.GuiWrapper;
import org.biojava.bio.structure.align.helper.AlignTools;
import org.biojava.bio.structure.align.helper.JointFragments;
import org.biojava.bio.structure.align.pairwise.AlignmentProgressListener;
import org.biojava.bio.structure.align.pairwise.AltAligComparator;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.align.pairwise.FragmentJoiner;
import org.biojava.bio.structure.align.pairwise.FragmentPair;
import org.biojava.bio.structure.io.PDBFileParser;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.jama.Matrix;


/**
 * Perform a pairwise protein structure superimposition.
 *
 * <p>
 * The algorithm is a distance matrix based, rigid body protein structure superimposition.
 * It is based on a variation of the PSC++ algorithm provided by Peter Lackner
 * (Peter.Lackner@sbg.ac.at, personal communication) .
 * </p>
 *
 *
 *
 * <h2>Example</h2>
 *  <pre>
 *  public void run(){

		// first load two example structures
		{@link InputStream} inStream1 = this.getClass().getResourceAsStream("/files/5pti.pdb");
		{@link InputStream} inStream2 = this.getClass().getResourceAsStream("/files/1tap.pdb");

		{@link Structure} structure1 = null;
		{@link Structure} structure2 = null;

		{@link PDBFileParser} pdbpars = new {@link PDBFileParser}();
		try {
			structure1 = pdbpars.parsePDBFile(inStream1) ;
			structure2 = pdbpars.parsePDBFile(inStream2);

		} catch ({@link IOException} e) {
			e.printStackTrace();
			return;
		}

		// calculate structure superimposition for two complete structures
		{@link StructurePairAligner} aligner = new {@link StructurePairAligner}();


		try {
			// align the full 2 structures with default parameters.
			// see StructurePairAligner for more options and how to align
			// any set of Atoms
			aligner.align(structure1,structure2);

			{@link AlternativeAlignment}[] aligs = aligner.getAlignments();
			{@link AlternativeAlignment} a = aligs[0];
			System.out.println(a);

			//display the alignment in Jmol

			// first get an artificial structure for the alignment
			{@link Structure} artificial = a.getAlignedStructure(structure1, structure2);


			// and then send it to Jmol (only will work if Jmol is in the Classpath)

			BiojavaJmol jmol = new BiojavaJmol();
			jmol.setTitle(artificial.getName());
			jmol.setStructure(artificial);

			// color the two structures


			jmol.evalString("select *; backbone 0.4; wireframe off; spacefill off; " +
					"select not protein and not solvent; spacefill on;");
			jmol.evalString("select *"+"/1 ; color red; model 1; ");


			// now color the equivalent residues ...

			String[] pdb1 = a.getPDBresnum1();
			for (String res : pdb1 ){
				jmol.evalString("select " + res + "/1 ; backbone 0.6; color white;");
			}

			jmol.evalString("select *"+"/2; color blue; model 2;");
			String[] pdb2 = a.getPDBresnum2();
			for (String res :pdb2 ){
				jmol.evalString("select " + res + "/2 ; backbone 0.6; color yellow;");
			}


			// now show both models again.
			jmol.evalString("model 0;");

		} catch ({@link StructureException} e){
			e.printStackTrace();
		}
	}
 *  </pre>
 *
 *
 *
 * @author Andreas Prlic
 * @author Peter Lackner
 * @since 1.4
 * @version %I% %G%
 */
public class StructurePairAligner {

	AlternativeAlignment[] alts;
	Matrix distanceMatrix;
	StrucAligParameters params;
	FragmentPair[]  fragPairs;

	List<AlignmentProgressListener> listeners = new ArrayList<AlignmentProgressListener>();

	boolean debug = false;

	public StructurePairAligner() {
		super();
		params = StrucAligParameters.getDefaultParameters();
		reset();
		alts = new AlternativeAlignment[0];
		distanceMatrix = new Matrix(0,0);
	}

	  public void addProgressListener(AlignmentProgressListener li){
	     listeners.add(li);
	  }

	  public void clearListeners(){
	     listeners.clear();
	  }


	/** example usage of this class
	 *
	 * @param args
	 */
	public static void main(String[] args){
		try {


			// UPDATE THE FOLLOWING LINES TO MATCH YOUR SETUP

			PDBFileReader pdbr = new PDBFileReader();
			pdbr.setPath("/Users/andreas/WORK/PDB/");


			//String pdb1 = "1crl";
			//String pdb2 = "1ede";

			String pdb1 = "1buz";
			String pdb2 = "1ali";
			String outputfile = "/tmp/alig_"+pdb1+"_"+pdb2+".pdb";

			// NO NEED TO DO CHANGE ANYTHING BELOW HERE...

			StructurePairAligner sc = new StructurePairAligner();


			// step1 : read molecules

			System.out.println("aligning " + pdb1 + " vs. " + pdb2);

			Structure s1 = pdbr.getStructureById(pdb1);
			Structure s2 = pdbr.getStructureById(pdb2);

			// step 2 : do the calculations
			sc.align(s1,s2);


			AlternativeAlignment[] aligs = sc.getAlignments();

			//cluster similar results together
			ClusterAltAligs.cluster(aligs);


			// print the result:
			// the AlternativeAlignment object gives also access to rotation matrices / shift vectors.
			for (int i=0 ; i< aligs.length; i ++){
				AlternativeAlignment aa = aligs[i];
				System.out.println(aa);
			}



			// convert AlternativeAlignemnt 1 to PDB file, so it can be opened with a viewer (e.g. Jmol, Rasmol)

			if ( aligs.length > 0) {
				AlternativeAlignment aa1 =aligs[0];
				String pdbstr = aa1.toPDB(s1,s2);

				System.out.println("writing alignment to " + outputfile);
				FileOutputStream out= new FileOutputStream(outputfile);
				PrintStream p =  new PrintStream( out );

				p.println (pdbstr);

				p.close();
				out.close();
			}


			// display the alignment in Jmol
			// only will work if Jmol is in the Classpath

			if ( aligs.length > 0) {

				if (! GuiWrapper.isGuiModuleInstalled()){
					System.err.println("Could not find structure-gui modules in classpath, please install first!");
					return;
				}

				//AlternativeAlignment aa1 =aligs[0];

				// first get an artificial structure for the alignment
				//Structure artificial = aa1.getAlignedStructure(s1, s2);


				// and then send it to Jmol (only will work if Jmol is in the Classpath)
				
				//GuiWrapper.display(afpChain, ca1, ca2, hetatms1, nucs1, hetatms2, nucs2);
			}


		} catch (Exception e){
			e.printStackTrace();
		}

	}



	private void reset(){
		alts = new AlternativeAlignment[0];
		distanceMatrix = new Matrix(0,0);
		fragPairs = new FragmentPair[0];

	}


	/** get the results of step 1 - the FragmentPairs used for seeding the alignment
	 * @return a FragmentPair[] array
	 */

	public FragmentPair[] getFragmentPairs() {
		return fragPairs;
	}



	public void setFragmentPairs(FragmentPair[] fragPairs) {
		this.fragPairs = fragPairs;
	}


	/** return the alternative alignments that can be found for the two structures
	 *
	 * @return AlternativeAlignment[] array
	 */
	public AlternativeAlignment[] getAlignments() {
		return alts;
	}

	/** return the difference of distance matrix between the two structures
	 *
	 * @return a Matrix
	 */
	public Matrix getDistMat(){
		return distanceMatrix;
	}

	/** get the parameters.
	 *
	 * @return the Parameters.
	 */
	public StrucAligParameters getParams() {
		return params;
	}

	/** set the parameters to be used for the algorithm
	 *
	 * @param params the Parameter object
	 */
	public void setParams(StrucAligParameters params) {
		this.params = params;
	}


	/** check if debug mode is set on
	 *
	 * @return debug flag
	 */
	public boolean isDebug() {
		return debug;
	}

	/** set the debug flag
	 *
	 * @param debug flag
	 */
	public void setDebug(boolean debug) {
		this.debug = debug;
	}




	/** Calculate the alignment between the two full structures with default parameters
	 *
	 * @param s1
	 * @param s2
	 * @throws StructureException
	 */
	public void align(Structure s1, Structure s2)
	throws StructureException {

		align(s1,s2,params);
	}

	
	
	

	/** Calculate the alignment between the two full structures with user provided parameters
	 *
	 * @param s1
	 * @param s2
	 * @param params
	 * @throws StructureException
	 */
	public void align(Structure s1, Structure s2, StrucAligParameters params)
	throws StructureException {
		// step 1 convert the structures to Atom Arrays


		Atom[] ca1 = getAlignmentAtoms(s1);
		Atom[] ca2 = getAlignmentAtoms(s2);

		notifyStartingAlignment(s1.getName(),ca1,s2.getName(),ca2);
		align(ca1,ca2,params);
	}


	/** Align two chains from the structures. Uses default parameters.
	 * 
	 * @param s1
	 * @param chainId1
	 * @param s2
	 * @param chainId2
	 */
	public void align(Structure s1, String chainId1, Structure s2, String chainId2) 
	throws StructureException{
		align(s1,chainId1,s2,chainId2, params);
	}
	
	/** Aligns two chains from the structures using user provided parameters.
	 * 
	 * @param s1
	 * @param chainId1
	 * @param s2
	 * @param chainId2
	 * @param params
	 * @throws StructureException
	 */
	public void align(Structure s1, String chainId1, Structure s2, String chainId2, StrucAligParameters params)
	throws StructureException{
		reset();
		this.params = params;
		
		Chain c1 = s1.getChainByPDB(chainId1);
		Chain c2 = s2.getChainByPDB(chainId2);
		
		Structure s3 = new StructureImpl();
		s3.addChain(c1);
		
		Structure s4 = new StructureImpl();
		s4.addChain(c2);
		
		Atom[] ca1 = getAlignmentAtoms(s3);
		Atom[] ca2 = getAlignmentAtoms(s4); 
				
		notifyStartingAlignment(s1.getName(),ca1,s2.getName(),ca2);
		align(ca1,ca2,params);
	}
	
	/** Returns the atoms that are being used for the alignment. (E.g. Calpha only, etc.)
	 * 
	 * @param s
	 * @return an array of Atoms objects
	 */
	public  Atom[] getAlignmentAtoms(Structure s){
		String[] atomNames = params.getUsedAtomNames();
		return StructureTools.getAtomArray(s,atomNames);
	}

	/** calculate the  protein structure superimposition, between two sets of atoms.
	 *
	 *
	 *
	 * @param ca1 set of Atoms of structure 1
	 * @param ca2 set of Atoms of structure 2
	 * @param params the parameters to use for the alignment
	 * @throws StructureException
	 */
	public void align(Atom[] ca1, Atom[] ca2, StrucAligParameters params)
	throws StructureException {


		reset();
		this.params = params;
		
		long timeStart = System.currentTimeMillis();

//		step 1 get all Diagonals of length X that are similar between both structures
		if ( debug ) {
			System.out.println(" length atoms1:" + ca1.length);
			System.out.println(" length atoms2:" + ca2.length);

			System.out.println("step 1 - get fragments with similar intramolecular distances ");
		}
		int k  = params.getDiagonalDistance();
		int k2 = params.getDiagonalDistance2();
		int fragmentLength = params.getFragmentLength();

		if ( ca1.length < (fragmentLength + 1) )  {
			throw new StructureException("structure 1 too short ("+ca1.length+"), can not align");
		}
		if ( ca2.length < (fragmentLength + 1) ){
		   throw new StructureException("structure 2 too short ("+ca2.length+"), can not align");
		}
		int rows = ca1.length - fragmentLength + 1;
		int cols = ca2.length - fragmentLength + 1;
		//System.out.println("rows "  + rows + " " + cols +
		//      " ca1 l " + ca1.length + " ca2 l " + ca2.length);
		distanceMatrix = new Matrix(rows,cols,0.0);

		double[] dist1 = AlignTools.getDiagonalAtK(ca1, k);

		double[] dist2 = AlignTools.getDiagonalAtK(ca2, k);
		double[] dist3 = new double[0];
		double[] dist4 = new double[0];
		if ( k2 > 0) {
			dist3 = AlignTools.getDiagonalAtK(ca1, k2);
			dist4 = AlignTools.getDiagonalAtK(ca2, k2);
		}

		double[][] utmp = new double[][] {{0,0,1}};
		//Matrix unitv = new Matrix(utmp);
		Atom unitvector = new AtomImpl();
		unitvector.setCoords(utmp[0]);

		List<FragmentPair> fragments = new ArrayList<FragmentPair>();

		for ( int i = 0 ; i< rows; i++){

			Atom[] catmp1  = AlignTools.getFragment( ca1,  i, fragmentLength);
			Atom   center1 = AlignTools.getCenter( ca1, i, fragmentLength);

			for ( int j = 0 ; j < cols ; j++){

				double rdd1 = AlignTools.rms_dk_diag(dist1,dist2,i,j,fragmentLength,k);
				double rdd2 = 0;
				if ( k2 > 0)
					rdd2 = AlignTools.rms_dk_diag(dist3,dist4,i,j,fragmentLength,k2);
				double rdd = rdd1 + rdd2;
				distanceMatrix.set(i,j,rdd);


				if ( rdd < params.getFragmentMiniDistance()) {
					FragmentPair f = new FragmentPair(fragmentLength,i,j);
					//System.out.println("i " + i + " " + j );
					try {

						Atom[] catmp2 = AlignTools.getFragment(ca2, j, fragmentLength);
						Atom  center2 = AlignTools.getCenter(ca2,j,fragmentLength);

						//System.out.println("c1 : " + center1 + " c2: " + center2);
						f.setCenter1(center1);
						f.setCenter2(center2);

						SVDSuperimposer svd = new SVDSuperimposer(catmp1,catmp2);
						Matrix rotmat = svd.getRotation();
						//rotmat.print(3,3);
						f.setRot(rotmat);

						Atom aunitv = (Atom)unitvector.clone();
						Calc.rotate(aunitv,rotmat);
						f.setUnitv(aunitv);

						boolean doNotAdd = false;
						if ( params.reduceInitialFragments()) {
							doNotAdd = FragmentJoiner.reduceFragments(fragments,f, distanceMatrix);

						}
						if ( doNotAdd)
							continue;

						fragments.add(f);

					} catch (StructureException e){
						e.printStackTrace();
						break;
					}
				}
			}
		}

		notifyFragmentListeners(fragments);

		FragmentPair[] fp = (FragmentPair[]) fragments.toArray(new FragmentPair[fragments.size()]);
		setFragmentPairs(fp);

		if (  debug )
			System.out.println( " got # fragment pairs: " + fp.length);

		if ( debug )
			System.out.println("step 2 - join fragments");

		// step 2 combine them to possible models
		FragmentJoiner joiner = new FragmentJoiner();


		JointFragments[] frags;

		if ( params.isJoinFast()) {
		   // apply the quick alignment procedure.
		   // less quality in alignments, better for DB searches...
		   frags =  joiner.approach_ap3(ca1,ca2,fp,params);

		   joiner.extendFragments(ca1,ca2,frags,params);

		} else if ( params.isJoinPlo()){
		// this approach by StrComPy (peter lackner):
           frags =  joiner.frag_pairwise_compat(fp,
                   params.getAngleDiff(),
                   params.getFragCompat(),
                   params.getMaxrefine());

		} else {

		   // my first implementation
			frags =  joiner.approach_ap3(
					ca1,ca2, fp, params);
		}

		notifyJointFragments(frags);

		if ( debug )
			System.out.println(" number joint fragments:"+frags.length);
		
		if ( debug )
			System.out.println("step 3 - refine alignments");

		List<AlternativeAlignment> aas = new ArrayList<AlternativeAlignment>();
		for ( int i = 0 ; i < frags.length;i++){
			JointFragments f = frags[i];
			AlternativeAlignment a = new AlternativeAlignment();
			//System.out.println(f.getRms());
			//a.setRms(f.getRms());
			a.apairs_from_idxlst(f);
			a.setAltAligNumber(i+1);
			a.setDistanceMatrix(distanceMatrix);

			try {
				if ( params.getMaxIter() > 0 ){

					a.refine(params,ca1,ca2);
				}
				else {

				      a.finish(params,ca1,ca2);


				}
			} catch (StructureException e){
				e.printStackTrace();
			}
			a.calcScores(ca1,ca2);
			// a.getRotationMatrix().print(3,3);
			aas.add(a);
		}


		// sort the alternative alignments
		Comparator<AlternativeAlignment> comp = new AltAligComparator();
		Collections.sort(aas,comp);
		Collections.reverse(aas);

		alts = (AlternativeAlignment[])aas.toArray(new AlternativeAlignment[aas.size()]);
		// do final numbering of alternative solutions
		int aanbr = 0;
		for ( int i = 0 ; i < alts.length; i++){
			AlternativeAlignment a = alts[i];
			aanbr++;
			a.setAltAligNumber(aanbr);
			//System.out.println(aanbr);
			//a.getRotationMatrix().print(3,3);
		}
		//System.out.println("calc done");

		if (debug){
		   long timeEnd = System.currentTimeMillis();
		   System.out.println("total calculation time: "+ (timeEnd-timeStart) + " ms.");
		}
	}

	private void notifyStartingAlignment(String name1, Atom[] ca1, String name2, Atom[] ca2){
	   for (AlignmentProgressListener li : listeners){
	      li.startingAlignment(name1, ca1, name2, ca2);
	   }
	}

	private void notifyFragmentListeners(List<FragmentPair> fragments){

	   for (AlignmentProgressListener li : listeners){
	      li.calculatedFragmentPairs(fragments);
	   }

	}

	private void notifyJointFragments(JointFragments[] fragments){
	   for (AlignmentProgressListener li : listeners){
	      li.jointFragments(fragments);
	   }
	}

}
